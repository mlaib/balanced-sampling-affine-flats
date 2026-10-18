"""Cube method with the flight continued to a vertex in the complete factorial frames of Section 5.1.

Uses the Python transcription of `cube` in BalancedSampling 2.1.1 (audit_cube.Cube, with the null vectors of
the reduced row-echelon form, as in the package) on the balancing matrix of mc_cube.R: the inclusion
probabilities and the indicators of levels 1, ..., q-1 of every variable.  For each draw the flight phase
is run once, and from its end point two landings are applied, as in audit_cube.py:

  bs         the landing of BalancedSampling (suppression of the last balancing variables);
  continued  the flight continued while a move preserving all balancing equations exists, followed by the
             same landing for any units still fractional (needed when the vertex reached is not integral).

With seven (2^7) or five (3^5) stratifications the constraint matrix need not be totally unimodular, so the
continued flight may end at a fractional vertex; the share of draws in which it ended at an integral one
is recorded.  Writes results/cube_factorial.csv (summary), results/cube_samples/transcription_*.txt.gz
(selected units of every draw and landing) and prints the summary.
Usage: python3 cube_factorial.py [B] [processes]
"""
import csv, gzip, os, sys, time
from multiprocessing import Pool
import numpy as np
import audit_cube as ac

OUT = os.path.dirname(os.path.abspath(__file__))
FRAMES = [(2, 7, (8, 16, 32)), (3, 5, (9, 27, 81))]
_D = {}


def frame(q, m):
    if (q, m) not in _D:
        X = np.loadtxt(f"{OUT}/data/fact_p{q}_m{m}_X.csv", delimiter=",")
        with open(f"{OUT}/data/fact_p{q}_m{m}_Y.csv") as fh:
            names = fh.readline().strip().split(",")
        Y = np.loadtxt(f"{OUT}/data/fact_p{q}_m{m}_Y.csv", delimiter=",", skiprows=1)
        _D[(q, m)] = (X, Y, names)
    return _D[(q, m)]


def level_counts(X, q, s):
    c = X[s].sum(axis=0).reshape(-1, q - 1)
    return np.c_[len(s) - c.sum(axis=1), c]          # one row per variable, levels 0, ..., q-1


def one_draw(args):
    q, m, k, seed = args
    X, Y, _ = frame(q, m)
    N = len(X)
    prob = np.full(N, k / N)
    c = ac.Cube(prob, np.c_[prob, X], np.random.default_rng(seed))
    c.flight(graph=False)
    out = []
    for name in ("bs", "continued"):
        d = c.copy()
        steps = d.continue_flight() if name == "continued" else 0
        frac = len(d.fractional())
        s = d.landing()
        cnt = level_counts(X, q, s)
        dev = cnt - k / q
        out.append({"q": q, "m": m, "k": k, "seed": seed, "landing": name, "continuation_steps": steps,
                    "fractional_before_landing": frac, "size": len(s), "exact": bool(np.abs(dev).max() < 1e-9),
                    "tM": float(q / k ** 2 * (dev ** 2).sum()), "means": Y[s].mean(axis=0), "sample": s})
    return out


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    os.makedirs(f"{OUT}/results/cube_samples", exist_ok=True)
    summary = []
    for q, m, ks in FRAMES:
        X, Y, names = frame(q, m)
        N, Ybar = len(X), Y.mean(axis=0)
        jobs = [(q, m, k, int(s)) for k in ks for s in np.random.SeedSequence(3000 + k).generate_state(B)]
        t0 = time.time()
        with Pool(procs) as pool:
            draws = [r for res in pool.imap(one_draw, jobs, chunksize=50) for r in res]
        print(f"{q}^{m}: {len(jobs)} draws in {time.time() - t0:.0f} s", file=sys.stderr, flush=True)
        with gzip.open(f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz", "wt") as fh:
            for r in draws:
                fh.write(f"{r['k']} {r['landing']} " + " ".join(map(str, r["sample"])) + "\n")
        for k in ks:
            srs = (1 / k - 1 / N) * Y.var(axis=0, ddof=1)
            for name in ("bs", "continued"):
                d = [r for r in draws if r["k"] == k and r["landing"] == name]
                err2 = np.array([(r["means"] - Ybar) ** 2 for r in d])
                ex = np.array([r["exact"] for r in d])
                tm = np.array([r["tM"] for r in d])
                integral = np.array([r["fractional_before_landing"] == 0 for r in d])
                for j, yname in enumerate(names):
                    summary.append({"q": q, "m": m, "k": k, "landing": name, "response": yname, "B": len(d),
                                    "deff": err2[:, j].mean() / srs[j], "deff_se": err2[:, j].std(ddof=1) / np.sqrt(len(d)) / srs[j],
                                    "share_exact": ex.mean(), "share_exact_se": np.sqrt(ex.mean() * (1 - ex.mean()) / len(d)),
                                    "share_integral_vertex": integral.mean() if name == "continued" else "",
                                    "tau": tm.mean() / (N / k - 1), "tau_se": tm.std(ddof=1) / np.sqrt(len(d)) / (N / k - 1),
                                    "size_min": min(r["size"] for r in d), "size_max": max(r["size"] for r in d),
                                    "bias_over_se": (np.array([r["means"][j] for r in d]).mean() - Ybar[j])
                                                    / (np.array([r["means"][j] for r in d]).std(ddof=1) / np.sqrt(len(d)))})
    with open(f"{OUT}/results/cube_factorial.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader(); w.writerows(summary)
    for r in summary:
        if r["response"] in ("additive", "higher-order"):
            print({a: (round(b, 4) if isinstance(b, float) else b) for a, b in r.items()})
