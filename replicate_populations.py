"""Variability of the factorial comparisons of Section 7.1 across generated populations (Appendix A).

The populations of Section 7.1 use seeds 7 (2^7) and 5 (3^5).  Here the same generator (study_factorial
.responses) is run with seeds 1000, ..., 1000 + P - 1 for each frame.  The designs do not depend on the
response, so they are held fixed: the uniform and minimax coset designs (the stored minimax mixtures), the
replicated design C(h-1, nu, q) and lattice sampling with the arrays of Table 2 are evaluated exactly for
every population; the cube method of BalancedSampling and the transcription with the flight continued to a
vertex are evaluated on their stored samples (20 000 draws per configuration, the same draws for every
population).  Writes results/replicated_populations.csv (quantiles of design effects and of ratios to the
minimax coset design) and results/replicated_populations.npz (all design effects).
Usage: python3 replicate_populations.py [P]
"""
import csv, gzip, itertools, json, math, os, sys
import numpy as np
from affsamp import Space, CosetDesign, admissible_directions, rref_mod
from study_factorial import coordinate_factors, responses

OUT = os.path.dirname(os.path.abspath(__file__))
FRAMES = [(2, 7, (8, 16, 32)), (3, 5, (9, 27, 81))]
QS = (0.05, 0.25, 0.5, 0.75, 0.95)


def lattice_design(sp, basis):
    """Lattice sampling with the linear array whose dual basis is `basis`, randomised by level permutations,
    i.e. the coset design with direction uniform on the images under coordinate scalings (q <= 3)."""
    p, m = sp.p, sp.m
    best = np.array(basis, dtype=np.int64)
    orbit = {}
    for a in itertools.product(range(1, p), repeat=m):
        ainv = np.array([pow(int(x), -1, p) for x in a])
        R_, _ = rref_mod((best * ainv) % p, p)
        orbit[R_.tobytes()] = R_
    return CosetDesign(sp, list(orbit.values()), r=1)


def read_samples(path, N, select=None):
    """Indicator matrices (draws x N) by sample size from a file of 'k [landing] units...' lines."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if select is not None:
                if f[1] != select:
                    continue
                k, units = int(f[0]), f[2:]
            else:
                k, units = int(f[0]), f[1:]
            out.setdefault(k, []).append(np.array(units, dtype=np.int64))
    mats = {}
    for k, lst in out.items():
        S = np.zeros((len(lst), N))
        for b, s in enumerate(lst):
            S[b, s] = 1.0 / k
        mats[k] = S
    return mats


if __name__ == "__main__":
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    mixtures = json.load(open(f"{OUT}/results/factorial_minimax_mixtures.json"))
    lattice = {(int(r["p"]), int(r["k"])): json.loads(r["dual_basis"]) for r in csv.DictReader(open(f"{OUT}/results/factorial_lattice.csv"))}
    rows, store = [], {}
    for q, m, ks in FRAMES:
        sp = Space(q, m)
        fac = coordinate_factors(m)
        N = sp.N
        pops = [responses(sp, 1000 + i) for i in range(P)]
        names = list(pops[0])
        Yall = np.array([[pop[n] for n in names] for pop in pops])            # P x responses x N
        power = np.array([[sp.fourier_power(y) for y in pr] for pr in Yall])  # P x responses x N
        tot = power[:, :, 1:].sum(axis=2)                                     # sigma_y^2
        cube = read_samples(f"{OUT}/results/cube_samples/factorial_p{q}_m{m}.txt.gz", N)
        tpath = f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz"
        cont = read_samples(tpath, N, select="continued") if os.path.exists(tpath) else {}
        for k in ks:
            h = round(math.log(k, q))
            A = admissible_directions(sp, fac, h)
            mix = next(x for x in mixtures if x["p"] == q and x["k"] == k)
            des = {"coset-uniform": CosetDesign(sp, A, r=1),
                   "coset-minimax": CosetDesign(sp, [np.array(s["dual_basis"]) for s in mix["support"]],
                                                [s["probability"] for s in mix["support"]], r=1),
                   "coset-replicated": CosetDesign(sp, admissible_directions(sp, fac, h - 1), r=q),
                   "lattice": lattice_design(sp, lattice[(q, k)])}
            wsrs = (N / k - 1) / (N - 1)
            deff = {}
            for name, d in des.items():
                assert d.k == k
                deff[name] = (power[:, :, 1:] @ d.weights()[1:]) / (wsrs * tot)
            ybar = Yall.mean(axis=2)                                              # P x responses
            vsrs = (1 / k - 1 / N) * Yall.var(axis=2, ddof=1)
            YY = Yall.reshape(-1, N).T                                           # N x (P responses)
            for name, mats in (("cube", cube), ("cube-continued", cont)):
                if k in mats:
                    S, acc = mats[k], 0.0
                    for c0 in range(0, len(S), 2000):                             # squared errors, in chunks of draws
                        acc = acc + ((S[c0:c0 + 2000] @ YY - ybar.reshape(-1)) ** 2).sum(axis=0)
                    deff[name] = (acc / len(S)).reshape(ybar.shape) / vsrs
            for name, v in deff.items():
                store[f"p{q}_k{k}_{name}"] = v
                for j, yname in enumerate(names):
                    x = v[:, j]
                    row = {"q": q, "m": m, "k": k, "response": yname, "design": name, "statistic": "deff", "P": P,
                           "mean": x.mean()}
                    row.update({f"q{int(100 * a):02d}": np.quantile(x, a) for a in QS})
                    rows.append(row)
                    if name != "coset-minimax":
                        r_ = x / deff["coset-minimax"][:, j]
                        row = {"q": q, "m": m, "k": k, "response": yname, "design": name, "statistic": "ratio_to_minimax",
                               "P": P, "mean": r_.mean(), "share_below_1": (r_ < 1).mean()}
                        row.update({f"q{int(100 * a):02d}": np.quantile(r_, a) for a in QS})
                        rows.append(row)
        print(f"{q}^{m} done", file=sys.stderr, flush=True)
    keys = ["q", "m", "k", "response", "design", "statistic", "P", "mean", "q05", "q25", "q50", "q75", "q95", "share_below_1"]
    with open(f"{OUT}/results/replicated_populations.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k_: (f"{v:.6g}" if isinstance(v, float) else v) for k_, v in r.items()})
    np.savez_compressed(f"{OUT}/results/replicated_populations.npz", **store)
    for r in rows:
        if r["statistic"] == "ratio_to_minimax" and r["design"] in ("lattice", "coset-uniform", "cube"):
            print(r["q"], r["k"], r["response"], r["design"], {a: round(r[a], 3) for a in ("q05", "q50", "q95")})
