"""Section 7.2: the cube method balanced on S4 on the elevation grid, with numerically sound balancing directions.

BalancedSampling 2.1.1 computes each balancing direction from a reduced row-echelon form whose pivots are found by
exact comparisons with zero (src/ReducedRowEchelonForm.cc), without a tolerance or partial pivoting.  With rows and
columns alone the elimination stays exact in floating point; with three or more line stratifications it does not,
rounding errors of order 1e-16 are taken as pivots, and the direction then violates the balancing equations,
including the sample size.  Each move is still a martingale move, so the inclusion probabilities are preserved,
but the size of the sample drifts (Section 7.2).  This script runs the transcription of `cube` in audit_cube.py
with one change: each direction is a null vector of the candidates' balancing matrix computed from a QR
factorisation (LAPACK), which satisfies the balancing equations to rounding error.  When that matrix has full rank
the direction is the package's direction up to scale and sign; otherwise both are valid null vectors.  The unit
order, the flight and landing phases (the landing drops the last balancing variables first, so the sample size is
kept), the random step and the integrality tolerance are those of the package.  For each draw it records the size,
the Horvitz-Thompson estimate (with the fixed denominator k), the balance on S2 and S4 and the main-effect weights,
and stores the selected units.  The function tolerant_rref (partial pivoting, entries at most 1e-9 times the largest
treated as zero) is the alternative used in the check of Appendix A; it gives the same sizes.
Usage: python3 cube_grid_fixed.py [B] [processes] [qr|rref]   Writes results/cube_grid_fixed.csv,
results/cube_grid_fixed_draws.csv and results/cube_samples/grid_cubeS4_fixed_k{61,183}.txt.gz.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import csv, gzip, sys, time
from multiprocessing import Pool
import numpy as np
import audit_cube as ac

OUT = os.path.dirname(os.path.abspath(__file__))
q = 61
N = q * q
TOLREL = 1e-9


def tolerant_rref(M):
    """Reduced row echelon form with partial pivoting; entries of absolute value at most TOLREL times the largest
    entry count as zero.  Same output conventions as audit_cube.rref (pivots exactly 1, pivot columns exactly 0)."""
    nrows, ncols = M.shape
    tol = TOLREL * max(np.abs(M).max(), 1.0)
    lead = 0
    for r in range(nrows):
        if ncols <= lead:
            return
        while True:
            col = np.abs(M[r:, lead])
            i = int(np.argmax(col))
            if col[i] > tol:
                i += r
                break
            M[r:, lead] = 0.0
            lead += 1
            if lead == ncols:
                return
        if i != r:
            M[[i, r]] = M[[r, i]]
        br = M[r]
        br[lead + 1:] /= br[lead]
        br[lead] = 1.0
        t = M[:, lead].copy()
        t[r] = 0.0
        M[:, lead] = 0.0
        M[r, lead] = 1.0
        rows = np.flatnonzero(t)
        if rows.size:
            M[rows, lead + 1:] -= np.outer(t[rows], br[lead + 1:])
        small = np.abs(M) <= tol
        M[small] = 0.0
        lead += 1


def qr_null_vector(M):
    """A unit null vector of M (nrows x nrows+1): the last column of Q in the complete QR factorisation of M'."""
    Q, _ = np.linalg.qr(M.T, mode="complete")
    return Q[:, -1].copy()


def use(method):
    if method == "qr":
        ac.rref = lambda M: None                 # the raw balancing matrix goes to qr_null_vector
        ac.null_vector = qr_null_vector
    else:
        ac.rref = tolerant_rref                  # the package's null vector, from a tolerant row reduction
        ac.null_vector = ORIGINAL_NULL_VECTOR


ORIGINAL_NULL_VECTOR = ac.null_vector
METHOD = sys.argv[3] if __name__ == "__main__" and len(sys.argv) > 3 else os.environ.get("CUBE_DIRECTIONS", "qr")
use(METHOD)
X = np.loadtxt(f"{OUT}/data/grid_X.csv", delimiter=",")          # rows, columns, anti-diagonals, diagonals: 4 x 60
y = np.loadtxt(f"{OUT}/data/grid_y.csv")
I, J = np.divmod(np.arange(N), q)
LINES = (I, J, (I + J) % q, (I - J) % q)


def tM(s, k, lines):
    return sum(q * (np.bincount(v[s], minlength=q) ** 2).sum() - len(s) ** 2 for v in lines) / k ** 2


def balanced(s, k, lines):
    return all(np.all(np.bincount(v[s], minlength=q) * q == k) for v in lines)


def one_draw(args):
    k, b, seed = args
    prob = np.full(N, k / N)
    x = np.c_[prob, X[:, :240]]
    t0 = time.time()
    d = ac.Cube(prob, x, np.random.default_rng(seed))
    d.flight(graph=False)
    s = d.landing()
    dev4 = max(np.abs(np.bincount(v[s], minlength=q) - k / q).max() for v in LINES)
    return {"k": k, "draw": b, "size": len(s), "ht": y[s].sum() / k, "mean": y[s].mean(),
            "balanced_S2": int(len(s) == k and balanced(s, k, LINES[:2])), "balanced_S4": int(len(s) == k and balanced(s, k, LINES)),
            "maxdev_S4": dev4, "tM_S2": tM(s, k, LINES[:2]), "tM_S4": tM(s, k, LINES), "seconds": time.time() - t0,
            "sample": " ".join(map(str, s))}


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    procs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    jobs = []
    for k in (61, 183):
        seeds = np.random.SeedSequence(5000 + k).spawn(B)
        jobs += [(k, b, seeds[b]) for b in range(B)]
    rows = []
    with Pool(procs) as pool:
        for n, r in enumerate(pool.imap_unordered(one_draw, jobs), 1):
            rows.append(r)
            if n % 20 == 0:
                print(f"{n}/{len(jobs)} draws, {np.mean([x['seconds'] for x in rows]):.1f} s per draw", flush=True)
    rows.sort(key=lambda r: (r["k"], r["draw"]))
    fields = [f for f in rows[0] if f != "sample"]
    with open(f"{OUT}/results/cube_grid_fixed_draws.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields + ["sample"])
        w.writeheader()
        for r in rows:
            w.writerow({a: (f"{v:.10g}" if isinstance(v, float) else v) for a, v in r.items()})
    summary = []
    Ybar = y.mean()
    for k in (61, 183):
        R = [r for r in rows if r["k"] == k]
        with gzip.open(f"{OUT}/results/cube_samples/grid_cubeS4_fixed_k{k}.txt.gz", "wt") as fh:
            for r in R:
                fh.write(f"{k} cubeS4fixed {r['sample']}\n")
        vsrs = (1 / k - 1 / N) * y.var(ddof=1)
        e_ht = np.array([(r["ht"] - Ybar) ** 2 for r in R])
        e_mean = np.array([(r["mean"] - Ybar) ** 2 for r in R])
        t4 = np.array([r["tM_S4"] for r in R]) / (N / k - 1)
        t2 = np.array([r["tM_S2"] for r in R]) / (N / k - 1)
        b2 = np.mean([r["balanced_S2"] for r in R]); b4 = np.mean([r["balanced_S4"] for r in R])
        summary.append({"k": k, "B": len(R), "size_exceptions": sum(r["size"] != k for r in R),
                        "deff_ht": e_ht.mean() / vsrs, "deff_ht_se": e_ht.std(ddof=1) / np.sqrt(len(R)) / vsrs,
                        "deff_mean": e_mean.mean() / vsrs, "deff_mean_se": e_mean.std(ddof=1) / np.sqrt(len(R)) / vsrs,
                        "share_exact_S2": b2, "share_exact_S4": b4, "mean_maxdev_S4": np.mean([r["maxdev_S4"] for r in R]),
                        "tau_S2": t2.mean(), "tau_S4": t4.mean(), "tau_S4_se": t4.std(ddof=1) / np.sqrt(len(R)),
                        "sec_per_draw": np.mean([r["seconds"] for r in R])})
        print(summary[-1], flush=True)
    with open(f"{OUT}/results/cube_grid_fixed.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        for r in summary:
            w.writerow({a: (f"{v:.6g}" if isinstance(v, float) else v) for a, v in r.items()})
