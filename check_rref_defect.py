"""Appendix A: the size exceptions of the cube method come from the row reduction of BalancedSampling 2.1.1.

The package finds the pivots of its reduced row-echelon form by exact comparisons with zero.  When the
elimination is not exact in floating point, rounding errors are taken as pivots and the computed direction
violates the balancing equations, the sample size included.  The transcription audit_cube.Cube reproduces the
package (validate_transcription.py); here it is run from the same seeds with the package's row reduction and with
the tolerant one of cube_grid_fixed.py (partial pivoting, entries at most 1e-9 times the largest treated as zero):
  - in the factorial frames of Section 7.1, flight and package landing, for every seed of cube_factorial.py whose
    stored draw has a size other than k, and for 300 other seeds of each configuration chosen at random;
  - on a 19 x 19 grid balanced on two, three and four line stratifications (20 draws per case).
check_size_drift.R counts the same for the package's cube and lcube on the 19 x 19 grid.
Usage: python3 check_rref_defect.py [processes]   Writes results/rref_defect_check.csv.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import csv, sys
from multiprocessing import Pool
import numpy as np
import audit_cube as ac
PACKAGE_RREF, PACKAGE_NULL = ac.rref, ac.null_vector       # before cube_grid_fixed replaces them
import cube_grid_fixed as cg

OUT = os.path.dirname(os.path.abspath(__file__))
_X = {}


def set_rref(name):
    ac.rref = PACKAGE_RREF if name == "package" else cg.tolerant_rref
    ac.null_vector = PACKAGE_NULL


def factorial_size(args):
    q, m, k, seed, name = args
    set_rref(name)
    if (q, m) not in _X:
        _X[(q, m)] = np.loadtxt(f"{OUT}/data/fact_p{q}_m{m}_X.csv", delimiter=",")
    X = _X[(q, m)]
    prob = np.full(len(X), k / len(X))
    c = ac.Cube(prob, np.c_[prob, X], np.random.default_rng(seed))
    c.flight(graph=False)
    return (q, m, k, name), len(c.landing())


def grid_size(args):
    q, nstrat, k, seed, name = args
    set_rref(name)
    N = q * q
    i, j = np.divmod(np.arange(N), q)
    lines = [i, j, (i + j) % q, (i - j) % q][:nstrat]
    ind = np.concatenate([np.stack([(v == b).astype(float) for b in range(1, q)], axis=1) for v in lines], axis=1)
    prob = np.full(N, k / N)
    c = ac.Cube(prob, np.c_[prob, ind], np.random.default_rng(seed))
    c.flight(graph=False)
    return (q, nstrat, k, name), len(c.landing())


def stored_sizes(q, m, ks, B=20000):
    """Sizes of the stored transcription draws with the package landing, by k and seed index (job order of
    cube_factorial.py: k, then seed, then the landings bs and continued)."""
    import gzip
    out = {k: np.zeros(B, dtype=int) for k in ks}
    with gzip.open(f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz", "rt") as fh:
        for n, line in enumerate(fh):
            k_i, rest = divmod(n, 2 * B)
            b, landing = divmod(rest, 2)
            f = line.split()
            assert int(f[0]) == ks[k_i] and f[1] == ("bs", "continued")[landing]
            if landing == 0:
                out[ks[k_i]][b] = len(f) - 2
    return out


if __name__ == "__main__":
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    rng = np.random.default_rng(2026)
    rows = []
    with Pool(procs) as pool:
        for q, m, ks in [(2, 7, (8, 16, 32)), (3, 5, (9, 27, 81))]:
            sizes = stored_sizes(q, m, ks)
            for k in ks:
                seeds = [int(x) for x in np.random.SeedSequence(3000 + k).generate_state(20000)]
                bad = np.flatnonzero(sizes[k] != k)
                control = rng.choice(np.setdiff1d(np.arange(20000), bad), size=300, replace=False)
                for group, idx in (("stored draws with size other than k", bad), ("300 other draws", control)):
                    for name in ("package", "tolerant"):
                        got = np.array([sz for _, sz in pool.imap(factorial_size, [(q, m, k, seeds[b], name) for b in idx], chunksize=20)])
                        if name == "package":
                            assert np.array_equal(got, sizes[k][idx]), (q, m, k, group)      # the transcription is reproduced
                        rows.append({"population": f"{q}^{m}", "stratifications": m, "k": k, "draws": group, "row_reduction": name,
                                     "n": len(idx), "size_other_than_k": int((got != k).sum()),
                                     "size_min": int(got.min()), "size_max": int(got.max())})
                        print(rows[-1], flush=True)
        for nstrat in (2, 3, 4):
            for k in (38, 76):
                for name in ("package", "tolerant"):
                    got = np.array([sz for _, sz in pool.imap(grid_size, [(19, nstrat, k, sd, name) for sd in range(20)])])
                    rows.append({"population": "19x19 grid", "stratifications": nstrat, "k": k, "draws": "seeds 0 to 19",
                                 "row_reduction": name, "n": 20, "size_other_than_k": int((got != k).sum()),
                                 "size_min": int(got.min()), "size_max": int(got.max())})
                    print(rows[-1], flush=True)
    with open(f"{OUT}/results/rref_defect_check.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
