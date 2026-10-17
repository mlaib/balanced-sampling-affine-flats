"""Ranks of the sample main-effects model used by the regression variance estimator (Section 5.3, Table 7).

varest_exact.py takes the rank of the sample design matrix X_s (an intercept and the indicators of levels
1..q-1 of each variable, 0/1 entries, unscaled) as the number of singular values above 1e-9 times the largest.
For every sample of every design evaluated with the regression estimator, this script
  - records the largest singular value counted as zero and the smallest counted as nonzero, relative to the
    largest, so that any threshold between the two gives the same ranks;
  - computes the rank exactly, by Gaussian elimination modulo the prime p = 2^31 - 1.  The rank over GF(p) never
    exceeds the rank over the rationals, and it equals it here: a nonzero minor of a 0/1 matrix of order n <= 11
    has absolute value at most (n+1)^{(n+1)/2}/2^n <= 1458 (Hadamard's bound), which is below p.
On the grid design balanced on rows the model is the one-way model by rows, and every sample meets every row,
so its rank is 61 in every sample; it is not recomputed.  Writes results/rank_check.json.
"""
import itertools, json, math, os
import numpy as np
from affsamp import Space, CosetDesign, admissible_directions, coset_labels
from study_factorial import coordinate_factors

OUT = os.path.dirname(os.path.abspath(__file__))
P = 2 ** 31 - 1


def rank_mod_p(X):
    """Ranks of a batch of integer matrices X (B x rows x cols) over GF(P), by fraction-free elimination."""
    A = X.astype(np.int64) % P
    B, n, c = A.shape
    used = np.zeros((B, n), dtype=bool)
    rank = np.zeros(B, dtype=np.int64)
    ar = np.arange(B)
    for j in range(c):
        cand = (A[:, :, j] != 0) & ~used
        has = cand.any(axis=1)
        piv = np.argmax(cand, axis=1)
        pr = A[ar, piv]                                   # pivot rows, B x c
        pv = pr[:, j]                                     # pivot values
        f = A[:, :, j].copy()                             # multipliers
        f[ar, piv] = 0
        # row_r <- pv * row_r - f_r * pivot_row (mod P) for the matrices with a pivot; scaling by pv != 0 keeps the rank
        newA = (pv[:, None, None] * A - f[:, :, None] * pr[:, None, :]) % P
        newA[ar, piv] = pr
        A = np.where(has[:, None, None], newA, A)
        used[ar[has], piv[has]] = True
        rank += has
    return rank


if __name__ == "__main__":
    out = []
    for (q, m, ks) in [(2, 7, (16, 32)), (3, 5, (27, 81))]:
        sp = Space(q, m)
        fac = coordinate_factors(m)
        lev = sp.vecs
        Xfull = np.ones((sp.N, 1))
        for i in range(m):
            for a in range(1, q):
                Xfull = np.c_[Xfull, (lev[:, i] == a).astype(float)]
        for k in ks:
            h = round(math.log(k, q))
            A = admissible_directions(sp, fac, h - 1)
            d = CosetDesign(sp, A, r=q)
            assert d.k == k
            combos = np.array(list(itertools.combinations(range(d.Q), q)))
            zero_max, nonzero_min, agree, ranks = 0.0, np.inf, True, {}
            for D in d.D:
                lab = coset_labels(sp, D)
                members = np.array([np.nonzero(lab == c)[0] for c in range(d.Q)])
                idx = members[combos].reshape(len(combos), -1)
                Xs = Xfull[idx]
                S = np.linalg.svd(Xs, compute_uv=False)
                rel = S / S[:, :1]
                keep = rel > 1e-9
                r_svd = keep.sum(axis=1)
                if (~keep).any():
                    zero_max = max(zero_max, float(rel[~keep].max()))
                nonzero_min = min(nonzero_min, float(rel[keep].min()))
                r_exact = rank_mod_p(Xs.astype(np.int64))
                agree &= bool((r_exact == r_svd).all())
                for v in r_exact:
                    ranks[int(v)] = ranks.get(int(v), 0) + 1
            row = {"frame": f"{q}^{m}", "k": k, "samples": int(sum(ranks.values())), "columns": int(Xfull.shape[1]),
                   "largest_relative_singular_value_counted_zero": zero_max,
                   "smallest_relative_singular_value_counted_nonzero": nonzero_min,
                   "exact_rank_equals_svd_rank_in_every_sample": agree,
                   "exact_rank_counts": {str(a): b for a, b in sorted(ranks.items())}}
            out.append(row)
            print(row, flush=True)
    json.dump(out, open(f"{OUT}/results/rank_check.json", "w"), indent=1)
