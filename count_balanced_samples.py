"""Exhaustive counts quoted in Sections 1 and 6 (Appendix A, exhaustive counts).

1. Example 2 of Tille and Rivest: the polytope {u in [0,1]^N : A u = b}, with A the sample-size row and
   the indicators of the three binary stratifications (rank 4) and b = A pi for pi = 1/2, for N = 8
   (one unit per cell of 2^3) and N = 16 (every cell duplicated; duplicates are distinct coordinates).
   Every vertex has a basis B of four columns with the other coordinates at 0 or 1: all bases and all
   0/1 settings are enumerated, the solutions with 0 <= u_B <= 1 kept and duplicates removed.  A vertex
   is integral when every coordinate is within 1e-9 of 0 or 1; the 0/1 solutions are also counted
   directly over all samples of the right size.
2. Section 6: samples of size q of the q x q grid balanced on rows, columns and anti-diagonals i+j.
   They are the permutations s of Z_q for which i -> i + s(i) is also a permutation (complete mappings),
   listed by backtracking; flats are s(i) = a i + b with a not in {0, -1}.  For the uniform design over
   them, the weights w(u) = E|eta_s(u)|^2 are computed from the two-dimensional discrete Fourier
   transform of each sample and compared with the equal-weight value 1/(q-2) on R_S.
Writes results/exhaustive_counts.json.
"""
import itertools, json, os
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))


def polytope_vertices(units, pi=0.5):
    n = len(units)
    A = np.vstack([np.ones(n), np.asarray(units, dtype=float).T])
    b = A @ np.full(n, pi)
    rank = np.linalg.matrix_rank(A)
    verts = set()
    for B in itertools.combinations(range(n), rank):
        AB = A[:, B]
        if abs(np.linalg.det(AB)) < 1e-9:
            continue
        rest = [i for i in range(n) if i not in B]
        Z = np.array(list(itertools.product((0.0, 1.0), repeat=len(rest))))
        uB = np.linalg.solve(AB, b[:, None] - A[:, rest] @ Z.T).T
        ok = np.all((uB > -1e-9) & (uB < 1 + 1e-9), axis=1)
        for z, ub in zip(Z[ok], uB[ok]):
            u = np.empty(n)
            u[rest] = z
            u[list(B)] = ub
            verts.add(tuple(np.round(u, 9)))
    integral = sum(all(abs(v) < 1e-9 or abs(v - 1) < 1e-9 for v in u) for u in verts)
    k = int(round(pi * n))
    target = np.asarray(units).sum(axis=0) * pi
    balanced = sum(1 for s in itertools.combinations(range(n), k)
                   if np.array_equal(np.asarray(units)[list(s)].sum(axis=0), target))
    return {"N": n, "rank": int(rank), "vertices": len(verts), "integral_vertices": int(integral),
            "balanced_samples_by_direct_count": int(balanced)}


def complete_mappings(q):
    sols, s, used_col, used_anti = [], [0] * q, [False] * q, [False] * q

    def rec(i):
        if i == q:
            sols.append(tuple(s))
            return
        for j in range(q):
            a = (i + j) % q
            if not used_col[j] and not used_anti[a]:
                used_col[j] = used_anti[a] = True
                s[i] = j
                rec(i + 1)
                used_col[j] = used_anti[a] = False
    rec(0)
    return sols


def grid_counts(q):
    sols = complete_mappings(q)
    flats = sum(1 for s in sols if any(all(s[i] == (a * i + s[0]) % q for i in range(q))
                                       for a in range(1, q - 1)))
    w = np.zeros((q, q))
    for s in sols:
        ind = np.zeros((q, q))
        ind[np.arange(q), list(s)] = 1.0
        w += np.abs(np.fft.fft2(ind) / q) ** 2
    w /= len(sols)
    u1, u2 = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    main = (u2 == 0) | (u1 == 0) | (u1 == u2)        # rows, columns and anti-diagonals i + j
    R = ~main
    vals, counts = np.unique(np.round(w[R], 6), return_counts=True)
    eq = 1.0 / (q - 2)
    return {"q": q, "balanced_samples": len(sols), "flats": flats, "residual_frequencies": int(R.sum()),
            "weights_on_R": {f"{v:.6f}": int(c) for v, c in zip(vals, counts)},
            "equal_weight_value": eq, "worst_case_ratio": float(w[R].max() / eq),
            "max_weight_on_main_effects": float(w[main & ~((u1 == 0) & (u2 == 0))].max())}


if __name__ == "__main__":
    cells = list(itertools.product((0, 1), repeat=3))
    res = {"example2": [polytope_vertices(cells), polytope_vertices([c for c in cells for _ in range(2)])],
           "grid_rows_columns_antidiagonals": [grid_counts(q) for q in (5, 7, 11)]}
    with open(f"{OUT}/results/exhaustive_counts.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
