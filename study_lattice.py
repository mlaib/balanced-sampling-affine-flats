"""Study 1 addition: lattice sampling (Owen-randomised linear orthogonal array of maximal strength).

For q in {2, 3}, independent uniform permutations of the categories of each variable are the
affine maps x -> a x + b, so lattice sampling with a linear array A = c + V is the coset design
whose direction is uniform on {diag(a) V : a in (F_q^*)^m}.  We pick, among the directions
admissible for the one-variable margins, one with maximal dual distance (strength + 1) and, among
those, the fewest words of minimal weight; remaining ties go to the first direction in the
enumeration order of reduced row-echelon bases.  The chosen dual basis is saved with the results.
"""
import csv
import itertools
import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import Space, CosetDesign, admissible_directions, main_effect_set, srs_variance, rref_mod
from study_factorial import coordinate_factors, responses

OUT = __file__.rsplit("/", 1)[0]
rows = []
for (p, m, ks, seed) in [(2, 7, [8, 16, 32], 7), (3, 5, [9, 27, 81], 5)]:
    sp = Space(p, m)
    fac = coordinate_factors(m)
    M = main_effect_set(sp, fac)
    R = ~M
    R[0] = False
    wt = (sp.vecs != 0).sum(axis=1)
    ys = responses(sp, seed)
    for k in ks:
        h = round(math.log(k, p))
        A = admissible_directions(sp, fac, h)
        best, key_best = None, None
        for D in A:
            codes = sp.span(D)
            w_ = wt[codes[codes > 0]]
            dmin = w_.min()
            key = (-dmin, int((w_ == dmin).sum()))
            if key_best is None or key < key_best:
                best, key_best = D, key
        # orbit under coordinate scalings (dual direction D -> D diag(a^-1))
        orbit = {}
        for a in itertools.product(range(1, p), repeat=m):
            ainv = np.array([pow(int(x), -1, p) for x in a])
            Dn = (best * ainv) % p
            R_, _ = rref_mod(Dn, p)
            orbit[R_.tobytes()] = R_
        dirs = list(orbit.values())
        d = CosetDesign(sp, dirs, r=1)
        flat = (sp.N / k - 1) / R.sum()
        wmax = d.weights()[R].max()
        strength = -key_best[0] - 1
        row = {"p": p, "m": m, "k": k, "strength": strength, "n_directions": len(dirs),
               "worst_case_over_flat": wmax / flat, "dual_basis": str(best.tolist())}
        for yname, y in ys.items():
            row["deff_" + yname] = d.variance(y) / srs_variance(sp.N, k, y)
        rows.append(row)
        print({kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in row.items()})
with open(f"{OUT}/results/factorial_lattice.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
