"""Table 1: minimax mixing distributions (linear program (8)) for complete factorial frames.

Balancing on all one-variable stratifications (order 1) or on all one- and two-variable
cross-classifications (order 2).  The LP is reduced by symmetry: under the monomial group
the coverage c_nu(u) of an orbit-symmetrised nu depends only on the Hamming weight of u.
"""
import csv
import itertools
import math
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import Space, subspaces, main_effect_set

OUT = __file__.rsplit("/", 1)[0]


def margin_factors(m, order):
    fac = []
    for o in range(1, order + 1):
        for T in itertools.combinations(range(m), o):
            B = np.zeros((o, m), dtype=np.int64)
            for a, i in enumerate(T):
                B[a, i] = 1
            fac.append(B)
    return fac


rows = []
for (p, m, order) in [(2, 5, 1), (2, 6, 1), (2, 7, 1), (2, 6, 2), (2, 7, 2), (3, 4, 1), (3, 5, 1), (3, 4, 2)]:
    sp = Space(p, m)
    M = main_effect_set(sp, margin_factors(m, order))
    wt = (sp.vecs != 0).sum(axis=1)
    nR = int((~M).sum() - 1)
    for h in range(order, m):
        types, nadm = {}, 0
        for D in subspaces(m, m - h, p):
            codes = sp.span(D)
            if M[codes].any():
                continue
            nadm += 1
            key = tuple(np.bincount(wt[codes], minlength=m + 1))
            types[key] = types.get(key, 0) + 1
        flat = (p ** (m - h) - 1) / nR
        row = {"q": p, "m": m, "balance": "one-way" if order == 1 else "one- and two-way",
               "k": p ** h, "n_admissible": nadm, "flat_bound": flat}
        if nadm == 0:
            row.update({"uniform_ratio": None, "minimax_ratio": None})
            rows.append(row)
            continue
        js = list(range(order + 1, m + 1))
        size = np.array([math.comb(m, j) * (p - 1) ** j for j in js], dtype=float)
        T = list(types)
        a = np.array([[t[j] for j in js] for t in T], dtype=float) / size
        cnt = np.array([types[t] for t in T], dtype=float)
        unif = (cnt[:, None] * a).sum(0) / cnt.sum()
        nT = len(T)
        res = linprog(np.r_[np.zeros(nT), 1.0], A_ub=np.c_[a.T, -np.ones(len(js))], b_ub=np.zeros(len(js)),
                      A_eq=np.r_[np.ones(nT), 0.0].reshape(1, -1), b_eq=[1.0],
                      bounds=[(0, None)] * nT + [(None, None)], method="highs")
        row.update({"uniform_ratio": unif.max() / flat, "minimax_ratio": res.fun / flat})
        rows.append(row)
with open(f"{OUT}/results/minimax_lp.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r)
