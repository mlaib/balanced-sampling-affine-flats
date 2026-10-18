"""Preprocessing cost of the three ways of drawing a coset design (Section 4, implementation).

For each configuration of Table 1 and Section 5.1: time and peak memory to enumerate the admissible
directions, to solve the reduced program (8) over weight-distribution classes, and, for Section 5.1,
to solve the unreduced program over all admissible directions; the acceptance rate of rejection
sampling of a uniform admissible direction; and the size of the stored minimax support.
Writes results/preprocessing_costs.csv.
"""
import csv, math, os, platform, sys, time, tracemalloc
import numpy as np
from scipy.optimize import linprog
from affsamp import Space, subspaces, main_effect_set, gaussian_binomial
from study_minimax_lp import margin_factors
from study_factorial import minimax_design

OUT = os.path.dirname(os.path.abspath(__file__))


def measure(f):
    tracemalloc.start()
    t0 = time.perf_counter()
    res = f()
    t = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    return res, t, peak / 2 ** 20


rows = []
for (p, m, order, full) in [(2, 7, 1, True), (2, 7, 2, False), (3, 5, 1, True), (3, 4, 2, False)]:
    sp = Space(p, m)
    fac = margin_factors(m, order)
    M = main_effect_set(sp, fac)
    wt = (sp.vecs != 0).sum(axis=1)
    R = ~M
    R[0] = False
    for h in range(order, m):
        k = p ** h

        def enumerate_admissible():
            return [D for D in subspaces(m, m - h, p) if not M[sp.span(D)].any()]
        A, t_enum, mem_enum = measure(enumerate_admissible)
        total = gaussian_binomial(m, m - h, p)
        row = {"q": p, "m": m, "balanced_on": "one-way" if order == 1 else "one- and two-way", "k": k,
               "n_subspaces": total, "n_admissible": len(A), "acceptance_rate": len(A) / total,
               "enumeration_s": t_enum, "enumeration_MB": mem_enum}
        if A:
            def reduced_lp():
                types = {}
                for D in A:
                    key = tuple(np.bincount(wt[sp.span(D)], minlength=m + 1))
                    types[key] = types.get(key, 0) + 1
                js = list(range(order + 1, m + 1))
                size = np.array([math.comb(m, j) * (p - 1) ** j for j in js], dtype=float)
                a = np.array([[t[j] for j in js] for t in types], dtype=float) / size
                nT = len(types)
                return linprog(np.r_[np.zeros(nT), 1.0], A_ub=np.c_[a.T, -np.ones(len(js))], b_ub=np.zeros(len(js)),
                               A_eq=np.r_[np.ones(nT), 0.0].reshape(1, -1), b_eq=[1.0],
                               bounds=[(0, None)] * nT + [(None, None)], method="highs")
            _, t_red, mem_red = measure(reduced_lp)
            row.update({"reduced_lp_s": t_red, "reduced_lp_MB": mem_red})
            if full and h >= 2 and k in (8, 16, 32, 9, 27, 81):
                (des, _), t_full, mem_full = measure(lambda: minimax_design(sp, A, R, r=1))
                row.update({"unreduced_lp_s": t_full, "unreduced_lp_MB": mem_full,
                            "minimax_support": int((des.nu > 0).sum())})
        rows.append(row)
        print({a: (round(b, 4) if isinstance(b, float) else b) for a, b in row.items()}, flush=True)

keys = []
for r in rows:
    keys += [k_ for k_ in r if k_ not in keys]
with open(f"{OUT}/results/preprocessing_costs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)
print("machine:", platform.processor() or platform.machine(), "| python", sys.version.split()[0])
