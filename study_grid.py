"""Study 2: 61 x 61 elevation grid (R 'volcano', central 61 rows), units = F_61^2.

Exact variances of translation-invariant designs (spectral formula) and an exact
evaluation of the between-coset variance estimator for r-coset designs.
Writes the balancing indicators for the Monte Carlo runs (mc_grid.R).
"""
import itertools
import json
import math
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import (Space, CosetDesign, admissible_directions, main_effect_set,
                     decomposition, srs_variance, adjusted_r2, coset_labels, subspaces)
from scipy import stats

OUT = __file__.rsplit("/", 1)[0]
q = 61
sp = Space(q, 2)
Y = np.loadtxt(f"{OUT}/data/volcano61.csv", delimiter=",")
y = Y.reshape(-1)  # code(i, j) = 61 i + j

ROWS, COLS = np.array([[1, 0]]), np.array([[0, 1]])
ANTI, DIAG = np.array([[1, 1]]), np.array([[1, q - 1]])
STRAT = {"S1 rows": [ROWS], "S2 rows+cols": [ROWS, COLS],
         "S3 +anti-diagonals": [ROWS, COLS, ANTI], "S4 +diagonals": [ROWS, COLS, ANTI, DIAG]}


def centred(a):
    a = a % q
    return np.where(a > q // 2, a - q, a)


def shortest_dual(D):
    """Length of the shortest nonzero vector of the dual line D (wrapped coordinates)."""
    u = D[0]
    c = np.arange(1, q)[:, None]
    v = centred(c * u)
    return float(np.sqrt((v ** 2).sum(axis=1)).min())


rows = []
designs = {}
srs = {k: srs_variance(sp.N, k, y) for k in (61, 183)}
for sname, S in STRAT.items():
    EM, sR, dM, dR = decomposition(sp, S, y)
    A = admissible_directions(sp, S, 1)
    for r in (1, 3):
        d = CosetDesign(sp, A, r=r)
        designs[(sname, "flat", r)] = d
        rows.append({"strat": sname, "design": "coset flat", "r": r, "k": d.k, "n_dir": len(A),
                     "var": d.variance(y), "deff": d.variance(y) / srs[d.k],
                     "adjR2": adjusted_r2(sp, S, y), "R2": EM / (EM + sR), "dM": dM})
# lattice-quality ("spectral test") mixtures among the directions admissible for S4
A4 = admissible_directions(sp, STRAT["S4 +diagonals"], 1)   # dual lines D = span{(1, s)}, s = 2, ..., 59
quality = np.array([shortest_dual(D) for D in A4])
order = np.argsort(-quality, kind="stable")   # ties broken by increasing slope s (enumeration order)
with open(f"{OUT}/results/grid_lattice_ranking.csv", "w") as fh:
    fh.write("rank,slope_s,dual_basis,shortest_dual_vector_length,shortest_dual_vector_length_sq\n")
    for rank, i in enumerate(order, 1):
        fh.write(f"{rank},{A4[i][0][1]},\"{A4[i][0].tolist()}\",{quality[i]:.6f},{round(quality[i] ** 2)}\n")
tie_sensitivity = []
for L in (8, 16):   # every choice among the directions tied at the cut-off
    cut = quality[order[L - 1]]
    above = [i for i in order if quality[i] > cut + 1e-9]
    tied = [i for i in order if abs(quality[i] - cut) < 1e-9]
    vals = [CosetDesign(sp, [A4[i] for i in above + list(c)], r=1).variance(y) / srs[61]
            for c in itertools.combinations(tied, L - len(above))]
    tie_sensitivity.append({"L": L, "tied_slopes": [int(A4[i][0][1]) for i in tied], "n_choices": len(vals),
                            "deff_min": min(vals), "deff_max": max(vals)})
for L in (4, 8, 16):
    sel = order[:L]
    for r in (1, 3):
        d = CosetDesign(sp, [A4[i] for i in sel], r=r)
        designs[("S4 +diagonals", f"best{L}", r)] = d
        rows.append({"strat": "S4 +diagonals", "design": f"coset best-{L} lattices", "r": r, "k": d.k,
                     "n_dir": L, "var": d.variance(y), "deff": d.variance(y) / srs[d.k],
                     "min_shortest_dual": float(quality[sel].min()),
                     "pi2_all_positive": bool(np.all(d.pair_probability()[1:] > 0))})
for k in (61, 183):
    rows.append({"strat": "-", "design": "SRS", "r": None, "k": k, "var": srs[k], "deff": 1.0})

# exact evaluation of the variance estimator for r-coset designs (support enumeration)
vest = []
for (sname, kind, r), d in designs.items():
    if r < 2:
        continue
    combos = np.array(list(itertools.combinations(range(d.Q), r)))
    tq = stats.t.ppf(0.975, r - 1)
    Ev = Ev2 = cover = zero = 0.0
    for D, pr in zip(d.D, d.nu):
        lab = coset_labels(sp, D)
        cm = np.bincount(lab, weights=y, minlength=d.Q) / np.bincount(lab, minlength=d.Q)
        c = cm[combos]
        ybar = c.mean(axis=1)
        v = (1 - r / d.Q) * c.var(axis=1, ddof=1) / r
        Ev += pr * v.mean()
        Ev2 += pr * (v ** 2).mean()
        cover += pr * np.mean(np.abs(ybar - y.mean()) <= tq * np.sqrt(v))   # an interval with v = 0 is the point ybar
        zero += pr * np.mean(v == 0)
    V = d.variance(y)
    vest.append({"strat": sname, "design": kind, "r": r, "k": d.k, "true_var": V, "E_v": Ev,
                 "rel_bias": Ev / V - 1, "cv_v": math.sqrt(Ev2 - Ev ** 2) / Ev,
                 "coverage_t95": cover, "share_v_zero": zero})

# balancing indicators for the Monte Carlo runs (drop level 0 of each factor)
ind = []
for a in (ROWS, COLS, ANTI, DIAG):
    lev = (sp.vecs @ a.T % q).reshape(-1)
    for l in range(1, q):
        ind.append((lev == l).astype(np.int8))
np.savetxt(f"{OUT}/data/grid_X.csv", np.array(ind).T, fmt="%d", delimiter=",")
np.savetxt(f"{OUT}/data/grid_y.csv", y, delimiter=",")

with open(f"{OUT}/results/grid_exact.json", "w") as fh:
    json.dump({"designs": rows, "variance_estimation": vest, "tie_sensitivity": tie_sensitivity,
               "slopes_quality_top16": [(A4[i].tolist(), float(quality[i])) for i in order[:16]]}, fh, indent=1)
# residual Fourier energy under the wrapped (periodic) indexing, and the seams of the block
E = (np.abs(np.fft.fft2(Y) / Y.size) ** 2).reshape(-1)          # |yhat(u)|^2, code(a, b) = 61 a + b
Rmask = ~main_effect_set(sp, STRAT["S4 +diagonals"])
Rmask[0] = False
fa, fb = np.divmod(np.arange(q * q), q)
len2 = centred(fa) ** 2 + centred(fb) ** 2
ER = E[Rmask].sum()
spec = {"residual_frequencies": int(Rmask.sum()), "residual_energy": float(ER),
        "low": [{"max_length_sq": t, "n_freq": int((Rmask & (len2 <= t)).sum()),
                 "energy_share": float(E[Rmask & (len2 <= t)].sum() / ER)} for t in (5, 10, 18, 25, 49)],
        "best_L_lines": [],
        "seams_mean_abs_diff_m": {"row0_vs_row60": float(np.abs(Y[0] - Y[-1]).mean()),
                                  "adjacent_rows": float(np.abs(np.diff(Y, axis=0)).mean()),
                                  "col0_vs_col60": float(np.abs(Y[:, 0] - Y[:, -1]).mean()),
                                  "adjacent_cols": float(np.abs(np.diff(Y, axis=1)).mean())}}
for L in (4, 8, 16):
    U = np.zeros(q * q, dtype=bool)
    for i in order[:L]:
        U[sp.span(A4[i])] = True
    U[0] = False
    spec["best_L_lines"].append({"L": L, "n_freq": int(U.sum()), "freq_share": float(U.sum() / Rmask.sum()),
                                 "energy_share": float(E[U].sum() / ER), "min_length_sq": int(len2[U].min())})
with open(f"{OUT}/results/grid_spectrum.json", "w") as fh:
    json.dump(spec, fh, indent=1)
for r_ in rows:
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r_.items()})
for r_ in vest:
    print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r_.items()})
