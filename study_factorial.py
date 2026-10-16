"""Study 1: full factorial populations F_p^m balanced on all one-way margins.

Exact design variances (spectral formula) for SRS and coset designs; writes the
population and balancing matrix for the Monte Carlo cube-method runs (mc_cube.R).
"""
import itertools
import json
import math
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import (Space, CosetDesign, admissible_directions, main_effect_set,
                     decomposition, srs_variance, adjusted_r2)

OUT = __file__.rsplit("/", 1)[0]


def coordinate_factors(m):
    fac = []
    for i in range(m):
        B = np.zeros((1, m), dtype=np.int64)
        B[0, i] = 1
        fac.append(B)
    return fac


def project_out(sp, mask, y):
    """Remove the Fourier components of y on the frequencies in mask."""
    Y = np.fft.fftn(y.reshape((sp.p,) * sp.m)).reshape(-1)
    Y[mask] = 0
    return np.fft.ifftn(Y.reshape((sp.p,) * sp.m)).real.reshape(-1)


def responses(sp, seed):
    """Four responses with different shares of additive (main-effect) variance."""
    rng = np.random.default_rng(seed)
    X = sp.vecs.astype(float) - (sp.p - 1) / 2  # centred level codes
    m = sp.m
    main = X @ (1.0 / np.arange(1, m + 1))
    inter = np.zeros(sp.N)
    for i, j in itertools.combinations(range(m), 2):
        inter += rng.normal() * X[:, i] * X[:, j]
    for i in range(m):  # quadratic main effects for p > 2 are additive as well
        if sp.p > 2:
            main += 0.3 * rng.normal() * (X[:, i] ** 2 - (X[:, i] ** 2).mean())
    noise = rng.normal(size=sp.N)
    M = main_effect_set(sp, coordinate_factors(m))
    inter = project_out(sp, M, inter)  # pure interaction: no main-effect component
    inter /= inter.std()
    main /= main.std()
    noise_res = project_out(sp, M, noise)
    return {
        "additive": main + 0.3 * noise,
        "mixed": main + 0.8 * inter + 0.3 * noise,
        "interaction": inter + 0.3 * noise_res,
        "noise": noise,
    }


def minimax_design(sp, A, Rmask, r=1):
    """LP: nu over admissible directions minimising max_{u in R} Pr(u in D)."""
    Rcodes = np.nonzero(Rmask)[0]
    pos = {c: i for i, c in enumerate(Rcodes)}
    C = np.zeros((len(Rcodes), len(A)))
    for j, D in enumerate(A):
        for c in sp.span(D):
            if c in pos:
                C[pos[c], j] = 1.0
    nA = len(A)
    res = linprog(np.r_[np.zeros(nA), 1.0],
                  A_ub=np.c_[C, -np.ones(len(Rcodes))], b_ub=np.zeros(len(Rcodes)),
                  A_eq=np.r_[np.ones(nA), 0.0].reshape(1, -1), b_eq=[1.0],
                  bounds=[(0, None)] * nA + [(None, None)], method="highs")
    nu = np.clip(res.x[:nA], 0, None)
    return CosetDesign(sp, A, nu, r=r), res.fun


def run(p, m, ks, seed):
    sp = Space(p, m)
    fac = coordinate_factors(m)
    M = main_effect_set(sp, fac)
    R = ~M
    R[0] = False
    ys = responses(sp, seed)
    rows = []
    designs_info = []
    mixtures = []
    for k in ks:
        h = round(math.log(k, p))
        A = admissible_directions(sp, fac, h)
        uni = CosetDesign(sp, A, r=1)
        mm, lpval = minimax_design(sp, A, R, r=1)
        mixtures.append({"p": p, "m": m, "k": k, "lp_value": lpval, "bound": (sp.N / k - 1) / int(R.sum()),
                         "support": [{"dual_basis": D.tolist(), "probability": float(w)}
                                     for D, w in zip(mm.D, mm.nu) if w > 0]})
        des = {"coset-uniform": uni, "coset-minimax": mm}
        if h - 1 >= 1:
            A2 = admissible_directions(sp, fac, h - 1)
            if A2 and p == 2:
                des["coset-uniform r=2"] = CosetDesign(sp, A2, r=2)
            elif A2 and p == 3:
                des["coset-uniform r=3"] = CosetDesign(sp, A2, r=3)
        flat = (sp.N / k - 1) / R.sum()
        for name, d in des.items():
            assert d.k == k
            w = d.weights()
            designs_info.append({"p": p, "m": m, "k": k, "design": name, "n_directions": int((d.nu > 0).sum()),
                                 "max_w_over_flat": float(w[R].max() / flat),
                                 "min_pair_prob_positive": bool(np.all(d.pair_probability()[1:] > 0))})
        for yname, y in ys.items():
            v_srs = srs_variance(sp.N, k, y)
            EM, sR, dM, dR = decomposition(sp, fac, y)
            r2 = EM / (EM + sR)
            row = {"p": p, "m": m, "N": sp.N, "k": k, "response": yname, "R2": r2,
                   "adjR2": adjusted_r2(sp, fac, y), "var_srs": v_srs}
            for name, d in des.items():
                row["deff_" + name] = d.variance(y) / v_srs
            row["deff_flat_bound"] = 1 - row["adjR2"]
            rows.append(row)
    # files for the Monte Carlo cube runs
    ind = []
    for f in fac:
        lev = (sp.vecs @ f.T % p).reshape(-1)
        for a in range(1, p):  # drop level 0 of each factor (the constant is added in R)
            ind.append((lev == a).astype(int))
    np.savetxt(f"{OUT}/data/fact_p{p}_m{m}_X.csv", np.array(ind).T, fmt="%d", delimiter=",")
    np.savetxt(f"{OUT}/data/fact_p{p}_m{m}_Y.csv", np.c_[[ys[n] for n in ys]].T, delimiter=",",
               header=",".join(ys), comments="")
    return rows, designs_info, mixtures


if __name__ == "__main__":
    all_rows, all_info, all_mix = [], [], []
    for (p, m, ks, seed) in [(2, 7, [8, 16, 32], 7), (3, 5, [9, 27, 81], 5)]:
        rws, info, mix = run(p, m, ks, seed)
        all_rows += rws
        all_info += info
        all_mix += mix
    with open(f"{OUT}/results/factorial_minimax_mixtures.json", "w") as fh:   # nu of the minimax coset designs
        json.dump(all_mix, fh, indent=1)
    import csv
    keys = sorted({k for r in all_rows for k in r}, key=lambda s: (not s[0].isupper() and s not in ("p", "m", "N", "k", "response"), s))
    with open(f"{OUT}/results/factorial_exact.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()) + [k for k in keys if k not in all_rows[0]])
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    with open(f"{OUT}/results/factorial_designs.json", "w") as fh:
        json.dump(all_info, fh, indent=1)
    for r in all_rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    for i in all_info:
        print(i)
