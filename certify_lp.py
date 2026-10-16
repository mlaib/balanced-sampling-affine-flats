"""Exact certificates for Table 1 and the archived configurations of program (8).

Attainment (LP value = bound): rebuild a mixing distribution over weight-distribution classes of
admissible directions in rational arithmetic and verify exactly that c_nu(u) = (Q-1)/|R_S| for every
residual frequency u.  nu is uniform within each class, which is invariant under the monomial group,
so c_nu(u) depends only on the Hamming weight of u.

Non-attainment (LP value > bound): verify in rational arithmetic a dual solution, a probability
distribution y on R_S that is constant on Hamming-weight classes, with y(D minus 0) >= z for every
admissible dual direction D, enumerated one by one.  For any nu, max_u c_nu(u) >= sum_u y_u c_nu(u)
= E_nu y(D minus 0) >= z, so no coset design C(h, nu, 1) has max_u c_nu(u) below z.  A rational primal
mixture with max_u c_nu(u) = z, checked frequency by frequency, shows that z is the minimax value.
"""
import csv, itertools, json, math, sys
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from affsamp import Space, subspaces, main_effect_set
from study_minimax_lp import margin_factors

OUT = __file__.rsplit("/", 1)[0]

def solve_exact(A, b):
    """Solve A x = b exactly (Fractions); A has full column rank. Returns x or None."""
    A = [row[:] for row in A]; b = b[:]
    n, p = len(A), len(A[0])
    M = [A[i] + [b[i]] for i in range(n)]
    piv_cols, r = [], 0
    for c in range(p):
        piv = next((i for i in range(r, n) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [vi - f * vr for vi, vr in zip(M[i], M[r])]
        piv_cols.append(c); r += 1
    if len(piv_cols) < p:
        return None
    if any(all(v == 0 for v in M[i][:p]) and M[i][p] != 0 for i in range(n)):
        return None  # inconsistent
    x = [Fraction(0)] * p
    for i, c in enumerate(piv_cols):
        x[c] = M[i][p]
    return x



def certify_above_bound(p, m, order, h, sp, M, wt, T, types, adm, js, size, a, res, bound):
    """Exact primal and dual solutions of the reduced program, verified on the unreduced problem."""
    nT = len(T)
    lam, t = res.x[:nT], res.fun
    beta = -np.asarray(res.ineqlin.marginals)
    A = [[Fraction(T[i][j], s) for j, s in zip(js, size)] for i in range(nT)]   # a_{T,j} exactly
    # primal: lambda on its support, t, with the tight weight constraints as equations
    S = [i for i in range(nT) if lam[i] > 1e-9]
    tight = [jj for jj in range(len(js)) if a[:, jj] @ lam > t - 1e-9]
    x = solve_exact([[A[i][jj] for i in S] + [Fraction(-1)] for jj in tight] + [[Fraction(1)] * len(S) + [Fraction(0)]],
                    [Fraction(0)] * len(tight) + [Fraction(1)])
    lam_ex = dict(zip(S, x[:-1])) if x else {i: Fraction(lam[i]).limit_denominator(10 ** 9) for i in S}
    # dual: beta on the weights where it is positive, z, with the tight class constraints as equations
    act = [jj for jj in range(len(js)) if beta[jj] > 1e-9]
    tightT = [i for i in range(nT) if a[i] @ beta < t + 1e-9]
    yv = solve_exact([[A[i][jj] for jj in act] + [Fraction(-1)] for i in tightT] + [[Fraction(1)] * len(act) + [Fraction(0)]],
                     [Fraction(0)] * len(tightT) + [Fraction(1)])
    beta_ex = dict(zip(act, yv[:-1])) if yv else {jj: Fraction(beta[jj]).limit_denominator(10 ** 9) for jj in act}
    # dual check on every admissible direction: sum_j beta_j A_j(D) / n_j >= z
    dual_ok = all(v >= 0 for v in beta_ex.values()) and sum(beta_ex.values()) == 1
    z_d = min(sum(beta_ex[jj] * Fraction(key[js[jj]], size[jj]) for jj in beta_ex) for key, _ in adm)
    # primal check on every residual frequency: c_nu(u) with nu uniform within each class
    cls = {T[i]: i for i in S}
    inc = {i: np.zeros(sp.N, dtype=np.int64) for i in S}
    for key, codes in adm:
        if key in cls:
            inc[cls[key]][codes] += 1
    R = np.flatnonzero(~M)
    R = R[R != 0]
    primal_ok = all(v >= 0 for v in lam_ex.values()) and sum(lam_ex.values()) == 1
    z_p = max(sum(lam_ex[i] * Fraction(int(inc[i][u]), types[T[i]]) for i in S) for u in R)
    return {"q": p, "m": m, "balanced_on": "one-way" if order == 1 else "one- and two-way", "k": p ** h,
            "n_admissible": len(adm), "bound": str(bound), "minimax_value": str(z_d), "ratio": float(z_d / bound),
            "lp_value_float": t, "dual_verified": bool(dual_ok), "primal_verified": bool(primal_ok),
            "value_exact": bool(z_p == z_d), "above_bound": bool(z_d > bound),
            "dual_mass_by_hamming_weight": {str(js[jj]): str(beta_ex.get(jj, Fraction(0))) for jj in range(len(js))},
            "primal_mixture": [{"weight_distribution": list(T[i]), "n_directions_in_class": types[T[i]],
                                "probability_of_class": str(lam_ex[i])} for i in S]}


certs, diag, nonatt = [], [], []
for (p, m, order) in [(2, 5, 1), (2, 6, 1), (2, 7, 1), (2, 6, 2), (2, 7, 2), (3, 4, 1), (3, 5, 1), (3, 4, 2)]:
    sp = Space(p, m)
    M = main_effect_set(sp, margin_factors(m, order))
    wt = (sp.vecs != 0).sum(axis=1)
    nR = int((~M).sum() - 1)
    for h in range(order, m):
        types, adm = {}, []
        for D in subspaces(m, m - h, p):
            codes = sp.span(D)
            if M[codes].any():
                continue
            key = tuple(int(v) for v in np.bincount(wt[codes], minlength=m + 1))
            types[key] = types.get(key, 0) + 1
            adm.append((key, codes))
        if not types:
            continue
        js = list(range(order + 1, m + 1))
        size = [math.comb(m, j) * (p - 1) ** j for j in js]
        T = list(types)
        a = np.array([[t[j] / s for j, s in zip(js, size)] for t in T])
        nT = len(T)
        res = linprog(np.r_[np.zeros(nT), 1.0], A_ub=np.c_[a.T, -np.ones(len(js))], b_ub=np.zeros(len(js)),
                      A_eq=np.r_[np.ones(nT), 0.0].reshape(1, -1), b_eq=[1.0],
                      bounds=[(0, None)] * nT + [(None, None)], method="highs")
        lam = res.x[:nT]
        bound = Fraction(p ** (m - h) - 1, nR)
        resid = max(abs(a.T @ lam - res.x[-1]).max(), abs(lam.sum() - 1))
        diag.append({"q": p, "m": m, "order": order, "k": p ** h, "lp_value": res.fun, "bound": float(bound),
                     "ratio": res.fun / float(bound), "max_residual_equal_weight": float(abs(a.T @ lam - float(bound)).max()),
                     "lp_status": res.status, "mip_gap_like": float(res.fun - float(bound))})
        if res.fun / float(bound) > 1 + 1e-9:
            nonatt.append(certify_above_bound(p, m, order, h, sp, M, wt, T, types, adm, js, size, a, res, bound))
            c = nonatt[-1]
            print(f"q={p} m={m} order={order} k={p**h}: value={c['minimax_value']} ratio={c['ratio']:.6f} "
                  f"dual verified={c['dual_verified']} primal verified={c['primal_verified']} exact={c['value_exact']}")
            continue
        # exact reconstruction on the LP support
        S = [i for i in range(nT) if lam[i] > 1e-10]
        Aex = [[Fraction(T[i][j], s) for i in S] for j, s in zip(js, size)] + [[Fraction(1)] * len(S)]
        bex = [bound] * len(js) + [Fraction(1)]
        x = solve_exact(Aex, bex)
        ok = x is not None and all(v >= 0 for v in x)
        if ok:  # exact verification of c_nu(u) for every residual weight class
            cov = [sum(x[t] * Fraction(T[S[t]][j], s) for t in range(len(S))) for j, s in zip(js, size)]
            ok = all(c == bound for c in cov) and sum(x) == 1
        certs.append({"q": p, "m": m, "balanced_on": "one-way" if order == 1 else "one- and two-way", "k": p ** h,
                      "bound": str(bound), "certified": bool(ok),
                      "mixture": [{"weight_distribution": list(T[S[t]]), "n_directions_in_class": types[T[S[t]]],
                                   "probability_of_class": str(x[t])} for t in range(len(S))] if x else None})
        print(f"q={p} m={m} order={order} k={p**h}: certified={ok}  classes={len(S)}  LP residual={resid:.1e}")
json.dump({"certificates": certs, "nonattainment_certificates": nonatt, "lp_diagnostics": diag,
           "solver": "scipy.optimize.linprog, method='highs', default tolerances (primal/dual feasibility 1e-7)"},
          open(f"{OUT}/results/lp_certificates.json", "w"), indent=1)
print("max equal-weight residual over certified configurations:",
      max(d["max_residual_equal_weight"] for d in diag if d["ratio"] <= 1 + 1e-9))
