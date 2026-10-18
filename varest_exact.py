"""Exact evaluation of two variance estimators for replicated coset designs (Section 4 and Supplementary S6, Table S7).

For a design that selects r cosets of one random direction, the support is enumerated (every direction in
the support of nu and every set of r cosets, all equally likely), and for each response we compute:

  v      the random-group estimator (1 - r/Q) s_c^2 / r of Proposition 7(b), with the interval
         ybar +- t_{0.975, r-1} v^{1/2};
  v_reg  the regression estimator (1/k - 1/N) RSS / (k - rank), where RSS is the residual sum of squares of
         the least-squares fit, on the sample, of the main-effects model of the balancing variables; by
         Corollary 3 it estimates the variance of an equal-weight design, (N/k - 1) sigma_R^2 / |R_S|; it is
         defined when k exceeds the rank of the model on the sample, and the interval uses t_{0.975, k-rank}.

Designs: in the factorial frames of Section 5.1, the replicated design C(h-1, nu, q) with nu uniform on
A_{h-1}(S) (the '$r=q$' column of Table S3); on the elevation grid, the r = 3 equal-weight design
balanced on rows (S1), for which the sample regression is the one-way analysis of variance by rows.
Writes results/varest_exact.csv.
"""
import csv, itertools, math, os
import numpy as np
from scipy import stats
from affsamp import Space, CosetDesign, admissible_directions, coset_labels
from study_factorial import coordinate_factors, responses

OUT = os.path.dirname(os.path.abspath(__file__))


def summarise(prob, ybar, Ybar, V, est, df, label, extra):
    """Moments and coverage of one estimator over the support (prob: probability of each sample)."""
    ok = np.isfinite(est)
    Ev = float((prob * np.where(ok, est, 0)).sum() / prob[ok].sum())
    Ev2 = float((prob * np.where(ok, est, 0) ** 2).sum() / prob[ok].sum())
    tq = stats.t.ppf(0.975, df)                     # df may vary from sample to sample (rank of the sample model)
    cover = float((prob * (np.abs(ybar - Ybar) <= tq * np.sqrt(np.where(ok, est, 0)))).sum())
    row = {"estimator": label, "true_var": V, "E_est": Ev, "rel_bias": Ev / V - 1,
           "cv": math.sqrt(max(Ev2 - Ev ** 2, 0.0)) / Ev, "coverage_t95": cover,
           "share_zero": float((prob * (est == 0)).sum()), "share_defined": float(prob[ok].sum())}
    row.update(extra)
    return row


def factorial_rows():
    rows = []
    for (q, m, ks, seed) in [(2, 7, (8, 16, 32), 7), (3, 5, (9, 27, 81), 5)]:
        sp = Space(q, m)
        fac = coordinate_factors(m)
        ys = responses(sp, seed)
        names = list(ys)
        Y = np.array([ys[n] for n in names])                       # responses x N
        lev = sp.vecs                                              # level of each unit on each variable
        Xfull = np.ones((sp.N, 1))
        for i in range(m):
            for a in range(1, q):
                Xfull = np.c_[Xfull, (lev[:, i] == a).astype(float)]
        npar = Xfull.shape[1]                                      # 1 + m(q-1)
        for k in ks:
            h = round(math.log(k, q))
            A = admissible_directions(sp, fac, h - 1)
            d = CosetDesign(sp, A, r=q)
            assert d.k == k
            r, Q = q, d.Q
            combos = np.array(list(itertools.combinations(range(Q), r)))
            V = np.array([d.variance(y) for y in Y])
            acc = {"ybar": [], "v": [], "vreg": [], "df": []}
            for D in d.D:
                lab = coset_labels(sp, D)
                members = np.array([np.nonzero(lab == c)[0] for c in range(Q)])   # Q x q^(h-1) units
                cm = Y[:, members].mean(axis=2)                                     # responses x Q
                c = cm[:, combos]                                                   # responses x ncomb x r
                acc["ybar"].append(c.mean(axis=2))
                acc["v"].append((1 - r / Q) * c.var(axis=2, ddof=1) / r)
                if k > npar:
                    idx = members[combos].reshape(len(combos), -1)                  # ncomb x k units
                    U, S, _ = np.linalg.svd(Xfull[idx], full_matrices=False)
                    keep = S > 1e-9 * S[:, :1]
                    ysamp = Y[:, idx]                                               # responses x ncomb x k
                    coef = np.einsum("bki,rbk->rbi", U, ysamp) * keep[None]
                    rss = (ysamp ** 2).sum(axis=2) - (coef ** 2).sum(axis=2)
                    df = k - keep.sum(axis=1)
                    acc["vreg"].append((1 / k - 1 / sp.N) * rss / df[None])
                    acc["df"].append(df)
            ybar = np.concatenate(acc["ybar"], axis=1)
            v = np.concatenate(acc["v"], axis=1)
            prob = np.full(ybar.shape[1], 1.0 / ybar.shape[1])
            df_all = np.concatenate(acc["df"]) if acc["df"] else None
            for j, yname in enumerate(names):
                base = {"population": f"{q}^{m}", "k": k, "r": r, "response": yname, "design": "C(h-1,nu,q), nu uniform",
                        "n_samples": len(prob)}
                rows.append(summarise(prob, ybar[j], Y[j].mean(), V[j], v[j], r - 1, "random groups",
                                      dict(base, df_min=r - 1, df_max=r - 1)))
                if df_all is not None:
                    vr = np.concatenate(acc["vreg"], axis=1)[j]
                    rows.append(summarise(prob, ybar[j], Y[j].mean(), V[j], vr, df_all, "regression",
                                          dict(base, df_min=int(df_all.min()), df_max=int(df_all.max()),
                                               share_full_rank=float((df_all == k - npar).mean()))))
        print(f"{q}^{m} done", flush=True)
    return rows


def grid_rows():
    """r = 3 equal-weight design balanced on rows: each coset is a line meeting every row once."""
    q, k, r = 61, 183, 3
    sp = Space(q, 2)
    y = np.loadtxt(f"{OUT}/data/volcano61.csv", delimiter=",").reshape(-1)
    A = admissible_directions(sp, [np.array([[1, 0]])], 1)
    d = CosetDesign(sp, A, r=r)
    V = d.variance(y)
    combos = np.array(list(itertools.combinations(range(q), r)))
    ybar_all, v_all, vr_all = [], [], []
    row_of = sp.vecs[:, 0]
    for D in d.D:
        lab = coset_labels(sp, D)
        Yc = np.zeros((q, q))                     # Yc[c, i]: value of the unit of coset c in row i
        Yc[lab, row_of] = y
        cm = Yc.mean(axis=1)
        SS = (Yc ** 2).sum(axis=1)
        G = Yc @ Yc.T
        c = cm[combos]
        ybar_all.append(c.mean(axis=1))
        v_all.append((1 - r / q) * c.var(axis=1, ddof=1) / r)
        a, b, e = combos.T
        rss = (2 / 3) * (SS[a] + SS[b] + SS[e]) - (2 / 3) * (G[a, b] + G[a, e] + G[b, e])
        vr_all.append((1 / k - 1 / sp.N) * rss / (k - q))
    ybar, v, vr = (np.concatenate(x) for x in (ybar_all, v_all, vr_all))
    prob = np.full(len(ybar), 1.0 / len(ybar))
    base = {"population": "61x61 grid", "k": k, "r": r, "response": "elevation", "design": "equal-weight, S1",
            "n_samples": len(prob)}
    return [summarise(prob, ybar, y.mean(), V, v, r - 1, "random groups", dict(base, df_min=2, df_max=2)),
            summarise(prob, ybar, y.mean(), V, vr, k - q, "regression", dict(base, df_min=k - q, df_max=k - q))]


if __name__ == "__main__":
    rows = factorial_rows() + grid_rows()
    keys = ["population", "k", "r", "response", "design", "estimator", "n_samples", "df_min", "df_max", "share_full_rank",
            "true_var", "E_est", "rel_bias", "cv", "coverage_t95", "share_zero", "share_defined"]
    with open(f"{OUT}/results/varest_exact.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r_ in rows:
            w.writerow({k_: (f"{x:.6g}" if isinstance(x, float) else x) for k_, x in r_.items()})
    for r_ in rows:
        print(r_["population"], r_["k"], r_["response"], r_["estimator"], "df", r_["df_min"],
              "bias %.4f cv %.2f cover %.3f zero %.4f" % (r_["rel_bias"], r_["cv"], r_["coverage_t95"], r_["share_zero"]))
