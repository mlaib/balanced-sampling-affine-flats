"""Factorial cube results with the Horvitz-Thompson mean (Section 7.1, Tables 2 and 3, Appendix A).

The cube method occasionally returns k - 1 or k + 1 units.  With inclusion probabilities k/N the
Horvitz-Thompson estimator of the mean is sum(y_s)/k, with the fixed denominator k; it is design-unbiased,
so its mean squared error is its variance, whereas the realised sample mean is not the Horvitz-Thompson
estimator on the draws whose size differs from k.  From the stored samples of the package (mc_cube.R) and
of the transcription with both landings (cube_factorial.py), this script computes, for every configuration
and response:
  - the number of size exceptions;
  - the design effect of the Horvitz-Thompson mean, and of the realised sample mean for comparison;
  - the decomposition of the squared error of the Horvitz-Thompson mean,
        sum(y_s)/k - Ybar = (|s|/k - 1) Ybar + sum(mu_s)/k + sum(rho_s)/k,
    into the parts due to the size, the main effects mu, the residual rho and the cross terms;
  - the main-effect weight t_M = sum_{u in M} |eta_s(u)|^2 with eta_s(u) = k^-1 sum_{x in s} chi_u(x), that
    is sum_W (q sum_b n_b^2 - |s|^2)/k^2, and its share tau = t_M/(N/k - 1), with Monte Carlo standard errors.
Writes results/factorial_ht.csv.
"""
import csv, gzip, os
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
FRAMES = [(2, 7, (8, 16, 32)), (3, 5, (9, 27, 81))]


def load(q, m):
    X = np.loadtxt(f"{OUT}/data/fact_p{q}_m{m}_X.csv", delimiter=",")
    with open(f"{OUT}/data/fact_p{q}_m{m}_Y.csv") as fh:
        names = fh.readline().strip().split(",")
    Y = np.loadtxt(f"{OUT}/data/fact_p{q}_m{m}_Y.csv", delimiter=",", skiprows=1)
    Z = np.c_[np.ones(len(X)), X]
    coef, *_ = np.linalg.lstsq(Z, Y, rcond=None)
    MU = Z @ coef - Y.mean(axis=0)               # main-effect part, centred
    RHO = Y - Y.mean(axis=0) - MU                # residual part
    L = X.reshape(len(X), m, q - 1)                                       # indicators of levels 1..q-1 per variable
    levels = np.concatenate([1 - L.sum(axis=2, keepdims=True), L], axis=2)  # units x variables x levels 0..q-1
    return names, Y, MU, RHO, levels


def samples(path, landing=None):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if landing is None:
                k, units = int(f[0]), f[1:]
            else:
                if f[1] != landing:
                    continue
                k, units = int(f[0]), f[2:]
            out.setdefault(k, []).append(np.array(units, dtype=np.int64))
    return out


if __name__ == "__main__":
    rows = []
    for q, m, ks in FRAMES:
        names, Y, MU, RHO, LEV = load(q, m)
        N, Ybar = len(Y), Y.mean(axis=0)
        sources = {"package": samples(f"{OUT}/results/cube_samples/factorial_p{q}_m{m}.txt.gz"),
                   "transcription, package landing": samples(f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz", "bs"),
                   "transcription, flight continued": samples(f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz", "continued")}
        for k in ks:
            vsrs = (1 / k - 1 / N) * Y.var(axis=0, ddof=1)
            T = N / k - 1
            for src, S in sources.items():
                draws = S[k]
                B = len(draws)
                size = np.array([len(s) for s in draws])
                ht = np.array([Y[s].sum(axis=0) / k for s in draws])
                mean = np.array([Y[s].mean(axis=0) for s in draws])
                emu = np.array([MU[s].sum(axis=0) / k for s in draws])
                erho = np.array([RHO[s].sum(axis=0) / k for s in draws])
                esize = (size / k - 1)[:, None] * Ybar[None, :]
                assert np.abs(esize + emu + erho - (ht - Ybar)).max() < 1e-9
                cnt = np.array([LEV[s].sum(axis=0) for s in draws])        # draws x m x q level counts
                tM = (q * (cnt ** 2).sum(axis=2) - size[:, None] ** 2).sum(axis=1) / k ** 2
                exact = (size == k) & np.all(np.abs(cnt - k / q) < 1e-9, axis=(1, 2))
                for j, yname in enumerate(names):
                    e_ht, e_mean = (ht[:, j] - Ybar[j]) ** 2, (mean[:, j] - Ybar[j]) ** 2
                    parts = {"size": esize[:, j] ** 2, "main": emu[:, j] ** 2, "resid": erho[:, j] ** 2,
                             "cross": 2 * (esize[:, j] * emu[:, j] + esize[:, j] * erho[:, j] + emu[:, j] * erho[:, j])}
                    row = {"q": q, "m": m, "k": k, "source": src, "response": yname.replace(".", "-"), "B": B,
                           "size_exceptions": int((size != k).sum()), "size_below": int((size < k).sum()),
                           "size_above": int((size > k).sum()),
                           "deff_ht": e_ht.mean() / vsrs[j], "deff_ht_se": e_ht.std(ddof=1) / np.sqrt(B) / vsrs[j],
                           "deff_mean": e_mean.mean() / vsrs[j], "deff_mean_se": e_mean.std(ddof=1) / np.sqrt(B) / vsrs[j],
                           "share_exact": exact.mean(), "share_exact_se": np.sqrt(exact.mean() * (1 - exact.mean()) / B),
                           "tau": tM.mean() / T, "tau_se": tM.std(ddof=1) / np.sqrt(B) / T}
                    for name, v in parts.items():
                        row[f"{name}_part"] = v.mean() / vsrs[j]
                        row[f"{name}_part_se"] = v.std(ddof=1) / np.sqrt(B) / vsrs[j]
                    rows.append(row)
    with open(f"{OUT}/results/factorial_ht.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({a: (f"{b:.6g}" if isinstance(b, float) else b) for a, b in r.items()})
    for r in rows:
        if r["response"] == "additive":
            print(r["q"], r["k"], f'{r["source"]:32s}', "B", r["B"], "exceptions", r["size_exceptions"],
                  "deff HT %.4f (%.4f) mean %.4f" % (r["deff_ht"], r["deff_ht_se"], r["deff_mean"]),
                  "size part %.2e" % r["size_part"], "main %.4f resid %.4f cross %.4f" % (r["main_part"], r["resid_part"], r["cross_part"]),
                  "tau %.5f (%.5f)" % (r["tau"], r["tau_se"]), "exact %.4f" % r["share_exact"])
