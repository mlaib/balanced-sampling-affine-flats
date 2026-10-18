"""Plug-in lower bounds of Proposition 3(a) for the Monte Carlo designs, with the correction for a random size.

For a design with inclusion probabilities k/N whose size |s| may differ from k, eta_s(u) = k^-1 sum_{x in s}
chi_u(x) has sum_{u != 0} E|eta_s(u)|^2 = N/k - 1 - var(|s|)/k^2, because eta_s(0) = |s|/k and E|s| = k.  For
responses with mean zero the error of the Horvitz-Thompson mean is sum_{u != 0} yhat(u) eta_s(u), and the proof
of Proposition 3(a) gives, relative to the equal-weight design (worst case sigma^2 T0/|R_S|, T0 = N/k - 1),

    1 + tau (C |R_S|/d_S - 1) - var(|s|)/(k^2 T0),        tau = t_M/T0,

the bound (9) of the manuscript.  With E|s| = k, var(|s|) = E(|s| - k)^2, so the bound is the mean over draws of
    Z = 1 + (C |R_S|/d_S - 1) t_M(s)/T0 - (|s| - k)^2/(k^2 T0),
and its Monte Carlo standard error is the standard deviation of Z divided by sqrt(B), which accounts for the
correlation of the two terms.  t_M(s) = sum_W (q^{n_W} sum_b n_b^2 - |s|^2)/k^2 over the stratifications W of S,
whose subspaces W^perp meet only in 0.  C = 1 throughout, as in Table 6.

Factorial frames: S = the one-variable margins; samples of the package (mc_cube.R) and of the transcription
with the flight continued (cube_factorial.py).  Grid: S = S4 (four line stratifications); samples of the
package (cube on S2 and S4, local cube on S4), of the cube method with directions by QR factorisation, and of
the audit (package landing and flight continued, S2), and of samplecube of the package sampling (S2 and S4).  The local pivotal method has fixed size and no stored
samples; its bound is 1 + 13.5 tau from results/mc_grid_v3.  Writes results/bound_centred.csv.
"""
import csv, glob, gzip, os
import numpy as np
from ht_factorial import FRAMES, load, samples

OUT = os.path.dirname(os.path.abspath(__file__))
C = 1.0


def grid_tm(draws, k, q=61):
    """Per-draw sizes and t_M for S4 (rows i, columns j, anti-diagonals i + j, diagonals i - j), units 61 i + j."""
    size = np.array([len(s) for s in draws])
    tm = np.zeros(len(draws))
    for b, s in enumerate(draws):
        i, j = s // q, s % q
        for lab in (i, j, (i + j) % q, (i - j) % q):
            n_b = np.bincount(lab, minlength=q)
            tm[b] += q * (n_b ** 2).sum() - len(s) ** 2
    return size, tm / k ** 2


def summary(frame, k, design, size, tm, N, d):
    T, R = N / k - 1, N - 1 - d
    a = C * R / d - 1
    corr = (size - k) ** 2 / (k ** 2 * T)
    Z = 1 + a * tm / T - corr
    B = len(Z)
    return {"frame": frame, "k": k, "design": design, "B": B, "size_exceptions": int((size != k).sum()),
            "size_min": int(size.min()), "size_max": int(size.max()), "mean_size": size.mean(),
            "var_size": float(((size - k) ** 2).mean()), "tau": tm.mean() / T, "tau_se": tm.std(ddof=1) / np.sqrt(B) / T,
            "factor": a, "correction": corr.mean(), "bound_uncorrected": 1 + a * tm.mean() / T,
            "bound": Z.mean(), "bound_se": Z.std(ddof=1) / np.sqrt(B)}


def read_grid(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            out.setdefault(int(f[0]), []).append(np.array(f[2:], dtype=np.int64))
    return out


if __name__ == "__main__":
    rows = []
    for q, m, ks in FRAMES:
        names, Y, MU, RHO, LEV = load(q, m)
        N, d = len(Y), m * (q - 1)
        src = {"cube (package)": samples(f"{OUT}/results/cube_samples/factorial_p{q}_m{m}.txt.gz"),
               "cube, flight continued": samples(f"{OUT}/results/cube_samples/transcription_p{q}_m{m}.txt.gz", "continued")}
        for impl, name in (("scube1", "samplecube, LP landing"), ("scube2", "samplecube, suppression")):
            path = f"{OUT}/results/cube_samples/factorial_{impl}_p{q}_m{m}.txt.gz"
            if os.path.exists(path):
                src[name] = samples(path)
        for k in ks:
            for name, S in src.items():
                draws = S[k]
                size = np.array([len(s) for s in draws])
                cnt = np.array([LEV[s].sum(axis=0) for s in draws])
                tm = (q * (cnt ** 2).sum(axis=2) - size[:, None] ** 2).sum(axis=1) / k ** 2
                rows.append(summary(f"{q}^{m}", k, name, size, tm, N, d))
    N, d = 61 ** 2, 4 * 60
    for k in (61, 183):
        for name, path in (("cube S2 (package)", f"grid_cubeS2_k{k}"), ("cube S4 (package)", f"grid_cubeS4_k{k}"),
                           ("local cube S4 (package)", f"grid_lcubeS4_k{k}"),
                           ("cube S4, QR directions", f"grid_cubeS4_fixed_k{k}")):
            draws = read_grid(f"{OUT}/results/cube_samples/{path}.txt.gz")[k]
            size, tm = grid_tm(draws, k)
            rows.append(summary("61x61", k, name, size, tm, N, d))
        for des, name in (("scubeS2", "samplecube S2"), ("scubeS4", "samplecube S4")):
            paths = sorted(glob.glob(f"{OUT}/results/cube_samples/grid_{des}_k{k}_c*.txt.gz"))
            if paths:                                       # samplecube of the package sampling, in chunks
                draws = [s for path in paths for s in read_grid(path).get(k, [])]
                size, tm = grid_tm(draws, k)
                rows.append(summary("61x61", k, name, size, tm, N, d))
        audit = [r for r in csv.DictReader(open(f"{OUT}/results/cube_audit_S2_draws.csv")) if int(r["k"]) == k]
        for landing, name in (("bs", "cube S2, transcription"), ("continued", "cube S2, flight continued")):
            draws = [np.array(r["sample"].split(), dtype=np.int64) for r in audit if r["landing"] == landing]
            size, tm = grid_tm(draws, k)
            rows.append(summary("61x61", k, name, size, tm, N, d))
        for name, f in (("local pivotal, planar", "lpm2"), ("local pivotal, wrapped", "lpm2torus")):
            r = next(csv.DictReader(open(f"{OUT}/results/mc_grid_v3/{f}_k{k}.csv")))
            T, a = N / k - 1, C * (N - 1 - d) / d - 1
            tau, tau_se = float(r["tau_S4"]), float(r["tau_S4_se"])
            rows.append({"frame": "61x61", "k": k, "design": name, "B": int(r["B"]), "size_exceptions": 0,
                         "size_min": k, "size_max": k, "mean_size": float(k), "var_size": 0.0, "tau": tau, "tau_se": tau_se,
                         "factor": a, "correction": 0.0, "bound_uncorrected": 1 + a * tau, "bound": 1 + a * tau,
                         "bound_se": a * tau_se})
    keys = list(rows[0])
    with open(f"{OUT}/results/bound_centred.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({a: (f"{b:.8g}" if isinstance(b, float) else b) for a, b in r.items()})
    for r in rows:
        print(f'{r["frame"]:6s} {r["k"]:4d} {r["design"]:28s} B {r["B"]:6d} exc {r["size_exceptions"]:4d} '
              f'tau {r["tau"]:.5f} ({r["tau_se"]:.5f}) corr {r["correction"]:.2e} '
              f'bound {r["bound_uncorrected"]:.4f} -> {r["bound"]:.4f} ({r["bound_se"]:.4f})')
