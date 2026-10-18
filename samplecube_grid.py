"""samplecube of the package sampling on the 61 x 61 grid (Section 5.2, Tables 3 and S5): summaries of the stored
draws of run_samplecube.sh, combined over chunks.  For S2 and S4 and k = 61, 183: the number of draws, the sizes,
the design effects of the Horvitz-Thompson mean (fixed denominator k) and of the realised sample mean with Monte
Carlo standard errors, the shares of samples exactly balanced on S2 and on S4 with binomial standard errors, the
mean largest deviation of a stratum count of S4 from its target, and the main-effect weight shares on S2 and S4.
Writes results/samplecube_grid.csv.
"""
import csv, glob, gzip, json, math, os
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
q = 61


def read(path):
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            out.append((int(f[0]), np.array(f[2:], dtype=np.int64)))
    return out


if __name__ == "__main__":
    y = np.loadtxt(f"{OUT}/data/grid_y.csv")
    N, Ybar = len(y), y.mean()
    srs = {r["k"]: r["var"] for r in json.load(open(f"{OUT}/results/grid_exact.json"))["designs"] if r["design"] == "SRS"}
    rows = []
    for des in ("scubeS2", "scubeS4"):
        for k in (61, 183):
            draws = [s for path in sorted(glob.glob(f"{OUT}/results/cube_samples/grid_{des}_k{k}_c*.txt.gz"))
                     for kk, s in read(path) if kk == k]
            if not draws:
                continue
            B = len(draws)
            size = np.array([len(s) for s in draws])
            ht = np.array([y[s].sum() / k for s in draws])
            mean = np.array([y[s].mean() for s in draws])
            e_ht, e_mean = (ht - Ybar) ** 2, (mean - Ybar) ** 2
            ex2, ex4, dev4, tm2, tm4 = [], [], [], [], []
            for s in draws:
                i, j = s // q, s % q
                labs = (i, j, (i + j) % q, (i - j) % q)
                cnt = [np.bincount(l, minlength=q) for l in labs]
                dev = [np.abs(c - k / q).max() for c in cnt]
                ex2.append(len(s) == k and max(dev[:2]) < 1e-9)
                ex4.append(len(s) == k and max(dev) < 1e-9)
                dev4.append(max(dev))
                t = [(q * (c ** 2).sum() - len(s) ** 2) / k ** 2 for c in cnt]
                tm2.append(sum(t[:2])); tm4.append(sum(t))
            T = N / k - 1
            p2, p4 = np.mean(ex2), np.mean(ex4)
            rows.append({"design": des, "k": k, "B": B, "size_exceptions": int((size != k).sum()),
                         "size_min": int(size.min()), "size_max": int(size.max()),
                         "deff_ht": e_ht.mean() / srs[k], "deff_ht_se": e_ht.std(ddof=1) / math.sqrt(B) / srs[k],
                         "deff_mean": e_mean.mean() / srs[k], "deff_mean_se": e_mean.std(ddof=1) / math.sqrt(B) / srs[k],
                         "share_exact_S2": p2, "share_exact_S2_se": math.sqrt(p2 * (1 - p2) / B),
                         "share_exact_S4": p4, "share_exact_S4_se": math.sqrt(p4 * (1 - p4) / B),
                         "mean_maxdev_S4": float(np.mean(dev4)),
                         "tau_S2": float(np.mean(tm2)) / T, "tau_S4": float(np.mean(tm4)) / T,
                         "tau_S4_se": float(np.std(tm4, ddof=1)) / math.sqrt(B) / T})
    with open(f"{OUT}/results/samplecube_grid.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({a: (f"{b:.6g}" if isinstance(b, float) else b) for a, b in r.items()})
    for r in rows:
        print(r["design"], r["k"], "B", r["B"], "size exc", r["size_exceptions"], "deff_ht %.4f (%.4f)" % (r["deff_ht"], r["deff_ht_se"]),
              "exact S2 %.3f S4 %.3f" % (r["share_exact_S2"], r["share_exact_S4"]), "tau_S4 %.4f" % r["tau_S4"], "maxdev %.2f" % r["mean_maxdev_S4"])
