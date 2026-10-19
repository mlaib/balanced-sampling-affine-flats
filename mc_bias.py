"""Mean error of the Horvitz-Thompson mean, with its Monte Carlo standard error (Section 5; Supplementary S2).

With inclusion probabilities k/N the Horvitz-Thompson estimator of the mean is sum(y_s)/k, with the fixed
denominator k, and it is design-unbiased also when a draw returns another size. From every stored Monte Carlo
sample (the cube methods in the factorial frames and on the grid, the local cube method and samplecube), this
script computes for each configuration and response the number of draws of another size, the mean error
e = mean(sum(y_s)/k) - Ybar, its standard error sd/sqrt(B) and z = e/se. The local pivotal draws are not stored
and have fixed size. Writes results/mc_bias.csv.
"""
import csv, glob, gzip, os
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))


def read(paths, landing=None):
    out = {}
    for path in paths:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.split()
                if landing is not None and f[1] != landing:
                    continue
                skip = 1 if f[1].isdigit() else 2                  # "k units" or "k label units"
                out.setdefault(int(f[0]), []).append(np.array(f[skip:], dtype=np.int64))
    return out


rows = []


def evaluate(population, source, draws, Y, names):
    for k, S in sorted(draws.items()):
        B = len(S)
        size = np.array([len(s) for s in S])
        ht = np.array([Y[s].sum(axis=0) / k for s in S])
        for j, name in enumerate(names):
            e = ht[:, j] - Y[:, j].mean()
            se = e.std(ddof=1) / np.sqrt(B)
            rows.append({"population": population, "source": source, "k": k, "response": name, "B": B,
                         "size_exceptions": int((size != k).sum()), "mean_error": e.mean(), "se": se,
                         "z": e.mean() / se})


for q, m in ((2, 7), (3, 5)):
    with open(f"{OUT}/data/fact_p{q}_m{m}_Y.csv") as fh:
        names = fh.readline().strip().split(",")
    Y = np.loadtxt(f"{OUT}/data/fact_p{q}_m{m}_Y.csv", delimiter=",", skiprows=1)
    cs = f"{OUT}/results/cube_samples"
    evaluate(f"{q}^{m}", "package", read([f"{cs}/factorial_p{q}_m{m}.txt.gz"]), Y, names)
    evaluate(f"{q}^{m}", "transcription, package landing", read([f"{cs}/transcription_p{q}_m{m}.txt.gz"], "bs"), Y, names)
    evaluate(f"{q}^{m}", "flight continued", read([f"{cs}/transcription_p{q}_m{m}.txt.gz"], "continued"), Y, names)
    for impl in ("scube1", "scube2", "bs_rev", "scube1_rev", "scube2_rev"):
        evaluate(f"{q}^{m}", impl, read([f"{cs}/factorial_{impl}_p{q}_m{m}.txt.gz"]), Y, names)

y = np.loadtxt(f"{OUT}/data/grid_y.csv")[:, None]
for des in ("cubeS2", "cubeS4", "cubeS4_fixed", "lcubeS4", "scubeS2", "scubeS4"):
    paths = sorted(glob.glob(f"{OUT}/results/cube_samples/grid_{des}_k*.txt.gz"))
    evaluate("grid", des, read(paths), y, ["elevation"])

with open(f"{OUT}/results/mc_bias.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow({a: (f"{b:.6g}" if isinstance(b, float) else b) for a, b in r.items()})
z = np.array([r["z"] for r in rows])
print(len(rows), "configurations and responses; |z| > 1.96 in", int((np.abs(z) > 1.96).sum()), "; max |z| = %.2f" % np.abs(z).max())
for r in rows:
    if abs(r["z"]) > 1.96:
        print(r["population"], r["source"], r["k"], r["response"], r["B"], r["size_exceptions"], "z = %.2f" % r["z"])
