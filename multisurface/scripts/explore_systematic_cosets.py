"""EXPLORATORY (not in PLAN.md): k = 183 with the 3 cosets chosen systematically, labels {c, c+20, c+40} mod 61,
c uniform, instead of by simple random sampling. Same best-16 directions. Exact weights:
for u = t d in V^perp (d the dual basis vector), w(u) = |(1/3) sum_j e(t j delta)|^2, averaged over nu."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, ".."))
from affsamp import Space, admissible_directions, main_effect_set
q, N, k, delta = 61, 3721, 183, 20
sp = Space(q, 2)
S4 = [np.array([[1, 0]]), np.array([[0, 1]]), np.array([[1, 1]]), np.array([[1, q - 1]])]
A4 = admissible_directions(sp, S4, 1)
def centred(a):
    a = a % q; return np.where(a > q // 2, a - q, a)
qual = np.array([np.sqrt((centred(np.arange(1, q)[:, None] * D[0]) ** 2).sum(1)).min() for D in A4])
best16 = [A4[i] for i in np.argsort(-qual, kind="stable")[:16]]
w = np.zeros(N)
for D in best16:
    for t in range(1, q):
        u = sp.code((t * D[0]) % q)
        w[u] += abs(np.exp(2j * np.pi * t * np.arange(3) * delta / q).mean()) ** 2 / len(best16)
assert abs(w.sum() - (N / k - 1)) < 1e-9          # Lemma 2(b): total weight N/k - 1
M = main_effect_set(sp, S4); R = ~M; R[0] = False
print("W over Y_1(S4) relative to equal weight:", round(w[R].max() / ((N / k - 1) / R.sum()), 3))
surf = [np.loadtxt(os.path.join(ROOT, "..", "data", "volcano61.csv"), delimiter=",")]
eb = np.load(os.path.join(ROOT, "data", "elevation_blocks.npz")); names = ["volcano"] + list(eb["names"])
surf += list(eb["blocks"])
Ys = np.array(surf, float); P = (np.abs(np.fft.fft2(Ys) / N) ** 2).reshape(len(Ys), N)
srs = (1 - k / N) * Ys.reshape(len(Ys), -1).var(axis=1, ddof=1) / k
deff_sys = P @ w / srs
import csv, collections
lpm = {}
for r in csv.DictReader(open(os.path.join(ROOT, "results", "deff_long.csv"))):
    if r["k"] == "183" and r["design"] in ("lpm2", "best16_S4") and r["source"] in ("elevation", "volcano"):
        lpm.setdefault(r["design"], []).append(float(r["deff"]))
lp, b16 = np.array(lpm["lpm2"]), np.array(lpm["best16_S4"])
e = np.arange(len(Ys)) > 0
for lab, d in (("best-16, random cosets", b16), ("best-16, systematic cosets", deff_sys)):
    ratio = d[e] / lp[e]
    print(f"{lab:28s} volcano {d[0]:.4f} | real blocks: wins vs lpm2 {np.mean(ratio < 1):.1%}, median ratio {np.median(ratio):.3f}")
