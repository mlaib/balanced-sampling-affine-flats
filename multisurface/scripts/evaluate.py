"""Design effects of every design of PLAN.md on every surface, and the verdict of the decision rule.

Exact designs: spectral formula with the coset weights of affsamp.py. Monte Carlo designs: the stored samples of
draw_mc.R, the Horvitz-Thompson mean with fixed denominator k, mean squared error about the known mean, standard
error of a mean of squared errors. Writes results/deff_long.csv, results/design_properties.json and
results/summary.md.
"""
import gzip
import json
import os
import sys

import numpy as np
from scipy import sparse

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
CODE = os.path.join(ROOT, "..")
sys.path.insert(0, CODE)
from affsamp import Space, CosetDesign, admissible_directions, main_effect_set  # noqa: E402

q, N = 61, 61 * 61
KS = (61, 183)
sp = Space(q, 2)
ROWS, COLS = np.array([[1, 0]]), np.array([[0, 1]])
ANTI, DIAG = np.array([[1, 1]]), np.array([[1, q - 1]])
S2, S4 = [ROWS, COLS], [ROWS, COLS, ANTI, DIAG]


def centred(a):
    a = a % q
    return np.where(a > q // 2, a - q, a)


# ---------------------------------------------------------------- surfaces (code 61 i + j = Y.reshape(-1))
surf, meta = [], []
vol = np.loadtxt(os.path.join(CODE, "data", "volcano61.csv"), delimiter=",")
surf.append(vol)
meta.append({"source": "volcano", "group": "volcano (paper)", "id": "volcano"})
eb = np.load(os.path.join(ROOT, "data", "elevation_blocks.npz"))
for b, n, p in zip(eb["blocks"], eb["names"], eb["pos"]):
    surf.append(b)
    meta.append({"source": "elevation", "group": str(n), "id": f"{n} {p[0]}{p[1]}"})
sf = np.load(os.path.join(ROOT, "data", "simulated_fields.npz"))
for t, (f, s) in enumerate(zip(sf["fields"], sf["settings"])):
    surf.append(f)
    meta.append({"source": "simulated", "group": f"nu={s[0]:g} range={s[1]:g} nugget={s[2]:g}", "id": f"sim{t}"})
Ys = np.array([s.reshape(-1) for s in surf], dtype=float)          # (n_surf, N)
n_surf = len(Ys)
Ybar = Ys.mean(axis=1)
S2y = Ys.var(axis=1, ddof=1)
srs = {k: (1 - k / N) * S2y / k for k in KS}
P = (np.abs(np.fft.fft2(np.array(surf, dtype=float)) / N) ** 2).reshape(n_surf, N)   # |yhat(u)|^2

# ---------------------------------------------------------------- exact coset designs
A2 = admissible_directions(sp, S2, 1)
A4 = admissible_directions(sp, S4, 1)


def shortest_dual(D):
    v = centred(np.arange(1, q)[:, None] * D[0])
    return float(np.sqrt((v ** 2).sum(axis=1)).min())


quality = np.array([shortest_dual(D) for D in A4])
order = np.argsort(-quality, kind="stable")
exact = {"eqw_S2": A2, "eqw_S4": A4}
for L in (4, 8, 16):
    exact[f"best{L}_S4"] = [A4[i] for i in order[:L]]
M4 = main_effect_set(sp, S4)
R4 = ~M4
R4[0] = False
d4, nR4 = int(M4.sum()), int(R4.sum())
props = {}
rows = []
for name, dirs in exact.items():
    for k in KS:
        d = CosetDesign(sp, dirs, r=k // q)
        w = d.weights()
        var = P @ w
        deff = var / srs[k]
        if k == 61:
            w_eq = (N / k - 1) / nR4                       # equal-weight value on R_S4
            props[name] = {"n_directions": len(dirs),
                           "W_Y1_S4_rel": float((w[M4].max() * 1 + w[R4].max()) / w_eq)}   # Prop 3(b), C = 1
        for i in range(n_surf):
            rows.append((i, k, name, deff[i], 0.0))
# the paper's numbers for the volcano block (Table 3)
paper = {"eqw_S2": 0.1621, "eqw_S4": 0.0878, "best4_S4": 0.0102, "best8_S4": 0.0113, "best16_S4": 0.0118}
check = {nm: float([r[3] for r in rows if r[0] == 0 and r[1] == 61 and r[2] == nm][0]) for nm in paper}
assert all(abs(check[nm] - paper[nm]) < 6e-5 for nm in paper), check

# ---------------------------------------------------------------- Monte Carlo designs
ii, jj = np.divmod(np.arange(N), q)
lines = [ii, jj, (ii + jj) % q, (ii - jj) % q]


def load(path):
    with gzip.open(path, "rt") as fh:
        return [np.array(l.split(), dtype=np.int64) for l in fh]


for des in ("lpm2", "lcubeS2", "cubeS2"):
    for k in KS:
        path = os.path.join(ROOT, "results", f"samples_{des}_k{k}.txt.gz")
        if not os.path.exists(path) or os.environ.get("EXACT_ONLY"):
            print("skipped", path)
            continue
        smp = load(path)
        B = len(smp)
        sizes = np.array([len(s) for s in smp])
        ptr = np.concatenate([[0], np.cumsum(sizes)])
        S = sparse.csr_matrix((np.ones(ptr[-1]), np.concatenate(smp), ptr), shape=(B, N))
        est = (S @ Ys.T) / k                                          # (B, n_surf), HT with fixed k
        se2 = (est - Ybar) ** 2
        mse, mse_se = se2.mean(axis=0), se2.std(axis=0, ddof=1) / np.sqrt(B)
        # balance and main-effect weight on S4 (t_M from stratum counts, Section S2)
        tM = np.zeros(B)
        bal2 = np.zeros(B, dtype=bool)
        for b, s in enumerate(smp):
            cnt = [np.bincount(v[s], minlength=q) for v in lines]
            tM[b] = sum(q * (c ** 2).sum() - len(s) ** 2 for c in cnt) / k ** 2
            bal2[b] = all(c.min() == c.max() for c in cnt[:2])
        tau = tM.mean() / (N / k - 1)
        props[f"{des}_k{k}"] = {"B": B, "size_exceptions": int((sizes != k).sum()),
                                "share_balanced_S2": float(bal2.mean()), "tau_S4": float(tau),
                                "W_Y1_S4_rel_lower_bound": float(1 + (nR4 / d4 - 1) * tau)}
        for i in range(n_surf):
            rows.append((i, k, des, mse[i] / srs[k][i], mse_se[i] / srs[k][i]))

# ---------------------------------------------------------------- output
with open(os.path.join(ROOT, "results", "deff_long.csv"), "w") as fh:
    fh.write("surface,source,group,k,design,deff,se\n")
    for i, k, des, v, s in rows:
        m = meta[i]
        fh.write(f"{m['id']},{m['source']},{m['group']},{k},{des},{v:.6g},{s:.3g}\n")
with open(os.path.join(ROOT, "results", "design_properties.json"), "w") as fh:
    json.dump(props, fh, indent=1)

tab = {}
for i, k, des, v, s in rows:
    tab[(i, k, des)] = v
src = np.array([m["source"] for m in meta])
grp = np.array([m["group"] for m in meta])
out = ["# Multi-surface test: results", "",
       f"Volcano check against Table 3 (k=61): {json.dumps({k: round(v, 4) for k, v in check.items()})}", ""]
for k in KS:
    have = all((0, k, d) in tab for d in ("lpm2", "lcubeS2"))
    if not have:
        out.append(f"k={k}: Monte Carlo samples missing")
        continue
    comp = np.array([min(tab[(i, k, "lpm2")], tab[(i, k, "lcubeS2")]) for i in range(n_surf)])
    b16 = np.array([tab[(i, k, "best16_S4")] for i in range(n_surf)])
    ratio = b16 / comp
    e = src == "elevation"
    win, med = float((ratio[e] < 1).mean()), float(np.median(ratio[e]))
    verdict = "holds" if win >= 0.8 and med <= 0.5 else ("partial" if win > 0.5 else "fails")
    out += [f"## k = {k}", "",
            f"Real elevation blocks ({e.sum()}): best-16 beats the better of lpm2 and lcube-S2 on {win:.1%}, "
            f"median ratio {med:.3f}, quartiles {np.quantile(ratio[e], 0.25):.3f} / {np.quantile(ratio[e], 0.75):.3f}. "
            f"**Verdict: {verdict}.**", "",
            f"Volcano: ratio {ratio[0]:.3f}.", "",
            "| Group | n | eqw S4 | best-4 | best-16 | lpm2 | lcube S2 | cube S2 | best-16 / competitor (median) | best-16 wins |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for g in list(dict.fromkeys(grp)):
        sel = np.nonzero(grp == g)[0]

        def med_(d):
            return np.median([tab[(i, k, d)] for i in sel]) if (sel[0], k, d) in tab else float("nan")
        out.append(f"| {g} | {len(sel)} | {med_('eqw_S4'):.4f} | {med_('best4_S4'):.4f} | {med_('best16_S4'):.4f} | "
                   f"{med_('lpm2'):.4f} | {med_('lcubeS2'):.4f} | {med_('cubeS2'):.4f} | "
                   f"{np.median(ratio[sel]):.3f} | {(ratio[sel] < 1).mean():.0%} |")
    out.append("")
out += ["## Design properties", "", "```", json.dumps(props, indent=1), "```"]
with open(os.path.join(ROOT, "results", "summary.md"), "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\n".join(out[:12]))

# ---------------------------------------------------------------- Table S8 of the supplementary material (v4)
if all((0, k, d) in tab for k in KS for d in ("lpm2", "lcubeS2")):
    def ratio_(i, k):
        return tab[(i, k, "best16_S4")] / min(tab[(i, k, "lpm2")], tab[(i, k, "lcubeS2")])

    def line(label, sel):
        med = lambda d, k=61: np.median([tab[(i, k, d)] for i in sel])
        rat = [np.median([ratio_(i, k) for i in sel]) for k in KS]
        win = [sum(ratio_(i, k) < 1 for i in sel) for k in KS]
        return (f"| {label} | {len(sel)} | {med('eqw_S4'):.3f} | {med('best16_S4'):.3f} | "
                f"{med('lpm2', 61):.3f} / {med('lpm2', 183):.3f} | {med('lcubeS2', 61):.3f} / {med('lcubeS2', 183):.3f} | "
                f"{rat[0]:.2f} / {rat[1]:.2f} | {win[0]} / {win[1]} |")
    s8 = ["| Terrain | Blocks | Equal weight, $\\mathcal S_4$ | Best 16, $\\mathcal S_4$ | Local pivotal, $k=61$ / $k=183$ | "
          "Local cube, $\\mathcal S_2$, $k=61$ / $k=183$ | Ratio, $k=61$ / $k=183$ | Best 16 smaller, $k=61$ / $k=183$ |",
          "|---|---|---|---|---|---|---|---|", line("Maunga Whau (Section 5.2)", [0])]
    for g in dict.fromkeys(grp[src == "elevation"]):
        s8.append(line(g, list(np.nonzero(grp == g)[0])))
    s8.append(line("All terrains", list(np.nonzero(src == "elevation")[0])))
    with open(os.path.join(ROOT, "results", "tableS8.md"), "w") as fh:
        fh.write("\n".join(s8) + "\n")
