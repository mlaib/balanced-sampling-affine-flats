"""Figure 1 (Section 5.2): coset samples on the 61 x 61 elevation grid.

(a) a flat whose direction is admissible for S4 but whose dual has the short vector (1, 2), the last-ranked
direction of results/grid_lattice_ranking.csv; (b) a flat with the first-ranked direction, dual (1, 8);
(c) three randomly selected cosets of that direction (k = 183). Each flat is {(i, j): i + s j = c mod 61} for the
dual basis (1, s); the elevations are drawn as contours every 15 m. Writes results/fig1_coset_samples.pdf and .png,
and a compact 4-level copy, fig1_coset_samples_web.png, for BrainLog.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
q = 61
Y = np.loadtxt(os.path.join(HERE, "data", "volcano61.csv"), delimiter=",")
with open(os.path.join(HERE, "results", "grid_lattice_ranking.csv")) as fh:
    slopes = [int(r["slope_s"]) for r in csv.DictReader(fh)]                 # in rank order
s_best, s_worst = slopes[0], slopes[-4]      # dual bases (1, 8), first, and (1, 2), first of the four of length sqrt 5
assert (s_best, s_worst) == (8, 2)

rng = np.random.default_rng(20261019)
c_one = int(rng.integers(q))
c_three = np.sort(rng.choice(q, 3, replace=False))
ii, jj = np.divmod(np.arange(q * q), q)


def flat(s, cs):
    keep = np.isin((ii + s * jj) % q, cs)
    return ii[keep], jj[keep]


panels = [("(a) Dual vector (1, 2), $k=61$", flat(s_worst, [c_one])),
          ("(b) Dual vector (1, 8), $k=61$", flat(s_best, [c_one])),
          ("(c) Dual vector (1, 8), $k=183$", flat(s_best, c_three))]
for _, (i, j) in panels:          # exactly balanced on rows, columns and both wrapped diagonals
    for lab in (i, j, (i + j) % q, (i - j) % q):
        cnt = np.bincount(lab, minlength=q)
        assert cnt.min() == cnt.max()

INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
CONTOUR = "#a9a8a3"                                                          # light, recessive elevation contours
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif", "font.size": 8, "axes.titlesize": 8,
                     "axes.edgecolor": INK2, "axes.linewidth": 0.5, "text.color": INK})
fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.35), constrained_layout=True)
for ax, (title, (i, j)) in zip(axes, panels):
    ax.contour(np.arange(q), np.arange(q), Y, levels=np.arange(100, 200, 15), colors=CONTOUR, linewidths=0.5)
    ax.set_xlim(-0.5, q - 0.5)
    ax.set_ylim(q - 0.5, -0.5)
    ax.scatter(j, i, s=7, c=INK, edgecolors=SURFACE, linewidths=0.4, zorder=3)
    ax.set_title(title, loc="left", color=INK)
    ax.set_xticks([0, 30, 60])
    ax.set_yticks([0, 30, 60])
    ax.tick_params(length=2, width=0.5, colors=INK2, labelsize=7)
    ax.set_xlabel("column $j$", color=INK2, fontsize=7)
    ax.set_aspect("equal")
axes[0].set_ylabel("row $i$", color=INK2, fontsize=7)
fig.savefig(os.path.join(HERE, "results", "fig1_coset_samples.pdf"), facecolor="white")
fig.savefig(os.path.join(HERE, "results", "fig1_coset_samples.png"), facecolor="white", dpi=220)
# a compact copy for BrainLog, whose figure upload takes the image inline: 4 fixed colours (2 bits), 90 dpi, no dithering
from PIL import Image
fig.savefig(os.path.join(HERE, "results", "fig1_tmp.png"), facecolor="white", dpi=90)
pal = Image.new("P", (1, 1))
pal.putpalette([255, 255, 255, 0xa9, 0xa8, 0xa3, 0x52, 0x51, 0x4e, 0x0b, 0x0b, 0x0b])
img = Image.open(os.path.join(HERE, "results", "fig1_tmp.png")).convert("RGB")
img.quantize(palette=pal, dither=Image.Dither.NONE).save(os.path.join(HERE, "results", "fig1_coset_samples_web.png"),
                                                         optimize=True, bits=2)
os.remove(os.path.join(HERE, "results", "fig1_tmp.png"))
print("cosets", c_one, c_three.tolist())
