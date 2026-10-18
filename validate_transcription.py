"""Transition-level check of audit_cube.py against cube() of BalancedSampling 2.1.1 (Supplementary material S4).

For each seed, R draws a sample with cube() after set.seed(seed), and saves the uniform numbers runif()
returns after the same set.seed(seed).  The transcription consumes those numbers in the package's order
(the shuffle of the units, one number per flight or landing step, one for the last unit) and its sample is
compared with the package's.  Balancing matrices: the rows and columns of Section 5.2 (flight with the
graph null vectors), and the one-variable indicators of the factorial frames of Section 5.1 (flight with
the row-echelon null vectors, as in cube_factorial.py).
Writes results/transcription_check.csv.  Usage: python3 validate_transcription.py [seeds per configuration]
"""
import csv, os, subprocess, sys, tempfile
import numpy as np
import audit_cube as ac

OUT = os.path.dirname(os.path.abspath(__file__))
R_CODE = r'''
.libPaths(c("~/R/library", .libPaths())); suppressMessages(library(BalancedSampling))
a <- commandArgs(TRUE); k <- as.integer(a[1]); seed <- as.integer(a[2]); out <- a[3]; xfile <- a[4]; ncol_x <- as.integer(a[5])
X <- as.matrix(read.csv(xfile, header = FALSE))[, 1:ncol_x, drop = FALSE]; N <- nrow(X); prob <- rep(k / N, N)
set.seed(seed); s <- cube(prob, cbind(prob, X))
set.seed(seed); U <- runif(50000)
writeLines(as.character(s), paste0(out, "_s.txt")); writeBin(U, paste0(out, "_u.bin"))
'''

CONFIGS = [("61x61 grid, rows and columns", "data/grid_X.csv", 120, (61, 183), True),
           ("2^7, one-variable margins", "data/fact_p2_m7_X.csv", 7, (8, 16, 32), False),
           ("3^5, one-variable margins", "data/fact_p3_m5_X.csv", 10, (9, 27, 81), False)]

if __name__ == "__main__":
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        rfile = os.path.join(tmp, "draw.R")
        open(rfile, "w").write(R_CODE)
        for label, xfile, ncol_x, ks, graph in CONFIGS:
            X = np.loadtxt(f"{OUT}/{xfile}", delimiter=",")[:, :ncol_x]
            N = len(X)
            for k in ks:
                prob = np.full(N, k / N)
                for seed in range(1, n_seeds + 1):
                    base = os.path.join(tmp, f"k{k}_s{seed}")
                    subprocess.run(["Rscript", rfile, str(k), str(seed), base, xfile, str(ncol_x)], cwd=OUT, check=True)
                    s_pkg = sorted(int(v) - 1 for v in open(base + "_s.txt").read().split())
                    stream = ac.UniformStream(np.fromfile(base + "_u.bin", dtype="<f8"))
                    c = ac.Cube(prob, np.c_[prob, X], None, stream=stream)
                    c.flight(graph=graph)
                    s_tr = sorted(int(v) for v in c.landing())
                    rows.append({"balancing": label, "k": k, "seed": seed, "package_size": len(s_pkg),
                                 "transcription_size": len(s_tr), "identical": s_pkg == s_tr, "uniforms_used": stream.i})
                    print(rows[-1], flush=True)
    with open(f"{OUT}/results/transcription_check.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    for label, *_ in CONFIGS:
        sub = [r for r in rows if r["balancing"] == label]
        print(label, ":", sum(r["identical"] for r in sub), "of", len(sub), "samples identical")
