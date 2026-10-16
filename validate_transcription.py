"""Transition-level check of audit_cube.py against cube() of BalancedSampling 2.1.1 (Appendix A).

For each seed, R draws a sample with cube() after set.seed(seed) on the rows-and-columns balancing
matrix of Section 7.2, and saves the uniform numbers runif() returns after the same set.seed(seed).
The transcription consumes those numbers in the package's order (the shuffle of the units, one number
per flight or landing step, one for the last unit) and its sample is compared with the package's.
Writes results/transcription_check.csv.  Usage: python3 validate_transcription.py [seeds per k]
"""
import csv, os, subprocess, sys, tempfile
import numpy as np
import audit_cube as ac

OUT = os.path.dirname(os.path.abspath(__file__))
R_CODE = r'''
.libPaths(c("~/R/library", .libPaths())); suppressMessages(library(BalancedSampling))
a <- commandArgs(TRUE); k <- as.integer(a[1]); seed <- as.integer(a[2]); out <- a[3]
X <- as.matrix(read.csv("data/grid_X.csv", header = FALSE))[, 1:120]; N <- nrow(X); prob <- rep(k / N, N)
set.seed(seed); s <- cube(prob, cbind(prob, X))
set.seed(seed); U <- runif(50000)
writeLines(as.character(s), paste0(out, "_s.txt")); writeBin(U, paste0(out, "_u.bin"))
'''

if __name__ == "__main__":
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    X, y = ac.data()
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        rfile = os.path.join(tmp, "draw.R")
        open(rfile, "w").write(R_CODE)
        for k in (61, 183):
            prob = np.full(len(y), k / len(y))
            for seed in range(1, n_seeds + 1):
                base = os.path.join(tmp, f"k{k}_s{seed}")
                subprocess.run(["Rscript", rfile, str(k), str(seed), base], cwd=OUT, check=True)
                s_pkg = sorted(int(v) - 1 for v in open(base + "_s.txt").read().split())
                stream = ac.UniformStream(np.fromfile(base + "_u.bin", dtype="<f8"))
                c = ac.Cube(prob, np.c_[prob, X], None, stream=stream)
                c.flight()
                s_tr = sorted(int(v) for v in c.landing())
                rows.append({"k": k, "seed": seed, "package_size": len(s_pkg), "transcription_size": len(s_tr),
                             "identical": s_pkg == s_tr, "uniforms_used": stream.i})
                print(rows[-1], flush=True)
    with open(f"{OUT}/results/transcription_check.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(sum(r["identical"] for r in rows), "of", len(rows), "samples identical")
