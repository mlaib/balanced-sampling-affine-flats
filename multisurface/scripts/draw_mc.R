# Monte Carlo samples on the 61 x 61 grid, stored once and evaluated on every surface (evaluate.py).
# Designs: lpm2 (planar spread), lcubeS2 (planar spread + balance on rows and columns), cubeS2 (balance on rows
# and columns). Units are coded 61 i + j (row i, column j), as in study_grid.py; prob is a balancing variable and
# one category of each stratification is dropped, as in mc_grid.R. Units are randomly permuted before each draw.
# Usage: Rscript draw_mc.R <design> <k> <B> <seed> <out.txt.gz>
# Output: one line per draw with the 0-based codes of the selected units.
.libPaths(c("~/R/library", .libPaths()))
suppressMessages(library(BalancedSampling))
a <- commandArgs(trailingOnly = TRUE)
des <- a[1]; k <- as.integer(a[2]); B <- as.integer(a[3]); seed <- as.integer(a[4]); out <- a[5]
set.seed(seed)
q <- 61L; N <- q * q
coords <- cbind(rep(0:(q - 1), each = q), rep(0:(q - 1), times = q))
ind <- function(v) sapply(1:(q - 1), function(l) as.integer(v == l))
XS2 <- cbind(ind(coords[, 1]), ind(coords[, 2]))
prob <- rep(k / N, N)
con <- gzfile(out, "w")
t0 <- proc.time()[3]
for (b in 1:B) {
  perm <- sample.int(N)
  s <- switch(des,
    lpm2 = perm[lpm2(prob[perm], coords[perm, ])],
    lcubeS2 = perm[lcube(prob[perm], coords[perm, ], cbind(prob, XS2)[perm, ])],
    cubeS2 = perm[cube(prob[perm], cbind(prob, XS2)[perm, ])])
  writeLines(paste(sort(s) - 1L, collapse = " "), con)
}
close(con)
cat(sprintf("%s k=%d B=%d seed=%d: %.3f s/draw\n", des, k, B, seed, (proc.time()[3] - t0) / B))
