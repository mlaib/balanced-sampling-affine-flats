# Monte Carlo designs on the 61 x 61 grid: cube (balancing on stratification indicators),
# local pivotal method (lpm2, spatial spread) and local cube (lcube, spread + balance).
# Usage: Rscript mc_grid.R <k1,k2> <B> <out.csv> [seed] [designs]
.libPaths(c("~/R/library", .libPaths()))
suppressMessages(library(BalancedSampling))
a <- commandArgs(trailingOnly = TRUE)
ks <- as.integer(strsplit(a[1], ",")[[1]]); B <- as.integer(a[2]); out <- a[3]
seed <- if (length(a) >= 4) as.integer(a[4]) else 1L
which_designs <- if (length(a) >= 5) strsplit(a[5], ",")[[1]] else c("cubeS2", "cubeS4", "lpm2", "lcubeS4")
set.seed(seed)
q <- 61L
X <- as.matrix(read.csv("data/grid_X.csv", header = FALSE))   # 4 factors x 60 indicators
y <- scan("data/grid_y.csv", quiet = TRUE)
N <- length(y); Ybar <- mean(y)
coords <- cbind(rep(0:(q - 1), each = q), rep(0:(q - 1), times = q))   # (row i, column j), code = 61 i + j
S2cols <- 1:120; S4cols <- 1:240
imbalance <- function(s, cols, k) {
  cnt <- colSums(X[s, cols, drop = FALSE]); cnt <- matrix(cnt, nrow = q - 1)
  cnt <- rbind(length(s) - colSums(cnt), cnt)
  max(abs(cnt - k / q))
}
res <- list()
for (k in ks) {
  prob <- rep(k / N, N)
  for (des in which_designs) {
    est <- numeric(B); dev2 <- dev4 <- numeric(B)
    t0 <- proc.time()[3]
    for (b in 1:B) {
      perm <- sample.int(N)
      s <- switch(des,
        cubeS2 = perm[cube(prob[perm], cbind(prob, X[, S2cols])[perm, ])],
        cubeS4 = perm[cube(prob[perm], cbind(prob, X[, S4cols])[perm, ])],
        lpm2 = perm[lpm2(prob[perm], coords[perm, ])],
        lcubeS4 = perm[lcube(prob[perm], coords[perm, ], cbind(prob, X[, S4cols])[perm, ])])
      est[b] <- mean(y[s])
      dev2[b] <- imbalance(s, S2cols, k); dev4[b] <- imbalance(s, S4cols, k)
    }
    el <- proc.time()[3] - t0
    se2 <- (est - Ybar)^2
    res[[length(res) + 1]] <- data.frame(k = k, design = des, B = B, mse = mean(se2),
      mse_se = sd(se2) / sqrt(B), prop_exact_S2 = mean(dev2 < 1e-9), prop_exact_S4 = mean(dev4 < 1e-9),
      mean_maxdev_S4 = mean(dev4), sec_per_sample = el / B)
    cat(sprintf("k=%d %s: mse=%.4f (se %.4f) exactS2=%.3f exactS4=%.3f %.3fs/sample\n", k, des,
                mean(se2), sd(se2) / sqrt(B), mean(dev2 < 1e-9), mean(dev4 < 1e-9), el / B))
    write.csv(do.call(rbind, res), out, row.names = FALSE)
  }
}
