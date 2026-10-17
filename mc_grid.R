# Monte Carlo designs on the 61 x 61 grid: cube (balancing on stratification indicators),
# local pivotal method (lpm2, spatial spread; lpm2torus with wrapped distances, through the embedding
# of the torus in R^4 by the angles 2 pi i/61 and 2 pi j/61) and local cube (lcube, spread + balance).
# Also records the main-effect weight t_M = sum over the main-effect frequencies of |eta_s(u)|^2, which for
# line stratifications is (q/k^2) times the sum of squared deviations of the stratum counts from k/q.
# Every draw also records its size, the Horvitz-Thompson estimate sum(y_s)/k (with the fixed denominator k) and
# t_M computed as sum over main-effect frequencies of |eta_s(u)|^2 with eta_s(u) = k^-1 sum_{x in s} chi_u(x),
# that is (q sum_b n_b^2 - |s|^2)/k^2 for each line stratification, which equals the expression above when |s| = k.
# Usage: Rscript mc_grid.R <k1,k2> <B> <out.csv> [seed] [designs] [samples.txt.gz]
.libPaths(c("~/R/library", .libPaths()))
suppressMessages(library(BalancedSampling))
a <- commandArgs(trailingOnly = TRUE)
ks <- as.integer(strsplit(a[1], ",")[[1]]); B <- as.integer(a[2]); out <- a[3]
seed <- if (length(a) >= 4) as.integer(a[4]) else 1L
which_designs <- if (length(a) >= 5) strsplit(a[5], ",")[[1]] else c("cubeS2", "cubeS4", "lpm2", "lcubeS4")
samples_out <- if (length(a) >= 6) a[6] else NA
if (!is.na(samples_out)) con <- gzfile(samples_out, "w")
set.seed(seed)
q <- 61L
X <- as.matrix(read.csv("data/grid_X.csv", header = FALSE))   # 4 factors x 60 indicators
y <- scan("data/grid_y.csv", quiet = TRUE)
N <- length(y); Ybar <- mean(y)
coords <- cbind(rep(0:(q - 1), each = q), rep(0:(q - 1), times = q))   # (row i, column j), code = 61 i + j
ang <- 2 * pi * coords / q
torus <- cbind(cos(ang[, 1]), sin(ang[, 1]), cos(ang[, 2]), sin(ang[, 2])) * q / (2 * pi)   # wrapped geometry
tm <- function(s, k) {   # main-effect weights on rows + columns (S2) and on all four line families (S4)
  i <- coords[s, 1]; j <- coords[s, 2]
  d <- sapply(list(i, j, (i + j) %% q, (i - j) %% q), function(v) q * sum(tabulate(v + 1, q)^2) - length(s)^2) / k^2
  c(sum(d[1:2]), sum(d))
}
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
    est <- numeric(B); dev2 <- dev4 <- numeric(B); tM <- matrix(0, B, 2); ht <- numeric(B); size <- integer(B)
    t0 <- proc.time()[3]
    for (b in 1:B) {
      perm <- sample.int(N)
      s <- switch(des,
        cubeS2 = perm[cube(prob[perm], cbind(prob, X[, S2cols])[perm, ])],
        cubeS4 = perm[cube(prob[perm], cbind(prob, X[, S4cols])[perm, ])],
        lpm2 = perm[lpm2(prob[perm], coords[perm, ])],
        lpm2torus = perm[lpm2(prob[perm], torus[perm, ])],
        lcubeS4 = perm[lcube(prob[perm], coords[perm, ], cbind(prob, X[, S4cols])[perm, ])])
      est[b] <- mean(y[s]); ht[b] <- sum(y[s]) / k; size[b] <- length(s)
      if (!is.na(samples_out)) writeLines(paste(c(k, des, sort(s) - 1L), collapse = " "), con)
      dev2[b] <- imbalance(s, S2cols, k); dev4[b] <- imbalance(s, S4cols, k)
      tM[b, ] <- tm(s, k)
    }
    el <- proc.time()[3] - t0
    se2 <- (est - Ybar)^2; seht <- (ht - Ybar)^2
    res[[length(res) + 1]] <- data.frame(k = k, design = des, B = B, mse = mean(se2),
      mse_se = sd(se2) / sqrt(B), prop_exact_S2 = mean(dev2 < 1e-9), prop_exact_S4 = mean(dev4 < 1e-9),
      mean_maxdev_S4 = mean(dev4), sec_per_sample = el / B,
      tM_S2 = mean(tM[, 1]), tM_S2_se = sd(tM[, 1]) / sqrt(B), tM_S4 = mean(tM[, 2]), tM_S4_se = sd(tM[, 2]) / sqrt(B),
      tau_S2 = mean(tM[, 1]) / (N / k - 1), tau_S4 = mean(tM[, 2]) / (N / k - 1),
      tau_S4_se = sd(tM[, 2]) / sqrt(B) / (N / k - 1), size_exceptions = sum(size != k),
      size_min = min(size), size_max = max(size), mse_ht = mean(seht), mse_ht_se = sd(seht) / sqrt(B))
    cat(sprintf("k=%d %s: mse=%.4f (se %.4f) exactS2=%.3f exactS4=%.3f %.3fs/sample\n", k, des,
                mean(se2), sd(se2) / sqrt(B), mean(dev2 < 1e-9), mean(dev4 < 1e-9), el / B))
    write.csv(do.call(rbind, res), out, row.names = FALSE)
  }
}
if (!is.na(samples_out)) close(con)
