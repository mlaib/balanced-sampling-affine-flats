# Monte Carlo comparison: cube method (equal probabilities) balancing on level indicators.
# Usage: Rscript mc_cube.R <X.csv> <Y.csv> <q> <nfactors> <k1,k2,...> <B> <out.csv> [seed]
# X.csv: N x (nfactors*(q-1)) 0/1 indicators (level 0 of each factor dropped).
.libPaths(c("~/R/library", .libPaths()))
suppressMessages(library(BalancedSampling))
a <- commandArgs(trailingOnly = TRUE)
X <- as.matrix(read.csv(a[1], header = FALSE))
Y <- as.matrix(read.csv(a[2], header = TRUE))
q <- as.integer(a[3]); nf <- as.integer(a[4])
ks <- as.integer(strsplit(a[5], ",")[[1]]); B <- as.integer(a[6]); out <- a[7]
seed <- if (length(a) >= 8) as.integer(a[8]) else 1L
set.seed(seed)
N <- nrow(X)
res <- list()
for (k in ks) {
  prob <- rep(k / N, N)
  xb <- cbind(prob, X)
  means <- matrix(NA_real_, B, ncol(Y)); exact <- logical(B); maxdev <- numeric(B)
  t0 <- proc.time()[3]
  for (b in 1:B) {
    perm <- sample.int(N)                     # random unit order before each draw
    s <- perm[cube(prob[perm], xb[perm, , drop = FALSE])]
    means[b, ] <- colMeans(Y[s, , drop = FALSE])
    # level counts for every factor (level 0 count = k - sum of the others)
    cnt <- colSums(X[s, , drop = FALSE])
    cnt <- matrix(cnt, nrow = q - 1)
    cnt <- rbind(length(s) - colSums(cnt), cnt)
    dev <- abs(cnt - k / q)
    exact[b] <- (length(s) == k) && all(dev < 1e-9)
    maxdev[b] <- max(dev)
  }
  el <- proc.time()[3] - t0
  ybar <- colMeans(Y)
  mse <- colMeans((means - matrix(ybar, B, ncol(Y), byrow = TRUE))^2)
  se <- apply((means - matrix(ybar, B, ncol(Y), byrow = TRUE))^2, 2, sd) / sqrt(B)
  for (j in seq_len(ncol(Y))) {
    res[[length(res) + 1]] <- data.frame(k = k, response = colnames(Y)[j], mse = mse[j], mse_se = se[j],
                                         prop_exact = mean(exact), mean_maxdev = mean(maxdev),
                                         B = B, sec_per_sample = el / B)
  }
  cat(sprintf("k=%d: exact balance in %.3f of %d samples; %.4f s/sample\n", k, mean(exact), B, el / B))
}
write.csv(do.call(rbind, res), out, row.names = FALSE)
