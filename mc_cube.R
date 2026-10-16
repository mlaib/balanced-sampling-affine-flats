# Monte Carlo comparison: cube method (equal probabilities) balancing on level indicators.
# Usage: Rscript mc_cube.R <X.csv> <Y.csv> <q> <nfactors> <k1,k2,...> <B> <out.csv> [seed] [samples.txt.gz]
# X.csv: N x (nfactors*(q-1)) 0/1 indicators (level 0 of each factor dropped).
# Also splits each squared error into the parts due to the main effects (zero for balanced samples),
# the residual and their cross term, y - Ybar = mu + rho with mu the population main-effects fit, and
# records the main-effect weight t_M = sum over main-effect frequencies of |eta_s(u)|^2, which for the
# one-variable stratifications is (q/k^2) times the sum of squared deviations of the level counts from k/q.
# With a ninth argument, writes the selected units of every draw (one line per draw) to that file.
.libPaths(c("~/R/library", .libPaths()))
suppressMessages(library(BalancedSampling))
a <- commandArgs(trailingOnly = TRUE)
X <- as.matrix(read.csv(a[1], header = FALSE))
Y <- as.matrix(read.csv(a[2], header = TRUE))
q <- as.integer(a[3]); nf <- as.integer(a[4])
ks <- as.integer(strsplit(a[5], ",")[[1]]); B <- as.integer(a[6]); out <- a[7]
seed <- if (length(a) >= 8) as.integer(a[8]) else 1L
samples_out <- if (length(a) >= 9) a[9] else NA
if (!is.na(samples_out)) { con <- gzfile(samples_out, "w") }
set.seed(seed)
N <- nrow(X)
Z <- cbind(1, X)
M <- apply(Y, 2, function(y) drop(Z %*% qr.coef(qr(Z), y)) - mean(y))   # main-effect part, centred
Rho <- sweep(Y, 2, colMeans(Y)) - M                                         # residual part
res <- list()
for (k in ks) {
  prob <- rep(k / N, N)
  xb <- cbind(prob, X)
  means <- matrix(NA_real_, B, ncol(Y)); exact <- logical(B); maxdev <- numeric(B); tM <- numeric(B)
  em <- er <- matrix(NA_real_, B, ncol(Y))
  t0 <- proc.time()[3]
  for (b in 1:B) {
    perm <- sample.int(N)                     # random unit order before each draw
    s <- perm[cube(prob[perm], xb[perm, , drop = FALSE])]
    means[b, ] <- colMeans(Y[s, , drop = FALSE])
    em[b, ] <- colMeans(M[s, , drop = FALSE]); er[b, ] <- colMeans(Rho[s, , drop = FALSE])
    # level counts for every factor (level 0 count = k - sum of the others)
    cnt <- colSums(X[s, , drop = FALSE])
    cnt <- matrix(cnt, nrow = q - 1)
    cnt <- rbind(length(s) - colSums(cnt), cnt)
    dev <- abs(cnt - k / q)
    exact[b] <- (length(s) == k) && all(dev < 1e-9)
    maxdev[b] <- max(dev)
    tM[b] <- q / k^2 * sum((cnt - k / q)^2)
    if (!is.na(samples_out)) writeLines(paste(c(k, sort(s) - 1L), collapse = " "), con)
  }
  el <- proc.time()[3] - t0
  ybar <- colMeans(Y)
  mse <- colMeans((means - matrix(ybar, B, ncol(Y), byrow = TRUE))^2)
  se <- apply((means - matrix(ybar, B, ncol(Y), byrow = TRUE))^2, 2, sd) / sqrt(B)
  err2 <- (means - matrix(ybar, B, ncol(Y), byrow = TRUE))^2
  stopifnot(max(abs((em + er)^2 - err2)) < 1e-9)   # the three parts add up to every squared error
  for (j in seq_len(ncol(Y))) {
    res[[length(res) + 1]] <- data.frame(k = k, response = colnames(Y)[j], mse = mse[j], mse_se = se[j],
                                         prop_exact = mean(exact), mean_maxdev = mean(maxdev),
                                         B = B, sec_per_sample = el / B,
                                         main = mean(em[, j]^2), main_se = sd(em[, j]^2) / sqrt(B),
                                         resid = mean(er[, j]^2), resid_se = sd(er[, j]^2) / sqrt(B),
                                         cross = mean(2 * em[, j] * er[, j]), cross_se = sd(2 * em[, j] * er[, j]) / sqrt(B),
                                         tM = mean(tM), tM_se = sd(tM) / sqrt(B), tau = mean(tM) / (N / k - 1))
  }
  cat(sprintf("k=%d: exact balance in %.3f of %d samples; %.4f s/sample\n", k, mean(exact), B, el / B))
}
if (!is.na(samples_out)) close(con)
write.csv(do.call(rbind, res), out, row.names = FALSE)
