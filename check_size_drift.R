# Appendix A: sizes of the samples returned by cube and lcube of BalancedSampling 2.1.1 on a 19 x 19 grid balanced
# on two, three and four line stratifications (rows, columns, anti-diagonals i+j, diagonals i-j, wrapped), with the
# inclusion probabilities as the first balancing variable, as in mc_grid.R.  With rows and columns the sizes are
# fixed; with three or four stratifications they drift, because the package's row reduction tests pivots for
# exact zero (check_rref_defect.py).  30 draws per case.  Writes results/size_drift_check.csv.
.libPaths(c("~/R/library", .libPaths()))
suppressMessages(library(BalancedSampling))
set.seed(2026)
q <- 19; N <- q * q; i <- (0:(N - 1)) %/% q; j <- (0:(N - 1)) %% q
ind <- function(v) sapply(1:(q - 1), function(b) as.numeric(v == b))
out <- NULL
for (k in c(38, 76, 114)) {
  prob <- rep(k / N, N)
  for (nb in 2:4) {
    X <- cbind(prob, ind(i), ind(j), if (nb >= 3) ind((i + j) %% q), if (nb >= 4) ind((i - j) %% q))
    for (method in c("cube", "lcube")) {
      sz <- replicate(30, length(if (method == "cube") cube(prob, X) else lcube(prob, cbind(i, j), X)))
      out <- rbind(out, data.frame(k = k, stratifications = nb, balancing_variables = ncol(X), method = method, draws = 30,
                                   size_other_than_k = sum(sz != k), size_min = min(sz), size_max = max(sz)))
      cat(sprintf("k=%d strata=%d %s: %d of 30 with size other than k, sizes %d-%d\n", k, nb, method, sum(sz != k), min(sz), max(sz)))
    }
  }
}
write.csv(out, "results/size_drift_check.csv", row.names = FALSE)
