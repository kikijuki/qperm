# Pure-R emulation of the maximum-likelihood quantum amplitude estimation (MLQAE)
# estimator. This is a simulation of the estimator's statistical behaviour, not a faster
# computation: it draws its measurement counts from the true amplitude (the corrected p
# that the classical engine already computed), so it demonstrates the 1/N versus
# 1/sqrt(N) query-complexity scaling rather than producing an independent or faster
# p-value. See the package vignette for the honest scoping.

# Maximum-likelihood fit of the Grover angle theta from measurement records.
# records: matrix with columns m (Grover power), good (successes), shots.
.qperm_mle_theta <- function(records, grid = 2000L) {
  th <- seq(0, pi / 2, length.out = grid + 1L)[-1L]
  best <- th[1]; bestll <- -Inf
  m <- records[, "m"]; good <- records[, "good"]; shots <- records[, "shots"]
  for (t in th) {
    s2 <- pmin(pmax(sin((2 * m + 1) * t)^2, 1e-12), 1 - 1e-12)
    ll <- sum(good * log(s2) + (shots - good) * log(1 - s2))
    if (ll > bestll) { bestll <- ll; best <- t }
  }
  best
}

# Emulate MLQAE for a true amplitude p_true. Returns the final estimate and a
# convergence trace (cumulative oracle queries versus absolute error) that grows as
# each Grover power is added to the schedule.
.qperm_emulate <- function(p_true, powers = c(0, 1, 2, 3, 4, 6, 8, 12, 16, 24),
                           shots = 100L) {
  p_true <- min(max(p_true, 0), 1)
  theta_true <- asin(sqrt(p_true))
  records <- matrix(nrow = 0L, ncol = 3L,
                    dimnames = list(NULL, c("m", "good", "shots")))
  queries <- 0
  q_vec <- numeric(length(powers)); err_vec <- numeric(length(powers))
  est <- p_true
  for (k in seq_along(powers)) {
    m <- powers[k]
    prob <- sin((2 * m + 1) * theta_true)^2
    good <- stats::rbinom(1L, shots, min(max(prob, 0), 1))
    records <- rbind(records, c(m = m, good = good, shots = shots))
    queries <- queries + shots * (2 * m + 1)
    est <- sin(.qperm_mle_theta(records))^2
    q_vec[k] <- queries
    err_vec[k] <- abs(est - p_true)
  }
  list(estimate = est,
       trace = data.frame(queries = q_vec, error = err_vec, power = powers),
       shots = shots)
}

# Classical Monte Carlo error at a set of query budgets, for the same amplitude.
# Each classical query is one Bernoulli(p_true) draw, so the RMSE falls as 1/sqrt(M).
.qperm_classical_trace <- function(p_true, budgets, reps = 40L) {
  p_true <- min(max(p_true, 0), 1)
  err <- vapply(budgets, function(M) {
    est <- stats::rbinom(reps, M, p_true) / M
    sqrt(mean((est - p_true)^2))
  }, numeric(1))
  data.frame(queries = budgets, error = err)
}

# A geometric budget schedule spanning roughly the quantum query range, for a fair
# same-axis comparison in the scaling plot.
.qperm_budgets <- function(qmax) {
  hi <- max(3, log10(max(qmax, 1000)))
  unique(round(10^seq(1.3, hi, by = 0.23)))
}
