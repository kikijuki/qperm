# Brute-force exact corrected p for the max |mean| statistic, computed independently of
# the package internals, so the test is a genuine cross-check.
brute_force_p <- function(Y, alternative = "two.sided") {
  n <- nrow(Y); N <- 2^n; nullmax <- numeric(N)
  orient <- switch(alternative, two.sided = abs, greater = identity,
                   less = function(x) -x)
  for (j in 0:(N - 1)) {
    s <- ifelse(bitwAnd(bitwShiftR(j, 0:(n - 1)), 1L) == 1L, -1, 1)
    nullmax[j + 1] <- max(orient(colMeans(s * Y)))
  }
  obs <- max(orient(colMeans(Y)))
  mean(nullmax >= obs)
}

test_that("classical max-statistic p matches brute-force enumeration", {
  data(qae_demo, package = "qperm")
  fit <- qperm_test(qae_demo, method = "classical")
  expect_equal(fit$p_corrected, brute_force_p(qae_demo), tolerance = 1e-12)
})

test_that("classical p matches brute force for one-sided alternatives", {
  data(qae_demo, package = "qperm")
  for (alt in c("greater", "less")) {
    fit <- qperm_test(qae_demo, method = "classical", alternative = alt)
    expect_equal(fit$p_corrected, brute_force_p(qae_demo, alt), tolerance = 1e-12)
  }
})

test_that("corrected p is a valid probability and peak is the map minimum", {
  data(qae_demo, package = "qperm")
  fit <- qperm_test(qae_demo, method = "classical")
  expect_true(fit$p_corrected >= 0 && fit$p_corrected <= 1)
  expect_equal(fit$p_corrected, min(fit$p_map))
})

test_that("quantum emulation recovers the exact corrected p", {
  data(qae_demo, package = "qperm")
  set.seed(1)
  fit <- qperm_test(qae_demo, method = "quantum", seed = 1, shots = 400L)
  expect_lt(abs(fit$q_estimate - fit$p_true), 0.03)
  expect_true(all(c("queries", "error") %in% names(fit$q_trace)))
})

test_that("resource estimate reflects the quadratic-in-both scaling", {
  r1 <- qperm_resources(20, 1000, 0.01)
  r2 <- qperm_resources(20, 4000, 0.01)   # 4x voxels
  # classical grows linearly in V, quantum as sqrt(V), so the ratio grows
  expect_gt(r2$speedup_factor, r1$speedup_factor)
  expect_true(r1$classical_queries > r1$quantum_queries)
})

test_that("cluster-extent finds a planted contiguous cluster", {
  set.seed(3)
  g <- 6
  coords <- as.matrix(expand.grid(x = 1:g, y = 1:g))
  Y <- matrix(rnorm(10 * nrow(coords)), 10)
  block <- which(coords[, 1] <= 2 & coords[, 2] <= 2)   # a 2x2 = 4-voxel block
  Y[, block] <- Y[, block] + 2
  fit <- qperm_test(Y, method = "classical", stat = "cluster",
                    threshold = 0.5, coords = coords, seed = 3)
  expect_gte(fit$observed_peak, length(block))
  expect_true(fit$p_corrected >= 0 && fit$p_corrected <= 1)
})

test_that("cluster extent requires coordinates", {
  data(qae_demo, package = "qperm")
  expect_error(qperm_test(qae_demo, stat = "cluster", threshold = 0.5),
               "coords")
})

test_that("bad input is rejected", {
  expect_error(qperm_test(matrix("a", 2, 2)), "numeric")
  expect_error(qperm_test(matrix(1, 1, 3)), "2 subjects")
})
