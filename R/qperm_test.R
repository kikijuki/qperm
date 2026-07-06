#' Sign-flip permutation test with a classical or quantum-emulated estimator
#'
#' Runs a one-sample / paired sign-flip permutation test for group neuroimaging data and
#' reports the family-wise-error (FWER) corrected p-value. The \code{"classical"} method
#' is the usable analysis: it enumerates the sign-flip null exactly when it is small enough
#' and otherwise samples it (Monte Carlo), producing a corrected p-value on your data
#' today. The \code{"quantum"} method offers the proposed quantum amplitude-estimation
#' route as an alternative: by default a faithful emulation of the maximum-likelihood
#' quantum amplitude estimation (MLQAE) estimator that reproduces the same corrected
#' p-value and shows its 1/N versus 1/sqrt(N) query-complexity scaling, and optionally, via
#' \code{backend = "qiskit"}, actual quantum circuits run on a statevector simulator for
#' small studies.
#'
#' The quantum method is a simulation, not a speedup on classical hardware: emulating the
#' sign-flip null is itself the classical computation, and the emulated measurement counts
#' are drawn from the exact amplitude the classical engine computes. It demonstrates the
#' estimator and its scaling; it does not give an independent or faster p-value. See
#' \code{\link{qperm_resources}} for the projected query advantage at larger scales, and
#' the package vignette for the full scoping.
#'
#' The per-voxel statistic is linear in the signs, \eqn{t_v(s) = (1/n) \sum_i s_i y_{iv}},
#' which is what makes the quantum oracle a constant-weighted adder. A variance-normalised
#' (t) statistic would break that linearity and is not used here.
#'
#' @param data A numeric subjects-by-voxels matrix (rows = subjects, columns = voxels), or
#'   the list returned by \code{\link{qperm_read_nifti}}.
#' @param method \code{"classical"} or \code{"quantum"}.
#' @param stat \code{"max"} (max over voxels) or \code{"cluster"} (cluster extent). Cluster
#'   extent needs voxel geometry, so supply \code{coords} or use NIfTI input.
#' @param alternative \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#' @param coords Optional voxels-by-d integer coordinate matrix for cluster extent. Taken
#'   from NIfTI input automatically when available.
#' @param threshold Cluster-forming threshold on the oriented per-voxel statistic. Required
#'   for \code{stat = "cluster"}.
#' @param exact If \code{TRUE} (default) enumerate the full sign-flip null when it has at
#'   most \code{exact_max} members; otherwise use Monte Carlo.
#' @param exact_max Largest exact null to enumerate. Default \code{2^20}.
#' @param n_perm Monte Carlo permutations when the exact null is too large. Default 5000.
#' @param backend For \code{method = "quantum"}: \code{"emulation"} (pure R) or
#'   \code{"qiskit"} (real circuits, small studies, needs reticulate + qiskit).
#' @param powers Grover-power schedule for the MLQAE readout.
#' @param shots Measurement shots per Grover power.
#' @param epsilon Target precision for the accompanying resource estimate.
#' @param seed Optional RNG seed for reproducibility.
#' @return An object of class \code{"qperm"}; see \code{print}, \code{summary}, and
#'   \code{plot} methods.
#' @examples
#' data(qae_demo)
#' # Classical corrected p-value on the bundled 10-subject synthetic data
#' qperm_test(qae_demo, method = "classical")
#' # Quantum-emulated estimator plus the scaling it would follow
#' fit <- qperm_test(qae_demo, method = "quantum", seed = 1)
#' fit
#' @seealso \code{\link{qperm_resources}}, \code{\link{qperm_read_nifti}}
#' @export
qperm_test <- function(data,
                       method = c("classical", "quantum"),
                       stat = c("max", "cluster"),
                       alternative = c("two.sided", "greater", "less"),
                       coords = NULL,
                       threshold = NULL,
                       exact = TRUE, exact_max = 2^20, n_perm = 5000L,
                       backend = c("emulation", "qiskit"),
                       powers = c(0, 1, 2, 3, 4, 6, 8, 12, 16, 24), shots = 100L,
                       epsilon = 0.005, seed = NULL) {
  method <- match.arg(method); stat <- match.arg(stat)
  alternative <- match.arg(alternative); backend <- match.arg(backend)
  if (!is.null(seed)) set.seed(seed)

  if (inherits(data, "qperm_nifti")) {
    if (is.null(coords)) coords <- data$coords
    data <- data$data
  }
  data <- .qperm_check_data(data)
  n <- nrow(data); V <- ncol(data)

  # ---- classical engine (always run; both methods and the p-value need it) ----
  if (stat == "max") {
    eng <- .qperm_maxnull(data, alternative, exact, exact_max, n_perm)
    obs_v <- .qperm_observed(data, alternative)
    observed_peak <- max(obs_v)
    p_peak <- .qperm_tail_p(eng$null, observed_peak, eng$exact)
    p_map <- .qperm_tail_p(eng$null, obs_v, eng$exact)
    clusters <- NULL
  } else {
    if (is.null(threshold))
      stop("stat = 'cluster' needs a cluster-forming `threshold`.", call. = FALSE)
    eng <- .qperm_clusternull(data, alternative, threshold, coords, exact, exact_max, n_perm)
    observed_peak <- eng$obs_size
    p_peak <- .qperm_tail_p(eng$null, observed_peak, eng$exact)
    p_map <- NULL
    obs_v <- NULL
    lab <- eng$obs_label
    if (any(lab > 0L)) {
      sizes <- tabulate(lab[lab > 0L])
      clusters <- data.frame(cluster = seq_along(sizes), size = sizes,
                             p_corrected = .qperm_tail_p(eng$null, sizes, eng$exact))
      clusters <- clusters[order(-clusters$size), ]
    } else clusters <- data.frame(cluster = integer(0), size = integer(0),
                                  p_corrected = numeric(0))
  }

  out <- list(
    method = method, stat = stat, alternative = alternative,
    n_subjects = n, n_voxels = V,
    engine_exact = eng$exact, n_eval = eng$n_eval,
    observed_peak = observed_peak, p_corrected = p_peak,
    p_map = p_map, clusters = clusters, threshold = threshold,
    null = eng$null, p_true = p_peak)

  # ---- quantum method ----
  if (method == "quantum") {
    if (backend == "qiskit") {
      if (stat != "max")
        stop("backend = 'qiskit' supports stat = 'max' only. The cluster-connectivity ",
             "oracle is an open problem; use backend = 'emulation' for cluster extent.",
             call. = FALSE)
      qk <- .qperm_qiskit(data, threshold = observed_peak, alternative = alternative,
                          powers = powers, shots = shots,
                          seed = if (is.null(seed)) 1234L else seed)
      out$q_estimate <- qk$estimate
      out$q_exact_fraction <- qk$exact_fraction
      out$q_trace <- .qperm_trace_from_counts(powers, qk$good, shots, p_peak)
    } else {
      em <- .qperm_emulate(p_peak, powers = powers, shots = shots)
      out$q_estimate <- em$estimate
      out$q_trace <- em$trace
    }
    budgets <- .qperm_budgets(max(out$q_trace$queries))
    out$c_trace <- .qperm_classical_trace(p_peak, budgets)
    out$q_backend <- backend
    out$resources <- qperm_resources(n, V, epsilon)
  }

  class(out) <- "qperm"
  out
}

# Build a convergence trace (cumulative queries vs error) from per-power success counts.
.qperm_trace_from_counts <- function(powers, goods, shots, p_true) {
  rec <- matrix(nrow = 0L, ncol = 3L, dimnames = list(NULL, c("m", "good", "shots")))
  q <- 0; qv <- numeric(length(powers)); ev <- numeric(length(powers))
  for (k in seq_along(powers)) {
    rec <- rbind(rec, c(m = powers[k], good = goods[k], shots = shots))
    q <- q + shots * (2 * powers[k] + 1)
    est <- sin(.qperm_mle_theta(rec))^2
    qv[k] <- q; ev[k] <- abs(est - p_true)
  }
  data.frame(queries = qv, error = ev, power = powers)
}

#' @export
print.qperm <- function(x, ...) {
  cat(sprintf("Sign-flip permutation test (%s, %s statistic, %s)\n",
              x$method, x$stat, x$alternative))
  cat(sprintf("  %d subjects, %d voxels; null %s over %s relabelings\n",
              x$n_subjects, x$n_voxels,
              if (x$engine_exact) "enumerated exactly" else "sampled (Monte Carlo)",
              format(x$n_eval, big.mark = ",")))
  if (x$stat == "max") {
    cat(sprintf("  observed peak statistic : %.4f\n", x$observed_peak))
    cat(sprintf("  FWER-corrected p (peak) : %.4g\n", x$p_corrected))
    if (!is.null(x$p_map))
      cat(sprintf("  voxels significant at .05: %d of %d\n",
                  sum(x$p_map <= 0.05), x$n_voxels))
  } else {
    cat(sprintf("  cluster-forming threshold: %.4f\n", x$threshold))
    cat(sprintf("  largest observed cluster : %d voxels\n", x$observed_peak))
    cat(sprintf("  FWER-corrected p (peak)  : %.4g\n", x$p_corrected))
    if (nrow(x$clusters) > 0)
      cat(sprintf("  clusters significant .05 : %d\n", sum(x$clusters$p_corrected <= 0.05)))
  }
  if (x$method == "quantum") {
    cat(sprintf("\n  quantum estimator (%s): p ~ %.4g  (exact %.4g)\n",
                x$q_backend, x$q_estimate, x$p_true))
    cat("  scaling: quantum error ~ 1/N vs classical ~ 1/sqrt(N); see plot().\n")
    cat(sprintf("  projected advantage at this size: ~%.3g x fewer evaluations (eps=%g)\n",
                x$resources$speedup_factor, x$resources$epsilon))
  }
  invisible(x)
}

#' @export
summary.qperm <- function(object, ...) {
  print(object)
  if (object$stat == "cluster" && nrow(object$clusters) > 0) {
    cat("\nObserved clusters (size, corrected p):\n")
    print(utils::head(object$clusters, 10), row.names = FALSE)
  }
  if (object$method == "quantum") {
    cat("\nMLQAE convergence (cumulative queries, error):\n")
    print(object$q_trace, row.names = FALSE)
    cat("\n"); print(object$resources)
  }
  invisible(object)
}

#' @export
plot.qperm <- function(x, ...) {
  if (x$method == "quantum" && !is.null(x$q_trace)) {
    qt <- x$q_trace; ct <- x$c_trace
     q_err <- pmax(qt$error, 1e-4); c_err <- pmax(ct$error, 1e-4)
    xr <- range(c(qt$queries, ct$queries))
    yr <- range(c(q_err, c_err))
    plot(qt$queries, q_err, log = "xy", type = "b", pch = 19, col = "#5b8fd6",
         xlim = xr, ylim = yr, xlab = "oracle / statistic evaluations",
         ylab = "estimation error",
         main = "Quantum 1/N vs classical 1/sqrt(N) scaling")
    graphics::lines(ct$queries, c_err, type = "b", pch = 19, col = "#dd6175")
    graphics::legend("bottomleft",
                     c("quantum MLQAE  ~ 1/N", "classical MC  ~ 1/sqrt(N)"),
                     col = c("#5b8fd6", "#dd6175"), pch = 19, lty = 1, bty = "n")
  } else {
    graphics::hist(x$null, breaks = 40, col = "#dbe4f0", border = "white",
                   main = "Sign-flip null of the test statistic",
                   xlab = "test statistic")
    graphics::abline(v = x$observed_peak, col = "#dd6175", lwd = 2)
    graphics::legend("topright", "observed", col = "#dd6175", lwd = 2, bty = "n")
  }
  invisible(x)
}
