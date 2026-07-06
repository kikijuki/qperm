#' Quantum versus classical resource estimate for a permutation test
#'
#' Reports order-of-magnitude query counts and a qubit count for the compound
#' quantum amplitude estimation algorithm at a given problem size, next to the classical
#' Monte Carlo cost. This is a forward-looking projection of the algorithm's
#' query-complexity, not a measurement of runtime on any real machine.
#'
#' The classical cost to resolve an amplitude to additive precision \code{epsilon} scales
#' as \eqn{O(\epsilon^{-2})}, and evaluating the max-over-voxels event touches all
#' \code{n_voxels} voxels, giving \eqn{O(\epsilon^{-2} V)}. The compound quantum
#' algorithm nests a Grover existence search over voxels, \eqn{O(\sqrt{V})}, inside
#' amplitude estimation over permutations, \eqn{O(\epsilon^{-1})}, giving
#' \eqn{O(\epsilon^{-1}\sqrt{V})}. Constants are approximate and the coherent nesting is
#' an open problem (see the vignette), so treat these as scaling, not timings.
#'
#' @param n_subjects Number of subjects (size of the sign-flip register).
#' @param n_voxels Number of voxels tested.
#' @param epsilon Target additive precision on the corrected p-value. Default 0.005.
#' @return An object of class \code{"qperm_resources"} with the query counts, their ratio,
#'   and an approximate qubit count.
#' @examples
#' qperm_resources(n_subjects = 20, n_voxels = 50000, epsilon = 0.005)
#' @export
qperm_resources <- function(n_subjects, n_voxels, epsilon = 0.005) {
  if (epsilon <= 0 || epsilon >= 1) stop("`epsilon` must be in (0, 1).", call. = FALSE)
  if (n_voxels < 1 || n_subjects < 1) stop("Sizes must be positive.", call. = FALSE)

  classical_queries <- n_voxels / epsilon^2
  quantum_queries   <- (pi / 4) * sqrt(n_voxels) / epsilon

  qubits <- n_subjects +                       # one qubit per subject (sign-flip register)
    ceiling(log2(max(n_voxels, 2))) +          # voxel index register for the inner Grover
    ceiling(log2(n_subjects + 1)) +            # accumulator for the constant-weighted adder
    2L                                         # comparator ancilla + label qubit

  out <- list(
    n_subjects = n_subjects, n_voxels = n_voxels, epsilon = epsilon,
    classical_queries = classical_queries,
    quantum_queries = quantum_queries,
    speedup_factor = classical_queries / quantum_queries,
    qubits = qubits)
  class(out) <- "qperm_resources"
  out
}

#' @export
print.qperm_resources <- function(x, ...) {
  cat("Quantum resource estimate (query complexity, constants approximate)\n")
  cat(sprintf("  problem size    : %d subjects, %s voxels, target precision eps = %g\n",
              x$n_subjects, format(x$n_voxels, big.mark = ","), x$epsilon))
  cat(sprintf("  classical cost  : ~%.3g oracle evaluations   O(eps^-2 * V)\n",
              x$classical_queries))
  cat(sprintf("  quantum cost    : ~%.3g oracle evaluations   O(eps^-1 * sqrt(V))\n",
              x$quantum_queries))
  cat(sprintf("  projected ratio : ~%.3g x fewer quantum evaluations\n", x$speedup_factor))
  cat(sprintf("  qubits (approx) : %d\n", x$qubits))
  cat("  note: a query-complexity projection, not a wall-clock speedup on today's hardware.\n")
  invisible(x)
}
