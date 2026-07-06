# Internal helpers: validation, sign-flip generation, tail probabilities.
# None of these are exported.

# Validate and coerce the input to a numeric subjects-by-voxels matrix.
.qperm_check_data <- function(data) {
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.matrix(data) || !is.numeric(data))
    stop("`data` must be a numeric subjects-by-voxels matrix ",
         "(rows = subjects, columns = voxels).", call. = FALSE)
  if (anyNA(data))
    stop("`data` contains NA. Remove or impute missing values first.", call. = FALSE)
  if (nrow(data) < 2L)
    stop("Need at least 2 subjects (rows).", call. = FALSE)
  if (nrow(data) > 30L)
    warning("With ", nrow(data), " subjects the exact sign-flip null has 2^",
            nrow(data), " members; the exact engine will fall back to Monte Carlo. ",
            "This is expected.", call. = FALSE)
  data
}

# Turn a per-voxel signed statistic into the one used by the chosen alternative.
# "two.sided" -> magnitude, "greater" -> as is, "less" -> negated.
.qperm_orient <- function(x, alternative) {
  switch(alternative,
    two.sided = abs(x),
    greater   = x,
    less      = -x,
    stop("Unknown alternative: ", alternative, call. = FALSE))
}

# Observed per-voxel statistic under the identity sign-flip: the column means,
# oriented for the alternative. Linear in the signs by construction, which is what
# makes the quantum oracle a constant-weighted adder.
.qperm_observed <- function(data, alternative) {
  .qperm_orient(colMeans(data), alternative)
}

# Build a block of sign-flip vectors as a K-by-n matrix of +/-1.
# When `index` is given, row r encodes integer index[r] in binary: bit i set -> subject
# i flipped to -1 (this matches the |s> basis-state convention used by the quantum code,
# where subject i is qubit i). When `index` is NULL, draw K random sign vectors.
.qperm_sign_block <- function(n, K, index = NULL) {
  if (is.null(index)) {
    matrix(sample(c(1, -1), n * K, replace = TRUE), nrow = K, ncol = n)
  } else {
    S <- matrix(1, nrow = K, ncol = n)
    for (i in seq_len(n)) {
      bit <- bitwAnd(bitwShiftR(index, i - 1L), 1L)
      S[bit == 1L, i] <- -1
    }
    S
  }
}

# A memory-safe chunk size so that a (chunk x n_voxels) statistic block stays near a
# target number of cells regardless of how many voxels there are.
.qperm_chunk_size <- function(n_voxels, target_cells = 4e6) {
  max(1L, as.integer(floor(target_cells / max(n_voxels, 1L))))
}

# Fast per-row maximum of a matrix without external packages.
.qperm_row_max <- function(M) {
  if (ncol(M) == 1L) return(as.numeric(M[, 1L]))
  M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
}

# Upper-tail probabilities P(null >= thr) for a vector of thresholds.
# `exact` uses the enumerated null directly; Monte Carlo uses the (1 + count)/(B + 1)
# convention so a p-value is never exactly zero.
.qperm_tail_p <- function(null, thresholds, exact) {
  s <- sort(null)
  N <- length(s)
  # number of null values >= thr = N - (number strictly < thr)
  ge <- N - findInterval(thresholds - .Machine$double.eps^0.5, s)
  if (exact) ge / N else (1 + ge) / (N + 1)
}
