# Internal engines that build the sign-flip null distribution of a test statistic.
# The max-over-voxels statistic and the cluster-extent statistic each get an engine.
# Both are classical: they enumerate the sign-flip null exactly when it is small enough,
# and fall back to Monte Carlo sampling otherwise.

# Decide whether the exact (full enumeration) engine is usable.
.qperm_use_exact <- function(n, exact, exact_max) {
  exact && n <= 30L && (2^n) <= exact_max
}

# ---- max-over-voxels statistic -------------------------------------------------------

# Returns list(null, exact, n_eval). `null` is the sign-flip distribution of
# max_v (oriented per-voxel statistic).
.qperm_maxnull <- function(data, alternative, exact, exact_max, n_perm) {
  n <- nrow(data); V <- ncol(data)
  use_exact <- .qperm_use_exact(n, exact, exact_max)
  N <- if (use_exact) 2^n else n_perm
  null <- numeric(N)
  chunk <- .qperm_chunk_size(V)
  done <- 0L
  while (done < N) {
    K <- min(chunk, N - done)
    idx <- if (use_exact) (done):(done + K - 1L) else NULL
    S <- .qperm_sign_block(n, K, index = idx)
    Tmat <- (S %*% data) / n
    Tmat <- .qperm_orient(Tmat, alternative)
    null[(done + 1L):(done + K)] <- .qperm_row_max(Tmat)
    done <- done + K
  }
  list(null = null, exact = use_exact, n_eval = N)
}

# ---- cluster-extent statistic --------------------------------------------------------

# Build a face-connectivity edge list from integer voxel coordinates (V x d).
# Two voxels are neighbours if they differ by 1 on exactly one axis.
.qperm_adjacency <- function(coords) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "integer"
  V <- nrow(coords); d <- ncol(coords)
  keys <- do.call(paste, c(as.data.frame(coords), sep = ","))
  from <- integer(0); to <- integer(0)
  for (a in seq_len(d)) {
    nb <- coords; nb[, a] <- nb[, a] + 1L
    nbkeys <- do.call(paste, c(as.data.frame(nb), sep = ","))
    j <- match(nbkeys, keys)
    hit <- which(!is.na(j))
    from <- c(from, hit); to <- c(to, j[hit])
  }
  cbind(from = from, to = to)
}

# Largest connected component among the suprathreshold voxels, via union-find.
.qperm_max_component <- function(supra, edges, V) {
  if (!any(supra)) return(0L)
  parent <- seq_len(V)
  find <- function(x) {
    while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
    x
  }
  if (nrow(edges) > 0L) {
    keep <- supra[edges[, 1L]] & supra[edges[, 2L]]
    e <- edges[keep, , drop = FALSE]
    for (r in seq_len(nrow(e))) {
      ra <- find(e[r, 1L]); rb <- find(e[r, 2L])
      if (ra != rb) parent[ra] <- rb
    }
  }
  roots <- vapply(which(supra), find, integer(1))
  max(tabulate(match(roots, unique(roots))))
}

# Full component labelling of a suprathreshold map (used for reporting observed clusters).
.qperm_label_components <- function(supra, edges, V) {
  label <- integer(V)
  if (!any(supra)) return(label)
  parent <- seq_len(V)
  find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  if (nrow(edges) > 0L) {
    keep <- supra[edges[, 1L]] & supra[edges[, 2L]]
    e <- edges[keep, , drop = FALSE]
    for (r in seq_len(nrow(e))) {
      ra <- find(e[r, 1L]); rb <- find(e[r, 2L]); if (ra != rb) parent[ra] <- rb
    }
  }
  idx <- which(supra)
  roots <- vapply(idx, find, integer(1))
  label[idx] <- match(roots, unique(roots))
  label
}

# Returns list(null, exact, n_eval, obs_size, obs_label, edges).
.qperm_clusternull <- function(data, alternative, threshold, coords, exact, exact_max, n_perm) {
  n <- nrow(data); V <- ncol(data)
  if (is.null(coords))
    stop("Cluster-extent statistics need voxel `coords` (a V x d integer matrix) or ",
         "NIfTI input. A plain data matrix has no spatial geometry.", call. = FALSE)
  if (nrow(coords) != V)
    stop("`coords` must have one row per voxel (", V, " rows).", call. = FALSE)
  edges <- .qperm_adjacency(coords)

  obs_v <- .qperm_observed(data, alternative)
  obs_label <- .qperm_label_components(obs_v >= threshold, edges, V)
  obs_size <- if (any(obs_label > 0L)) max(tabulate(obs_label[obs_label > 0L])) else 0L

  use_exact <- .qperm_use_exact(n, exact, exact_max)
  N <- if (use_exact) 2^n else n_perm
  null <- numeric(N)
  chunk <- .qperm_chunk_size(V)
  done <- 0L
  while (done < N) {
    K <- min(chunk, N - done)
    idx <- if (use_exact) (done):(done + K - 1L) else NULL
    S <- .qperm_sign_block(n, K, index = idx)
    Tmat <- .qperm_orient((S %*% data) / n, alternative)
    for (r in seq_len(K)) {
      null[done + r] <- .qperm_max_component(Tmat[r, ] >= threshold, edges, V)
    }
    done <- done + K
  }
  list(null = null, exact = use_exact, n_eval = N,
       obs_size = obs_size, obs_label = obs_label, edges = edges)
}
