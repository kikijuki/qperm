#' Synthetic one-sample contrast dataset
#'
#' A small synthetic group dataset of subject-level contrast values, 10 subjects by 6
#' voxels, used by the interactive demo in the original project. It is small enough that
#' the entire sign-flip null (2^10 = 1024 relabelings) can be enumerated exactly, which
#' lets the quantum estimate be checked against exact ground truth.
#'
#' @format A numeric matrix with 10 rows (subjects) and 6 columns (voxels).
#' @source Synthetic benchmark data from the "Quantum-accelerated permutation inference in
#'   fMRI" project.
"qae_demo"
