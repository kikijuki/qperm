#' Read subject contrast images into a qperm input
#'
#' Reads one NIfTI contrast map per subject and returns a subjects-by-voxels matrix with
#' voxel coordinates, ready to pass to \code{\link{qperm_test}}. This gives cluster-extent
#' tests real spatial geometry, and lets researchers point the test at image files from an
#' SPM or FSL first-level pipeline.
#'
#' @param files Character vector of NIfTI paths, one contrast image per subject.
#' @param mask Optional path to a NIfTI mask, or a mask array. If omitted, voxels that are
#'   finite and non-zero across all subjects are kept.
#' @return A list of class \code{"qperm_nifti"} with \code{data} (subjects x voxels),
#'   \code{coords} (voxels x 3 integer grid coordinates), \code{dims}, and
#'   \code{mask_index}. Pass it directly as the \code{data} argument of
#'   \code{\link{qperm_test}}.
#' @examples
#' \dontrun{
#' img <- qperm_read_nifti(list.files("con_maps", full.names = TRUE))
#' qperm_test(img, method = "classical", stat = "cluster", threshold = 2.3)
#' }
#' @export
qperm_read_nifti <- function(files, mask = NULL) {
  if (!requireNamespace("RNifti", quietly = TRUE))
    stop("qperm_read_nifti needs the RNifti package. install.packages('RNifti').",
         call. = FALSE)
  files <- as.character(files)
  if (length(files) < 2L)
    stop("Provide one contrast image per subject (at least 2 files).", call. = FALSE)

  imgs <- lapply(files, RNifti::readNifti)
  dims <- dim(imgs[[1L]])
  ok <- vapply(imgs, function(im) identical(dim(im), dims), logical(1))
  if (!all(ok)) stop("All images must share the same dimensions.", call. = FALSE)

  nvox_all <- prod(dims)
  Y_all <- matrix(0, nrow = length(imgs), ncol = nvox_all)
  for (i in seq_along(imgs)) Y_all[i, ] <- as.numeric(imgs[[i]])

  if (is.null(mask)) {
    keep <- colSums(!is.finite(Y_all) | Y_all == 0) == 0L
  } else {
    if (is.character(mask)) mask <- RNifti::readNifti(mask)
    if (!identical(dim(as.array(mask)), dims))
      stop("mask dimensions must match the images.", call. = FALSE)
    keep <- as.numeric(mask) != 0
  }
  keep <- which(keep)
  if (length(keep) < 1L) stop("Mask leaves no voxels.", call. = FALSE)

  out <- list(
    data = Y_all[, keep, drop = FALSE],
    coords = arrayInd(keep, .dim = dims),
    dims = dims,
    mask_index = keep)
  class(out) <- "qperm_nifti"
  out
}
