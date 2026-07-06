#' qperm: quantum-accelerated sign-flip permutation inference for neuroimaging
#'
#' qperm runs one-sample / paired sign-flip permutation tests for group neuroimaging data
#' and reports family-wise-error corrected p-values, with two methods for estimating the
#' corrected p. The classical method is a usable permutation test on your data. The quantum
#' method offers the proposed quantum amplitude-estimation route as an alternative,
#' emulated in pure R by default and optionally run as real circuits through Qiskit for
#' small studies. The quantum path is a simulation and a query-complexity demonstration,
#' not a runtime speedup on classical hardware.
#'
#' Start with \code{\link{qperm_test}}. See \code{\link{qperm_resources}} for the projected
#' query advantage at scale and \code{\link{qperm_read_nifti}} for image input.
#'
#' @keywords internal
"_PACKAGE"
