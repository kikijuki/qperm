# Bridge to the real Qiskit backend (inst/python/qperm_qiskit.py) via reticulate.
# Optional: only reached when the user asks for backend = "qiskit". Kept out of Imports so
# the package installs and runs fully without Python.

.qperm_qiskit <- function(data, threshold, alternative, powers, shots, seed = 1234L) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("backend = 'qiskit' needs the reticulate package plus a Python with qiskit. ",
         "Install reticulate, then in Python: pip install qiskit. ",
         "Or use backend = 'emulation'.", call. = FALSE)
  if (nrow(data) > 16L)
    stop("The qiskit backend uses one qubit per subject and enumerates 2^n oracle ",
         "phases, so it is limited to <= 16 subjects (ideally <= 12). Use ",
         "backend = 'emulation' for larger studies.", call. = FALSE)

  pyfile <- system.file("python", "qperm_qiskit.py", package = "qperm")
  if (!nzchar(pyfile)) pyfile <- file.path("inst", "python", "qperm_qiskit.py")
  if (!file.exists(pyfile))
    stop("Could not locate qperm_qiskit.py.", call. = FALSE)

  mod <- reticulate::import_from_path("qperm_qiskit", path = dirname(pyfile))
  res <- mod$estimate_p(Y = data, threshold = threshold, alternative = alternative,
                        powers = as.integer(powers), shots = as.integer(shots),
                        seed = as.integer(seed))
  list(estimate = as.numeric(res$estimate),
       exact_fraction = as.numeric(res$exact_fraction),
       good = as.integer(unlist(res$good)),
       shots = as.integer(res$shots))
}
