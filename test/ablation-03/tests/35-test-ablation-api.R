# The public ablation API is intentionally a single representation workflow.
source(file.path("R", "ablation.R"))

stopifnot(!"experiment" %in% names(formals(ablation)))
stopifnot(identical(
  names(formals(ablation)),
  c("object", "data", "metadata", "output.dir", "params", "seed",
    "verbose", "step", "input", "cache.root")
))

message("Ablation public API contract passed.")
