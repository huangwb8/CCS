# Path resolution, prepared-cache I/O and provenance; no scientific computation.
env_path <- c("00.Environment.R", "test/ablation-03/00.Environment.R")
env_path <- env_path[file.exists(env_path)][1L]
if (is.na(env_path)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(env_path, local = TRUE)
source(.ablation03_path("02-ablation03-representation_functions.R"), local = TRUE)
.wf_root <- .ablation03_dir
.wf_path <- function(...) .ablation03_path(...)
.wf_output <- function(...) .ablation03_path("tmp", ...)
.wf_validate <- function(stage_dir) {
  .ae_validate_stage_receipt(stage_dir)
  receipt <- readRDS(file.path(stage_dir, "stage-receipt.rds"))
  upstream <- names(receipt$hashes)
  upstream <- upstream[basename(upstream) == "stage-receipt.rds"]
  for (path in upstream) .wf_validate(dirname(path))
  invisible(TRUE)
}
.wf_read <- function(stage, filename) {
  stage_dir <- .wf_output(stage)
  if (!file.exists(file.path(stage_dir, filename))) {
    stop("Missing cached product: ", file.path(stage_dir, filename),
      "; run its preparation/analysis script first (see README).", call. = FALSE)
  }
  .wf_validate(stage_dir)
  readRDS(file.path(stage_dir, filename))
}
.wf_receipt <- function(stage, script, inputs = character(), outputs) {
  .ae_write_stage_receipt(.wf_output(stage),
    inputs = c(.wf_path(script), .wf_path("00-workflow_functions.R"),
      .ablation03_repo_path("R", "ablation.R"),
      .ablation03_path("02-ablation03-representation_functions.R"), inputs),
    outputs = outputs)
}
