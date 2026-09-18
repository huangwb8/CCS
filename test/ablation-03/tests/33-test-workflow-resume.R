# Verify force/resume routing and corrupted-RDS fallback without real analysis data.
source("test/ablation-03/templates/workflow_helpers.R", local = TRUE)

Sys.unsetenv(c("BENSZ_FORCE_STEP", "BENSZ_RESUME_FROM"))
stopifnot(identical(
  bensz_run_decision("01.00.00", list(valid = TRUE)),
  "hit"
))
stopifnot(identical(
  bensz_run_decision("01.00.00", list(valid = FALSE)),
  "run"
))

Sys.setenv(BENSZ_FORCE_STEP = "01.00.00")
stopifnot(identical(
  bensz_run_decision("01.00.00", list(valid = TRUE)),
  "run"
))
Sys.unsetenv("BENSZ_FORCE_STEP")

Sys.setenv(BENSZ_RESUME_FROM = "03.00.00")
resume_error <- tryCatch(
  bensz_run_decision("02.00.00", list(valid = FALSE)),
  error = identity
)
stopifnot(inherits(resume_error, "error"))
stopifnot(identical(
  bensz_run_decision("02.00.00", list(valid = TRUE)),
  "hit"
))
Sys.unsetenv("BENSZ_RESUME_FROM")

cache_path <- tempfile("workflow-safe-rds-", fileext = ".rds")
.wf_atomic_save_rds(list(status = "complete", value = 1L), cache_path)
validator <- function(value) identical(value$status, "complete") && identical(value$value, 1L)
stopifnot(identical(.wf_read_rds_safe(cache_path, validator), list(status = "complete", value = 1L)))
writeLines("corrupted", cache_path)
stopifnot(is.null(.wf_read_rds_safe(cache_path, validator)))
unlink(cache_path, force = TRUE)

# Formal stage products are unreadable while the stage state is running, even
# when an old result file is still present on disk.
running_stage <- tempfile("workflow-running-stage-")
dir.create(running_stage)
.wf_atomic_save_rds(list(value = "stale-formal-result"), file.path(running_stage, "result.rds"))
.wf_atomic_save_rds(
  list(status = "running", stage = "ablation-experiment"),
  file.path(running_stage, "run-state.rds")
)
original_wf_output <- .wf_output
.wf_output <- function(...) running_stage
running_error <- tryCatch(
  .wf_read("ablation-experiment", "result.rds"),
  error = identity
)
.wf_output <- original_wf_output
stopifnot(
  inherits(running_error, "error"),
  grepl("not complete", conditionMessage(running_error), fixed = TRUE)
)
unlink(running_stage, recursive = TRUE, force = TRUE)

cat("workflow resume tests passed\n")
