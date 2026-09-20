# Verify force/resume routing and corrupted-RDS fallback without real analysis data.
old_cache_root <- Sys.getenv("CCS_ABLATION_CACHE_ROOT", unset = NA_character_)
old_mode <- Sys.getenv("CCS_ABLATION_MODE", unset = NA_character_)
old_test_entry <- Sys.getenv("CCS_ABLATION_ALLOW_TEST_ENTRY", unset = NA_character_)
workflow_test_root <- tempfile("ablation-workflow-test-")
Sys.setenv(
  CCS_ABLATION_CACHE_ROOT = workflow_test_root,
  CCS_ABLATION_MODE = "formal",
  CCS_ABLATION_ALLOW_TEST_ENTRY = "1"
)
source("test/ablation-03/scripts/helpers/workflow_helpers.R", local = TRUE)
stopifnot(startsWith(
  normalizePath(.wf_product_dir("01-data"), winslash = "/", mustWork = FALSE),
  paste0(normalizePath(.ablation03_dir, winslash = "/", mustWork = TRUE), "/products/main/")
))

Sys.unsetenv(c("BENSZ_FORCE_STEP", "BENSZ_RESUME_FROM"))
stopifnot(identical(
  bensz_run_decision("01.01.00", list(valid = TRUE)),
  "hit"
))
stopifnot(identical(
  bensz_run_decision("01.01.00", list(valid = FALSE)),
  "run"
))

Sys.setenv(BENSZ_FORCE_STEP = "01.01.00")
stopifnot(identical(
  bensz_run_decision("01.01.00", list(valid = TRUE)),
  "run"
))
Sys.unsetenv("BENSZ_FORCE_STEP")

Sys.setenv(BENSZ_RESUME_FROM = "01.03.00")
resume_error <- tryCatch(
  bensz_run_decision("01.02.00", list(valid = FALSE)),
  error = identity
)
stopifnot(inherits(resume_error, "error"))
stopifnot(identical(
  bensz_run_decision("01.02.00", list(valid = TRUE)),
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

failed_stage <- file.path(workflow_test_root, "ablation-experiment")
dir.create(failed_stage, recursive = TRUE)
.wf_atomic_save_rds(
  list(
    schema_version = 2L, status = "running", stage = "ablation-experiment",
    run_id = "test-run", mode = "formal", pid = Sys.getpid(),
    hostname = unname(Sys.info()[["nodename"]])
  ),
  file.path(failed_stage, "run-state.rds")
)
.wf_fail("ablation-experiment", "synthetic workflow failure")
failed_state <- readRDS(file.path(failed_stage, "run-state.rds"))
stopifnot(
  identical(failed_state$status, "failed"),
  grepl("synthetic workflow failure", failed_state$error_summary, fixed = TRUE)
)

# Lightweight runs keep even their reviewable checkpoint products under the
# explicit cache root and must not mutate the formal project products.
project_product_hashes <- function() {
  paths <- list.files(
    "test/ablation-03/products/main",
    all.files = TRUE,
    full.names = TRUE,
    no.. = TRUE,
    recursive = TRUE
  )
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  stats::setNames(unname(tools::md5sum(paths)), paths)
}
formal_products_before <- project_product_hashes()
lightweight_root <- tempfile("ablation-lightweight-products-")
dir.create(lightweight_root)
Sys.setenv(
  CCS_ABLATION_CACHE_ROOT = lightweight_root,
  CCS_ABLATION_MODE = "lightweight",
  CCS_ABLATION_ALLOW_TEST_ENTRY = "1"
)
lightweight_env <- new.env(parent = globalenv())
sys.source(
  "test/ablation-03/scripts/helpers/workflow_helpers.R",
  envir = lightweight_env
)
lightweight_product_dir <- lightweight_env$.wf_product_dir("ablation-experiment")
stopifnot(startsWith(
  normalizePath(lightweight_product_dir, winslash = "/", mustWork = FALSE),
  paste0(normalizePath(lightweight_root, winslash = "/", mustWork = TRUE), "/products/main/")
))
lightweight_env$.wf_begin("ablation-experiment", list(test = TRUE))
lightweight_output <- file.path(lightweight_root, "ablation-experiment", "synthetic.rds")
saveRDS(list(value = 1L), lightweight_output)
lightweight_env$.wf_receipt(
  "ablation-experiment",
  lightweight_env$.wf_stage_units[["ablation-experiment"]],
  outputs = lightweight_output
)
stopifnot(
  file.exists(file.path(lightweight_product_dir, "SUCCESS")),
  identical(formal_products_before, project_product_hashes())
)
unlink(lightweight_root, recursive = TRUE, force = TRUE)

unlink(workflow_test_root, recursive = TRUE, force = TRUE)
if (is.na(old_cache_root)) Sys.unsetenv("CCS_ABLATION_CACHE_ROOT") else {
  Sys.setenv(CCS_ABLATION_CACHE_ROOT = old_cache_root)
}
if (is.na(old_mode)) Sys.unsetenv("CCS_ABLATION_MODE") else {
  Sys.setenv(CCS_ABLATION_MODE = old_mode)
}
if (is.na(old_test_entry)) Sys.unsetenv("CCS_ABLATION_ALLOW_TEST_ENTRY") else {
  Sys.setenv(CCS_ABLATION_ALLOW_TEST_ENTRY = old_test_entry)
}

cat("workflow resume tests passed\n")
