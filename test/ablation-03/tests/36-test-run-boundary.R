# Verify the pure-R entry boundary, explicit paths and stale-lock recovery.
source("test/ablation-03/scripts/helpers/run_identity_helpers.R", local = TRUE)

root <- tempfile("ablation-run-boundary-")
dir.create(root)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

formal <- .ablation03_bind_cache_identity(root, "formal")
stopifnot(
  identical(formal$mode, "formal"),
  identical(formal$cache_root, normalizePath(root, winslash = "/", mustWork = TRUE))
)
same <- .ablation03_bind_cache_identity(root, "formal")
stopifnot(identical(formal, same))

conflict <- tryCatch(
  .ablation03_bind_cache_identity(root, "lightweight"),
  error = identity
)
stopifnot(
  inherits(conflict, "error"),
  grepl("identity conflicts", conditionMessage(conflict), fixed = TRUE)
)

legacy_root <- tempfile("ablation-legacy-boundary-")
dir.create(legacy_root)
on.exit(unlink(legacy_root, recursive = TRUE, force = TRUE), add = TRUE)
legacy_root <- normalizePath(legacy_root, winslash = "/", mustWork = TRUE)
saveRDS(
  list(
    schema_version = 1L,
    analysis = "ablation-03",
    mode = "formal",
    cache_root = legacy_root,
    store = file.path(legacy_root, "targets")
  ),
  file.path(legacy_root, ".ablation03-root.rds")
)
migrated <- .ablation03_bind_cache_identity(legacy_root, "formal")
stopifnot(
  identical(migrated$schema_version, 2L),
  identical(migrated$runner, "numbered-r-scripts")
)

entry <- normalizePath(
  "test/ablation-03/run-ablation-03.R", winslash = "/", mustWork = TRUE
)
lock <- .ablation03_acquire_entry_lock(root, "formal", entry, "formal-test-1")
stopifnot(file.exists(file.path(lock$lock_path, "owner.rds")))
active <- tryCatch(
  .ablation03_acquire_entry_lock(root, "formal", entry, "formal-test-2"),
  error = identity
)
stopifnot(inherits(active, "error"), grepl("active", conditionMessage(active), fixed = TRUE))
.ablation03_release_entry_lock(lock)

stale_lock <- file.path(root, ".ablation-entry-lock")
dir.create(stale_lock)
saveRDS(
  list(
    run_id = "stale-run",
    mode = "formal",
    hostname = unname(Sys.info()[["nodename"]]),
    pid = 2147483647L
  ),
  file.path(stale_lock, "owner.rds")
)
recovered <- .ablation03_acquire_entry_lock(root, "formal", entry, "formal-test-3")
stopifnot(
  dir.exists(file.path(root, "stale-locks")),
  length(list.files(file.path(root, "stale-locks"), pattern = "^entry-")) == 1L
)
.ablation03_release_entry_lock(recovered)

metadata <- c(formal, list(run_id = "formal-test-4", entry_script = entry))
record <- .ablation03_write_run_record(metadata)
stopifnot(file.exists(record), identical(readRDS(record)$run_id, "formal-test-4"))

bad_id <- tryCatch(
  .ablation03_write_run_record(within(metadata, run_id <- "../escape")),
  error = identity
)
stopifnot(inherits(bad_id, "error"))

entry_lines <- readLines(entry, warn = FALSE, encoding = "UTF-8")
stage_ids <- c("01.01.00", "01.02.00", "01.03.00", "02.01.00", "02.02.00", "02.03.00")
stage_lines <- vapply(stage_ids, function(stage) {
  matches <- grep(paste0("`", stage, "` ="), entry_lines, fixed = TRUE)
  if (length(matches) != 1L) NA_integer_ else matches
}, integer(1))
stopifnot(!anyNA(stage_lines), identical(stage_lines, sort(stage_lines)))
stopifnot(any(grepl("--cache-root is required", entry_lines, fixed = TRUE)))
stopifnot(
  any(grepl("CCS_ABLATION_STAGE_ID", entry_lines, fixed = TRUE)),
  any(grepl("CCS_ABLATION_PROJECT_DIR", entry_lines, fixed = TRUE)),
  any(grepl("Sys.setlocale('LC_CTYPE', 'Chinese_China.utf8')", entry_lines, fixed = TRUE)),
  any(grepl("encoding = 'UTF-8'", entry_lines, fixed = TRUE))
)

legacy <- c(
  "test/ablation-03/scripts/run-fresh-analysis.ps1",
  "test/ablation-03/scripts/run-targets-renv.ps1",
  "test/ablation-03/scripts/helpers/run-boundary.ps1"
)
stopifnot(!any(file.exists(legacy)))

benchmark_lines <- readLines(
  "test/ablation-03/scripts/benchmark-learning-curve.R",
  warn = FALSE,
  encoding = "UTF-8"
)
stopifnot(
  any(grepl("--cache-root and --output-dir are required", benchmark_lines, fixed = TRUE)),
  any(grepl("--output-dir must be outside --cache-root", benchmark_lines, fixed = TRUE))
)

cat("run boundary tests passed\n")
