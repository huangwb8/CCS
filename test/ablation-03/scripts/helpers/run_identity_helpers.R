# Explicit cache identity, entry locking and run records for ablation-03.

.ablation03_normalize_path <- function(path, must_work = FALSE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

.ablation03_validate_mode <- function(mode) {
  allowed <- c("formal", "lightweight", "benchmark")
  if (length(mode) != 1L || is.na(mode) || !mode %in% allowed) {
    stop(
      "ablation-03: mode must be one of ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  mode
}

.ablation03_validate_run_id <- function(run_id) {
  if (length(run_id) != 1L || is.na(run_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", run_id)) {
    stop("ablation-03: run_id contains unsupported path characters.", call. = FALSE)
  }
  run_id
}

.ablation03_pid_is_alive <- function(pid) {
  if (length(pid) != 1L || is.na(pid) || !is.finite(pid) || pid < 1L) return(FALSE)
  pid <- as.integer(pid)
  if (identical(.Platform$OS.type, "windows")) {
    output <- tryCatch(
      suppressWarnings(system2(
        "tasklist", c("/FI", shQuote(paste("PID eq", pid)), "/NH"),
        stdout = TRUE, stderr = FALSE
      )),
      error = function(error) character()
    )
    return(any(grepl(paste0("\\b", pid, "\\b"), output)))
  }
  status <- tryCatch(
    system2("kill", c("-0", as.character(pid)), stdout = FALSE, stderr = FALSE),
    error = function(error) 1L
  )
  identical(status, 0L)
}

.ablation03_atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(value, temporary, version = 3)
  if (!file.rename(temporary, path)) {
    stop("ablation-03: cannot commit RDS file: ", path, call. = FALSE)
  }
  invisible(path)
}

.ablation03_bind_cache_identity <- function(cache_root, mode) {
  mode <- .ablation03_validate_mode(mode)
  cache_root <- .ablation03_normalize_path(cache_root)
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cache_root)) {
    stop("ablation-03: cannot create cache root: ", cache_root, call. = FALSE)
  }
  expected <- list(
    schema_version = 2L,
    analysis = "ablation-03",
    mode = mode,
    cache_root = cache_root,
    runner = "numbered-r-scripts"
  )
  path <- file.path(cache_root, ".ablation03-root.rds")
  if (file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = function(error) NULL)
    comparable <- c("schema_version", "analysis", "mode", "cache_root", "runner")
    legacy <- is.list(existing) && identical(existing$schema_version, 1L) &&
      identical(existing$analysis, expected$analysis) &&
      identical(existing$mode, expected$mode) &&
      identical(existing$cache_root, expected$cache_root)
    if (legacy) {
      .ablation03_atomic_save_rds(expected, path)
      return(invisible(expected))
    }
    if (!is.list(existing) || !identical(existing[comparable], expected[comparable])) {
      stop(
        "ablation-03: cache root identity conflicts with mode/store: ", cache_root,
        call. = FALSE
      )
    }
    return(invisible(existing))
  }
  .ablation03_atomic_save_rds(expected, path)
  invisible(expected)
}

.ablation03_entry_lock_path <- function(cache_root) {
  file.path(cache_root, ".ablation-entry-lock")
}

.ablation03_archive_stale_lock <- function(lock_path, cache_root) {
  stale_parent <- file.path(cache_root, "stale-locks")
  dir.create(stale_parent, recursive = TRUE, showWarnings = FALSE)
  stale_path <- file.path(
    stale_parent,
    paste0("entry-", format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid())
  )
  if (!file.rename(lock_path, stale_path)) {
    stop("ablation-03: cannot archive stale entry lock: ", lock_path, call. = FALSE)
  }
  invisible(stale_path)
}

.ablation03_acquire_entry_lock <- function(cache_root, mode, entry_script, run_id = NULL) {
  mode <- .ablation03_validate_mode(mode)
  cache_root <- .ablation03_normalize_path(cache_root)
  entry_script <- .ablation03_normalize_path(entry_script, must_work = TRUE)
  .ablation03_bind_cache_identity(cache_root, mode)
  if (is.null(run_id) || !nzchar(run_id)) {
    run_id <- paste0(mode, "-", format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid())
  }
  run_id <- .ablation03_validate_run_id(run_id)
  lock_path <- .ablation03_entry_lock_path(cache_root)
  owner_path <- file.path(lock_path, "owner.rds")
  if (dir.exists(lock_path)) {
    owner <- tryCatch(readRDS(owner_path), error = function(error) NULL)
    local_host <- unname(Sys.info()[["nodename"]])
    active <- is.list(owner) &&
      (!identical(owner$hostname, local_host) || .ablation03_pid_is_alive(owner$pid))
    if (active) {
      stop(
        "ablation-03: cache root is active in run ", owner$run_id, ": ", cache_root,
        call. = FALSE
      )
    }
    .ablation03_archive_stale_lock(lock_path, cache_root)
  }
  if (!dir.create(lock_path, showWarnings = FALSE)) {
    stop("ablation-03: cache root was locked concurrently: ", cache_root, call. = FALSE)
  }
  owner <- list(
    schema_version = 1L,
    run_id = run_id,
    mode = mode,
    cache_root = cache_root,
    entry_script = entry_script,
    hostname = unname(Sys.info()[["nodename"]]),
    pid = Sys.getpid(),
    started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  committed <- FALSE
  on.exit({
    if (!committed && dir.exists(lock_path)) unlink(lock_path, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  .ablation03_atomic_save_rds(owner, owner_path)
  committed <- TRUE
  owner$lock_path <- lock_path
  owner
}

.ablation03_release_entry_lock <- function(lock) {
  if (!is.list(lock) || is.null(lock$lock_path) || is.null(lock$run_id)) {
    stop("ablation-03: invalid entry lock object.", call. = FALSE)
  }
  owner_path <- file.path(lock$lock_path, "owner.rds")
  owner <- tryCatch(readRDS(owner_path), error = function(error) NULL)
  if (is.list(owner) && !identical(owner$run_id, lock$run_id)) {
    stop("ablation-03: refusing to release an entry lock owned by another run.", call. = FALSE)
  }
  if (dir.exists(lock$lock_path) &&
      unlink(lock$lock_path, recursive = TRUE, force = TRUE) != 0L) {
    stop("ablation-03: cannot release entry lock: ", lock$lock_path, call. = FALSE)
  }
  invisible(TRUE)
}

.ablation03_write_run_record <- function(metadata) {
  if (!is.list(metadata) || is.null(metadata$cache_root) || is.null(metadata$run_id)) {
    stop("ablation-03: invalid run metadata.", call. = FALSE)
  }
  metadata$run_id <- .ablation03_validate_run_id(metadata$run_id)
  run_dir <- file.path(metadata$cache_root, "runs", metadata$run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  .ablation03_atomic_save_rds(metadata, file.path(run_dir, "runtime.rds"))
}
