# targets helpers for ablation-03
#
# This file is intentionally limited to project wiring. Scientific calculations
# belong to the installed CCS package; no helper here sources the repository's
# R/ directory or owns a cache/checkpoint protocol.

.ablation03_cache_root <- function() {
  default <- if (identical(.Platform$OS.type, "windows")) {
    "D:/cache/ccs/_ablation-03"
  } else {
    ""
  }
  root <- Sys.getenv("CCS_ABLATION_CACHE_ROOT", unset = default)
  if (!nzchar(root)) {
    stop(
      "Set CCS_ABLATION_CACHE_ROOT before running ablation-03 on this platform.",
      call. = FALSE
    )
  }
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

.ablation03_store <- function(cache_root = .ablation03_cache_root()) {
  file.path(cache_root, "targets")
}

.ablation03_observability_config <- function(cache_root = .ablation03_cache_root()) {
  root <- file.path(cache_root, "logs", "targets-crew")
  worker_log_dir <- file.path(root, "workers")
  dir.create(worker_log_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(worker_log_dir)) {
    stop("ablation-03: cannot create crew worker log directory: ", worker_log_dir,
      call. = FALSE)
  }
  workers <- suppressWarnings(as.integer(Sys.getenv(
    "CCS_ABLATION_TARGET_WORKERS", unset = "2"
  )))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) workers <- 1L
  list(
    schema_version = 1L,
    root = normalizePath(root, winslash = "/", mustWork = FALSE),
    worker_log_dir = normalizePath(worker_log_dir, winslash = "/", mustWork = FALSE),
    main_log = file.path(root, "main-process.log"),
    metrics_interval = 1,
    workers = workers,
    backend = "crew_controller_local"
  )
}

.ablation03_crew_controller <- function(observability) {
  if (!requireNamespace("crew", quietly = TRUE)) {
    stop(
      "ablation-03 requires the crew package for targets parallel execution.",
      call. = FALSE
    )
  }
  crew::crew_controller_local(
    name = "ablation03-local",
    workers = as.integer(observability$workers),
    seconds_idle = 60,
    options_local = crew::crew_options_local(
      log_directory = observability$worker_log_dir
    ),
    options_metrics = crew::crew_options_metrics(
      path = observability$worker_log_dir,
      seconds_interval = observability$metrics_interval
    )
  )
}

.ablation03_read_resource_metrics <- function(runtime_config) {
  if (!requireNamespace("autometric", quietly = TRUE)) {
    stop("ablation-03 requires autometric to read worker resource metrics.", call. = FALSE)
  }
  observability <- runtime_config$observability
  paths <- c(
    observability$main_log,
    list.files(observability$worker_log_dir, full.names = TRUE)
  )
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  if (length(paths) == 0L) return(data.frame())
  records <- lapply(paths, function(path) {
    tryCatch(autometric::log_read(path), error = function(error) NULL)
  })
  records <- records[!vapply(records, is.null, logical(1L))]
  if (length(records) == 0L) data.frame() else do.call(rbind, records)
}

.ablation03_worker_health <- function(resource_metrics) {
  if (!is.data.frame(resource_metrics) || nrow(resource_metrics) == 0L) {
    return(data.frame())
  }
  fields <- intersect(c("name", "pid", "status", "phase"), names(resource_metrics))
  resource_metrics[, fields, drop = FALSE]
}

.ablation03_assert_writable <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("ablation-03: cannot create cache directory: ", path, call. = FALSE)
  }
  probe <- tempfile("write-probe-", tmpdir = path)
  ok <- file.create(probe)
  if (isTRUE(ok)) unlink(probe, force = TRUE)
  if (!isTRUE(ok)) {
    stop("ablation-03: cache directory is not writable: ", path, call. = FALSE)
  }
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

.ablation03_runtime_config <- function() {
  if (!requireNamespace("CCS", quietly = TRUE)) {
    stop(
      "ablation-03 requires an installed CCS package; do not use pkgload::load_all() here.",
      call. = FALSE
    )
  }
  cache_root <- .ablation03_cache_root()
  store <- .ablation03_store(cache_root)
  .ablation03_assert_writable(cache_root)
  .ablation03_assert_writable(store)
  run_id <- Sys.getenv("CCS_ABLATION_RUN_ID", unset = "manual")
  package_path <- tryCatch(find.package("CCS"), error = function(e) NA_character_)
  description <- tryCatch(utils::packageDescription("CCS"), error = function(e) NULL)
  required_version <- "0.8.3"
  package_version <- if (is.null(description)) NA_character_ else {
    as.character(description$Version)
  }
  if (is.na(package_version) ||
      utils::compareVersion(package_version, required_version) != 0) {
    stop(
      "ablation-03 requires the installed CCS package version ",
      required_version, "; found ", package_version, ".",
      call. = FALSE
    )
  }
  metadata <- list(
    schema_version = 1L,
    run_id = as.character(run_id),
    cache_root = normalizePath(cache_root, winslash = "/", mustWork = FALSE),
    store = normalizePath(store, winslash = "/", mustWork = FALSE),
    R = R.version.string,
    package = list(
      name = "CCS",
      version = package_version,
      path = package_path,
      git_commit = Sys.getenv("CCS_GIT_COMMIT", unset = NA_character_),
      build_id = Sys.getenv("CCS_BUILD_ID", unset = NA_character_)
    ),
    api = "ablation"
  )
  metadata$project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  lockfile <- file.path(metadata$project_dir, "renv.lock")
  if (!file.exists(lockfile) || !requireNamespace("renv", quietly = TRUE)) {
    stop("ablation-03 requires an initialized and activated renv project.", call. = FALSE)
  }
  metadata$renv <- list(
    project = metadata$project_dir,
    lockfile = normalizePath(lockfile, winslash = "/", mustWork = TRUE),
    lockfile_md5 = unname(tools::md5sum(lockfile)),
    activated_project = Sys.getenv("RENV_PROJECT", unset = NA_character_),
    version = as.character(utils::packageVersion("renv"))
  )
  metadata$seed <- as.integer(Sys.getenv("CCS_ABLATION_SEED", unset = "20260727"))
  metadata$input_rds <- Sys.getenv("CCS_ABLATION_INPUT_RDS", unset = NA_character_)
  metadata
}

.ablation03_read_inputs <- function(runtime_config) {
  input_path <- runtime_config$input_rds
  if (is.na(input_path) || !nzchar(input_path) || !file.exists(input_path)) {
    stop(
      "Set CCS_ABLATION_INPUT_RDS to an RDS containing object, data and metadata.",
      call. = FALSE
    )
  }
  inputs <- readRDS(input_path)
  if (!is.list(inputs)) {
    stop("ablation-03 input RDS must contain a named list.", call. = FALSE)
  }
  # The staged CCS 0.8.3 workflow uses object/data/metadata.  Preserve
  # compatibility with the existing 01-data product, which names the same
  # values resCCS_ablation/data_all/ablation_metadata.
  aliases <- list(
    object = c("object", "resCCS_ablation"),
    data = c("data", "data_all"),
    metadata = c("metadata", "ablation_metadata")
  )
  for (field in names(aliases)) {
    candidates <- aliases[[field]]
    present <- candidates[candidates %in% names(inputs)]
    present <- present[!vapply(present, function(name) is.null(inputs[[name]]), logical(1))]
    if (length(present) > 0L) inputs[[field]] <- inputs[[present[1L]]]
  }
  if (is.null(inputs$biology_inputs)) {
    biology_path <- file.path(
      runtime_config$cache_root,
      "01-biology",
      "expression-anchor-cache.rds"
    )
    if (file.exists(biology_path)) inputs$biology_inputs <- readRDS(biology_path)
  }
  if (is.null(inputs$structural_inputs)) {
    structural_path <- file.path(
      runtime_config$cache_root,
      "01-representations",
      "structural-inputs.rds"
    )
    if (file.exists(structural_path)) inputs$structural_inputs <- readRDS(structural_path)
  }
  required <- names(aliases)
  missing <- required[vapply(required, function(field) is.null(inputs[[field]]), logical(1))]
  if (length(missing) > 0L) {
    stop("ablation-03 input RDS is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!methods::is(inputs$object, "CCS")) {
    stop("ablation-03 input$object must be a CCS object.", call. = FALSE)
  }
  inputs
}

.ablation03_representation_params <- function() {
  path <- Sys.getenv("CCS_ABLATION_PARAMS_RDS", unset = "")
  if (!nzchar(path)) return(list())
  if (!file.exists(path)) stop("CCS_ABLATION_PARAMS_RDS does not exist: ", path, call. = FALSE)
  params <- readRDS(path)
  if (!is.list(params)) stop("CCS_ABLATION_PARAMS_RDS must contain a list.", call. = FALSE)
  params
}

.ablation03_require_optional_input <- function(inputs, field) {
  value <- inputs[[field]]
  if (is.null(value)) {
    stop(
      "ablation-03 input RDS must provide `", field,
      "` before the corresponding target can run.",
      call. = FALSE
    )
  }
  value
}

.ablation03_split_jobs <- function(jobs, metrics) {
  if (!is.data.frame(jobs) || nrow(jobs) == 0L) return(list())
  if (!is.data.frame(metrics) || nrow(metrics) == 0L) {
    return(lapply(seq_len(nrow(jobs)), function(i) data.frame()))
  }
  fractions <- sort(unique(metrics$requested_fraction))
  lapply(seq_len(nrow(jobs)), function(i) {
    job <- jobs[i, , drop = FALSE]
    requested_fraction <- fractions[job$fraction_index]
    keep <- metrics$requested_fraction == requested_fraction &
      metrics$repeat_id == job$repeat_id &
      metrics$representation == job$representation
    metrics[keep, , drop = FALSE]
  })
}

.ablation03_passthrough_result <- function(value, label) {
  if (is.null(value)) stop("ablation-03 missing ", label, " input.", call. = FALSE)
  list(status = "provided", label = label, value = value)
}
