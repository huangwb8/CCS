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

.ablation03_render_report <- function(
    report_file,
    dependency = NULL,
    cache_root = .ablation03_cache_root()) {
  if (!is.null(dependency)) invisible(dependency)
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("ablation-03 requires rmarkdown to render HTML reports.", call. = FALSE)
  }
  report_path <- file.path(getwd(), report_file)
  if (!file.exists(report_path)) {
    stop("ablation-03 report is missing: ", report_path, call. = FALSE)
  }
  output_file <- sub("\\.Rmd$", ".html", basename(report_path), ignore.case = TRUE)
  preview_root <- file.path(cache_root, "logs", "report-previews")
  dir.create(preview_root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(preview_root)) {
    stop(
      "ablation-03: cannot create report preview directory: ",
      preview_root,
      call. = FALSE
    )
  }
  old_cache_root <- Sys.getenv("CCS_ABLATION_CACHE_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(old_cache_root)) {
      Sys.unsetenv("CCS_ABLATION_CACHE_ROOT")
    } else {
      Sys.setenv(CCS_ABLATION_CACHE_ROOT = old_cache_root)
    }
  }, add = TRUE)
  Sys.setenv(CCS_ABLATION_CACHE_ROOT = cache_root)
  rendered <- rmarkdown::render(
    input = normalizePath(report_path, winslash = "/", mustWork = TRUE),
    output_file = output_file,
    output_dir = normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    knit_root_dir = normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    params = list(plot_run_dir = preview_root),
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )
  if (!file.exists(rendered)) {
    stop("ablation-03 report produced no HTML: ", rendered, call. = FALSE)
  }
  normalizePath(rendered, winslash = "/", mustWork = TRUE)
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
  run_id <- Sys.getenv("CCS_ABLATION_RUN_ID", unset = "targets")
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
  if (is.na(metadata$input_rds) || !nzchar(metadata$input_rds)) {
    stop("Set CCS_ABLATION_INPUT_RDS to the formal or test input bundle.", call. = FALSE)
  }
  metadata$input_rds <- normalizePath(metadata$input_rds, winslash = "/", mustWork = TRUE)
  metadata$analysis_contract <- "same-targets-same-parameters-input-only-differs"
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
  # Normalize the historical stage names once at the input boundary. Every
  # downstream target consumes this same contract for formal and test data.
  inputs$resCCS_ablation <- inputs$object
  inputs$resCCS_full <- if (is.null(inputs$resCCS_full)) inputs$object else inputs$resCCS_full
  inputs$data_all <- inputs$data
  inputs$ablation_metadata <- inputs$metadata
  inputs$n_cores <- suppressWarnings(as.integer(Sys.getenv("CCS_ABLATION_CORES", unset = "1")))
  if (is.na(inputs$n_cores) || inputs$n_cores < 1L) inputs$n_cores <- 1L
  if (is.null(inputs$filtered_cohorts)) {
    inputs$filtered_cohorts <- as.character(inputs$object@Data$filtered.cohort)
  }
  inputs$full_d1 <- inputs$resCCS_full@Data$Probability$d1
  if (is.null(inputs$tissue_resolution_audit)) {
    inputs$tissue_resolution_audit <- data.frame()
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

# targets is the only workflow owner.  Numbered R files remain scientific
# implementations, but are evaluated as target commands in the current R
# process; they do not start child R sessions or own cache/lock/receipt state.
.ablation03_target_stage <- function(stage_file, dependency = NULL, cache_root = .ablation03_cache_root()) {
  if (!is.null(dependency)) invisible(dependency)
  stage_path <- file.path(getwd(), stage_file)
  if (!file.exists(stage_path)) {
    stop("ablation-03 target stage is missing: ", stage_path, call. = FALSE)
  }
  old <- Sys.getenv(c(
    "CCS_ABLATION_TARGETS", "CCS_ABLATION_CACHE_ROOT", "CCS_ABLATION_MODE",
    "CCS_ABLATION_RUN_ID", "CCS_ABLATION_ALLOW_TEST_ENTRY"
  ), unset = NA_character_)
  on.exit({
    names(old) <- c(
      "CCS_ABLATION_TARGETS", "CCS_ABLATION_CACHE_ROOT", "CCS_ABLATION_MODE",
      "CCS_ABLATION_RUN_ID", "CCS_ABLATION_ALLOW_TEST_ENTRY"
    )
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else {
        do.call(Sys.setenv, stats::setNames(list(old[[name]]), name))
      }
    }
  }, add = TRUE)
  Sys.setenv(
    CCS_ABLATION_TARGETS = "1",
    CCS_ABLATION_CACHE_ROOT = cache_root,
    CCS_ABLATION_MODE = "formal",
    CCS_ABLATION_RUN_ID = "targets",
    CCS_ABLATION_ALLOW_TEST_ENTRY = "1"
  )
  stage_env <- new.env(parent = globalenv())
  # `sys.source()` in R 4.3.1 has no `encoding` argument. Use `source()`
  # with an explicit UTF-8 encoding while preserving the isolated stage
  # environment used by the targets command.
  source(stage_path, local = stage_env, encoding = "UTF-8")
  .ablation03_stage_artifacts(stage_file, cache_root)
}

.ablation03_stage_artifacts <- function(stage_file, cache_root) {
  stage_id <- sub("^([0-9]{2}\\.[0-9]{2}\\.[0-9]{2}).*$", "\\1", basename(stage_file))
  directory <- switch(
    stage_id,
    `01.01.00` = "01-data",
    `01.02.00` = "01-representations",
    `01.03.00` = "01-biology",
    `02.01.00` = "ablation-experiment",
    `02.02.00` = "ablation-biology",
    `02.03.00` = "ablation-structural-reproducibility",
    stop("Unknown ablation-03 target stage: ", stage_file, call. = FALSE)
  )
  directory <- file.path(cache_root, directory)
  if (!dir.exists(directory)) stop("Target stage produced no output directory: ", directory, call. = FALSE)
  list(
    stage = stage_id,
    cache_root = normalizePath(cache_root, winslash = "/", mustWork = FALSE),
    directory = normalizePath(directory, winslash = "/", mustWork = TRUE),
    files = sort(list.files(directory, full.names = TRUE, recursive = FALSE))
  )
}

.ablation03_prepare_data_target <- function(runtime_config) {
  output_path <- file.path(runtime_config$cache_root, '01-data', 'inputs.rds')
  if (identical(tolower(normalizePath(runtime_config$input_rds, winslash = '/')),
      tolower(normalizePath(output_path, winslash = '/', mustWork = FALSE)))) {
    stop('ablation-03 input must not be its generated 01-data/inputs.rds output.', call. = FALSE)
  }
  inputs <- .ablation03_read_inputs(runtime_config)
  root <- file.path(runtime_config$cache_root, "01-data")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(inputs, file.path(root, "inputs.rds"), version = 3)
  profile <- inputs$data_profile
  if (is.null(profile)) {
    profile_helpers <- new.env(parent = globalenv())
    source(file.path(getwd(), '01.01.00. 数据准备_functions.R'),
      local = profile_helpers, encoding = 'UTF-8')
    resolved_index <- profile_helpers$.atd_cohort_index(inputs$data)
    original_index <- resolved_index
    cohort_tissues <- unique(inputs$metadata[, c('cohort', 'bank_tissue')])
    if (anyDuplicated(cohort_tissues$cohort)) {
      stop('ablation-03: cohort-to-bank-tissue mapping is ambiguous.', call. = FALSE)
    }
    original_index$tissue <- cohort_tissues$bank_tissue[
      match(original_index$cohort, cohort_tissues$cohort)
    ]
    if (anyNA(original_index$tissue)) {
      stop('ablation-03: bank tissue is missing from input metadata.', call. = FALSE)
    }
    original_index$cohort_key <- paste(
      original_index$tissue, original_index$cohort, sep = '/'
    )
    query_keys <- unique(inputs$metadata[
      inputs$metadata$analysis_set == 'external_query',
      c('bank_cohort_key', 'cohort_key'), drop = FALSE
    ])
    profile <- profile_helpers$.atd_build_data_profile(
      raw_cohort_index = original_index,
      resolved_cohort_index = resolved_index,
      metadata = inputs$metadata,
      tissue_resolution_audit = inputs$tissue_resolution_audit,
      filtered_model_cohorts = query_keys$bank_cohort_key,
      filtered_cohorts = query_keys$cohort_key
    )
  }
  saveRDS(profile, file.path(root, "data-profile.rds"), version = 3)
  list(input = inputs, artifacts = list(
    stage = "01.01.00",
    cache_root = normalizePath(runtime_config$cache_root, winslash = "/", mustWork = FALSE),
    directory = normalizePath(root, winslash = "/", mustWork = TRUE),
    files = sort(list.files(root, full.names = TRUE, recursive = FALSE))
  ))
}

.ablation03_materialize_optional_inputs <- function(data_target, runtime_config) {
  inputs <- data_target$input
  if (!is.null(inputs$biology_inputs)) {
    dir.create(file.path(runtime_config$cache_root, "01-biology"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(inputs$biology_inputs,
      file.path(runtime_config$cache_root, "01-biology", "expression-anchor-cache.rds"), version = 3)
    if (!is.null(inputs$biology_inputs$structural_anchor_cache)) {
      saveRDS(inputs$biology_inputs$structural_anchor_cache,
        file.path(runtime_config$cache_root, "01-biology", "structural-anchor-cache.rds"), version = 3)
    }
  }
  if (!is.null(inputs$structural_inputs)) {
    dir.create(file.path(runtime_config$cache_root, "01-representations"), recursive = TRUE, showWarnings = FALSE)
    saveRDS(inputs$structural_inputs,
      file.path(runtime_config$cache_root, "01-representations", "structural-inputs.rds"), version = 3)
    if (!is.null(inputs$structural_inputs$structural_anchor_cache)) {
      dir.create(file.path(runtime_config$cache_root, "01-biology"), recursive = TRUE, showWarnings = FALSE)
      saveRDS(inputs$structural_inputs$structural_anchor_cache,
        file.path(runtime_config$cache_root, "01-biology", "structural-anchor-cache.rds"), version = 3)
    }
  }
  invisible(data_target)
}

.ablation03_biology_target <- function(data_target, representation_target, runtime_config) {
  invisible(data_target)
  inputs <- data_target$input
  if (!is.null(inputs$biology_inputs)) {
    .ablation03_materialize_optional_inputs(data_target, runtime_config)
    if (!file.exists(file.path(runtime_config$cache_root, "01-biology", "structural-anchor-cache.rds"))) {
      stop(
        "biology_inputs must include structural_anchor_cache when targets bypasses external expression preparation.",
        call. = FALSE
      )
    }
    return(.ablation03_stage_artifacts("01.03.00. biology-inputs.R", runtime_config$cache_root))
  }
  .ablation03_target_stage(
    "01.03.00. 生物输入准备.R",
    dependency = representation_target,
    cache_root = runtime_config$cache_root
  )
}

.ablation03_read_artifact <- function(artifacts, filename, required = TRUE) {
  path <- file.path(artifacts$directory, filename)
  if (!file.exists(path)) {
    if (required) stop("Missing target artifact: ", path, call. = FALSE)
    return(NULL)
  }
  readRDS(path)
}
