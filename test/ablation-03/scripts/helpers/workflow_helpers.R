# Path resolution, external-cache I/O and provenance; no scientific computation.
env_path <- c("00.Environment.R", "test/ablation-03/00.Environment.R")
env_path <- env_path[file.exists(env_path)][1L]
if (is.na(env_path)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(env_path, local = TRUE)
source(.ablation03_path("scripts", "helpers", "checkpoint_helpers.R"), local = TRUE)
source(.ablation03_path("scripts", "helpers", "stage_receipt_helpers.R"), local = TRUE)
source(.ablation03_path("scripts", "helpers", "run_identity_helpers.R"), local = TRUE)
.ablation03_ccs_version <- as.character(utils::packageVersion("CCS"))
if (!identical(.ablation03_ccs_version, "0.8.3")) {
  stop(
    "ablation-03 requires installed CCS 0.8.3; found ", .ablation03_ccs_version, ".",
    call. = FALSE
  )
}
.ablation03_ccs_namespace <- asNamespace("CCS")
.ablation03_ccs_symbols <- grep(
  "^\\.ablation_", ls(.ablation03_ccs_namespace, all.names = TRUE), value = TRUE
)
for (.ablation03_symbol in .ablation03_ccs_symbols) {
  assign(
    .ablation03_symbol,
    get(.ablation03_symbol, envir = .ablation03_ccs_namespace, inherits = FALSE),
    envir = environment()
  )
}
rm(.ablation03_symbol)
.ablation03_ccs_description <- file.path(find.package("CCS"), "DESCRIPTION")
.wf_root <- .ablation03_dir
.wf_path <- function(...) .ablation03_path(...)
.wf_output <- function(...) file.path(.ablation03_cache_root, ...)
.wf_lock_dir <- file.path(.ablation03_cache_root, ".workflow-lock")
.wf_mode <- Sys.getenv("CCS_ABLATION_MODE", unset = "formal")
if (!.wf_mode %in% c("formal", "lightweight")) {
  stop("CCS_ABLATION_MODE must be formal or lightweight.", call. = FALSE)
}
.wf_product_root <- if (identical(.wf_mode, "formal")) {
  .ablation03_dir
} else {
  .ablation03_cache_root
}
.ablation03_bind_cache_identity(.ablation03_cache_root, .wf_mode)
.wf_entry_lock <- file.path(.ablation03_cache_root, ".ablation-entry-lock", "owner.rds")
if (!identical(Sys.getenv("CCS_ABLATION_ALLOW_TEST_ENTRY"), "1")) {
  .wf_entry_owner <- tryCatch(readRDS(.wf_entry_lock), error = function(error) NULL)
  .wf_run_id <- Sys.getenv("CCS_ABLATION_RUN_ID", unset = "")
  if (!is.list(.wf_entry_owner) || !nzchar(.wf_run_id) ||
      !identical(.wf_entry_owner$run_id, .wf_run_id) ||
      !identical(.wf_entry_owner$mode, .wf_mode)) {
    stop(
      "ablation-03 numbered stages must run through run-ablation-03.R.",
      call. = FALSE
    )
  }
}
.wf_stage_units <- c(
  `01-data` = "01.01.00. 数据准备",
  `01-representations` = "01.02.00. 表示输入准备",
  `01-biology` = "01.03.00. 生物输入准备",
  `ablation-experiment` = "02.01.00. 表示分析",
  `ablation-biology` = "02.02.00. 生物锚点分析",
  `ablation-structural-reproducibility` = "02.03.00. 结构复现分析"
)
.wf_product_dir <- function(stage) {
  unit_stem <- unname(.wf_stage_units[[stage]])
  if (is.null(unit_stem)) stop("Unknown workflow stage: ", stage, call. = FALSE)
  bensz_product_dir("main", unit_stem, root = .wf_product_root)
}
.wf_product_status <- function(stage) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(.wf_product_root)
  bensz_checkpoint_status(.wf_product_dir(stage))
}
.wf_atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(value, temporary, version = 3)
  invisible(readRDS(temporary))
  backup <- paste0(path, ".bak-", Sys.getpid())
  if (file.exists(path) && !file.rename(path, backup)) {
    stop("Cannot stage workflow state for replacement: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("Cannot commit workflow state: ", path, call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  invisible(path)
}
.wf_identity <- function(parameters = list()) {
  digest::digest(parameters, algo = "md5")
}
.wf_unit_id <- function(stage) {
  unit_stem <- unname(.wf_stage_units[[stage]])
  if (is.null(unit_stem)) stop("Unknown workflow stage: ", stage, call. = FALSE)
  sub("^([0-9]{2}\\.[0-9]{2}\\.[0-9]{2}).*$", "\\1", unit_stem)
}
.wf_read_rds_safe <- function(path, validator = NULL) {
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(readRDS(path), error = function(error) NULL)
  if (is.null(value)) return(NULL)
  if (!is.null(validator) && !isTRUE(validator(value))) return(NULL)
  value
}
.wf_cache_hit <- function(stage, parameters = list()) {
  stage_dir <- .wf_output(stage)
  state_path <- file.path(stage_dir, "run-state.rds")
  expected_identity <- .wf_identity(parameters)
  valid <- tryCatch({
    .wf_validate(stage_dir)
    state <- readRDS(state_path)
    if (!identical(state$cache_identity, expected_identity)) {
      stop("stage parameter identity changed", call. = FALSE)
    }
    product_status <- .wf_product_status(stage)
    if (!isTRUE(product_status$valid)) stop(product_status$reason, call. = FALSE)
    TRUE
  }, error = function(error) FALSE)
  decision <- bensz_run_decision(
    .wf_unit_id(stage),
    list(valid = isTRUE(valid), reason = if (isTRUE(valid)) "cache hit" else "invalid cache")
  )
  if (identical(decision, "hit")) {
    message("[CACHE HIT] ", unname(.wf_stage_units[[stage]]))
    return(TRUE)
  }
  FALSE
}
.wf_begin <- function(stage, parameters = list()) {
  stage_dir <- .wf_output(stage)
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(.ablation03_cache_root, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(.wf_lock_dir)) {
    owner <- .wf_read_rds_safe(file.path(.wf_lock_dir, "owner.rds"))
    local_host <- unname(Sys.info()[["nodename"]])
    active <- is.list(owner) &&
      (!identical(owner$hostname, local_host) || .ablation03_pid_is_alive(owner$pid))
    if (active) {
      stop(
        "Another ablation-03 stage owns the external-cache lock: ", .wf_lock_dir,
        " (run ", owner$run_id, ").",
        call. = FALSE
      )
    }
    stale_dir <- file.path(
      .ablation03_cache_root, "stale-locks",
      paste0("stage-", format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid())
    )
    dir.create(dirname(stale_dir), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(.wf_lock_dir, stale_dir)) {
      stop("Cannot archive stale workflow lock: ", .wf_lock_dir, call. = FALSE)
    }
  }
  if (!dir.create(.wf_lock_dir, showWarnings = FALSE)) {
    stop(
      "Another ablation-03 stage owns the external-cache lock: ", .wf_lock_dir,
      ". Do not run stages concurrently. If a prior R process crashed, confirm it ",
      "has stopped and remove this lock directory manually.",
      call. = FALSE
    )
  }
  lock_ready <- FALSE
  on.exit({
    if (!lock_ready && dir.exists(.wf_lock_dir)) {
      unlink(.wf_lock_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  .wf_atomic_save_rds(
    list(
      schema_version = 2L,
      stage = stage,
      run_id = Sys.getenv("CCS_ABLATION_RUN_ID", unset = "manual"),
      mode = .wf_mode,
      pid = Sys.getpid(),
      hostname = unname(Sys.info()[["nodename"]]),
      started = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    file.path(.wf_lock_dir, "owner.rds")
  )
  receipt <- file.path(stage_dir, "stage-receipt.rds")
  if (file.exists(receipt) && unlink(receipt) != 0L) {
    stop("Cannot invalidate old stage receipt: ", receipt, call. = FALSE)
  }
  success <- file.path(.wf_product_dir(stage), "SUCCESS")
  if (file.exists(success) && unlink(success) != 0L) {
    stop("Cannot invalidate old product marker: ", success, call. = FALSE)
  }
  .wf_atomic_save_rds(
    list(
      schema_version = 2L,
      status = "running",
      stage = stage,
      cache_identity = .wf_identity(parameters),
      run_id = Sys.getenv("CCS_ABLATION_RUN_ID", unset = "manual"),
      mode = .wf_mode,
      pid = Sys.getpid(),
      hostname = unname(Sys.info()[["nodename"]]),
      started = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    file.path(stage_dir, "run-state.rds")
  )
  options(
    bensz.ablation03.previous_error = getOption("error"),
    bensz.ablation03.active_stage = stage,
    error = function() {
      previous <- getOption("bensz.ablation03.previous_error")
      options(error = previous)
      active <- getOption("bensz.ablation03.active_stage")
      if (is.character(active) && length(active) == 1L && nzchar(active)) {
        .wf_fail(active, geterrmessage())
      }
      if (is.function(previous)) previous()
    }
  )
  lock_ready <- TRUE
  invisible(stage_dir)
}
.wf_fail <- function(stage, message) {
  state_path <- file.path(.wf_output(stage), "run-state.rds")
  state <- .wf_read_rds_safe(state_path)
  if (!is.list(state)) state <- list(schema_version = 2L, stage = stage)
  state$status <- "failed"
  state$updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  state$failed <- state$updated
  state$error_summary <- substr(gsub("[\r\n\t]+", " ", as.character(message)), 1L, 1000L)
  .wf_atomic_save_rds(state, state_path)
  if (dir.exists(.wf_lock_dir)) unlink(.wf_lock_dir, recursive = TRUE, force = TRUE)
  invisible(state)
}
.wf_validate <- function(stage_dir) {
  state_path <- file.path(stage_dir, "run-state.rds")
  if (!file.exists(state_path)) stop("Missing run state: ", stage_dir, call. = FALSE)
  state <- tryCatch(readRDS(state_path), error = function(error) NULL)
  if (!is.list(state) || !identical(state$status, "complete")) {
    stop("Workflow stage is not complete: ", stage_dir, call. = FALSE)
  }
  .ae_validate_stage_receipt(stage_dir)
  receipt_path <- file.path(stage_dir, "stage-receipt.rds")
  receipt_md5 <- unname(tools::md5sum(receipt_path))
  if (!is.character(state$receipt_md5) || length(state$receipt_md5) != 1L ||
      !identical(state$receipt_md5, receipt_md5)) {
    stop("Workflow state/receipt identity mismatch: ", stage_dir, call. = FALSE)
  }
  stage <- state$stage
  if (!is.character(stage) || length(stage) != 1L || is.null(.wf_stage_units[[stage]])) {
    stop("Workflow state has an unknown stage: ", stage_dir, call. = FALSE)
  }
  product_status <- .wf_product_status(stage)
  if (!isTRUE(product_status$valid)) {
    stop("Invalid project checkpoint for ", stage, ": ", product_status$reason, call. = FALSE)
  }
  product <- readRDS(file.path(.wf_product_dir(stage), "main.rds"))
  if (!is.list(product) || !identical(product$stage, stage) ||
      !identical(product$receipt_md5, receipt_md5)) {
    stop("Project checkpoint/stage receipt mismatch: ", stage, call. = FALSE)
  }
  receipt <- readRDS(receipt_path)
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
.wf_receipt <- function(stage, unit_stem, inputs = character(), outputs) {
  expected_stem <- unname(.wf_stage_units[[stage]])
  if (!identical(unit_stem, expected_stem)) {
    stop("Workflow stage/unit mismatch: ", stage, " / ", unit_stem, call. = FALSE)
  }
  script <- paste0(unit_stem, ".R")
  stage_dir <- .wf_output(stage)
  .ae_write_stage_receipt(.wf_output(stage),
    inputs = c(
      .wf_path(script),
      .wf_path("scripts", "helpers", "workflow_helpers.R"),
      .wf_path("scripts", "helpers", "checkpoint_helpers.R"),
      .wf_path("scripts", "helpers", "stage_receipt_helpers.R"),
      inputs
    ),
    outputs = outputs)
  receipt_path <- file.path(stage_dir, "stage-receipt.rds")
  receipt <- readRDS(receipt_path)
  running_state <- readRDS(file.path(stage_dir, "run-state.rds"))

  # Keep only a lightweight, reviewable product in the repository. Large RDS
  # files remain under CCS_ABLATION_CACHE_ROOT.
  relative_outputs <- substring(
    normalizePath(outputs, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(stage_dir, winslash = "/", mustWork = TRUE)) + 2L
  )
  product <- list(
    stage = stage,
    cache_root_env = "CCS_ABLATION_CACHE_ROOT",
    output_files = relative_outputs,
    receipt_md5 = unname(tools::md5sum(receipt_path))
  )
  cache_identity <- bensz_hash_value(list(unit_stem = unit_stem, receipt = receipt))
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(.wf_product_root)
  product_dir <- .wf_product_dir(stage)
  summary_lines <- c(
    paste0("# ", unit_stem), "",
    paste0("- External cache stage: ", stage),
    paste0("- Output files: ", length(relative_outputs)),
    paste0("- Stage receipt MD5: ", product$receipt_md5),
    if (identical(.wf_mode, "formal")) {
      "- Large objects remain under CCS_ABLATION_CACHE_ROOT."
    } else {
      "- Lightweight product and large objects remain under CCS_ABLATION_CACHE_ROOT."
    }
  )
  bensz_write_checkpoint(
    object = product,
    product_dir = product_dir,
    unit_id = sub("^([0-9]{2}\\.[0-9]{2}\\.[0-9]{2}).*$", "\\1", unit_stem),
    cache_identity = cache_identity,
    summary_lines = summary_lines,
    extra_metadata = list(external_cache_stage = stage),
    output_contract_version = 1L
  )
  product_status <- bensz_checkpoint_status(product_dir)
  if (!isTRUE(product_status$valid)) {
    stop("Cannot commit completed workflow stage: ", product_status$reason, call. = FALSE)
  }
  .wf_atomic_save_rds(
    list(
      schema_version = 2L,
      status = "complete",
      stage = stage,
      cache_identity = running_state$cache_identity,
      run_id = running_state$run_id,
      mode = running_state$mode,
      pid = running_state$pid,
      hostname = running_state$hostname,
      started = running_state$started,
      updated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      receipt_md5 = unname(tools::md5sum(receipt_path))
    ),
    file.path(stage_dir, "run-state.rds")
  )
  if (dir.exists(.wf_lock_dir) && unlink(.wf_lock_dir, recursive = TRUE, force = TRUE) != 0L) {
    stop("Stage completed but the workflow lock could not be released: ", .wf_lock_dir, call. = FALSE)
  }
  options(
    error = getOption("bensz.ablation03.previous_error"),
    bensz.ablation03.active_stage = NULL,
    bensz.ablation03.previous_error = NULL
  )
  invisible(product)
}
