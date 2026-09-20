#!/usr/bin/env Rscript

# The single supported ablation-03 runner. It delegates scientific work to the
# numbered R scripts in fresh R sessions and requires an explicit cache root.

args <- commandArgs(trailingOnly = TRUE)

.entry_usage <- function() {
  paste(
    "Usage:",
    "Rscript --vanilla test/ablation-03/run-ablation-03.R",
    "--cache-root <path> [--profile formal|lightweight]",
    "[--cores <n>] [--workers <n>] [--memory-gb <n>] [--run-id <id>]",
    "[--through <01.01.00|01.02.00|01.03.00|02.01.00|02.02.00|02.03.00>]"
  )
}

.entry_parse_args <- function(args) {
  value_flags <- c(
    "--cache-root", "--profile", "--cores", "--workers", "--memory-gb",
    "--run-id", "--through"
  )
  if (length(args) == 1L && args[[1L]] %in% c("-h", "--help")) {
    cat(.entry_usage(), "\n")
    quit(save = "no", status = 0L)
  }
  unknown <- args[grepl("^-", args) & !args %in% value_flags]
  if (length(unknown) > 0L) {
    stop("Unknown option(s): ", paste(unknown, collapse = ", "), "\n", .entry_usage(), call. = FALSE)
  }
  parsed <- list()
  i <- 1L
  while (i <= length(args)) {
    flag <- args[[i]]
    if (!flag %in% value_flags || i == length(args)) {
      stop("Missing value for ", flag, ".\n", .entry_usage(), call. = FALSE)
    }
    parsed[[sub("^--", "", flag)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  parsed
}

.entry_integer <- function(value, label) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(label, " must be a positive integer.", call. = FALSE)
  }
  parsed
}

.entry_number <- function(value, label) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || is.na(parsed) || !is.finite(parsed) || parsed <= 0) {
    stop(label, " must be a positive number.", call. = FALSE)
  }
  parsed
}

.entry_script_file <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- all_args[grepl("^--file=", all_args)]
  if (length(file_arg) == 0L) {
    stop("ablation-03: cannot resolve the runner path.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
}

parsed <- .entry_parse_args(args)
if (is.null(parsed[["cache-root"]]) || !nzchar(parsed[["cache-root"]])) {
  stop("--cache-root is required; no default cache path is permitted.\n", .entry_usage(), call. = FALSE)
}
profile <- parsed$profile
if (is.null(profile)) profile <- "formal"
if (!profile %in% c("formal", "lightweight")) {
  stop("--profile must be formal or lightweight.", call. = FALSE)
}
cores <- .entry_integer(if (is.null(parsed$cores)) "1" else parsed$cores, "--cores")
workers <- .entry_integer(if (is.null(parsed$workers)) "1" else parsed$workers, "--workers")
memory_gb <- .entry_number(
  if (is.null(parsed[["memory-gb"]])) "16" else parsed[["memory-gb"]],
  "--memory-gb"
)
if (workers > cores) stop("--workers cannot exceed --cores.", call. = FALSE)

entry_script <- .entry_script_file()
project_dir <- normalizePath(dirname(entry_script), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(project_dir, "..", ".."), winslash = "/", mustWork = TRUE)
identity_helper <- file.path(project_dir, "scripts", "helpers", "run_identity_helpers.R")
source(identity_helper, local = TRUE)

cache_root <- normalizePath(parsed[["cache-root"]], winslash = "/", mustWork = FALSE)
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(cache_root)) stop("Cannot create --cache-root: ", cache_root, call. = FALSE)

stages <- c(
  `01.01.00` = "01.01.00. 数据准备.R",
  `01.02.00` = "01.02.00. 表示输入准备.R",
  `01.03.00` = "01.03.00. 生物输入准备.R",
  `02.01.00` = "02.01.00. 表示分析.R",
  `02.02.00` = "02.02.00. 生物锚点分析.R",
  `02.03.00` = "02.03.00. 结构复现分析.R"
)
through <- parsed$through
if (!is.null(through)) {
  if (!through %in% names(stages)) {
    stop("--through must name a numbered analysis stage.", call. = FALSE)
  }
  stages <- stages[seq_len(match(through, names(stages)))]
}

package_version <- tryCatch(
  as.character(utils::packageVersion("CCS")),
  error = function(error) NA_character_
)
if (is.na(package_version) || !identical(package_version, "0.8.3")) {
  stop("ablation-03 requires installed CCS 0.8.3; found ", package_version, ".", call. = FALSE)
}

lock <- .ablation03_acquire_entry_lock(
  cache_root = cache_root,
  mode = profile,
  entry_script = entry_script,
  run_id = parsed[["run-id"]]
)
lock_released <- FALSE
on.exit({
  if (!lock_released) .ablation03_release_entry_lock(lock)
}, add = TRUE)

Sys.setenv(
  CCS_ABLATION_CACHE_ROOT = cache_root,
  CCS_ABLATION_MODE = profile,
  CCS_ABLATION_RUN_ID = lock$run_id,
  CCS_ABLATION_ENTRY_SCRIPT = entry_script,
  CCS_ABLATION_CORES = as.character(cores),
  CCS_ABLATION_WORKERS = as.character(workers),
  CCS_ABLATION_MEMORY_GB = as.character(memory_gb)
)
if (identical(.Platform$OS.type, "windows")) {
  Sys.setenv(LANG = "Chinese_China.utf8", LC_CTYPE = "Chinese_China.utf8")
}
if (identical(profile, "lightweight")) {
  defaults <- c(
    CCS_ABLATION_MAX_REFERENCE_SAMPLES = "1000",
    CCS_ABLATION_MAX_QUERY_SAMPLES = "400",
    CCS_ABLATION_GEOMETRY_SAMPLES = "500",
    CCS_ABLATION_DISTANCE_PAIRS = "1000",
    CCS_ABLATION_LEARNING_FRACTIONS = "0.5,1",
    CCS_ABLATION_REPEATS = "1",
    CCS_ABLATION_NROUNDS = "2",
    CCS_ABLATION_SCALING_ENABLED = "false",
    CCS_ABLATION_NULL_CONTROLS = "false",
    CCS_ABLATION_DECODER_ENABLED = "false",
    CCS_ABLATION_STRUCTURAL_MIN_ENTITY_N = "2",
    CCS_ABLATION_STRUCTURAL_BOOTSTRAP = "20",
    CCS_ABLATION_MATCHED_REPEATS = "2"
  )
  for (name in names(defaults)) {
    if (!nzchar(Sys.getenv(name, unset = ""))) {
      do.call(Sys.setenv, as.list(stats::setNames(defaults[[name]], name)))
    }
  }
}

metadata <- list(
  schema_version = 1L,
  run_id = lock$run_id,
  status = "running",
  mode = profile,
  cache_root = cache_root,
  entry_script = entry_script,
  started_at = lock$started_at,
  pid = Sys.getpid(),
  hostname = unname(Sys.info()[["nodename"]]),
  R = R.version.string,
  package = list(
    name = "CCS",
    version = package_version,
    path = find.package("CCS"),
    ablation_source_md5 = unname(tools::md5sum(file.path(repo_root, "R", "ablation.R")))
  ),
  resources = list(cores = cores, workers = workers, memory_gb = memory_gb),
  stages = unname(stages),
  completed_stages = character()
)
.ablation03_write_run_record(metadata)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)
rscript <- file.path(R.home("bin"), if (identical(.Platform$OS.type, "windows")) "Rscript.exe" else "Rscript")
stage_expression <- paste(
  "if (identical(.Platform$OS.type, 'windows'))",
  "invisible(Sys.setlocale('LC_CTYPE', 'Chinese_China.utf8'));",
  "stage_id <- Sys.getenv('CCS_ABLATION_STAGE_ID');",
  "stage_dir <- Sys.getenv('CCS_ABLATION_PROJECT_DIR');",
  "pattern <- paste0('^', gsub('[.]', '[.]', stage_id), '[.] .*[.]R$');",
  "scripts <- list.files(stage_dir, pattern = pattern, full.names = TRUE);",
  "scripts <- scripts[!grepl('_functions[.]R$', scripts)];",
  "if (length(scripts) != 1L) stop('Cannot resolve one numbered stage script: ', stage_id);",
  "source(scripts[[1L]], encoding = 'UTF-8')"
)
Sys.setenv(CCS_ABLATION_PROJECT_DIR = project_dir)

message("ablation-03 run: ", lock$run_id)
message("profile: ", profile)
message("cache.root: ", cache_root)
message("resources: ", workers, " worker(s), ", cores, " total core(s), ", memory_gb, " GB")

failed <- TRUE
on.exit({
  if (failed) {
    metadata$status <- "failed"
    metadata$finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    try(.ablation03_write_run_record(metadata), silent = TRUE)
  }
}, add = TRUE)

for (stage_id in names(stages)) {
  stage <- unname(stages[[stage_id]])
  stage_path <- normalizePath(
    file.path("test", "ablation-03", stage),
    winslash = "/",
    mustWork = TRUE
  )
  message("\n>>> ", stage_path)
  Sys.setenv(CCS_ABLATION_STAGE_ID = stage_id)
  status <- system2(rscript, c("--vanilla", "-e", shQuote(stage_expression)))
  if (!identical(status, 0L)) {
    metadata$failed_stage <- unname(stage)
    metadata$exit_status <- status
    stop("ablation-03 stage failed with exit status ", status, ": ", stage, call. = FALSE)
  }
  metadata$completed_stages <- c(metadata$completed_stages, unname(stage))
  .ablation03_write_run_record(metadata)
}

metadata$status <- "complete"
metadata$finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
.ablation03_write_run_record(metadata)
failed <- FALSE
.ablation03_release_entry_lock(lock)
lock_released <- TRUE
message("\nablation-03 completed: ", cache_root)
