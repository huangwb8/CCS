#!/usr/bin/env Rscript

# Optional, read-only benchmark for an existing ablation-03 cache. Benchmark
# outputs must be written outside the scientific cache root.

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(NULL)
  args[[index + 1L]]
}
usage <- paste(
  "Usage: Rscript --vanilla test/ablation-03/scripts/benchmark-learning-curve.R",
  "--cache-root <existing-cache> --output-dir <benchmark-output> [--workers <n>]"
)
if (any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
cache_root <- value_after("--cache-root")
output_dir <- value_after("--output-dir")
if (is.null(cache_root) || !nzchar(cache_root) || is.null(output_dir) || !nzchar(output_dir)) {
  stop("--cache-root and --output-dir are required; no default paths are permitted.\n", usage, call. = FALSE)
}
workers_arg <- value_after("--workers")
workers <- suppressWarnings(as.integer(if (is.null(workers_arg)) "2" else workers_arg))
if (length(workers) != 1L || is.na(workers) || workers < 1L) {
  stop("--workers must be a positive integer.", call. = FALSE)
}

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[grepl("^--file=", all_args)]
if (length(file_arg) == 0L) stop("Cannot resolve benchmark script path.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(project_dir, "..", ".."), winslash = "/", mustWork = TRUE)
cache_root <- normalizePath(cache_root, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
cache_prefix <- paste0(tolower(cache_root), "/")
if (identical(tolower(output_dir), tolower(cache_root)) ||
    startsWith(paste0(tolower(output_dir), "/"), cache_prefix)) {
  stop("--output-dir must be outside --cache-root so benchmark data cannot alter formal cache.", call. = FALSE)
}

if (!dir.exists(file.path(cache_root, "targets"))) {
  stop("--cache-root does not contain a targets store.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Cannot create --output-dir: ", output_dir, call. = FALSE)
benchmark_lock <- file.path(output_dir, ".benchmark-lock")
if (!dir.create(benchmark_lock, showWarnings = FALSE)) {
  stop("Benchmark output is already locked: ", output_dir, call. = FALSE)
}
on.exit(unlink(benchmark_lock, recursive = TRUE, force = TRUE), add = TRUE)

Sys.setenv(
  CCS_ABLATION_CACHE_ROOT = cache_root,
  CCS_ABLATION_TARGETS = "1"
)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)
source(file.path(project_dir, "scripts", "helpers", "workflow_helpers.R"), local = TRUE)
source(.ablation03_path("02.01.00. 表示分析_functions.R"), local = TRUE)

bundle <- .wf_read("01-representations", "representation-inputs.rds")
prepared <- bundle$analysis$prepared
config <- .ae_apply_runtime_config(bundle$config, bundle$analysis)
view <- prepared$query_views[["cancer_readout"]]
query_rows <- match(view$metadata$sample_id, prepared$query_metadata$sample_id)

fraction_text <- Sys.getenv(
  "CCS_ABLATION_BENCHMARK_FRACTIONS",
  unset = as.character(config$validation$learning_fractions[[1L]])
)
fractions <- as.numeric(strsplit(fraction_text, ",", fixed = TRUE)[[1L]])
if (length(fractions) < 1L || any(!is.finite(fractions)) ||
    any(!fractions %in% config$validation$learning_fractions)) {
  stop("Benchmark fractions must be a comma-separated subset of learning_fractions.")
}
repeats <- as.integer(Sys.getenv("CCS_ABLATION_BENCHMARK_REPEATS", unset = "1"))
if (!is.finite(repeats) || repeats < 1L || repeats > config$validation$repeats) {
  stop("CCS_ABLATION_BENCHMARK_REPEATS is outside the configured repeat range.")
}
total_threads <- as.integer(Sys.getenv(
  "CCS_ABLATION_CORES",
  unset = as.character(config$validation$numCores * config$validation$workers)
))
workers <- max(1L, min(workers, total_threads))

representations <- list(
  `Direct-GSClassifier` = list(
    train = prepared$reference_direct,
    test = prepared$query_direct[query_rows, , drop = FALSE],
    blocks = NULL
  ),
  `Cohort-d1` = list(
    train = prepared$reference_d1,
    test = prepared$query_d1[query_rows, , drop = FALSE],
    blocks = prepared$selected_blocks
  )
)
run_once <- function(worker_count) {
  threads <- max(1L, floor(total_threads / worker_count))
  timing <- system.time({
    value <- .ablation_learning_curve(
      representations = representations,
      train_metadata = prepared$reference_metadata,
      test_metadata = view$metadata,
      label_column = config$anchors$primary,
      fractions = fractions,
      repeats = repeats,
      lambda = config$validation$lambda,
      inner_folds = config$validation$inner_folds,
      nrounds = config$validation$nrounds,
      numCores = threads,
      workers = worker_count,
      seed = bundle$seed + 30000L
    )
  })
  list(value = value, timing = timing, threads = threads)
}

serial <- run_once(1L)
parallel <- run_once(workers)
comparison <- all.equal(
  serial$value$metrics,
  parallel$value$metrics,
  tolerance = 1e-12,
  check.attributes = TRUE
)
if (!isTRUE(comparison)) {
  stop("Serial and parallel benchmark results differ: ", paste(comparison, collapse = "; "))
}

report <- data.frame(
  mode = c("serial", "parallel"),
  workers = c(1L, workers),
  threads_per_worker = c(serial$threads, parallel$threads),
  total_thread_budget = total_threads,
  fractions = paste(fractions, collapse = ","),
  repeats = repeats,
  elapsed_seconds = c(serial$timing[["elapsed"]], parallel$timing[["elapsed"]]),
  speedup_vs_serial = c(1, serial$timing[["elapsed"]] / parallel$timing[["elapsed"]]),
  result_equivalent = TRUE,
  r_version = R.version.string,
  ccs_version = as.character(utils::packageVersion("CCS")),
  xgboost_version = as.character(utils::packageVersion("xgboost")),
  stringsAsFactors = FALSE
)
.ablation_atomic_write_csv(report, file.path(output_dir, "learning-curve-benchmark.csv"))
print(report)
