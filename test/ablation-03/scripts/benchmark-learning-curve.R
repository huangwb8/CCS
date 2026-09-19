# Benchmark serial versus bounded-PSOCK learning-curve execution on persisted
# ablation-03 inputs. This tool does not alter scientific result products.
bootstrap <- c(
  file.path("scripts", "helpers", "workflow_helpers.R"),
  file.path("test", "ablation-03", "scripts", "helpers", "workflow_helpers.R")
)
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
source(.ablation03_repo_path("R", "ablation.R"))
source(.ablation03_path("02.01.00. 表示分析_functions.R"))

if (!dir.create(.wf_lock_dir, showWarnings = FALSE)) {
  stop(
    "The ablation-03 external cache is locked. Do not benchmark concurrently; ",
    "if the owner crashed, confirm that process has stopped and remove the lock manually."
  )
}
on.exit(unlink(.wf_lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
.wf_atomic_save_rds(
  list(stage = "learning-curve-benchmark", pid = Sys.getpid()),
  file.path(.wf_lock_dir, "owner.rds")
)

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
workers <- as.integer(Sys.getenv("CCS_ABLATION_BENCHMARK_WORKERS", unset = "2"))
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
  xgboost_version = as.character(utils::packageVersion("xgboost")),
  stringsAsFactors = FALSE
)
output_dir <- file.path(.wf_output("ablation-experiment"), "performance")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
.ablation_atomic_write_csv(report, file.path(output_dir, "learning-curve-benchmark.csv"))
print(report)
