#' Stream every legal cutpoint combination and retain bounded Top N results
#'
#' @description Internal adapter for the native discrete normCCS search. All
#'   retained AUCs are checked against the existing R batch scorer, and a
#'   deterministic subset is additionally checked with the original ROC code.
#' @param df Validated patient-by-state input.
#' @param ranked Fixed ordered subtype table.
#' @param k Number of nonempty grades, between 2 and 8.
#' @param method Existing Raw or OptimizeRank ROC semantics.
#' @param started Elapsed-time origin of the calling normalization search.
#' @param verbose Whether to report progress.
#' @param objective Existing normCCS selection objective.
#' @param top_n Maximum number of retained candidates.
#' @param num_cores Requested native threads, or NULL for available cores minus one.
#' @return Raw and Summary candidate tables and explicit search diagnostics.
#' @examples
#' if (FALSE) {
#'   subtypeNorm(df, 5, norm_search_method = "stream_exact", norm_exact_top_n = 2000)
#' }
#' @useDynLib CCS, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @noRd
subtype_norm_exact <- function(df, ranked, k, method, started, verbose,
                               objective = "median", top_n = 2000L,
                               num_cores = NULL) {
  # Test
  if (FALSE) {
    df <- expand.grid(Cohort = c("A", "B"), Subtype = 1:6, response = 0:1)
    df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
    df$tumor_type <- "Lung"
    subtypeNorm(df, 5, numCores = 1, norm_search_method = "stream_exact")
  }
  norm_validate_exact_top_n(top_n)
  objective <- norm_validate_objective(objective)
  m <- nrow(ranked)
  if (k < 2L || k > 8L || k > m || m > 512L) {
    stop("stream_exact supports 2 <= K <= min(8, M) and at most 512 observed states.")
  }
  if (!is.null(num_cores) && num_cores > 1024L) {
    stop("numCores must not exceed 1024 for stream_exact.")
  }
  requested_workers <- if (is.null(num_cores)) {
    as.integer(parallelly::availableCores(omit = 1L))
  } else as.integer(num_cores)
  requested_workers <- max(1L, requested_workers)
  total <- choose(m - 1L, k - 1L)
  retained <- as.integer(min(top_n, total))
  nc <- length(unique(df$Cohort)) + 1L
  # Include count tensors, per-thread heaps, interval strings and returned rows.
  estimate <- 64 * nc * m^2 + requested_workers * retained * 128 +
    retained * (nc * 256 + m * 32) + m^3 * 16
  if (!is.finite(estimate) || estimate > 2 * 1024^3) {
    stop("stream_exact exceeds the interval/Top-N memory budget; reduce norm_exact_top_n or numCores.")
  }
  if (verbose) luckyBase::LuckyVerbose(
    paste0("subtypeNorm: stream_exact, K = ", k, ", retaining up to ", retained, " candidates..."),
    levels = 2)
  preprocessing_started <- proc.time()[["elapsed"]]
  model <- norm_prepare_counts(df, ranked)
  model$ranked_labels <- enc2utf8(as.character(ranked$Subtype))
  preprocessing_seconds <- proc.time()[["elapsed"]] - preprocessing_started
  native <- norm_exact_cpp(model, as.integer(k), method, objective, retained, requested_workers)
  if (!isTRUE(native$complete) || native$visited != total ||
      native$total_candidates != total) {
    stop("stream_exact did not cover the complete search space; no optimal result was returned.")
  }
  cuts <- native$cuts
  if (!is.matrix(cuts) || nrow(cuts) != retained || ncol(cuts) != k - 1L ||
      anyDuplicated(as.data.frame(cuts))) {
    stop("stream_exact returned an invalid candidate set.")
  }
  verification_started <- proc.time()[["elapsed"]]
  auc <- norm_auc_batch(cuts, model, method)
  native_auc <- norm_exact_auc_cpp(cuts, model, method)
  scores <- norm_objective_score(auc, model$cohorts, objective)
  means <- norm_objective_score(auc, model$cohorts, "mean")
  tolerance <- 1e-10
  delta <- max(abs(auc - native_auc), abs(scores - native$scores),
               abs(means - native$means))
  if (any(!is.finite(auc)) || !is.finite(delta) || delta > tolerance) {
    stop("Native stream_exact scores disagree with the R reference; results were not returned.")
  }
  candidates <- lapply(seq_len(nrow(cuts)), function(i) as.integer(cuts[i, ]))
  infos <- vapply(candidates, norm_info, character(1), ranked = ranked)
  # K2 has additional binary performance columns that must come from its
  # original evaluator. Higher K has only ROCAUC, already verified for all rows.
  checked <- if (k == 2L) seq_len(retained) else
    unique(as.integer(round(seq(1L, retained, length.out = min(5L, retained)))))
  reference <- norm_reference_candidates(df, ranked, candidates[checked], model, method)
  expected <- as.vector(t(auc[checked, , drop = FALSE]))
  original_delta <- max(abs(reference$ROCAUC - expected))
  if (any(!is.finite(reference$ROCAUC)) || !is.finite(original_delta) || original_delta > tolerance) {
    stop("stream_exact scores disagree with the original ROC evaluator; results were not returned.")
  }
  if (k == 2L) {
    raw <- reference
  } else {
    raw <- data.frame(
      Cohort = rep(model$cohorts, times = retained),
      Info = rep(infos, each = nc),
      nSample = rep(vapply(model$rows, length, integer(1)), times = retained),
      ROCAUC = as.vector(t(auc)), stringsAsFactors = FALSE)
  }
  if (objective == "median") {
    summary <- data.frame(Info = infos, medianROCAUC = scores, stringsAsFactors = FALSE)
    summary <- summary[order(-summary$medianROCAUC, enc2utf8(summary$Info), method = "radix"), , drop = FALSE]
  } else {
    summary <- norm_objective_summary(raw, objective)
    summary <- summary[order(-summary$objectiveScore, -summary$cohortMeanROCAUC,
                             enc2utf8(summary$Info), method = "radix"), , drop = FALSE]
  }
  rownames(summary) <- NULL
  winner <- match(summary$Info[1L], infos)
  verification_seconds <- proc.time()[["elapsed"]] - verification_started
  diagnostics <- list(
    method = "stream_exact", objective = objective,
    roc_method = method, score_semantics = if (method == "Raw" || k == 2L) "legacy_raw_auto" else "legacy_optimize_rank",
    result_policy = "top_n_only", constraint_mode = "unrestricted",
    optimality_status = "optimal_proven", exact_global_optimum = TRUE,
    score_reference_tolerance = tolerance,
    tie_break = if (objective == "median") "score_desc_then_Info_UTF8_bytes" else
      "score_desc_then_cohort_mean_desc_then_Info_UTF8_bytes",
    M = m, K = k, top_n = as.integer(top_n), returned_candidates = retained,
    n_possible_candidates = total, n_candidates_evaluated = native$visited,
    selected_cutpoints = candidates[[winner]], ranked_subtypes = as.character(ranked$Subtype),
    cohorts = model$cohorts, requested_workers = requested_workers,
    workers = native$actual_threads, openmp = native$openmp,
    n_batch_reference_values = length(auc),
    n_original_roc_verifications = nrow(reference),
    reference_max_error = max(delta, original_delta),
    interval_storage_bytes = as.numeric(utils::object.size(model)),
    preprocessing_seconds = unname(preprocessing_seconds),
    search_seconds = native$search_elapsed_seconds,
    verification_seconds = unname(verification_seconds),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started))
  if (objective != "median") {
    settings <- norm_objective_settings(objective)
    diagnostics <- c(diagnostics, settings[setdiff(names(settings), "norm_objective")])
  }
  list(Raw = raw, Summary = summary, SearchDiagnostics = diagnostics)
}

norm_validate_exact_top_n <- function(top_n) {
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n) || top_n > .Machine$integer.max) {
    stop("norm_exact_top_n must be a positive integer no larger than .Machine$integer.max.")
  }
  invisible(top_n)
}
