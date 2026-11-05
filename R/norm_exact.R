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
#' @param min_grade_n Minimum pooled training records per grade.
#' @param regularization Nonnegative multiplier of mean inverse-root grade size.
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
                               num_cores = NULL, min_grade_n = 1L,
                               regularization = 0) {
  # Test
  if (FALSE) {
    df <- expand.grid(Cohort = c("A", "B"), Subtype = 1:6, response = 0:1)
    df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
    df$tumor_type <- "Lung"
    subtypeNorm(df, 5, numCores = 1, norm_search_method = "stream_exact")
  }
  norm_validate_exact_top_n(top_n)
  norm_validate_regularization(min_grade_n, regularization)
  regularized <- min_grade_n > 1L || regularization > 0
  objective <- norm_validate_objective(objective)
  m <- nrow(ranked)
  if (!is.null(num_cores) && num_cores > 1024L) {
    stop("numCores must not exceed 1024 for stream_exact.")
  }
  requested_workers <- if (is.null(num_cores)) {
    as.integer(parallelly::availableCores(omit = 1L))
  } else as.integer(num_cores)
  requested_workers <- max(1L, requested_workers)
  nc <- length(unique(df$Cohort)) + 1L
  cost <- norm_exact_search_cost(m, k, nc, method, requested_workers)
  norm_report_exact_cost(cost, verbose)
  if (k < 2L || k > 8L || k > m || m > 512L) {
    stop("stream_exact supports 2 <= K <= min(8, M) and at most 512 observed states.")
  }
  total <- cost$n_candidates
  retained <- as.integer(min(top_n, total))
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
  native <- norm_exact_cpp(model, as.integer(k), method, objective, retained,
                           requested_workers, min_grade_n, regularization)
  if (!isTRUE(native$complete) || native$visited != total ||
      native$total_candidates != total) {
    stop("stream_exact did not cover the complete search space; no optimal result was returned.")
  }
  retained <- as.integer(min(top_n, native$feasible_candidates))
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
  support <- if (regularized) norm_exact_support(cuts, df, ranked) else NULL
  selection_scores <- if (regularized) scores - regularization * support$penalty else scores
  tolerance <- 1e-10
  # Selection scores can be very large in magnitude for a large lambda;
  # use relative error there while preserving absolute AUC verification.
  selection_delta <- abs(selection_scores - native$scores) / pmax(1, abs(selection_scores))
  delta <- max(abs(auc - native_auc), selection_delta,
               abs(means - native$means))
  if (regularized && any(support$sizes < min_grade_n)) {
    stop("Native stream_exact returned a mapping below norm_min_grade_n.")
  }
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
  if (regularized) {
    indices <- match(summary$Info, infos)
    summary$supportPenalty <- support$penalty[indices]
    summary$regularizationPenalty <- regularization * support$penalty[indices]
    summary$regularizedScore <- selection_scores[indices]
    summary$minGradeN <- apply(support$sizes[indices, , drop = FALSE], 1L, min)
    summary$gradeSampleSizes <- lapply(indices, function(i) support$sizes[i, ])
    ordering <- if (objective == "median") {
      order(-summary$regularizedScore, enc2utf8(summary$Info), method = "radix")
    } else {
      order(-summary$regularizedScore, -summary$cohortMeanROCAUC,
            enc2utf8(summary$Info), method = "radix")
    }
    summary <- summary[ordering, , drop = FALSE]
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
    estimated_search_seconds = cost$seconds,
    runtime_estimate_basis = cost$basis,
    preprocessing_seconds = unname(preprocessing_seconds),
    search_seconds = native$search_elapsed_seconds,
    verification_seconds = unname(verification_seconds),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started))
  if (regularized) {
    diagnostics$constraint_mode <- "pooled_training_grade_support"
    diagnostics$norm_min_grade_n <- min_grade_n
    diagnostics$norm_regularization <- regularization
    diagnostics$selection_score <- "objective_minus_lambda_mean_inverse_sqrt_grade_n"
    diagnostics$optimality_scope <- "fixed_order_K_support_feasible_regularized_objective"
    diagnostics$n_feasible_candidates <- native$feasible_candidates
    diagnostics$n_rejected_candidates <- total - native$feasible_candidates
    diagnostics$support_unit <- "input_records_pooled_once"
    diagnostics$selected_grade_sizes <- summary$gradeSampleSizes[[1L]]
    diagnostics$tie_break <- if (objective == "median") {
      "regularized_score_desc_then_Info_UTF8_bytes"
    } else "regularized_score_desc_then_cohort_mean_desc_then_Info_UTF8_bytes"
  }
  if (objective != "median") {
    settings <- norm_objective_settings(objective)
    diagnostics <- c(diagnostics, settings[setdiff(names(settings), "norm_objective")])
  }
  list(Raw = raw, Summary = summary, SearchDiagnostics = diagnostics)
}

# Independent record-based support calculation, deliberately separate from
# native pooled-prefix counts. Columns follow the high-to-low state order.
norm_exact_support <- function(cuts, df, ranked) {
  # Test
  if (FALSE) {
    norm_exact_support(matrix(c(1L, 3L), 1L),
      data.frame(Subtype = c("a", "a", "b", "c", "d")),
      data.frame(Subtype = letters[1:4]))
  }
  state_n <- tabulate(match(as.character(df$Subtype), as.character(ranked$Subtype)),
                      nbins = nrow(ranked))
  sizes <- t(vapply(seq_len(nrow(cuts)), function(i) {
    endpoints <- c(0L, cuts[i, ], nrow(ranked))
    vapply(seq_len(length(endpoints) - 1L), function(g) {
      sum(state_n[seq.int(endpoints[g] + 1L, endpoints[g + 1L])])
    }, numeric(1))
  }, numeric(ncol(cuts) + 1L)))
  colnames(sizes) <- as.character(seq.int(ncol(cuts), 0L))
  list(sizes = sizes, penalty = rowMeans(1 / sqrt(sizes)))
}

norm_validate_regularization <- function(min_grade_n, regularization,
                                         search_method = "stream_exact") {
  if (!is.numeric(min_grade_n) || length(min_grade_n) != 1L ||
      !is.finite(min_grade_n) || min_grade_n < 1 ||
      min_grade_n != floor(min_grade_n) || min_grade_n > .Machine$integer.max) {
    stop("norm_min_grade_n must be a positive integer no larger than .Machine$integer.max.")
  }
  if (!is.numeric(regularization) || length(regularization) != 1L ||
      !is.finite(regularization) || regularization < 0) {
    stop("norm_regularization must be a finite nonnegative number.")
  }
  if (search_method != "stream_exact" && (min_grade_n != 1L || regularization != 0)) {
    stop("norm_min_grade_n and norm_regularization require norm_search_method = 'stream_exact'.")
  }
  invisible(NULL)
}

# Reference-hardware estimate, not a benchmark of the invoking machine or a CI.
# Calibration: tmp/normCCS-stream-exact-088/REPORT.md, i9-14900K, 16 threads,
# 21 scoring coordinates (20 real cohorts plus All), median objective, Top 2000.
norm_exact_search_cost <- function(m, k, nc, method, workers) {
  # Test
  if (FALSE) {
    norm_exact_search_cost(107L, 8L, 21L, "Raw", 16L)
  }
  total <- choose(m - 1L, k - 1L)
  ranked <- method == "OptimizeRank" && k > 2L
  per_candidate <- if (ranked) 7.570 / choose(106, 5) else
    409.035 / choose(106, 7)
  # Raw reuses unchanged prefixes; OptimizeRank recomputes K group scores.
  rank_scale <- if (ranked) k / 6 else 1
  thread_scale <- 16 / min(workers, 16L)
  seconds <- total * per_candidate * (nc / 21) * thread_scale * rank_scale
  list(n_candidates = total, seconds = seconds, M = m, K = k,
       method = method, workers = workers,
       explosion = total >= 1e8 || seconds >= 60,
       basis = "reference_i9_14900K_16_threads_C21_Top2000_median")
}

norm_format_search_time <- function(seconds) {
  if (!is.finite(seconds)) return("beyond a finite estimate")
  if (seconds < 1) return("less than 1 second")
  units <- c("seconds", "minutes", "hours", "days", "years")
  scales <- c(1, 60, 3600, 86400, 365.25 * 86400)
  unit <- max(which(seconds >= scales))
  paste(formatC(seconds / scales[unit], format = "f", digits = 1), units[unit])
}

norm_report_exact_cost <- function(cost, verbose) {
  if (!verbose) return(invisible(NULL))
  count <- format(cost$n_candidates, scientific = cost$n_candidates >= 1e15,
                  big.mark = ",", trim = TRUE, digits = 7)
  luckyBase::LuckyVerbose(
    paste0("subtypeNorm: stream_exact search size: M = ", cost$M,
           ", K = ", cost$K, ", choose(M - 1, K - 1) = ", count,
           " combinations; estimated native search time: ",
           norm_format_search_time(cost$seconds), " (", cost$method,
           ", ", cost$workers, " requested threads)."), levels = 2)
  if (cost$explosion) luckyBase::LuckyVerbose(
    paste0("subtypeNorm: WARNING: combinatorial explosion. Every combination ",
           "must be evaluated; reducing norm_exact_top_n only reduces retained ",
           "results. Consider a smaller K or norm_search_method = 'dp' ",
           "(approximate for the public objective)."), levels = 2)
  luckyBase::LuckyVerbose(
    paste0("subtypeNorm: Rough reference-hardware estimate, not a 95% interval ",
           "or a local benchmark: i9-14900K, 16 OpenMP threads, 20 cohorts + All, ",
           "median objective, Top 2000; scaled by cohort count and requested ",
           "threads (speedup capped at 16), with a K/6 factor for OptimizeRank. ",
           "Assumes OpenMP parallelism; hardware, objective, ties and Top N can ",
           "change runtime. Excludes preprocessing and verification. ",
           "K > 8 is an unsupported extrapolation."), levels = 2)
  invisible(NULL)
}

norm_validate_exact_top_n <- function(top_n) {
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n) || top_n > .Machine$integer.max) {
    stop("norm_exact_top_n must be a positive integer no larger than .Machine$integer.max.")
  }
  invisible(top_n)
}
