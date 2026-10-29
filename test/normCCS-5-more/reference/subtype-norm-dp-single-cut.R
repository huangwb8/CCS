# Experimental first iteration: single-boundary refinement only.
# Uses the same validated CCS 0.8.5 count/ROC kernels; no global optimum claim.
subtype_norm_dp_single_cut <- function(df, ranked, k, method, started, verbose) {
    if (FALSE) {
        model <- norm_prepare_counts(df, ranked)
        norm_solve_dp(model$loss[1L, , ], k)
    }
    if (verbose) 
        luckyBase::LuckyVerbose("subtypeNorm: DP starts and median-AUC refinement...", levels = 2)
    model <- norm_prepare_counts(df, ranked)
    nc <- length(model$cohorts)
    weights <- rbind(rep(1/nc, nc), c(0, rep(1/(nc - 1L), nc - 1L)), diag(nc))
    weights <- rbind(weights, -weights)
    loss_matrix <- matrix(model$loss, nrow = nc)
    starts <- lapply(seq_len(nrow(weights)), function(i) {
        cost <- apply(model$loss * weights[i, ], c(2, 3), sum)
        norm_solve_dp(cost, k)$cuts
    })
    starts <- unique(c(starts, list(as.integer(floor(seq_len(k - 1L) * model$m/k)))))
    score <- function(cuts) apply(norm_auc_batch(cuts, model, method), 1, stats::median)
    start_scores <- score(do.call(rbind, starts))
    selected <- order(-start_scores, seq_along(starts))[seq_len(min(6L, length(starts)))]
    candidates <- starts[selected]
    evaluations <- length(starts)
    refinements <- vector("list", length(selected))
    for (i in seq_along(selected)) {
        cuts <- starts[[selected[i]]]
        best <- start_scores[selected[i]]
        converged <- FALSE
        pair_moves <- 0L
        for (pass in seq_len(10L)) {
            improved <- FALSE
            for (j in seq_along(cuts)) {
                lower <- if (j == 1L) 
                  1L
                else cuts[j - 1L] + 1L
                upper <- if (j == length(cuts)) 
                  model$m - 1L
                else cuts[j + 1L] - 1L
                proposals <- matrix(rep(cuts, each = upper - lower + 1L), ncol = length(cuts))
                proposals[, j] <- lower:upper
                values <- score(proposals)
                evaluations <- evaluations + length(values)
                winner <- which.max(values)
                if (values[winner] > best + 1e-12) {
                  cuts <- as.integer(proposals[winner, ])
                  best <- values[winner]
                  improved <- TRUE
                }
            }
            NULL
            if (!improved) {
                converged <- TRUE
                break
            }
        }
        candidates[[length(candidates) + 1L]] <- cuts
        refinements[[i]] <- list(initial_auc = start_scores[selected[i]], final_auc = best, passes = pass, 
            converged = converged, pair_moves = pair_moves)
    }
    candidates <- unique(candidates)
    infos <- vapply(candidates, norm_info, character(1), ranked = ranked)
    raw <- lapply(seq_along(candidates), function(i) {
        cuts <- candidates[[i]]
        group <- k - 1L - findInterval(match(as.character(df$Subtype), as.character(ranked$Subtype)) - 
            1L, cuts)
        do.call(rbind, lapply(seq_along(model$rows), function(c) {
            ids <- model$rows[[c]]
            data.frame(Cohort = model$cohorts[c], Info = infos[i], nSample = length(ids), get_perform_markers(df$response[ids], 
                group[ids], k, method), stringsAsFactors = FALSE)
        }))
    })
    raw <- do.call(rbind, raw)
    final_scores <- vapply(split(raw$ROCAUC, raw$Info), stats::median, numeric(1), na.rm = TRUE)
    fast_auc <- norm_auc_batch(do.call(rbind, candidates), model, method)
    reference_auc <- matrix(raw$ROCAUC, nrow = length(candidates), byrow = TRUE)
    if (any(!is.finite(fast_auc)) || any(!is.finite(reference_auc)) || any(abs(reference_auc - fast_auc) > 
        1e-10)) {
        stop("Fast normCCS scoring disagrees with the reference ROC evaluator; results were not returned.")
    }
    summary <- data.frame(Info = names(final_scores), medianROCAUC = unname(final_scores), stringsAsFactors = FALSE)
    summary <- summary[order(-summary$medianROCAUC, summary$Info), , drop = FALSE]
    rownames(summary) <- NULL
    winner <- match(summary$Info[1L], infos)
    list(Raw = raw, Summary = summary, SearchDiagnostics = list(method = "dp", objective = "legacy_median_auc_including_All", 
        roc_method = method, refinement_start_limit = 6L, refinement_pass_limit = 10L, refinement_neighborhood = "single_boundary", 
        refinement_pair_batch_size = 512L, surrogate = "fixed_direction_weighted_mean_raw_auc", exact_global_optimum = FALSE, 
        dp_subproblems_exact = TRUE, M = model$m, K = k, n_dp_starts = nrow(weights), n_unique_starts = length(starts), 
        n_unique_candidates = length(candidates), n_scored_proposals = evaluations, n_reference_verifications = length(reference_auc), 
        selected_cutpoints = candidates[[winner]], ranked_subtypes = as.character(ranked$Subtype), best_initial_median_auc = max(start_scores), 
        refinements = refinements, cohort_only_median_auc = stats::median(raw$ROCAUC[raw$Info == summary$Info[1L] & 
            raw$Cohort != "All"]), workers = 1L, interval_storage_bytes = as.numeric(object.size(model)), 
        elapsed_seconds = unname(proc.time()[["elapsed"]] - started)))
}
