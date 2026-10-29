# Internal normCCS search kernels. States are ordered from high to low RR/NRR.
# No RNG, persistent cache, native code, or parallel backend is used here.

norm_validate_input <- function(df, k, method, cores) {
  required <- c("Cohort", "Subtype", "response", "SampleIDs", "tumor_type")
  if (!is.data.frame(df) || !all(required %in% names(df)) || !nrow(df)) {
    stop("df must be a nonempty data frame with Cohort, Subtype, response, SampleIDs and tumor_type.")
  }
  if (anyNA(df[required]) || !all(as.character(df$response) %in% c("0", "1"))) {
    stop("Required columns must be complete; response must contain only 0 and 1.")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2 || k != floor(k)) {
    stop("n_norm_subtype must be a single integer >= 2.")
  }
  if (length(method) != 1L || !method %in% c("Raw", "OptimizeRank")) {
    stop("norm_roc_method must be 'Raw' or 'OptimizeRank'.")
  }
  if (!is.null(cores) && (!is.numeric(cores) || length(cores) != 1L ||
      !is.finite(cores) || cores < 1 || cores != floor(cores))) {
    stop("numCores must be NULL or a positive integer.")
  }
  if (any(vapply(split(df$response, df$Cohort, drop = TRUE),
                 function(x) length(unique(x)) != 2L, logical(1)))) {
    stop("Every cohort must contain both response classes (0 and 1).")
  }
  if (any(as.character(df$Cohort) == "All")) {
    stop("Cohort name 'All' is reserved for the pooled summary.")
  }
}

# Prefix counts and interval loss relative to keeping all M distinct scores.
# L[a,b] = U[a,b] - P[a,b]*N[a,b]/2; negative loss is legal.
norm_prepare_counts <- function(df, ranked) {
  m <- nrow(ranked)
  cohorts <- sort(unique(as.character(df$Cohort)))
  state <- match(as.character(df$Subtype), as.character(ranked$Subtype))
  rows <- c(list(seq_len(nrow(df))), lapply(cohorts, function(x) which(df$Cohort == x)))
  names(rows) <- c("All", cohorts)
  # pROC respects factor level order, even when '0' is the case level.
  roc_levels <- levels(droplevels(as.factor(df$response)))
  p <- n <- matrix(0, length(rows), m)
  loss <- first <- array(0, c(length(rows), m, m))
  median_states <- matrix(0L, length(rows), 4L)
  for (c in seq_along(rows)) {
    ids <- rows[[c]]
    positive <- as.character(df$response[ids]) == roc_levels[2L]
    p[c, ] <- tabulate(state[ids[positive]], nbins = m)
    n[c, ] <- tabulate(state[ids[!positive]], nbins = m)
    cp <- c(0, cumsum(p[c, ]))
    cn <- c(0, cumsum(n[c, ]))
    pair_prefix <- c(0, cumsum(p[c, ] * (cn[-1L] - n[c, ] / 2)))
    first_state <- match(seq_len(m), state[ids])
    first_state[is.na(first_state)] <- Inf
    for (a in seq_len(m)) {
      b <- a:m
      pg <- cp[b + 1L] - cp[a]
      ng <- cn[b + 1L] - cn[a]
      loss[c, a, b] <- (pg * cn[b + 1L] -
        (pair_prefix[b + 1L] - pair_prefix[a]) - pg * ng / 2) / (sum(p[c, ]) * sum(n[c, ]))
      first[c, a, b] <- cummin(first_state[b])
    }
    for (j in 1:2) {
      counts <- if (j == 1L) p[c, ] else n[c, ]
      total <- sum(counts)
      median_states[c, (2L * j - 1L):(2L * j)] <-
        vapply(c(ceiling(total / 2), floor(total / 2) + 1),
               function(x) which(cumsum(counts) >= x)[1L], integer(1))
    }
  }
  prefix <- function(x) t(apply(cbind(0, x), 1, cumsum))
  base_auc <- rowSums(p * (rowSums(n) - t(apply(n, 1, cumsum)) + n / 2)) /
    (rowSums(p) * rowSums(n))
  list(m = m, cohorts = names(rows), rows = rows, p = p, n = n,
       roc_levels = roc_levels,
       cp = prefix(p), cn = prefix(n), loss = loss, first = first,
       median_states = median_states, base_auc = base_auc)
}

# Exact O(K*M^2) solution for a weighted fixed-direction mean AUC.
norm_solve_dp <- function(cost, k) {
  m <- nrow(cost)
  value <- matrix(Inf, k, m)
  parent <- matrix(0L, k, m)
  value[1L, ] <- cost[1L, ]
  for (g in 2:k) {
    for (b in g:m) {
      ends <- (g - 1L):(b - 1L)
      choices <- value[g - 1L, ends] + cost[cbind(ends + 1L, b)]
      best <- which.min(choices)
      value[g, b] <- choices[best]
      parent[g, b] <- ends[best]
    }
  }
  cuts <- integer(k - 1L)
  b <- m
  for (g in k:2) {
    cuts[g - 1L] <- parent[g, b]
    b <- cuts[g - 1L]
  }
  list(cuts = cuts, loss = value[k, m])
}

norm_median_bin <- function(counts) {
  cumulative <- counts
  if (ncol(counts) > 1L) for (j in 2:ncol(counts)) {
    cumulative[, j] <- cumulative[, j - 1L] + counts[, j]
  }
  total <- rowSums(counts)
  (rowSums(cumulative < ceiling(total / 2)) +
     rowSums(cumulative < (floor(total / 2) + 1))) / 2 + 1
}

# Batched legacy AUC, including pROC's median-based automatic direction.
# OptimizeRank preserves first-observation tie order in one-vs-rest AUCs.
norm_auc_batch <- function(cuts, model, method = "Raw") {
  if (is.null(dim(cuts))) cuts <- matrix(cuts, nrow = 1L)
  nr <- nrow(cuts)
  k <- ncol(cuts) + 1L
  m <- model$m
  left <- cbind(1L, cuts + 1L)
  right <- cbind(cuts, m)
  interval <- left + (right - 1L) * m
  result <- matrix(0, nr, length(model$cohorts), dimnames = list(NULL, model$cohorts))
  for (c in seq_along(model$cohorts)) {
    if (method == "Raw" || k == 2L) {
      auc <- model$base_auc[c] - rowSums(matrix(model$loss[c, , ][as.vector(interval)], nr, k))
      median_group <- matrix(vapply(model$median_states[c, ],
                                    function(s) rowSums(cuts < s), numeric(nr)), nr, 4L)
      reverse <- rowSums(median_group[, 3:4, drop = FALSE]) <
        rowSums(median_group[, 1:2, drop = FALSE])
    } else {
      p <- matrix(model$cp[c, right + 1L] - model$cp[c, left], nr, k)
      n <- matrix(model$cn[c, right + 1L] - model$cn[c, left], nr, k)
      pt <- rowSums(p)
      nt <- rowSums(n)
      # Match pROC's trapezoid arithmetic, not just its algebraic value:
      # last-bit differences change OptimizeRank's ordering of nearly tied AUCs.
      sensitivity <- p / pt
      specificity <- (nt - n) / nt
      ovr <- ((1 + sensitivity) / 2) * specificity +
        (sensitivity / 2) * (1 - specificity)
      # A binary indicator has median 0, 0.5 or 1 depending on its count.
      binary_median <- function(x, total) (x >= total / 2) / 2 + (x > total / 2) / 2
      flip <- binary_median(n, nt) > binary_median(p, pt)
      reverse_sensitivity <- (pt - p) / pt
      reverse_specificity <- n / nt
      reverse_auc <- ((1 + reverse_sensitivity) / 2) * reverse_specificity +
        (reverse_sensitivity / 2) * (1 - reverse_specificity)
      ovr[flip] <- reverse_auc[flip]
      first <- matrix(model$first[c, , ][as.vector(interval)], nr, k)
      # Absent groups never occur in unique(predictor) in the reference.
      ovr[p + n == 0] <- Inf
      ordering <- order(rep(seq_len(nr), k), as.vector(ovr), as.vector(first))
      p <- matrix(p[ordering], nr, k, byrow = TRUE)
      n <- matrix(n[ordering], nr, k, byrow = TRUE)
      cumulative_n <- n
      if (k > 1L) for (j in 2:k) cumulative_n[, j] <- cumulative_n[, j - 1L] + n[, j]
      auc <- rowSums(p * (cumulative_n - n / 2)) / (pt * nt)
      reverse <- norm_median_bin(n) > norm_median_bin(p)
    }
    result[, c] <- ifelse(reverse, 1 - auc, auc)
  }
  result
}

norm_info <- function(cuts, ranked) {
  k <- length(cuts) + 1L
  groups <- k - 1L - findInterval(seq_len(nrow(ranked)) - 1L, cuts)
  paste(vapply(0:(k - 1L), function(g) {
    paste0(g, "=", paste(ranked$Subtype[groups == g], collapse = ","))
  }, character(1)), collapse = "; ")
}

subtype_norm_dp <- function(df, ranked, k, method, started, verbose, objective = "median") {
  # Test
  if (FALSE) {
    # The synthetic generator and complete runnable examples live in
    # test/normCCS-5-more/; no patient data are needed.
    model <- norm_prepare_counts(df, ranked)
    norm_solve_dp(model$loss[1L, , ], k)
  }
  objective <- norm_validate_objective(objective)
  if (verbose) luckyBase::LuckyVerbose(paste0("subtypeNorm: DP starts and ", objective, " refinement..."), levels = 2)
  model <- norm_prepare_counts(df, ranked)
  nc <- length(model$cohorts)
  weights <- rbind(rep(1 / nc, nc), c(0, rep(1 / (nc - 1L), nc - 1L)), diag(nc))
  # Also seed reverse directions because legacy pROC auto-direction can flip.
  weights <- rbind(weights, -weights)
  loss_matrix <- matrix(model$loss, nrow = nc)
  starts <- lapply(seq_len(nrow(weights)), function(i) {
    cost <- matrix(as.vector(crossprod(weights[i, ], loss_matrix)), model$m, model$m)
    norm_solve_dp(cost, k)$cuts
  })
  starts <- unique(c(starts, list(as.integer(floor(seq_len(k - 1L) * model$m / k)))))
  score <- function(cuts) norm_objective_score(norm_auc_batch(cuts, model, method), model$cohorts, objective)
  candidate_order <- function(values, cuts) {
    if (objective == "median") return(order(-values, seq_along(values)))
    auc <- norm_auc_batch(cuts, model, method)
    means <- norm_objective_score(auc, model$cohorts, "mean")
    infos <- apply(cuts, 1L, norm_info, ranked = ranked)
    order(-values, -means, infos)
  }
  proposal_winner <- function(values, proposals) {
    if (objective == "median") return(which.max(values))
    ties <- which(values == max(values))
    if (length(ties) == 1L) return(ties)
    ties[candidate_order(values[ties], proposals[ties, , drop = FALSE])[1L]]
  }
  start_scores <- score(do.call(rbind, starts))
  selected <- candidate_order(start_scores, do.call(rbind, starts))[seq_len(min(6L, length(starts)))]
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
        lower <- if (j == 1L) 1L else cuts[j - 1L] + 1L
        upper <- if (j == length(cuts)) model$m - 1L else cuts[j + 1L] - 1L
        proposals <- matrix(rep(cuts, each = upper - lower + 1L), ncol = length(cuts))
        proposals[, j] <- lower:upper
        values <- score(proposals)
        evaluations <- evaluations + length(values)
        winner <- proposal_winner(values, proposals)
        if (values[winner] > best + 1e-12) {
          cuts <- as.integer(proposals[winner, ])
          best <- values[winner]
          improved <- TRUE
        }
      }
      # A single-boundary optimum can still improve by moving two neighboring
      # boundaries together. Search that neighborhood only after singles stall.
      if (!improved && length(cuts) > 1L) {
        for (j in seq_len(length(cuts) - 1L)) {
          lower <- if (j == 1L) 1L else cuts[j - 1L] + 1L
          upper <- if (j + 1L == length(cuts)) model$m - 1L else cuts[j + 2L] - 1L
          if (upper - lower < 2L) next
          first <- rep(lower:(upper - 1L), times = (upper - lower):1L)
          second <- unlist(lapply(lower:(upper - 1L), function(a) (a + 1L):upper), use.names = FALSE)
          pair_best <- best
          pair_cuts <- cuts
          for (offset in seq.int(1L, length(first), by = 512L)) {
            ids <- offset:min(offset + 511L, length(first))
            proposals <- matrix(rep(cuts, each = length(ids)), ncol = length(cuts))
            proposals[, j] <- first[ids]
            proposals[, j + 1L] <- second[ids]
            values <- score(proposals)
            evaluations <- evaluations + length(values)
            winner <- proposal_winner(values, proposals)
            if (values[winner] > pair_best + 1e-12) {
              pair_cuts <- as.integer(proposals[winner, ])
              pair_best <- values[winner]
            } else if (objective != "median" && values[winner] == pair_best &&
                       pair_best > best + 1e-12) {
              tied <- rbind(pair_cuts, proposals[winner, ])
              preferred <- proposal_winner(c(pair_best, values[winner]), tied)
              pair_cuts <- as.integer(tied[preferred, ])
            }
          }
          if (pair_best > best + 1e-12) {
            cuts <- pair_cuts
            best <- pair_best
            improved <- TRUE
            pair_moves <- pair_moves + 1L
          }
        }
      }
      if (!improved) {
        converged <- TRUE
        break
      }
    }
    candidates[[length(candidates) + 1L]] <- cuts
    refinements[[i]] <- list(initial_auc = start_scores[selected[i]], final_auc = best,
                             passes = pass, converged = converged, pair_moves = pair_moves)
  }
  candidates <- unique(candidates)
  infos <- vapply(candidates, norm_info, character(1), ranked = ranked)
  raw <- lapply(seq_along(candidates), function(i) {
    cuts <- candidates[[i]]
    group <- k - 1L - findInterval(match(as.character(df$Subtype), as.character(ranked$Subtype)) - 1L, cuts)
    do.call(rbind, lapply(seq_along(model$rows), function(c) {
      ids <- model$rows[[c]]
      data.frame(Cohort = model$cohorts[c], Info = infos[i], nSample = length(ids),
                 get_perform_markers(df$response[ids], group[ids], k, method),
                 stringsAsFactors = FALSE)
    }))
  })
  raw <- do.call(rbind, raw)
  final_scores <- vapply(split(raw$ROCAUC, raw$Info), stats::median, numeric(1), na.rm = TRUE)
  fast_auc <- norm_auc_batch(do.call(rbind, candidates), model, method)
  reference_auc <- matrix(raw$ROCAUC, nrow = length(candidates), byrow = TRUE)
  if (any(!is.finite(fast_auc)) || any(!is.finite(reference_auc)) ||
      any(abs(reference_auc - fast_auc) > 1e-10)) {
    stop("Fast normCCS scoring disagrees with the reference ROC evaluator; results were not returned.")
  }
  if (any(abs(norm_objective_score(fast_auc, model$cohorts, objective) -
              norm_objective_score(reference_auc, model$cohorts, objective)) > 1e-10)) {
    stop("Fast normCCS objective disagrees with reference ROC objective; results were not returned.")
  }
  summary <- data.frame(Info = names(final_scores), medianROCAUC = unname(final_scores), stringsAsFactors = FALSE)
  summary <- summary[order(-summary$medianROCAUC, summary$Info), , drop = FALSE]
  rownames(summary) <- NULL
  if (objective != "median") summary <- norm_objective_summary(raw, objective)
  winner <- match(summary$Info[1L], infos)
  result <- list(Raw = raw, Summary = summary, SearchDiagnostics = list(
    method = "dp", objective = "legacy_median_auc_including_All",
    roc_method = method, refinement_start_limit = 6L, refinement_pass_limit = 10L,
    refinement_neighborhood = "single_and_adjacent_pair_boundaries",
    refinement_pair_batch_size = 512L,
    surrogate = "fixed_direction_weighted_mean_raw_auc",
    exact_global_optimum = FALSE, dp_subproblems_exact = TRUE,
    M = model$m, K = k, n_dp_starts = nrow(weights), n_unique_starts = length(starts),
    n_unique_candidates = length(candidates), n_scored_proposals = evaluations,
    n_reference_verifications = length(reference_auc),
    selected_cutpoints = candidates[[winner]], ranked_subtypes = as.character(ranked$Subtype),
    best_initial_median_auc = max(start_scores), refinements = refinements,
    cohort_only_median_auc = stats::median(raw$ROCAUC[raw$Info == summary$Info[1L] & raw$Cohort != "All"]),
    workers = 1L, interval_storage_bytes = as.numeric(utils::object.size(model)),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started)))
  if (objective != "median") {
    diagnostics <- result$SearchDiagnostics
    diagnostics$objective <- objective
    diagnostics$objective_scope <- "cohorts_only"
    diagnostics$n_real_cohorts <- nc - 1L
    diagnostics$best_initial_median_auc <- NULL
    diagnostics$best_initial_objective_score <- max(start_scores)
    diagnostics$best_final_objective_score <- summary$objectiveScore[1L]
    diagnostics$refinements <- lapply(refinements, function(x) {
      names(x)[1:2] <- c("initial_objective_score", "final_objective_score")
      x
    })
    settings <- norm_objective_settings(objective)
    diagnostics <- c(diagnostics, settings[setdiff(names(settings), "norm_objective")])
    result$SearchDiagnostics <- diagnostics
  }
  result
}
