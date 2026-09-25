# Additional inference from frozen ablation-03 artifacts. This file is sourced
# by the existing targets graph; it does not create a second execution entry.

.asi_percentile <- function(x, min_valid = 100L) {
  x <- x[is.finite(x)]
  if (length(x) < min_valid) return(c(NA_real_, NA_real_))
  unname(stats::quantile(x, c(0.025, 0.975), names = FALSE))
}

.asi_decoder_metric <- function(truth, prediction, type) {
  if (identical(type, "gene_pair")) {
    probability <- pmin(1, pmax(0, prediction))
    if (!all(c(0, 1) %in% truth)) {
      return(c(balanced_accuracy = NA_real_,
        brier = mean((probability - truth)^2)))
    }
    called <- as.integer(probability >= 0.5)
    return(c(
      balanced_accuracy = mean(c(mean(called[truth == 0] == 0),
        mean(called[truth == 1] == 1))),
      brier = mean((probability - truth)^2)
    ))
  }
  correlation <- if (length(truth) >= 2L &&
      is.finite(stats::sd(truth)) && is.finite(stats::sd(prediction)) &&
      stats::sd(truth) > 0 && stats::sd(prediction) > 0) {
    suppressWarnings(stats::cor(truth, prediction, method = "spearman"))
  } else NA_real_
  if (identical(type, "single_bin")) {
    return(c(spearman = correlation, mae = mean(abs(prediction - truth))))
  }
  c(spearman = correlation, rmse = sqrt(mean((prediction - truth)^2)))
}

.asi_decoder <- function(prepared, decoder, reference_limit,
    reference_seed, n_boot = 2000L, seed = 20260925L) {
  prediction <- decoder$prediction
  ids <- rownames(prediction)
  query_rows <- match(ids, rownames(prepared$query_direct))
  if (anyNA(query_rows) || !identical(colnames(prediction), colnames(prepared$query_direct))) {
    stop("Decoder predictions and frozen query truth do not align.", call. = FALSE)
  }
  metadata <- prepared$query_metadata[
    match(ids, prepared$query_metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id)) stop("Decoder query cohort identity is missing.", call. = FALSE)
  truth <- prepared$query_direct[query_rows, , drop = FALSE]
  types <- prepared$feature_manifest$feature_manifest$feature_type[
    match(colnames(prediction), prepared$feature_manifest$feature_manifest$feature)]
  if (anyNA(types)) stop("Decoder feature types do not align.", call. = FALSE)
  reference_metadata <- getFromNamespace(".ablation_limit_metadata", "CCS")(
    prepared$reference_metadata, reference_limit, reference_seed)
  reference_rows <- match(reference_metadata$sample_id,
    prepared$reference_metadata$sample_id)
  if (length(reference_rows) != decoder$reference_sample_count || anyNA(reference_rows)) {
    stop("Decoder baseline reference subset does not match archived fit.", call. = FALSE)
  }
  reference_mean <- colMeans(prepared$reference_direct[reference_rows, , drop = FALSE])
  cohorts <- sort(unique(as.character(metadata$cohort_key)))
  measures <- list(gene_pair = c("balanced_accuracy", "brier"),
    single_bin = c("spearman", "mae"), set_pair = c("spearman", "rmse"))
  rows <- vector("list", length(cohorts) * ncol(truth))
  index <- 0L
  for (cohort in cohorts) {
    selected <- which(metadata$cohort_key == cohort)
    for (j in seq_len(ncol(truth))) {
      metrics <- .asi_decoder_metric(truth[selected, j], prediction[selected, j], types[j])
      baseline_loss <- rep(NA_real_, length(metrics))
      names(baseline_loss) <- names(metrics)
      if (identical(types[j], "gene_pair")) {
        baseline_loss["brier"] <- mean((reference_mean[j] - truth[selected, j])^2)
      } else if (identical(types[j], "single_bin")) {
        baseline_loss["mae"] <- mean(abs(reference_mean[j] - truth[selected, j]))
      } else {
        baseline_loss["rmse"] <- sqrt(mean((reference_mean[j] - truth[selected, j])^2))
      }
      index <- index + 1L
      rows[[index]] <- data.frame(cohort = cohort, feature_type = types[j],
        feature = colnames(truth)[j], metric = names(metrics),
        value = as.numeric(metrics), baseline_loss = as.numeric(baseline_loss),
        sample_count = length(selected))
    }
  }
  per_feature <- do.call(rbind, rows)
  per_cohort <- stats::aggregate(value ~ cohort + feature_type + metric,
    data = per_feature[is.finite(per_feature$value), , drop = FALSE], FUN = mean)
  counts <- stats::aggregate(sample_count ~ cohort, data = per_feature,
    FUN = function(x) x[1L])
  per_cohort <- merge(per_cohort, counts, by = "cohort", sort = FALSE)
  loss_rows <- per_feature[is.finite(per_feature$baseline_loss) &
    is.finite(per_feature$value), , drop = FALSE]
  loss_rows$improvement <- loss_rows$baseline_loss - loss_rows$value
  loss_cohort <- stats::aggregate(cbind(baseline_loss, improvement) ~
    cohort + feature_type + metric, data = loss_rows, FUN = mean)
  keys <- expand.grid(feature_type = names(measures),
    metric = unique(unlist(measures, use.names = FALSE)), stringsAsFactors = FALSE)
  keys <- keys[vapply(seq_len(nrow(keys)), function(i) {
    keys$metric[i] %in% measures[[keys$feature_type[i]]]
  }, logical(1)), , drop = FALSE]
  set.seed(seed)
  draws <- replicate(n_boot, sample(cohorts, length(cohorts), replace = TRUE),
    simplify = FALSE)
  result <- lapply(seq_len(nrow(keys)), function(i) {
    part <- per_cohort[per_cohort$feature_type == keys$feature_type[i] &
      per_cohort$metric == keys$metric[i], , drop = FALSE]
    value <- setNames(part$value, part$cohort)
    boot <- vapply(draws, function(sampled) {
      sampled_values <- value[sampled]
      if (sum(is.finite(sampled_values)) < 3L) return(NA_real_)
      mean(sampled_values[is.finite(sampled_values)])
    }, numeric(1))
    ci <- .asi_percentile(boot)
    data.frame(feature_type = keys$feature_type[i], metric = keys$metric[i],
      sample_weighted_estimate = decoder$summary[[keys$metric[i]]][
        match(keys$feature_type[i], decoder$summary$feature_type)],
      cohort_equal_estimate = mean(part$value), ci_low = ci[1L], ci_high = ci[2L],
      n_cohort = nrow(part), n_sample = sum(part$sample_count),
      resamples = n_boot, valid_resamples = sum(is.finite(boot)), seed = seed,
      unit = "external_query_cohort", method = "cohort_percentile_bootstrap",
      condition = "frozen_decoder_and_query_atlas",
      status = if (all(is.finite(ci))) "estimable" else "not_estimable",
      reason = if (all(is.finite(ci))) "" else "fewer_than_100_valid_resamples",
      stringsAsFactors = FALSE)
  })
  loss_groups <- split(loss_cohort,
    interaction(loss_cohort$feature_type, loss_cohort$metric, drop = TRUE))
  comparisons <- lapply(seq_along(loss_groups), function(i) {
    part <- loss_groups[[i]]
    delta <- part$improvement
    observed <- mean(delta)
    set.seed(seed + 100L + i)
    bootstrap <- replicate(n_boot, mean(sample(delta, length(delta), replace = TRUE)))
    ci <- .asi_percentile(bootstrap)
    if (length(delta) < 3L) {
      p_value <- NA_real_
      p_method <- "not_estimable"
    } else if (length(delta) <= 18L) {
      signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(delta))))
      null <- as.vector(signs %*% delta / length(delta))
      p_value <- mean(abs(null) >= abs(observed))
      p_method <- "exact_cohort_sign_flip"
    } else {
      null <- replicate(n_boot,
        mean(sample(c(-1, 1), length(delta), replace = TRUE) * delta))
      p_value <- (sum(abs(null) >= abs(observed)) + 1) / (n_boot + 1)
      p_method <- "monte_carlo_cohort_sign_flip"
    }
    data.frame(feature_type = part$feature_type[1L], metric = part$metric[1L],
      estimate_baseline_minus_decoder = observed,
      ci_low = ci[1L], ci_high = ci[2L], p_value = p_value,
      p_method = p_method, n_cohort = nrow(part),
      n_query = sum(metadata$cohort_key %in% part$cohort),
      resamples = n_boot, seed = seed + 100L + i,
      unit = "external_query_cohort", baseline = "reference_feature_mean_or_prevalence",
      multiplicity_family = "three_predeclared_decoder_loss_comparisons",
      stringsAsFactors = FALSE)
  })
  comparisons <- do.call(rbind, comparisons)
  comparisons$p_value_adj <- NA_real_
  valid_p <- is.finite(comparisons$p_value)
  comparisons$p_value_adj[valid_p] <- stats::p.adjust(
    comparisons$p_value[valid_p], method = "holm", n = 3L)
  list(summary = do.call(rbind, result), per_cohort = per_cohort,
    baseline_comparisons = comparisons,
    input_key = prepared$input_key, prediction_md5 = digest::digest(prediction, algo = "md5"))
}

.asi_local <- function(data, cohort_column, delta_column, endpoint,
    min_n = 20L, n_boot = 2000L, seed = 20260925L) {
  cohorts <- sort(unique(as.character(data[[cohort_column]])))
  do.call(rbind, lapply(seq_along(cohorts), function(i) {
    part <- data[data[[cohort_column]] == cohorts[i], delta_column]
    part <- part[is.finite(part)]
    ci <- c(NA_real_, NA_real_)
    if (length(part) >= min_n) {
      set.seed(seed + i - 1L)
      ci <- .asi_percentile(replicate(n_boot,
        mean(sample(part, length(part), replace = TRUE))))
    }
    data.frame(endpoint = endpoint, cohort = cohorts[i], estimate = if (length(part)) mean(part) else NA_real_,
      ci_low = ci[1L], ci_high = ci[2L], n_query = length(part),
      resamples = if (length(part) >= min_n) n_boot else 0L,
      valid_resamples = if (all(is.finite(ci))) n_boot else 0L,
      seed = seed + i - 1L, unit = "query_within_cohort",
      method = "within_cohort_query_percentile_bootstrap",
      status = if (all(is.finite(ci))) "estimable" else "not_estimable",
      reason = if (all(is.finite(ci))) "" else "below_minimum_query_count",
      stringsAsFactors = FALSE)
  }))
}

.asi_local_readout <- function(predictions, n_boot = 2000L, seed = 20260925L) {
  direct <- predictions[predictions$representation == "Direct-GSClassifier",
    c("sample_id", "cohort", "true_label", "predicted_label")]
  d1 <- predictions[predictions$representation == "Cohort-d1",
    c("sample_id", "cohort", "true_label", "predicted_label")]
  paired <- merge(direct, d1, by = c("sample_id", "cohort", "true_label"),
    suffixes = c("_direct", "_d1"))
  if (anyDuplicated(paired$sample_id)) stop("Readout query pairs are duplicated.", call. = FALSE)
  paired$delta <- as.numeric(paired$predicted_label_d1 == paired$true_label) -
    as.numeric(paired$predicted_label_direct == paired$true_label)
  .asi_local(paired, "cohort", "delta", "readout_accuracy_d1_minus_direct",
    n_boot = n_boot, seed = seed)
}

.asi_local_biology <- function(per_query, n_boot = 2000L, seed = 20260925L) {
  direct <- per_query[per_query$representation == "Direct-GSClassifier",
    c("anchor", "query_sample", "query_cohort", "utility")]
  d1 <- per_query[per_query$representation == "Cohort-d1",
    c("anchor", "query_sample", "query_cohort", "utility")]
  paired <- merge(direct, d1, by = c("anchor", "query_sample", "query_cohort"),
    suffixes = c("_direct", "_d1"))
  if (anyDuplicated(paired[c("anchor", "query_sample")])) {
    stop("Biology anchor query pairs are duplicated.", call. = FALSE)
  }
  paired$delta <- paired$utility_d1 - paired$utility_direct
  anchors <- sort(unique(as.character(paired$anchor)))
  result <- lapply(seq_along(anchors), function(i) {
    part <- paired[paired$anchor == anchors[i], , drop = FALSE]
    local <- .asi_local(part, "query_cohort", "delta",
      paste0("anchor_", anchors[i], "_utility_d1_minus_direct"),
      n_boot = n_boot, seed = seed + 100L * i)
    local$anchor <- anchors[i]
    local
  })
  do.call(rbind, result)
}

.asi_node_stats <- function(pairs, n_boot = 2000L, seed = 20260925L) {
  nodes <- sort(unique(c(as.character(pairs$cohort_a), as.character(pairs$cohort_b))))
  weight_a <- match(pairs$cohort_a, nodes)
  weight_b <- match(pairs$cohort_b, nodes)
  value <- pairs$delta_d1_minus_direct
  weighted_median <- function(x, w) {
    ordered <- order(x)
    x <- x[ordered]
    w <- w[ordered]
    x[which(cumsum(w) >= sum(w) / 2)[1L]]
  }
  set.seed(seed)
  draws <- replicate(n_boot, {
    multiplicity <- tabulate(sample.int(length(nodes), length(nodes), replace = TRUE),
      nbins = length(nodes))
    weights <- multiplicity[weight_a] * multiplicity[weight_b]
    if (sum(weights) == 0) return(c(NA_real_, NA_real_, NA_real_))
    c(mean = stats::weighted.mean(value, weights),
      median = weighted_median(value[weights > 0], weights[weights > 0]),
      fraction_nonnegative = stats::weighted.mean(value >= 0, weights))
  })
  list(draws = draws, valid = sum(is.finite(draws[1L, ])), nodes = nodes,
    index_a = weight_a, index_b = weight_b)
}

.asi_node_jackknife <- function(value, index_a, index_b, n_nodes) {
  estimates <- vapply(seq_len(n_nodes), function(i) {
    keep <- index_a != i & index_b != i
    if (sum(keep) < 3L) return(NA_real_)
    mean(value[keep])
  }, numeric(1))
  if (anyNA(estimates)) return(NA_real_)
  sqrt((n_nodes - 1) / n_nodes * sum((estimates - mean(estimates))^2))
}

.asi_node_coverage <- function(index_a, index_b, n_nodes, seed,
    n_simulation = 300L) {
  set.seed(seed)
  coverage <- vapply(c(0.1, 1), function(noise_sd) {
    covered <- replicate(n_simulation, {
      node_effect <- stats::rnorm(n_nodes)
      outcome <- node_effect[index_a] + node_effect[index_b] +
        stats::rnorm(length(index_a), sd = noise_sd)
      se <- .asi_node_jackknife(outcome, index_a, index_b, n_nodes)
      is.finite(se) && abs(mean(outcome)) <=
        stats::qt(0.975, df = n_nodes - 1L) * se
    })
    mean(covered)
  }, numeric(1))
  min(coverage)
}

.asi_structural <- function(pairs, n_boot = 2000L, seed = 20260925L) {
  directions <- unique(as.character(pairs$direction))
  rows <- list()
  index <- 0L
  for (direction in directions) {
    directional <- pairs[pairs$direction == direction, , drop = FALSE]
    for (scope in c("all_cohort_pairs", "same_cancer_type")) {
      part <- if (identical(scope, "same_cancer_type")) {
        directional[directional$same_cancer_type, , drop = FALSE]
      } else directional
      if (!nrow(part)) next
      index <- index + 1L
      boot <- .asi_node_stats(part, n_boot, seed + index - 1L)
      n_nodes <- length(boot$nodes)
      degree <- tabulate(c(boot$index_a, boot$index_b), nbins = n_nodes)
      mean_ci <- .asi_percentile(boot$draws[1L, ])
      median_ci <- .asi_percentile(boot$draws[2L, ])
      fraction_ci <- .asi_percentile(boot$draws[3L, ])
      status <- "estimable"
      reason <- ""
      p_value <- NA_real_
      coverage <- NA_real_
      if (n_nodes < 12L || min(degree) < 3L ||
          boot$valid < ceiling(0.9 * n_boot)) {
        status <- "not_estimable"
        reason <- "sparse_node_network_or_low_valid_bootstrap_rate"
      } else {
        coverage <- .asi_node_coverage(boot$index_a, boot$index_b,
          n_nodes, seed + 1000L + index)
        if (!is.finite(coverage) || coverage < 0.9 || coverage > 0.99) {
          status <- "not_estimable"
          reason <- "node_jackknife_failed_null_simulation_coverage"
        } else {
          se <- .asi_node_jackknife(part$delta_d1_minus_direct,
            boot$index_a, boot$index_b, n_nodes)
          if (!is.finite(se) || se <= 0) {
            status <- "not_estimable"
            reason <- "node_jackknife_variance_unavailable"
          } else {
            p_value <- 2 * stats::pt(-abs(mean(part$delta_d1_minus_direct) / se),
              df = n_nodes - 1L)
          }
        }
      }
      rows[[index]] <- data.frame(direction = direction, scope = scope,
        mean_delta = mean(part$delta_d1_minus_direct),
        mean_ci_low = mean_ci[1L], mean_ci_high = mean_ci[2L],
        median_delta = stats::median(part$delta_d1_minus_direct),
        median_ci_low = median_ci[1L], median_ci_high = median_ci[2L],
        fraction_d1_ge_direct = mean(part$delta_d1_minus_direct >= 0),
        fraction_ci_low = fraction_ci[1L], fraction_ci_high = fraction_ci[2L],
        p_value = p_value, n_cohort = n_nodes, n_pairs = nrow(part),
        min_node_degree = min(degree), resamples = n_boot,
        valid_resamples = boot$valid, null_simulation_coverage = coverage,
        seed = seed + index - 1L, unit = "shared_cohort_node",
        ci_method = "cohort_node_percentile_bootstrap",
        p_method = if (is.finite(p_value)) "node_jackknife_t" else "not_estimable",
        status = status, reason = reason, stringsAsFactors = FALSE)
    }
  }
  result <- do.call(rbind, rows)
  result$p_value_adj <- NA_real_
  estimable <- is.finite(result$p_value)
  if (any(estimable)) result$p_value_adj[estimable] <-
    stats::p.adjust(result$p_value[estimable], method = "holm", n = 4L)
  result$multiplicity_family <- "four_direction_scope_mean_deltas"
  result
}

.asi_matched_bank <- function(repeats) {
  groups <- split(repeats, interaction(repeats$direction, repeats$scope, drop = TRUE))
  do.call(rbind, lapply(groups, function(part) {
    values <- part$mean_delta[is.finite(part$mean_delta)]
    ci <- if (length(values) >= 5L) {
      unname(stats::quantile(values, c(0.025, 0.975)))
    } else c(NA_real_, NA_real_)
    data.frame(direction = part$direction[1L], scope = part$scope[1L],
      design_median = stats::median(values), design_q025 = ci[1L],
      design_q975 = ci[2L], n_design = length(values),
      unit = "matched_bank_design_repeat", method = "empirical_design_distribution",
      condition = "frozen_patients_and_module_pool", p_value = NA_real_,
      reason = "overlapping_modules_are_not_independent_new_cohorts",
      stringsAsFactors = FALSE)
  }))
}

.ablation03_statistical_inference <- function(representation_inputs,
    representation_analysis, biology_analysis, structural_analysis,
    cache_root) {
  prepared <- readRDS(file.path(representation_inputs$directory,
    "representation-inputs.rds"))$analysis$prepared
  representation_dir <- representation_analysis$directory
  structural_dir <- structural_analysis$directory
  retrieval <- readRDS(file.path(representation_dir, "retrieval.rds"))
  readout <- readRDS(file.path(representation_dir, "readout.rds"))
  decoder <- readRDS(file.path(representation_dir, "tradeoffs.rds"))$decoder
  manifest <- readRDS(file.path(representation_dir, "manifest.rds"))
  pairs <- utils::read.csv(file.path(structural_dir,
    "structural_directional_pair_comparisons.csv"))
  repeats <- utils::read.csv(file.path(structural_dir,
    "structural_matched_bank_repeats.csv"))
  biology_queries <- readRDS(file.path(biology_analysis$directory,
    "anchor_per_query_utility.rds"))
  retrieval_pairs <- retrieval$paired[retrieval$paired$k == 30L, , drop = FALSE]
  result <- list(
    decoder = .asi_decoder(prepared, decoder,
      reference_limit = manifest$config$tradeoffs$decoder_max_reference_samples,
      reference_seed = manifest$seed + 40000L),
    retrieval_local = .asi_local(retrieval_pairs, "cohort",
      "delta_top_k_label_rate", "retrieval_top30_d1_minus_direct"),
    readout_local = .asi_local_readout(readout$predictions),
    biology_local = .asi_local_biology(biology_queries),
    structural = .asi_structural(pairs),
    matched_bank = .asi_matched_bank(repeats),
    input_key = prepared$input_key,
    source_md5 = unname(tools::md5sum(c(
      file.path(representation_dir, "tradeoffs.rds"),
      file.path(biology_analysis$directory, "anchor_per_query_utility.rds"),
      file.path(structural_dir, "structural_directional_pair_comparisons.csv"),
      file.path(structural_dir, "structural_matched_bank_repeats.csv"))))
  )
  output <- file.path(cache_root, "statistical-inference", "summary.rds")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output)
  output
}

.ablation03_learning_query_inference <- function(representation_inputs,
    representation_analysis, cache_root, n_boot = 2000L,
    bootstrap_seed = 20260926L) {
  bundle <- readRDS(file.path(representation_inputs$directory,
    "representation-inputs.rds"))
  prepared <- bundle$analysis$prepared
  result_dir <- representation_analysis$directory
  manifest <- readRDS(file.path(result_dir, "manifest.rds"))
  original <- readRDS(file.path(result_dir, "learning_curve.rds"))$metrics
  fractions <- manifest$config$validation$learning_fractions
  fraction_index <- match(1, fractions)
  if (is.na(fraction_index)) stop("Learning curve has no 100% design.", call. = FALSE)
  repeat_id <- 1L
  fit_seed <- manifest$seed + 30000L + fraction_index * 1000L + repeat_id
  test_metadata <- prepared$query_views[["cancer_readout"]]$metadata
  test_rows <- match(test_metadata$sample_id, prepared$query_metadata$sample_id)
  if (anyNA(test_rows)) stop("Learning query identity does not align.", call. = FALSE)
  fit_fun <- getFromNamespace(".ablation_linear_readout", "CCS")
  arms <- list(
    `Direct-GSClassifier` = list(train = prepared$reference_direct,
      test = prepared$query_direct[test_rows, , drop = FALSE], blocks = NULL),
    `Cohort-d1` = list(train = prepared$reference_d1,
      test = prepared$query_d1[test_rows, , drop = FALSE],
      blocks = prepared$selected_blocks)
  )
  fits <- lapply(names(arms), function(name) {
    arm <- arms[[name]]
    fit <- fit_fun(train = arm$train, test = arm$test,
      train_metadata = prepared$reference_metadata,
      test_metadata = test_metadata, label_column = manifest$anchor,
      lambda = manifest$config$validation$lambda,
      inner_folds = manifest$config$validation$inner_folds,
      nrounds = manifest$config$validation$nrounds,
      numCores = manifest$config$validation$numCores,
      seed = fit_seed, blocks = arm$blocks)
    baseline <- original[original$requested_fraction == 1 &
      original$repeat_id == repeat_id & original$representation == name, ]
    if (nrow(baseline) != 1L ||
        !isTRUE(all.equal(fit$selected_lambda, baseline$selected_lambda,
          tolerance = 1e-10)) ||
        !isTRUE(all.equal(fit$overall$accuracy, baseline$accuracy,
          tolerance = 1e-10)) ||
        !isTRUE(all.equal(fit$overall$balanced_accuracy,
          baseline$balanced_accuracy, tolerance = 1e-10)) ||
        !isTRUE(all.equal(fit$overall$macro_auroc,
          baseline$macro_auroc, tolerance = 1e-10))) {
      stop("Recomputed 100% learning score differs from the archived arm: ",
        name, call. = FALSE)
    }
    fit$predictions
  })
  names(fits) <- names(arms)
  direct <- fits[["Direct-GSClassifier"]][,
    c("sample_id", "cohort", "true_label", "predicted_label")]
  d1 <- fits[["Cohort-d1"]][,
    c("sample_id", "cohort", "true_label", "predicted_label")]
  paired <- merge(direct, d1, by = c("sample_id", "cohort", "true_label"),
    suffixes = c("_direct", "_d1"))
  if (nrow(paired) != nrow(test_metadata) || anyDuplicated(paired$sample_id)) {
    stop("100% learning predictions are not complete query pairs.", call. = FALSE)
  }
  paired$delta <- as.numeric(paired$predicted_label_d1 == paired$true_label) -
    as.numeric(paired$predicted_label_direct == paired$true_label)
  cohort <- stats::aggregate(delta ~ cohort, data = paired, FUN = mean)
  set.seed(bootstrap_seed)
  sampled <- replicate(n_boot,
    mean(sample(cohort$delta, nrow(cohort), replace = TRUE)))
  ci <- .asi_percentile(sampled)
  result <- list(summary = data.frame(
    requested_fraction = 1, repeat_id = repeat_id,
    training_seed = fit_seed, bootstrap_seed = bootstrap_seed,
    estimate = mean(cohort$delta), ci_low = ci[1L], ci_high = ci[2L],
    n_cohort = nrow(cohort), n_query = nrow(paired),
    resamples = n_boot, valid_resamples = sum(is.finite(sampled)),
    unit = "external_query_cohort", method = "cohort_percentile_bootstrap",
    condition = "fixed_100_percent_training_design_and_seed",
    design_level_ci = NA_real_, design_level_p_value = NA_real_,
    score_recalculation = "matched_archived_100_percent_arm",
    stringsAsFactors = FALSE), cohort = cohort,
    input_key = prepared$input_key)
  output <- file.path(cache_root, "statistical-inference", "learning-query.rds")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output)
  output
}

.asi_weighted_rank <- function(value, weight) {
  ordering <- order(value)
  sorted <- value[ordering]
  groups <- cumsum(c(TRUE, diff(sorted) != 0))
  group_weight <- as.numeric(rowsum(weight[ordering], groups, reorder = FALSE))
  group_rank <- cumsum(group_weight) - group_weight / 2 + 0.5
  result <- numeric(length(value))
  result[ordering] <- group_rank[groups]
  result
}

.asi_geometry_bootstrap <- function(direct, d1, cohort, first, second,
    direct_distance, d1_distance, k, n_boot = 200L, seed = 20260925L) {
  cohort_ids <- match(cohort, sort(unique(cohort)))
  n_cohort <- max(cohort_ids)
  n <- nrow(direct)
  samples <- matrix(NA_real_, n_boot, 3L, dimnames = list(NULL,
    c("linear_cka", "distance_spearman", "knn_jaccard")))
  cka <- function(x, y, weight) {
    size <- sum(weight)
    x_sum <- colSums(x * weight)
    y_sum <- colSums(y * weight)
    x_weighted <- x * sqrt(weight)
    y_weighted <- y * sqrt(weight)
    xy <- crossprod(x_weighted, y_weighted) - tcrossprod(x_sum, y_sum) / size
    xx <- crossprod(x_weighted) - tcrossprod(x_sum) / size
    yy <- crossprod(y_weighted) - tcrossprod(y_sum) / size
    denominator <- sqrt(sum(xx^2) * sum(yy^2))
    if (!is.finite(denominator) || denominator == 0) return(NA_real_)
    sum(xy^2) / denominator
  }
  pair_spearman <- function(weight) {
    pair_weight <- weight[first] * weight[second]
    keep <- pair_weight > 0
    x <- direct_distance[keep]
    y <- d1_distance[keep]
    w <- pair_weight[keep]
    if (length(x) < 3L) return(NA_real_)
    x_rank <- .asi_weighted_rank(x, w)
    y_rank <- .asi_weighted_rank(y, w)
    x_rank <- x_rank - weighted.mean(x_rank, w)
    y_rank <- y_rank - weighted.mean(y_rank, w)
    denominator <- sqrt(sum(w * x_rank^2) * sum(w * y_rank^2))
    if (denominator == 0) return(NA_real_)
    sum(w * x_rank * y_rank) / denominator
  }
  neighbor_jaccard <- function(weight) {
    keep <- which(weight > 0)
    if (length(keep) <= k) return(NA_real_)
    x_neighbors <- getFromNamespace(".ablation_knn", "CCS")(
      direct[keep, , drop = FALSE], k)
    y_neighbors <- getFromNamespace(".ablation_knn", "CCS")(
      d1[keep, , drop = FALSE], k)
    agreement <- vapply(seq_along(keep), function(i) {
      length(intersect(x_neighbors[i, ], y_neighbors[i, ])) /
        length(union(x_neighbors[i, ], y_neighbors[i, ]))
    }, numeric(1))
    weighted.mean(agreement, weight[keep])
  }
  for (i in seq_len(n_boot)) {
    set.seed(seed + i - 1L)
    multiplicity <- tabulate(sample.int(n_cohort, n_cohort, replace = TRUE),
      nbins = n_cohort)
    weight <- multiplicity[cohort_ids]
    samples[i, ] <- c(cka(direct, d1, weight), pair_spearman(weight),
      neighbor_jaccard(weight))
  }
  metrics <- colnames(samples)
  intervals <- t(vapply(seq_along(metrics), function(i) {
    .asi_percentile(samples[, i])
  }, numeric(2)))
  valid <- colSums(is.finite(samples))
  data.frame(endpoint = metrics, ci_low = intervals[, 1L],
    ci_high = intervals[, 2L], p_value = NA_real_, p_value_adj = NA_real_,
    n_cohort = n_cohort, n_sample = n, resamples = n_boot,
    valid_resamples = as.integer(valid), seed = seed,
    unit = "reference_cohort", method = "cohort_percentile_bootstrap",
    condition = "frozen_model_and_reference_atlas_fixed_pairs",
    status = ifelse(valid >= 100L, "estimable", "not_estimable"),
    reason = ifelse(valid >= 100L, NA_character_, "too_few_valid_resamples"),
    stringsAsFactors = FALSE)
}

.asi_geometry_leave_one <- function(prepared, manifest, leave_limit = Inf,
    n_boot = 200L, bootstrap_seed = 20260925L) {
  config <- manifest$config$geometry
  reference_metadata <- prepared$reference_metadata
  sample_rows <- if (nrow(reference_metadata) > config$geometry_samples) {
    getFromNamespace(".ablation_stratified_sample", "CCS")(
      reference_metadata, config$geometry_samples, manifest$seed)
  } else seq_len(nrow(reference_metadata))
  direct_transform <- getFromNamespace(".ablation_scale_train_apply", "CCS")(
    prepared$reference_direct, prepared$query_direct)
  d1_transform <- getFromNamespace(".ablation_module_balanced_transform", "CCS")(
    prepared$reference_d1, prepared$query_d1, prepared$selected_blocks)
  direct <- direct_transform$train[sample_rows, , drop = FALSE]
  d1 <- d1_transform$reference[sample_rows, , drop = FALSE]
  cohort <- as.character(reference_metadata$cohort_key[sample_rows])
  all_cohorts <- sort(unique(cohort))
  cohorts <- head(all_cohorts, leave_limit)
  n <- nrow(direct)
  k <- min(max(config$k), n - 1L)
  max_omitted <- max(tabulate(match(cohort, all_cohorts),
    nbins = length(all_cohorts)))
  neighbor_count <- min(n - 1L, k + max_omitted)
  direct_neighbors <- getFromNamespace(".ablation_knn", "CCS")(
    direct, neighbor_count)
  d1_neighbors <- getFromNamespace(".ablation_knn", "CCS")(
    d1, neighbor_count)
  direct_sums <- colSums(direct)
  d1_sums <- colSums(d1)
  direct_cross <- crossprod(direct)
  d1_cross <- crossprod(d1)
  between_cross <- crossprod(direct, d1)
  n_pairs <- min(as.integer(config$distance_pairs), n * (n - 1) / 2)
  set.seed(manifest$seed)
  first <- sample.int(n, n_pairs, replace = TRUE)
  second <- sample.int(n, n_pairs, replace = TRUE)
  same <- first == second
  while (any(same)) {
    second[same] <- sample.int(n, sum(same), replace = TRUE)
    same <- first == second
  }
  direct_distance <- sqrt(rowSums((direct[first, , drop = FALSE] -
    direct[second, , drop = FALSE])^2))
  d1_distance <- sqrt(rowSums((d1[first, , drop = FALSE] -
    d1[second, , drop = FALSE])^2))
  cka_from_moments <- function(x_cross, y_cross, xy_cross, x_sum, y_sum, size) {
    xy_centered <- xy_cross - tcrossprod(x_sum, y_sum) / size
    x_centered <- x_cross - tcrossprod(x_sum) / size
    y_centered <- y_cross - tcrossprod(y_sum) / size
    denominator <- sqrt(sum(x_centered^2) * sum(y_centered^2))
    if (!is.finite(denominator) || denominator == 0) return(NA_real_)
    sum(xy_centered^2) / denominator
  }
  jaccard_after_omission <- function(keep, omit) {
    mean(vapply(keep, function(i) {
      x <- head(direct_neighbors[i, !direct_neighbors[i, ] %in% omit], k)
      y <- head(d1_neighbors[i, !d1_neighbors[i, ] %in% omit], k)
      length(intersect(x, y)) / length(union(x, y))
    }, numeric(1)))
  }
  rows <- lapply(seq_along(cohorts), function(i) {
    omit <- which(cohort == cohorts[i])
    keep <- which(cohort != cohorts[i])
    x <- direct[omit, , drop = FALSE]
    y <- d1[omit, , drop = FALSE]
    pair_keep <- cohort[first] != cohorts[i] & cohort[second] != cohorts[i]
    data.frame(omitted_cohort = cohorts[i], n_sample = length(keep),
      linear_cka = cka_from_moments(
        direct_cross - crossprod(x), d1_cross - crossprod(y),
        between_cross - crossprod(x, y),
        direct_sums - colSums(x), d1_sums - colSums(y), length(keep)),
      distance_spearman = suppressWarnings(stats::cor(
        direct_distance[pair_keep], d1_distance[pair_keep],
        method = "spearman")),
      knn_jaccard = jaccard_after_omission(keep, omit),
      unit = "reference_cohort_omission",
      method = "leave_one_cohort_fixed_pairs_sensitivity",
      condition = "frozen_model_and_reference_atlas",
      stringsAsFactors = FALSE)
  })
  baseline <- c(
    linear_cka = cka_from_moments(direct_cross, d1_cross,
      between_cross, direct_sums, d1_sums, n),
    distance_spearman = suppressWarnings(stats::cor(direct_distance,
      d1_distance, method = "spearman")),
    knn_jaccard = jaccard_after_omission(seq_len(n), integer())
  )
  inference <- .asi_geometry_bootstrap(direct, d1, cohort, first, second,
    direct_distance, d1_distance, k, n_boot, bootstrap_seed)
  list(values = do.call(rbind, rows), inference = inference,
    baseline = baseline,
    n_cohort = length(all_cohorts),
    n_sample = length(sample_rows), input_key = prepared$input_key,
    seed = manifest$seed)
}

.ablation03_geometry_sensitivity <- function(representation_inputs,
    representation_analysis, cache_root) {
  prepared <- readRDS(file.path(representation_inputs$directory,
    "representation-inputs.rds"))$analysis$prepared
  manifest <- readRDS(file.path(representation_analysis$directory, "manifest.rds"))
  result <- .asi_geometry_leave_one(prepared, manifest)
  native <- readRDS(file.path(representation_analysis$directory,
    "native_geometry.rds"))$metrics
  expected <- native$metric_value[match(names(result$baseline), native$metric_name)]
  if (anyNA(expected) || any(abs(result$baseline - expected) > 1e-8)) {
    stop("Geometry resampling baseline differs from the frozen native estimate.",
      call. = FALSE)
  }
  output <- file.path(cache_root, "statistical-inference", "geometry-leave-one.rds")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output)
  output
}
