#!/usr/bin/env Rscript

# Purpose: Lock nested grouped readout and paired learning-curve behavior.
# Input: Synthetic reference cohorts plus two independent external query cohorts.
# Parameters: Two lambdas, two grouped folds, and short deterministic linear fits.
# Output: Reproducible predictions, fold audit, and paired learning-curve metrics.

luckyBase::Plus.library(c("xgboost", "digest"))
source(file.path("R", "ablation.R"))

stopifnot(exists(".ablation_linear_readout", mode = "function"))
stopifnot(exists(".ablation_learning_curve", mode = "function"))

set.seed(601)
reference_cohort <- rep(paste0("R", seq_len(8)), each = 12)
reference_label <- rep(rep(c("A", "B"), each = 12), 4)
reference_signal <- ifelse(reference_label == "A", -1, 1)
reference <- cbind(
  signal = reference_signal + stats::rnorm(length(reference_signal), 0, 0.3),
  cohort_noise = rep(stats::rnorm(8, 0, 2), each = 12),
  noise = stats::rnorm(length(reference_signal))
)
rownames(reference) <- sprintf("RS%03d", seq_len(nrow(reference)))
reference_metadata <- data.frame(
  sample_id = rownames(reference),
  cohort = reference_cohort,
  tissue = "T",
  cancer_type = reference_label,
  stringsAsFactors = FALSE
)

query_label <- rep(c("A", "B"), each = 15)
query <- cbind(
  signal = ifelse(query_label == "A", -1, 1) + stats::rnorm(30, 0, 0.3),
  cohort_noise = rep(c(-3, 3), each = 15),
  noise = stats::rnorm(30)
)
rownames(query) <- sprintf("QS%03d", seq_len(nrow(query)))
query_metadata <- data.frame(
  sample_id = rownames(query),
  cohort = rep(c("Q1", "Q2"), each = 15),
  tissue = "T",
  cancer_type = query_label,
  d1_provenance = "external_frozen",
  stringsAsFactors = FALSE
)

fit <- .ablation_linear_readout(
  train = reference,
  test = query,
  train_metadata = reference_metadata,
  test_metadata = query_metadata,
  label_column = "cancer_type",
  lambda = c(0.1, 1),
  inner_folds = 2L,
  nrounds = 15L,
  numCores = 1L,
  seed = 602L
)
stopifnot(fit$status == "complete")
stopifnot(fit$selected_lambda %in% c(0.1, 1))
stopifnot(nrow(fit$inner_cv) == 2L)
stopifnot(nrow(fit$predictions) == nrow(query))
stopifnot(all(fit$predictions$sample_id == rownames(query)))
stopifnot(fit$overall$balanced_accuracy > 0.9)
stopifnot(nrow(fit$by_cohort) == 2L)
stopifnot(identical(
  fit$predictions$predicted_label,
  .ablation_linear_readout(
    train = reference,
    test = query,
    train_metadata = reference_metadata,
    test_metadata = query_metadata,
    label_column = "cancer_type",
    lambda = c(0.1, 1),
    inner_folds = 2L,
    nrounds = 15L,
    numCores = 1L,
    seed = 602L
  )$predictions$predicted_label
))

# Cohort-d1 is represented here by a cleaner one-column signal; both arms must
# share cohort subsets, seeds, folds, lambda grid, and external query samples.
d1_reference <- reference[, "signal", drop = FALSE]
d1_query <- query[, "signal", drop = FALSE]
curve <- .ablation_learning_curve(
  representations = list(
    `Direct-GSClassifier` = list(
      train = reference,
      test = query,
      blocks = NULL
    ),
    `Cohort-d1` = list(
      train = d1_reference,
      test = d1_query,
      blocks = list(M1 = 1L)
    )
  ),
  train_metadata = reference_metadata,
  test_metadata = query_metadata,
  label_column = "cancer_type",
  fractions = c(0.5, 1),
  repeats = 2L,
  lambda = c(0.1, 1),
  inner_folds = 2L,
  nrounds = 10L,
  numCores = 1L,
  seed = 603L
)
stopifnot(curve$status == "complete")
stopifnot(nrow(curve$metrics) == 8L)
stopifnot(nrow(curve$paired) == 4L)
stopifnot(all(curve$paired$test_sample_hash == curve$test_sample_hash))
stopifnot(all(curve$paired$cohort_subset_hash_direct ==
  curve$paired$cohort_subset_hash_d1))
stopifnot(all(curve$paired$delta_balanced_accuracy >= -1))
stopifnot(all(curve$paired$delta_balanced_accuracy <= 1))

message("12-test-ablation-readout: all tests passed")
