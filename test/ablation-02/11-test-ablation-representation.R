#!/usr/bin/env Rscript

# Purpose: Lock the representation-only ablation contracts before implementation.
# Input: Small deterministic matrices and two tiny frozen XGBoost ensembles.
# Parameters: Fixed seeds, exact neighbors, and label-balanced grouped folds.
# Output: Explicit stopifnot failures for any contract regression.

luckyBase::Plus.library(c("xgboost", "digest", "dbscan"))
source(file.path("R", "ablation.R"))

expect_function <- function(name) {
  stopifnot(exists(name, mode = "function", inherits = TRUE))
}

expect_function(".ablation_predict_module_from_direct")
expect_function(".ablation_module_balanced_transform")
expect_function(".ablation_query_reference_retrieval")
expect_function(".ablation_null_perm_eligibility")
expect_function(".ablation_evidence_level")
expect_function(".ablation_validate_neighbor_search")

# Should reproduce the frozen ensemble median when Direct features are precomputed.
set.seed(11)
direct <- matrix(
  stats::rnorm(80),
  nrow = 40,
  dimnames = list(sprintf("S%02d", seq_len(40)), c("f1", "f2"))
)
label <- as.integer(direct[, "f1"] > direct[, "f2"])

fit_binary <- function(target, seed) {
  set.seed(seed)
  xgboost::xgboost(
    data = direct,
    label = target,
    booster = "gblinear",
    objective = "binary:logistic",
    nrounds = 8L,
    eta = 0.2,
    lambda = 0.1,
    alpha = 0,
    nthread = 1L,
    verbose = 0
  )
}

repeat_one <- list(
  `1` = list(bst = fit_binary(label, 101L), genes = colnames(direct)),
  `2` = list(bst = fit_binary(1L - label, 102L), genes = colnames(direct))
)
repeat_two <- list(
  `1` = list(bst = fit_binary(label, 201L), genes = colnames(direct)),
  `2` = list(bst = fit_binary(1L - label, 202L), genes = colnames(direct))
)
model <- list(Model = list(repeat_one, repeat_two))

predicted <- .ablation_predict_module_from_direct(
  direct = direct,
  model = model,
  module_id = "T1|C1"
)
expected <- cbind(
  `T1|C1|1` = apply(cbind(
    predict(repeat_one[["1"]]$bst, direct),
    predict(repeat_two[["1"]]$bst, direct)
  ), 1, stats::median),
  `T1|C1|2` = apply(cbind(
    predict(repeat_one[["2"]]$bst, direct),
    predict(repeat_two[["2"]]$bst, direct)
  ), 1, stats::median)
)
stopifnot(identical(dim(predicted), dim(expected)))
stopifnot(max(abs(predicted - expected)) < 1e-12)

# Should give each d1 module equal total weight regardless of block width.
reference <- rbind(
  R1 = c(0, 1, 0),
  R2 = c(1, 0, 2),
  R3 = c(2, 1, 4),
  R4 = c(3, 0, 6)
)
query <- rbind(Q1 = c(1.5, 0.5, 3))
colnames(reference) <- colnames(query) <- c("M1|A", "M1|B", "M2|A")
blocks <- list(M1 = 1:2, M2 = 3L)
balanced <- .ablation_module_balanced_transform(reference, query, blocks)

z_reference <- sweep(
  sweep(reference, 2, balanced$center, "-"),
  2,
  balanced$scale,
  "/"
)
z_query <- sweep(
  sweep(query, 2, balanced$center, "-"),
  2,
  balanced$scale,
  "/"
)
manual_sq <- mean(c(
  mean((z_reference[1, blocks$M1] - z_query[1, blocks$M1])^2),
  mean((z_reference[1, blocks$M2] - z_query[1, blocks$M2])^2)
))
actual_sq <- sum((balanced$reference[1, ] - balanced$query[1, ])^2)
stopifnot(abs(manual_sq - actual_sq) < 1e-12)
stopifnot(identical(balanced$distance, "module_balanced_standardized_euclidean"))

# Should search only from query samples into the independent reference atlas.
reference_embedding <- rbind(
  A1 = c(0, 0),
  A2 = c(0.1, 0),
  B1 = c(10, 10),
  B2 = c(10.1, 10)
)
query_embedding <- rbind(QA = c(0.05, 0), QB = c(10.05, 10))
reference_metadata <- data.frame(
  sample_id = rownames(reference_embedding),
  cohort = c("RA1", "RA2", "RB1", "RB2"),
  cancer_type = c("A", "A", "B", "B"),
  assay_type = c("rna", "array", "rna", "array"),
  stringsAsFactors = FALSE
)
query_metadata <- data.frame(
  sample_id = rownames(query_embedding),
  cohort = c("QA", "QB"),
  cancer_type = c("A", "B"),
  assay_type = c("rna", "rna"),
  d1_provenance = "external_frozen",
  stringsAsFactors = FALSE
)
retrieval <- .ablation_query_reference_retrieval(
  reference = reference_embedding,
  query = query_embedding,
  reference_metadata = reference_metadata,
  query_metadata = query_metadata,
  label_column = "cancer_type",
  technical_columns = "assay_type",
  k = c(1L, 2L),
  search = "exact",
  seed = 301L
)
stopifnot(all(retrieval$per_sample$top1_label_match == 1))
stopifnot(all(retrieval$per_sample$mrr == 1))
stopifnot(!any(retrieval$neighbors$reference_cohort == retrieval$neighbors$query_cohort))
stopifnot(identical(sort(unique(retrieval$per_sample$k)), c(1L, 2L)))
stopifnot(all(c(
  "assay_type_match_rate",
  "assay_type_expected_rate",
  "assay_type_match_excess"
) %in% colnames(retrieval$per_sample)))
stopifnot(all(is.finite(retrieval$per_sample$assay_type_match_excess)))

# Should quantify approximate-search recall against exact query-to-reference kNN.
search_validation <- .ablation_validate_neighbor_search(
  reference = reference_embedding,
  query = query_embedding,
  reference_metadata = reference_metadata,
  query_metadata = query_metadata,
  label_column = "cancer_type",
  k = 2L,
  query_samples = 2L,
  n_trees = 100L,
  search_k = 1000L,
  seed = 302L
)
stopifnot(search_validation$recall == 1)
stopifnot(search_validation$query_sample_count == 2L)

# Should reject Null-Perm when the anchor is constant inside every query cohort.
constant_metadata <- data.frame(
  cohort = rep(c("C1", "C2"), each = 3),
  anchor = rep(c("A", "B"), each = 3),
  stringsAsFactors = FALSE
)
variable_metadata <- rbind(
  constant_metadata,
  data.frame(cohort = "C3", anchor = c("A", "B", "A"))
)
stopifnot(
  identical(
    .ablation_null_perm_eligibility(constant_metadata, "anchor")$status,
    "not_eligible"
  )
)
stopifnot(
  identical(
    .ablation_null_perm_eligibility(variable_metadata, "anchor")$status,
    "eligible"
  )
)

# Should reserve confirmatory evidence for independent anchors on external d1.
stopifnot(
  identical(
    .ablation_evidence_level(query_metadata, anchor_role = "independent")$level,
    "confirmatory"
  )
)
query_metadata$d1_provenance[1] <- "unknown"
stopifnot(
  identical(
    .ablation_evidence_level(query_metadata, anchor_role = "independent")$level,
    "descriptive"
  )
)
stopifnot(
  identical(
    .ablation_evidence_level(query_metadata, anchor_role = "bank_aligned")$level,
    "descriptive"
  )
)

# Should distribute label-bearing cohorts across finite grouped folds when feasible.
cohort <- paste0("C", seq_len(6))
anchor <- rep(c("A", "B"), each = 3)
fold <- .ablation_grouped_folds(
  cohort = cohort,
  n_folds = 3L,
  seed = 401L,
  label = anchor
)
coverage <- table(fold, anchor)
stopifnot(all(coverage > 0))
stopifnot(identical(fold, .ablation_grouped_folds(
  cohort = cohort,
  n_folds = 3L,
  seed = 401L,
  label = anchor
)))

message("11-test-ablation-representation: all tests passed")
