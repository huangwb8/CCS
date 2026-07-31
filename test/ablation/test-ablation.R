library(luckyBase)
Plus.library(c("CCS", "digest"))

source(file.path("test", "ablation", "ablation.R"))

make_test_fixture <- function() {
  set.seed(101)
  sample_ids <- sprintf("S%03d", seq_len(48))
  cohort <- rep(paste0("C", 1:4), each = 12)
  tissue <- rep(c("T1", "T1", "T2", "T2"), each = 12)
  biology <- rep(rep(c("B1", "B2"), each = 6), 4)

  expr <- matrix(
    stats::rnorm(4 * length(sample_ids)),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), sample_ids)
  )
  expr["g1", biology == "B1"] <- expr["g1", biology == "B1"] + 1
  expr["g4", biology == "B2"] <- expr["g4", biology == "B2"] + 1

  data <- list(
    T1 = list(
      C1 = list(expr = expr[, cohort == "C1", drop = FALSE], subtype = biology[cohort == "C1"]),
      C2 = list(expr = expr[, cohort == "C2", drop = FALSE], subtype = biology[cohort == "C2"])
    ),
    T2 = list(
      C3 = list(expr = expr[, cohort == "C3", drop = FALSE], subtype = biology[cohort == "C3"]),
      C4 = list(expr = expr[, cohort == "C4", drop = FALSE], subtype = biology[cohort == "C4"])
    )
  )

  tsp_features <- apply(utils::combn(rownames(expr), 2), 2, paste, collapse = ":")
  tsp <- vapply(strsplit(tsp_features, ":", fixed = TRUE), function(pair) {
    as.integer(expr[pair[1], ] >= expr[pair[2], ])
  }, integer(length(sample_ids)))
  rownames(tsp) <- sample_ids

  module_ids <- c("T1|M1", "T1|M2", "T1|M3", "T2|M4", "T2|M5", "T2|M6")
  d1 <- do.call(cbind, lapply(seq_along(module_ids), function(i) {
    score <- stats::plogis(tsp[, ((i - 1) %% ncol(tsp)) + 1] + stats::rnorm(nrow(tsp), 0, 0.1))
    block <- cbind(`1` = score, `2` = 1 - score)
    colnames(block) <- paste(module_ids[i], colnames(block), sep = "|")
    block
  }))
  rownames(d1) <- sample_ids

  make_model <- function() {
    class_model <- list(
      bst = NULL,
      breakVec = c(0, 0.5, 1),
      genes = c(rownames(expr), tsp_features)
    )
    list(
      Repeat = list(),
      Model = list(list(`1` = class_model, `2` = class_model))
    )
  }
  models <- list(
    T1 = setNames(rep(list(make_model()), 3), paste0("M", 1:3)),
    T2 = setNames(rep(list(make_model()), 3), paste0("M", 4:6))
  )

  object <- methods::new(
    "CCS",
    Repeat = list(
      method = "GSClassifier",
      geneSet = list(A = c("g1", "g2"), B = c("g3", "g4")),
      geneAnnotation = data.frame(ENSEMBL = rownames(expr)),
      geneid = "ensembl",
      seed = 101,
      model.dir = ""
    ),
    Model = models,
    Data = list(
      Probability = list(d1 = d1, d2 = matrix(0, 48, 4), d3 = matrix(0, 48, 2)),
      CCS = setNames(rep(1:2, length.out = 48), sample_ids),
      CancerType = setNames(tissue, sample_ids)
    )
  )

  metadata <- data.frame(
    sample_id = sample_ids,
    cohort = cohort,
    tissue = tissue,
    biology = biology,
    stringsAsFactors = FALSE
  )

  list(object = object, data = data, metadata = metadata, tsp_features = tsp_features)
}

fixture <- make_test_fixture()
fixture$metadata$pathway_label <- rep(c("P1", "P2"), length.out = nrow(fixture$metadata))

# Block manifest reflects the d1 naming convention and keeps complete blocks.
module_manifest <- .ablation_module_manifest(fixture$object)
stopifnot(
  nrow(module_manifest$modules) == 6,
  all(module_manifest$modules$block_width == 2),
  length(module_manifest$blocks) == 6
)

# Only pairwise TSP features actually consumed by frozen models enter Direct.
tsp_features <- .ablation_extract_tsp_features(fixture$object, module_manifest)
stopifnot(setequal(tsp_features, fixture$tsp_features))

# Input preparation aligns samples once and constructs a binary TSP matrix.
prepared <- .ablation_prepare_input(
  object = fixture$object,
  data = fixture$data,
  metadata = fixture$metadata,
  max_samples = Inf,
  seed = 303
)
stopifnot(
  identical(rownames(prepared$tsp), prepared$metadata$sample_id),
  identical(rownames(prepared$d1), prepared$metadata$sample_id),
  all(prepared$tsp %in% c(0L, 1L)),
  ncol(prepared$tsp) == length(fixture$tsp_features),
  identical(prepared$metadata$pathway_label, fixture$metadata$pathway_label)
)
probe_metadata <- prepared$metadata
probe_metadata$probe_label <- ifelse(
  probe_metadata$cohort == "C1",
  "single_cohort_label",
  probe_metadata$biology
)
eligible_probe <- .ablation_eligible_probe_labels(probe_metadata, "probe_label")
stopifnot(all(is.na(eligible_probe[probe_metadata$cohort == "C1"])))

# Grouped folds never split a cohort and are reproducible.
folds_a <- .ablation_grouped_folds(prepared$metadata$cohort, n_folds = 2, seed = 404)
folds_b <- .ablation_grouped_folds(prepared$metadata$cohort, n_folds = 2, seed = 404)
stopifnot(identical(folds_a, folds_b))
stopifnot(all(vapply(split(folds_a, prepared$metadata$cohort), function(x) {
  length(unique(x)) == 1
}, logical(1))))

# Null-Perm preserves every block's rows within cohort while breaking alignment.
permuted_a <- .ablation_permute_blocks(
  prepared$d1,
  prepared$module_manifest$blocks,
  prepared$metadata$cohort,
  seed = 505
)
permuted_b <- .ablation_permute_blocks(
  prepared$d1,
  prepared$module_manifest$blocks,
  prepared$metadata$cohort,
  seed = 505
)
stopifnot(identical(permuted_a, permuted_b))
for (block in prepared$module_manifest$blocks) {
  for (cohort_i in unique(prepared$metadata$cohort)) {
    rows <- prepared$metadata$cohort == cohort_i
    raw_rows <- apply(prepared$d1[rows, block, drop = FALSE], 1, paste, collapse = "|")
    perm_rows <- apply(permuted_a[rows, block, drop = FALSE], 1, paste, collapse = "|")
    stopifnot(identical(sort(unname(raw_rows)), sort(unname(perm_rows))))
  }
}
singleton_permutation <- .ablation_permute_blocks(
  prepared$d1[1:3, , drop = FALSE],
  prepared$module_manifest$blocks,
  cohort = c("shared", "singleton", "shared"),
  seed = 506
)
stopifnot(identical(
  singleton_permutation[2, , drop = FALSE],
  prepared$d1[2, , drop = FALSE]
))

# Null-RP and core geometry metrics are deterministic.
rp_a <- .ablation_random_projection(prepared$tsp, rank_q = 3, seed = 606)
rp_b <- .ablation_random_projection(prepared$tsp, rank_q = 3, seed = 606)
stopifnot(identical(rp_a, rp_b), identical(dim(rp_a), c(48L, 3L)))
stopifnot(abs(.ablation_linear_cka(rp_a, rp_a) - 1) < 1e-10)
stopifnot(abs(.ablation_knn_jaccard(rp_a, rp_a, k = 3) - 1) < 1e-10)

probe_train <- seq(1, nrow(prepared$tsp), by = 2)
probe_test <- setdiff(seq_len(nrow(prepared$tsp)), probe_train)
probe_metrics <- .ablation_probe(
  train = prepared$tsp[probe_train, , drop = FALSE],
  test = prepared$tsp[probe_test, , drop = FALSE],
  train_label = prepared$metadata$biology[probe_train],
  test_label = prepared$metadata$biology[probe_test],
  seed = 607,
  config = .ablation_merge_lists(
    .ablation_default_params(607),
    list(probe_nrounds = 5, numCores = 1)
  )
)
stopifnot(all(is.finite(probe_metrics)))

# The main function returns all four experiment-one groups and an audit trail.
result <- ablation(
  object = fixture$object,
  data = fixture$data,
  metadata = fixture$metadata,
  experiment = "cohort",
  output.dir = file.path("test", "ablation", "tmp", "synthetic"),
  params = list(
    rank = 3,
    k = 3,
    n_folds = 2,
    bootstrap = 10,
    rp_seeds = c(701, 702),
    permutation_seeds = c(801, 802),
    mechanism_samples = 24,
    probe = FALSE,
    cover = TRUE
  ),
  seed = 909,
  verbose = FALSE
)
stopifnot(inherits(result, "CCSAblation"))
stopifnot(setequal(
  unique(result$experiments$cohort$metrics$group_id),
  c("Direct", "Cohort", "Null-RP", "Null-Perm")
))
stopifnot(
  is.data.frame(result$experiments$cohort$contrasts),
  all(c("Cohort-Direct", "Cohort-Null-RP", "Cohort-Null-Perm") %in%
    unique(result$experiments$cohort$contrasts$contrast))
)
stopifnot(
  is.list(result$experiments$cohort$rank_sensitivity),
  is.data.frame(result$experiments$cohort$rank_sensitivity$metrics),
  nrow(result$experiments$cohort$rank_sensitivity$metrics) > 0
)
stopifnot(all(c(
  "experiment_id", "group_id", "sample_manifest_hash", "feature_manifest_hash",
  "fold_id", "rank_q", "k", "distance", "metric_name", "metric_value"
) %in% colnames(result$audit)))

gate_warning <- FALSE
invisible(withCallingHandlers(
  .ablation_gate_one(result$experiments$cohort, result$config$gate1),
  warning = function(condition) {
    gate_warning <<- TRUE
    invokeRestart("muffleWarning")
  }
))
stopifnot(!gate_warning)

# Scaling includes both the d1 saturation curve and representative d3/metaCCS
# stability; tissue-first returns paired one-stage/two-stage results.
layered_warnings <- character()
layered_result <- withCallingHandlers(ablation(
  object = fixture$object,
  data = fixture$data,
  metadata = fixture$metadata,
  experiment = c("scaling", "tissue_first"),
  output.dir = file.path("test", "ablation", "tmp", "synthetic-layered"),
  params = list(
    rank = 3,
    rank_sensitivity = 3,
    k = 3,
    n_folds = 2,
    bootstrap = 5,
    rp_seeds = 1001,
    permutation_seeds = 1002,
    mechanism_samples = 24,
    probe = FALSE,
    scaling_counts = c(2, 4, 6),
    scaling_sequences = 2,
    scaling_embedding_counts = c(2, 6),
    scaling_embedding_sequences = 1,
    scaling_embedding_seeds = c(1101, 1102),
    tissue_seeds = c(1201, 1202),
    fidelity_samples = 48,
    dr = list(
      method = "UWOT",
      dimension = c(3, 2),
      n_neighbors = 5,
      min_dist = 0.1,
      spread = 1,
      set_op_mix_ratio = 1
    ),
    cluster = list(eps = 0.5, minPts = 3),
    gate1 = list(enforce = FALSE),
    cover = TRUE
  ),
  seed = 1301,
  verbose = FALSE
), warning = function(condition) {
  layered_warnings <<- c(layered_warnings, conditionMessage(condition))
  invokeRestart("muffleWarning")
})
stopifnot(!any(grepl(
  "n_components > number of columns",
  layered_warnings,
  fixed = TRUE
)))
stopifnot(setequal(
  unique(layered_result$experiments$scaling$metrics$module_count),
  c(2, 4, 6)
))
stopifnot(
  is.data.frame(layered_result$experiments$scaling$summary$saturation),
  nrow(layered_result$experiments$scaling$summary$saturation) > 0
)
stopifnot(all(
  c("macro_auroc", "balanced_accuracy") %in%
    unique(layered_result$experiments$scaling$metrics$metric_name)
))
stopifnot(
  is.list(layered_result$experiments$scaling$embedding),
  is.data.frame(layered_result$experiments$scaling$embedding$metrics),
  nrow(layered_result$experiments$scaling$embedding$metrics) > 0
)
stopifnot(setequal(
  unique(layered_result$experiments$tissue_first$metrics$group_id),
  c("Two-stage", "One-stage")
))
stopifnot(
  "cluster_size_entropy" %in%
    unique(layered_result$experiments$tissue_first$metrics$metric_name),
  is.data.frame(layered_result$experiments$tissue_first$stratified),
  nrow(layered_result$experiments$tissue_first$stratified) > 0
)
stopifnot(
  is.data.frame(layered_result$experiments$tissue_first$contrasts),
  all(layered_result$experiments$tissue_first$contrasts$contrast ==
    "Two-stage-One-stage")
)
stopifnot(nrow(layered_result$experiments$tissue_first$stability) > 0)

message("All ablation tests passed.")
