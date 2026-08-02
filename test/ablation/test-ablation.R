library(luckyBase)
Plus.library(c("CCS", "digest"))

source(file.path("R", "ablation.R"))

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
  direct_features <- c(rownames(expr), tsp_features, "s1s2")

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
      bst = list(feature_names = direct_features),
      breakVec = c(0, 0.5, 1),
      genes = direct_features
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

  list(
    object = object,
    data = data,
    metadata = metadata,
    expr = expr,
    direct_features = direct_features,
    tsp_features = tsp_features
  )
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

# Direct uses the complete frozen GSClassifier feature contract; the legacy
# TSP extractor remains restricted to ordinary gene-pair features.
direct_features <- .ablation_extract_direct_features(fixture$object, module_manifest)
tsp_features <- .ablation_extract_tsp_features(fixture$object, module_manifest)
stopifnot(
  setequal(direct_features, fixture$direct_features),
  setequal(tsp_features, fixture$tsp_features)
)

# Input preparation must reproduce GSClassifier's complete feature builder.
prepared <- .ablation_prepare_input(
  object = fixture$object,
  data = fixture$data,
  metadata = fixture$metadata,
  max_samples = Inf,
  seed = 303
)
expected_direct <- GSClassifier:::trainDataProc_X(
  Xmat = fixture$expr,
  geneSet = fixture$object@Repeat$geneSet,
  breakVec = c(0, 0.5, 1)
)$dat$Xbin[, fixture$direct_features, drop = FALSE]
matched_expr <- GSClassifier::geneMatch(
  fixture$expr,
  fixture$object@Repeat$geneAnnotation,
  fixture$object@Repeat$geneid,
  "fix"
)$Subset
expected_prediction_input <- GSClassifier:::dataProc(
  X = matched_expr,
  mods = fixture$object@Model$T1$M1$Model[[1]][[1]],
  geneSet = fixture$object@Repeat$geneSet
)
stopifnot(
  identical(rownames(prepared$direct), prepared$metadata$sample_id),
  identical(rownames(prepared$tsp), prepared$metadata$sample_id),
  identical(rownames(prepared$d1), prepared$metadata$sample_id),
  isTRUE(all.equal(
    prepared$direct[, fixture$direct_features, drop = FALSE],
    expected_direct[prepared$metadata$sample_id, , drop = FALSE],
    check.attributes = TRUE
  )),
  isTRUE(all.equal(
    prepared$direct[, colnames(expected_prediction_input), drop = FALSE],
    expected_prediction_input[prepared$metadata$sample_id, , drop = FALSE],
    check.attributes = TRUE
  )),
  all(prepared$tsp %in% c(0L, 1L)),
  ncol(prepared$direct) == length(fixture$direct_features),
  ncol(prepared$tsp) == length(fixture$tsp_features),
  identical(
    as.integer(factor(
      prepared$feature_manifest$feature_manifest$feature_type,
      levels = c("single_bin", "gene_pair", "set_pair")
    ) |> table()),
    c(4L, 6L, 1L)
  ),
  identical(prepared$metadata$pathway_label, fixture$metadata$pathway_label)
)

# Cohort matrices are joined by gene union, matching CCS::getResData().
uneven_data <- fixture$data
uneven_samples <- colnames(uneven_data$T1$C2$expr)
uneven_data$T1$C2$expr <- uneven_data$T1$C2$expr[-4, , drop = FALSE]
uneven <- .ablation_flatten_expression(uneven_data)
stopifnot(
  "g4" %in% rownames(uneven$expr),
  all(is.na(uneven$expr["g4", uneven_samples])),
  all(is.finite(uneven$expr["g4", colnames(fixture$data$T1$C1$expr)]))
)

# One Direct feature name must have one frozen binning definition.
heterogeneous <- fixture$object
heterogeneous@Model$T1$M1$Model[[1]][[1]]$breakVec <- c(0, 0.25, 0.5, 0.75, 1)
stopifnot(inherits(
  try(.ablation_frozen_feature_manifest(heterogeneous, module_manifest), silent = TRUE),
  "try-error"
))
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
rp_a <- .ablation_random_projection(prepared$direct, rank_q = 3, seed = 606)
rp_b <- .ablation_random_projection(prepared$direct, rank_q = 3, seed = 606)
stopifnot(identical(rp_a, rp_b), identical(dim(rp_a), c(48L, 3L)))
stopifnot(abs(.ablation_linear_cka(rp_a, rp_a) - 1) < 1e-10)
stopifnot(abs(.ablation_knn_jaccard(rp_a, rp_a, k = 3) - 1) < 1e-10)

probe_train <- seq(1, nrow(prepared$direct), by = 2)
probe_test <- setdiff(seq_len(nrow(prepared$direct)), probe_train)
probe_metrics <- .ablation_probe(
  train = prepared$direct[probe_train, , drop = FALSE],
  test = prepared$direct[probe_test, , drop = FALSE],
  train_label = prepared$metadata$biology[probe_train],
  test_label = prepared$metadata$biology[probe_test],
  seed = 607,
  config = .ablation_resolve_config(
    607,
    list(general = list(probe_nrounds = 5, numCores = 1))
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
    general = list(
      rank = 3,
      k = 3,
      n_folds = 2,
      bootstrap = 10,
      probe = FALSE,
      cover = TRUE
    ),
    cohort = list(
      rp_seeds = c(701, 702),
      permutation_seeds = c(801, 802),
      mechanism_samples = 24
    )
  ),
  seed = 909,
  verbose = FALSE
)
stopifnot(inherits(result, "CCSAblation"))
stopifnot(
  result$manifest$version == 2L,
  identical(result$manifest$gsclassifier_feature_builder, "trainDataProc_X"),
  result$manifest$direct_feature_count == length(fixture$direct_features),
  result$manifest$tsp_feature_count == length(fixture$tsp_features)
)
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
  .ablation_gate_one(result$experiments$cohort, result$config$scaling$gate),
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
    general = list(
      rank = 3,
      k = 3,
      n_folds = 2,
      bootstrap = 5,
      probe = FALSE,
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
      cover = TRUE
    ),
    cohort = list(
      rank_sensitivity = 3,
      rp_seeds = 1001,
      permutation_seeds = 1002,
      mechanism_samples = 24
    ),
    scaling = list(
      counts = c(2, 4, 6),
      sequences = 2,
      embedding_counts = c(2, 6),
      embedding_sequences = 1,
      embedding_seeds = c(1101, 1102),
      gate = list(enforce = FALSE)
    ),
    tissue_first = list(seeds = c(1201, 1202))
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

# Experiment 4 uses tissue-specific frozen GSClassifier support and keeps
# both arms paired through two-stage reduction and DBSCAN.
metaccs_fixture <- make_test_fixture()
metaccs_features <- list(
  T1 = c("g1", "g2", "g1:g2", "s1s2"),
  T2 = c("g3", "g4", "g3:g4", "s1s2")
)
for (tissue_i in names(metaccs_fixture$object@Model)) {
  for (cohort_i in names(metaccs_fixture$object@Model[[tissue_i]])) {
    model_i <- metaccs_fixture$object@Model[[tissue_i]][[cohort_i]]
    model_i$Model <- lapply(model_i$Model, function(repeat_i) {
      lapply(repeat_i, function(class_i) {
        class_i$genes <- metaccs_features[[tissue_i]]
        class_i$bst$feature_names <- metaccs_features[[tissue_i]]
        class_i
      })
    })
    metaccs_fixture$object@Model[[tissue_i]][[cohort_i]] <- model_i
  }
}

metaccs_prepared <- .ablation_prepare_input(
  object = metaccs_fixture$object,
  data = metaccs_fixture$data,
  metadata = metaccs_fixture$metadata,
  max_samples = Inf,
  seed = 1401
)
stopifnot(
  setequal(names(metaccs_prepared$feature_manifest$tissue_features), c("T1", "T2")),
  setequal(metaccs_prepared$feature_manifest$tissue_features$T1, metaccs_features$T1),
  setequal(metaccs_prepared$feature_manifest$tissue_features$T2, metaccs_features$T2),
  ncol(metaccs_prepared$direct) == 7,
  ncol(metaccs_prepared$tsp) == 2
)

direct_blocks <- .ablation_direct_tissue_blocks(
  metaccs_prepared$direct,
  metaccs_prepared$feature_manifest$tissue_features
)
stopifnot(
  identical(names(direct_blocks), c("T1", "T2")),
  identical(colnames(direct_blocks$T1), metaccs_features$T1),
  identical(colnames(direct_blocks$T2), metaccs_features$T2)
)

paired_a <- .ablation_paired_two_stage_embeddings(
  direct_blocks = direct_blocks,
  cohort_d1 = metaccs_prepared$d1,
  dr_config = list(
    method = "UWOT",
    dimension = c(3, 2),
    n_neighbors = 5,
    min_dist = 0.1,
    spread = 1,
    set_op_mix_ratio = 1,
    metric = "euclidean",
    n_threads = 1
  ),
  seed = 1501
)
changed_d1 <- metaccs_prepared$d1
changed_d1[, 1] <- rev(changed_d1[, 1])
paired_b <- .ablation_paired_two_stage_embeddings(
  direct_blocks = direct_blocks,
  cohort_d1 = changed_d1,
  dr_config = list(
    method = "UWOT",
    dimension = c(3, 2),
    n_neighbors = 5,
    min_dist = 0.1,
    spread = 1,
    set_op_mix_ratio = 1,
    metric = "euclidean",
    n_threads = 1
  ),
  seed = 1501
)
stopifnot(
  identical(
    paired_a$embeddings$`Direct-GSClassifier`,
    paired_b$embeddings$`Direct-GSClassifier`
  ),
  !identical(paired_a$embeddings$`Cohort-d1`, paired_b$embeddings$`Cohort-d1`),
  identical(
    rownames(paired_a$embeddings$`Direct-GSClassifier`),
    rownames(paired_a$embeddings$`Cohort-d1`)
  ),
  all(paired_a$parameter_audit$direct_dimension ==
    paired_a$parameter_audit$cohort_dimension),
  all(paired_a$parameter_audit$common_dimension <= 3)
)

# Cluster-biology metrics must be invariant to arbitrary raw cluster labels.
cluster_truth <- rep(c("B1", "B2"), each = 4)
cluster_a <- c(1, 1, 1, 0, 2, 2, 3, 3)
cluster_b <- c(7, 7, 7, 0, 4, 4, 9, 9)
biology_a <- .ablation_cluster_biology(cluster_a, cluster_truth)
biology_b <- .ablation_cluster_biology(cluster_b, cluster_truth)
stopifnot(isTRUE(all.equal(biology_a, biology_b, tolerance = 1e-12)))

metaccs_result <- ablation(
  object = metaccs_fixture$object,
  data = metaccs_fixture$data,
  metadata = metaccs_fixture$metadata,
  experiment = "metaccs",
  output.dir = file.path("test", "ablation", "tmp", "synthetic-metaccs"),
  params = list(
    general = list(
      rank = 3,
      k = 3,
      bootstrap = 5,
      fidelity_samples = 48,
      dr = list(
        method = "UWOT",
        dimension = c(3, 2),
        n_neighbors = 5,
        min_dist = 0.1,
        spread = 1,
        set_op_mix_ratio = 1,
        metric = "euclidean",
        n_threads = 1
      ),
      cluster = list(eps = 0.5, minPts = 3),
      cover = TRUE
    ),
    cohort = list(
      rp_seeds = 1601,
      permutation_seeds = 1602
    ),
    metaccs = list(
      resample_seeds = 1701,
      umap_seeds = c(1801, 1802),
      subsample_fraction = 1,
      parameter_mode = "fixed",
      direct_feature_mode = "tissue_model_union",
      retain_assignments = TRUE
    )
  ),
  seed = 1901,
  verbose = FALSE
)
metaccs <- metaccs_result$experiments$metaccs
stopifnot(
  setequal(unique(metaccs$metrics$group_id), c("Direct-GSClassifier", "Cohort-d1")),
  all(c(
    "cluster_biology_ari", "cluster_biology_nmi", "weighted_cluster_purity",
    "non_noise_coverage"
  ) %in% unique(metaccs$metrics$metric_name)),
  all(metaccs$contrasts$contrast == "Cohort-d1-Direct-GSClassifier"),
  nrow(metaccs$assignments) == 48 * 2 * 2,
  nrow(metaccs$cross_arm_agreement) == 2,
  nrow(metaccs$stability) == 2,
  is.data.frame(metaccs$stability_contrasts),
  nrow(metaccs$parameter_manifest) == 1,
  all(c(
    "resample_id", "umap_seed", "dr_param_set_id", "cluster_param_set_id",
    "direct_feature_hash", "input_hash"
  ) %in% colnames(metaccs$audit)),
  file.exists(file.path(
    "test", "ablation", "tmp", "synthetic-metaccs", "experiment-04-metaccs.rds"
  ))
)

# Grid mode expands the shared parameter budget without duplicating public arguments.
grid_config <- .ablation_merge_lists(
  .ablation_default_params(2001),
  list(
    metaccs = list(
      parameter_mode = "grid",
      dr_grid = list(
        list(n_neighbors = 4, min_dist = 0.1),
        list(n_neighbors = 5, min_dist = 0.2)
      ),
      cluster_grid = list(
        list(eps = 0.4, minPts = 3),
        list(eps = 0.5, minPts = 4)
      )
    )
  )
)
grid_manifest <- .ablation_metaccs_parameter_manifest(grid_config)
stopifnot(
  length(unique(grid_manifest$dr_param_set_id)) == 2,
  length(unique(grid_manifest$cluster_param_set_id)) == 2,
  nrow(grid_manifest) == 4
)
invalid_grid <- grid_config
invalid_grid$metaccs$dr_grid <- list(list(n_neighbors = 1))
stopifnot(inherits(
  try(.ablation_metaccs_parameter_manifest(invalid_grid), silent = TRUE),
  "try-error"
))
stopifnot(
  is.na(.ablation_ari(1, 1)),
  is.na(.ablation_cluster_jaccard(c(0, 0), c(0, 0)))
)

message("All ablation tests passed.")
