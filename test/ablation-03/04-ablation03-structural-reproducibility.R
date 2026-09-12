#!/usr/bin/env Rscript
# Purpose: Evaluate cross-cohort reproducibility of biological-state geometry.
# Input: Frozen Direct/d1 inputs plus the existing expression-anchor cache.
# Parameters: Within-cohort tail states, minimum entity size, and node bootstrap.
# Output: Complete entity audits, cohort-pair similarities, matrices, and summary.

options(stringsAsFactors = FALSE)
env_candidates <- c(
  file.path(getwd(), "00.Environment.R"),
  file.path(getwd(), "test", "ablation-03", "00.Environment.R")
)
env_path <- env_candidates[file.exists(env_candidates)][1L]
if (is.na(env_path)) {
  stop("ablation-03: run from the project directory or repository root.",
    call. = FALSE
  )
}
source(env_path)
source(.ablation03_path("01-ablation03-test-data.R"))
source(.ablation03_path("02-ablation03-experiment_functions.R"))
source(.ablation03_repo_path("R", "ablation.R"))
source(.ablation03_path(
  "04-ablation03-structural-reproducibility_functions.R"
))

# Step 1: Freeze the analysis contract before inspecting structural results.
seed <- 20260912L
tail_fraction <- 1 / 3
min_entity_n <- 8L
min_shared_entities <- 8L
n_boot <- 2000L
output_dir <- .ablation03_path("tmp", "ablation-structural-reproducibility")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- readRDS(.ablation03_path(
  "tmp", "ablation-experiment", "manifest.rds"
))
retrieval <- readRDS(.ablation03_path(
  "tmp", "ablation-experiment", "retrieval.rds"
))
cache <- readRDS(.ablation03_path(
  "tmp", "ablation-biology", "expression-anchor-cache.rds"
))
if (!identical(cache$schema_version, 1L) ||
    !identical(cache$status, "complete")) {
  stop("structural reproducibility: unsupported biology cache schema.",
    call. = FALSE
  )
}
cache_neighbours <- retrieval$neighbors[
  retrieval$neighbors$neighbor_rank <= 15L,
  ,
  drop = FALSE
]
target_ids <- sort(unique(c(
  as.character(cache_neighbours$query_sample),
  as.character(cache_neighbours$reference_sample)
)))
expected_sample_hash <- digest::digest(
  paste(target_ids, collapse = "\n"),
  algo = "md5",
  serialize = FALSE
)
if (!identical(cache$sample_key_hash, expected_sample_hash)) {
  stop("structural reproducibility: biology cache sample hash mismatch.",
    call. = FALSE
  )
}

# Step 2: Reconstruct the exact transformed query representations used by 03.
ablation_params <- .ae_ablation_params(filtered_cohorts, n_cores)
config <- .ablation_resolve_representation_config(seed, ablation_params)
analysis_external_cohorts <- intersect(
  manifest$external_cohorts,
  names(cache$cohorts)
)
analysis_query_data <- partition_data(
  data_all,
  analysis_external_cohorts,
  keep = TRUE
)
analysis_data <- reference_data
for (tissue in names(analysis_query_data)) {
  for (cohort in names(analysis_query_data[[tissue]])) {
    analysis_data[[tissue]][[cohort]] <- analysis_query_data[[tissue]][[cohort]]
  }
}
analysis <- .ablation_prepare_representation_analysis(
  object = resCCS_ablation,
  data = analysis_data,
  metadata = ablation_metadata,
  config = config,
  output.dir = output_dir,
  seed = seed,
  verbose = TRUE
)
prepared <- analysis$prepared
direct_scaled <- .ablation_scale_train_apply(
  prepared$reference_direct,
  prepared$query_direct
)
d1_scaled <- .ablation_module_balanced_transform(
  prepared$reference_d1,
  prepared$query_d1,
  prepared$selected_blocks
)
representations <- list(
  `Direct-GSClassifier` = direct_scaled$test,
  `Cohort-d1` = d1_scaled$query
)

# Step 3: Define shared biological entities independently of either arm.
scores <- .asr_compute_anchor_scores(cache, manifest$external_cohorts)
states <- .asr_assign_anchor_states(
  scores,
  tail_fraction = tail_fraction,
  min_entity_n = min_entity_n
)
query_metadata <- prepared$query_metadata
query_metadata$cohort_key <- paste(
  query_metadata$tissue,
  query_metadata$cohort,
  sep = "/"
)
query_metadata$cancer_type <- as.character(query_metadata$biology)
query_ids <- rownames(representations[[1L]])
anchor_names <- sort(unique(states$anchor))
score_coverage <- vapply(anchor_names, function(anchor) {
  all(query_ids %in% scores$sample_id[scores$anchor == anchor])
}, logical(1L))
if (!all(score_coverage)) {
  stop(
    "structural reproducibility: anchor states do not cover every query sample: ",
    paste(anchor_names[!score_coverage], collapse = ", "),
    call. = FALSE
  )
}

direct_geometry <- .asr_build_cohort_geometries(
  representations[["Direct-GSClassifier"]],
  query_metadata,
  states
)
d1_geometry <- .asr_build_cohort_geometries(
  representations[["Cohort-d1"]],
  query_metadata,
  states
)
cohort_metadata <- unique(query_metadata[, c(
  "cohort_key", "cancer_type", "assay_type", "source_system"
)])

# Step 4: Compare matched distance-matrix entries for every valid cohort pair.
pair_comparisons <- .asr_compare_cohort_pairs(
  direct_geometry$geometries,
  d1_geometry$geometries,
  cohort_metadata,
  min_shared_entities = min_shared_entities
)
summary_statistics <- .asr_summarize_pairs(
  pair_comparisons,
  n_boot = n_boot,
  seed = seed
)
cohorts <- sort(unique(c(
  pair_comparisons$cohort_a,
  pair_comparisons$cohort_b
)))
similarity_matrices <- list(
  `Direct-GSClassifier` = .asr_similarity_matrix(
    pair_comparisons,
    "direct_similarity",
    cohorts
  ),
  `Cohort-d1` = .asr_similarity_matrix(
    pair_comparisons,
    "d1_similarity",
    cohorts
  )
)

entity_audit <- unique(states[, c(
  "cohort_key", "anchor", "state", "entity", "cohort_sample_count",
  "entity_sample_count", "gene_count"
)])
entity_audit <- entity_audit[order(
  entity_audit$cohort_key,
  entity_audit$anchor,
  entity_audit$state
), ]

# Step 5: Save complete numeric products; interpretation remains in the Rmd.
utils::write.csv(
  entity_audit,
  file.path(output_dir, "structural_entity_audit.csv"),
  row.names = FALSE
)
utils::write.csv(
  pair_comparisons,
  file.path(output_dir, "structural_pair_comparisons.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_statistics,
  file.path(output_dir, "structural_summary.csv"),
  row.names = FALSE
)
saveRDS(
  similarity_matrices,
  file.path(output_dir, "structural_similarity_matrices.rds")
)
saveRDS(
  list(
    schema_version = 1L,
    method = list(
      biological_entity = "within_cohort_anchor_low_high_tail",
      anchors = anchor_names,
      tail_fraction = tail_fraction,
      min_entity_n = min_entity_n,
      distance = "euclidean_in_existing_transformed_representation",
      structural_similarity = "spearman_upper_triangle",
      min_shared_entities = min_shared_entities,
      uncertainty = "cohort_node_bootstrap",
      bootstrap = n_boot,
      seed = seed
    ),
    provenance = list(
      representation_cache_key = prepared$cache_key,
      biology_cache_source_md5 = cache$source$md5,
      biology_cache_sample_key_hash = cache$sample_key_hash
    ),
    scores = scores,
    states = states,
    entity_audit = entity_audit,
    direct_prototypes = direct_geometry$prototypes,
    d1_prototypes = d1_geometry$prototypes,
    direct_geometries = direct_geometry$geometries,
    d1_geometries = d1_geometry$geometries,
    pair_comparisons = pair_comparisons,
    summary = summary_statistics,
    similarity_matrices = similarity_matrices
  ),
  file.path(output_dir, "ablation03-structural-reproducibility.rds")
)

cat(sprintf(
  paste0(
    "structural reproducibility: cohorts=%d pairs=%d entities=%d ",
    "mean_delta=%+.4f output=%s\n"
  ),
  length(cohorts),
  nrow(pair_comparisons),
  length(unique(states$entity)),
  summary_statistics$mean_delta[
    summary_statistics$scope == "all_cohort_pairs"
  ],
  output_dir
))
