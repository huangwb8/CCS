# Evaluate the original reciprocal structural design using prepared inputs only.
options(stringsAsFactors = FALSE, device = function(...) grDevices::pdf(file = NULL))
bootstrap <- c(file.path("scripts", "helpers", "workflow_helpers.R"),
  file.path("test", "ablation-03", "scripts", "helpers", "workflow_helpers.R"))
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
source(.ablation03_path("02.03.00. 结构复现分析_functions.R"))
# Step 1: Freeze the reciprocal validation contract before inspecting results.
seed <- 20260912L
tail_fraction <- 1 / 3
min_entity_n <- 8L
min_shared_entities <- 8L
n_boot <- 2000L
min_entity_n <- suppressWarnings(as.integer(Sys.getenv(
  "CCS_ABLATION_STRUCTURAL_MIN_ENTITY_N", unset = as.character(min_entity_n)
)))
n_boot <- suppressWarnings(as.integer(Sys.getenv(
  "CCS_ABLATION_STRUCTURAL_BOOTSTRAP", unset = as.character(n_boot)
)))
if (length(min_entity_n) != 1L || is.na(min_entity_n) || min_entity_n < 1L) {
  stop("CCS_ABLATION_STRUCTURAL_MIN_ENTITY_N must be a positive integer.", call. = FALSE)
}
if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 0L) {
  stop("CCS_ABLATION_STRUCTURAL_BOOTSTRAP must be a non-negative integer.", call. = FALSE)
}
matched_repeats <- suppressWarnings(as.integer(Sys.getenv(
  "CCS_ABLATION_MATCHED_REPEATS", unset = "20"
)))
if (length(matched_repeats) != 1L || is.na(matched_repeats) || matched_repeats < 1L) {
  stop("CCS_ABLATION_MATCHED_REPEATS must be a positive integer.", call. = FALSE)
}
stage_parameters <- list(
  seed = seed,
  tail_fraction = tail_fraction,
  min_entity_n = min_entity_n,
  min_shared_entities = min_shared_entities,
  n_boot = n_boot,
  matched_repeats = matched_repeats
)
if (.wf_cache_hit("ablation-structural-reproducibility", stage_parameters)) {
  quit(save = "no", status = 0L)
}
.wf_begin("ablation-structural-reproducibility", stage_parameters)
output_dir <- .wf_output("ablation-structural-reproducibility")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bundle <- .wf_read("01-representations", "structural-inputs.rds")
prepared <- bundle$analysis$prepared
full_d1 <- bundle$structural$full_d1
representation_bundle <- .wf_read("01-representations", "representation-inputs.rds")
manifest <- .wf_read("ablation-experiment", "manifest.rds")
biology_cache <- .wf_read("01-biology", "expression-anchor-cache.rds")
anchor_cache <- .wf_read("01-biology", "structural-anchor-cache.rds")
anchor_cache_key <- anchor_cache$cache_key
.asr_assert_representation_contract(
  structural_prepared = prepared,
  representation_prepared = representation_bundle$analysis$prepared,
  representation_manifest = manifest
)
full_module_manifest <- bundle$structural$full_module_manifest
module_table <- .asr_resolve_module_table(
  full_module_manifest,
  bundle$structural$tissue_resolution_audit,
  manifest$external_cohorts
)
reference_module_ids <- module_table$module_id[module_table$bank_role == "reference"]
external_module_ids <- module_table$module_id[module_table$bank_role == "external"]
if (length(intersect(reference_module_ids, external_module_ids)) > 0L ||
    length(c(reference_module_ids, external_module_ids)) != nrow(module_table)) {
  stop("structural reproducibility: module banks are not a disjoint partition.", call. = FALSE)
}

reference_ids <- Reduce(intersect, list(
  rownames(prepared$reference_direct), rownames(full_d1)
))
external_ids <- Reduce(intersect, list(
  rownames(prepared$query_direct), rownames(full_d1)
))
if (length(intersect(reference_ids, external_ids)) > 0L) {
  stop("structural reproducibility: reference and external samples overlap.", call. = FALSE)
}

# Step 4: Fit every scale on the module-bank side and score only its target side.
forward_base <- .asr_prepare_direction_base(
  prepared$reference_direct,
  prepared$query_direct,
  prepared$reference_metadata,
  prepared$query_metadata,
  anchor_cache,
  reference_ids,
  external_ids,
  tail_fraction,
  min_entity_n
)
reverse_base <- .asr_prepare_direction_base(
  prepared$query_direct,
  prepared$reference_direct,
  prepared$query_metadata,
  prepared$reference_metadata,
  anchor_cache,
  external_ids,
  reference_ids,
  tail_fraction,
  min_entity_n
)

forward_bank <- .asr_subset_bank(
  full_d1, full_d1, full_module_manifest$blocks, reference_module_ids
)
reverse_bank <- .asr_subset_bank(
  full_d1, full_d1, full_module_manifest$blocks, external_module_ids
)
forward_result <- .asr_evaluate_direction(
  forward_base,
  forward_bank$fit,
  forward_bank$target,
  forward_bank$blocks,
  "reference_bank_to_external_targets",
  "reference",
  "external",
  min_shared_entities,
  n_boot,
  seed
)
reverse_result <- .asr_evaluate_direction(
  reverse_base,
  reverse_bank$fit,
  reverse_bank$target,
  reverse_bank$blocks,
  "external_bank_to_reference_targets",
  "external",
  "reference",
  min_shared_entities,
  n_boot,
  seed + 100L
)

directional_pairs <- rbind(
  forward_result$pair_comparisons,
  reverse_result$pair_comparisons
)
directional_summary <- rbind(forward_result$summary, reverse_result$summary)
directional_matrices <- list(
  reference_bank_to_external_targets = forward_result$similarity_matrices,
  external_bank_to_reference_targets = reverse_result$similarity_matrices
)
directional_entity_audit <- rbind(
  transform(
    unique(forward_base$states[, c(
      "cohort_key", "anchor", "state", "entity", "cohort_sample_count",
      "entity_sample_count", "gene_count"
    )]),
    direction = "reference_bank_to_external_targets"
  ),
  transform(
    unique(reverse_base$states[, c(
      "cohort_key", "anchor", "state", "entity", "cohort_sample_count",
      "entity_sample_count", "gene_count"
    )]),
    direction = "external_bank_to_reference_targets"
  )
)
directional_entity_audit <- directional_entity_audit[, c(
  "direction", setdiff(names(directional_entity_audit), "direction")
)]

# Step 5: Repeat a shared-tissue, equal-size design to expose bank sensitivity.
matched_design <- .asr_matched_bank_design(
  module_table,
  n_repeats = matched_repeats,
  seed = seed + 1000L
)
matched_rows <- vector("list", matched_repeats * 2L)
matched_index <- 1L
for (repeat_id in seq_len(matched_repeats)) {
  for (bank_role in c("reference", "external")) {
    selected_ids <- matched_design$module_id[
      matched_design$repeat_id == repeat_id &
        matched_design$bank_role == bank_role
    ]
    bank <- .asr_subset_bank(
      full_d1, full_d1, full_module_manifest$blocks, selected_ids
    )
    if (bank_role == "reference") {
      result <- .asr_evaluate_direction(
        forward_base, bank$fit, bank$target, bank$blocks,
        "reference_bank_to_external_targets", "reference", "external",
        min_shared_entities, 0L, seed + 2000L + repeat_id
      )
    } else {
      result <- .asr_evaluate_direction(
        reverse_base, bank$fit, bank$target, bank$blocks,
        "external_bank_to_reference_targets", "external", "reference",
        min_shared_entities, 0L, seed + 3000L + repeat_id
      )
    }
    result$summary$repeat_id <- repeat_id
    result$summary$sensitivity_type <- "matched_bank_design"
    matched_rows[[matched_index]] <- result$summary
    matched_index <- matched_index + 1L
  }
}
matched_repeats_summary <- do.call(rbind, matched_rows)
matched_groups <- split(
  matched_repeats_summary,
  interaction(
    matched_repeats_summary$direction,
    matched_repeats_summary$scope,
    drop = TRUE,
    lex.order = TRUE
  )
)
matched_summary <- do.call(rbind, lapply(matched_groups, function(part) {
  data.frame(
    direction = part$direction[1L],
    bank_role = part$bank_role[1L],
    target_role = part$target_role[1L],
    scope = part$scope[1L],
    repeat_count = length(unique(part$repeat_id)),
    bank_module_count = unique(part$bank_module_count)[1L],
    mean_delta_median = stats::median(part$mean_delta),
    mean_delta_min = min(part$mean_delta),
    mean_delta_max = max(part$mean_delta),
    mean_direct_similarity_median = stats::median(part$mean_direct_similarity),
    mean_d1_similarity_median = stats::median(part$mean_d1_similarity),
    fraction_d1_ge_direct_median = stats::median(part$fraction_d1_ge_direct),
    range_interpretation = "design_sensitivity_not_confidence_interval",
    stringsAsFactors = FALSE
  )
}))
rownames(matched_summary) <- NULL

sample_audit <- data.frame(
  direction = c(
    "reference_bank_to_external_targets",
    "external_bank_to_reference_targets"
  ),
  bank_role = c("reference", "external"),
  target_role = c("external", "reference"),
  fit_sample_count = c(
    length(forward_base$fit_sample_ids), length(reverse_base$fit_sample_ids)
  ),
  target_sample_count = c(
    length(forward_base$target_sample_ids), length(reverse_base$target_sample_ids)
  ),
  fit_target_overlap = c(
    length(intersect(forward_base$fit_sample_ids, forward_base$target_sample_ids)),
    length(intersect(reverse_base$fit_sample_ids, reverse_base$target_sample_ids))
  ),
  stringsAsFactors = FALSE
)

# Step 6: Save complete products; the established files remain the forward arm.
forward_entity_audit <- directional_entity_audit[
  directional_entity_audit$direction ==
    "reference_bank_to_external_targets", -1L, drop = FALSE
]
utils::write.csv(
  forward_entity_audit,
  file.path(output_dir, "structural_entity_audit.csv"),
  row.names = FALSE
)
utils::write.csv(
  forward_result$pair_comparisons,
  file.path(output_dir, "structural_pair_comparisons.csv"),
  row.names = FALSE
)
utils::write.csv(
  forward_result$summary,
  file.path(output_dir, "structural_summary.csv"),
  row.names = FALSE
)
saveRDS(
  forward_result$similarity_matrices,
  file.path(output_dir, "structural_similarity_matrices.rds")
)
utils::write.csv(
  directional_entity_audit,
  file.path(output_dir, "structural_directional_entity_audit.csv"),
  row.names = FALSE
)
utils::write.csv(
  directional_pairs,
  file.path(output_dir, "structural_directional_pair_comparisons.csv"),
  row.names = FALSE
)
utils::write.csv(
  directional_summary,
  file.path(output_dir, "structural_directional_summary.csv"),
  row.names = FALSE
)
saveRDS(
  directional_matrices,
  file.path(output_dir, "structural_directional_similarity_matrices.rds")
)
utils::write.csv(
  module_table,
  file.path(output_dir, "structural_module_bank_audit.csv"),
  row.names = FALSE
)
utils::write.csv(
  sample_audit,
  file.path(output_dir, "structural_directional_sample_audit.csv"),
  row.names = FALSE
)
utils::write.csv(
  matched_design,
  file.path(output_dir, "structural_matched_bank_design.csv"),
  row.names = FALSE
)
utils::write.csv(
  matched_repeats_summary,
  file.path(output_dir, "structural_matched_bank_repeats.csv"),
  row.names = FALSE
)
utils::write.csv(
  matched_summary,
  file.path(output_dir, "structural_matched_bank_summary.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    schema_version = 2L,
    method = list(
      biological_entity = "within_target_cohort_anchor_low_high_tail",
      anchors = names(biology_cache$anchors),
      tail_fraction = tail_fraction,
      min_entity_n = min_entity_n,
      distance = "euclidean_in_fixed_bank_scaled_representation",
      structural_similarity = "spearman_upper_triangle",
      min_shared_entities = min_shared_entities,
      uncertainty = "cohort_node_bootstrap",
      bootstrap = n_boot,
      matched_repeats = matched_repeats,
      matched_range = "design_sensitivity_not_confidence_interval",
      seed = seed
    ),
    provenance = list(
      representation_cache_key = prepared$cache_key,
      biology_cache_source_md5 = biology_cache$source$md5,
      structural_anchor_cache_key = anchor_cache_key,
      module_manifest_hash = digest::digest(module_table, algo = "md5"),
      full_d1_hash = digest::digest(full_d1, algo = "md5")
    ),
    module_table = module_table,
    sample_audit = sample_audit,
    directions = list(
      reference_bank_to_external_targets = list(
        scores = forward_base$scores,
        states = forward_base$states,
        result = forward_result
      ),
      external_bank_to_reference_targets = list(
        scores = reverse_base$scores,
        states = reverse_base$states,
        result = reverse_result
      )
    ),
    directional_entity_audit = directional_entity_audit,
    directional_pair_comparisons = directional_pairs,
    directional_summary = directional_summary,
    directional_similarity_matrices = directional_matrices,
    matched = list(
      design = matched_design,
      repeats = matched_repeats_summary,
      summary = matched_summary
    ),
    # Schema-v1 compatibility fields: the established forward arm.
    scores = forward_base$scores,
    states = forward_base$states,
    entity_audit = forward_entity_audit,
    direct_prototypes = forward_result$direct_prototypes,
    d1_prototypes = forward_result$d1_prototypes,
    pair_comparisons = forward_result$pair_comparisons,
    summary = forward_result$summary,
    similarity_matrices = forward_result$similarity_matrices
  ),
  file.path(output_dir, "ablation03-structural-reproducibility.rds")
)

forward_all <- forward_result$summary[
  forward_result$summary$scope == "all_cohort_pairs", , drop = FALSE
]
.wf_receipt("ablation-structural-reproducibility", "02.03.00. 结构复现分析",
  inputs = c(.wf_output("ablation-experiment", "stage-receipt.rds"),
    .wf_output("01-biology", "stage-receipt.rds"),
    .wf_output("01-representations", "stage-receipt.rds"),
    .wf_path("scripts", "helpers", "workflow_helpers.R"),
    .wf_path("02.03.00. 结构复现分析.R"),
    .ablation03_path("02.03.00. 结构复现分析_functions.R"),
    .ablation03_ccs_description),
  outputs = list.files(output_dir, pattern = "^(structural_|ablation03-structural).*\\.(csv|rds)$", full.names = TRUE))
reverse_all <- reverse_result$summary[
  reverse_result$summary$scope == "all_cohort_pairs", , drop = FALSE
]
cat(sprintf(
  paste0(
    "structural reproducibility: reference_modules=%d external_modules=%d ",
    "forward_pairs=%d forward_delta=%+.4f reverse_pairs=%d ",
    "reverse_delta=%+.4f matched=%d output=%s\n"
  ),
  length(reference_module_ids),
  length(external_module_ids),
  forward_all$cohort_pair_count,
  forward_all$mean_delta,
  reverse_all$cohort_pair_count,
  reverse_all$mean_delta,
  unique(matched_summary$bank_module_count),
  output_dir
))
