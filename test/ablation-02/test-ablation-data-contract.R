#!/usr/bin/env Rscript

# Purpose: Verify that the real-data loader is self-contained and audit-ready.
# Input: The PADv20240810 expression RDS, filtered resCCS, and BatchInfo workbook.
# Parameters: Defaults declared by ablation-test-data.R; environment overrides are allowed.
# Output: Assertions covering sample identity, provenance, labels, and external cohorts.

source(file.path("test", "ablation-02", "ablation-test-data.R"))

required_objects <- c(
  "data_all",
  "raw_cohort_index",
  "resolved_cohort_index",
  "tissue_mapping",
  "tissue_resolution_audit",
  "ablation_data_profile",
  "resCCS",
  "batch_info",
  "cohort_summary",
  "ablation_metadata",
  "external_data",
  "reference_data",
  "filtered_cohorts",
  "filtered_model_cohorts"
)
stopifnot(all(vapply(required_objects, function(name) {
  exists(name, envir = .GlobalEnv, inherits = FALSE)
}, logical(1))))

required_columns <- c(
  "sample_id",
  "cohort",
  "tissue",
  "bank_tissue",
  "bank_cohort_key",
  "cancer_type",
  "assay_type",
  "platform_id",
  "source_system",
  "duplicate_sample_id_global",
  "analysis_set",
  "d1_provenance",
  "anchor_role"
)
stopifnot(all(required_columns %in% colnames(ablation_metadata)))
stopifnot(nrow(batch_info) == 39112L)
stopifnot(nrow(cohort_summary) == 193L)
stopifnot(nrow(raw_cohort_index) == 193L)
stopifnot(nrow(resolved_cohort_index) == 193L)
stopifnot(sum(raw_cohort_index$sample_count) == 39112L)
stopifnot(sum(resolved_cohort_index$sample_count) == 39112L)
stopifnot(sum(raw_cohort_index$tissue == "Undefined") == 45L)
stopifnot(!"Undefined" %in% names(data_all))
stopifnot(!any(ablation_metadata$tissue == "Undefined"))
stopifnot(nrow(tissue_resolution_audit) == 45L)
stopifnot(sum(tissue_resolution_audit$sample_count) == 8626L)
stopifnot(setequal(
  paste(tissue_resolution_audit$resolved_tissue, tissue_resolution_audit$cohort, sep = "/"),
  paste(tissue_mapping$tissue, tissue_mapping$cohort_id, sep = "/")
))
resolved_mapped_index <- resolved_cohort_index[
  resolved_cohort_index$cohort %in% tissue_mapping$cohort_id,
  ,
  drop = FALSE
]
stopifnot(setequal(
  resolved_mapped_index$cohort_key,
  paste(tissue_mapping$tissue, tissue_mapping$cohort_id, sep = "/")
))
stopifnot(length(filtered_cohorts) == 43L)
stopifnot(length(filtered_model_cohorts) == 43L)
stopifnot(sum(grepl("^Undefined/", filtered_model_cohorts)) == 24L)
stopifnot(!any(grepl("^Undefined/", filtered_cohorts)))
stopifnot(sum(ablation_metadata$analysis_set == "external_query") == 14083L)
stopifnot(sum(ablation_metadata$analysis_set == "reference_atlas") == 25029L)
stopifnot(all(
  ablation_metadata$d1_provenance[
    ablation_metadata$analysis_set == "external_query"
  ] == "external_frozen"
))
stopifnot(all(
  ablation_metadata$d1_provenance[
    ablation_metadata$analysis_set == "reference_atlas"
  ] == "in_sample"
))
stopifnot(all(ablation_metadata$anchor_role == "independent"))

d1_ids <- rownames(resCCS@Data$Probability$d1)
external_ids <- ablation_metadata$sample_id[
  ablation_metadata$analysis_set == "external_query"
]
stopifnot(length(intersect(d1_ids, external_ids)) == 0L)
stopifnot(sum(ablation_metadata$duplicate_sample_id_global) == 542L)
stopifnot(length(unique(
  ablation_metadata$sample_id[ablation_metadata$duplicate_sample_id_global]
)) == 271L)

external_keys <- unlist(lapply(names(external_data), function(tissue) {
  paste(tissue, names(external_data[[tissue]]), sep = "/")
}), use.names = FALSE)
reference_keys <- unlist(lapply(names(reference_data), function(tissue) {
  paste(tissue, names(reference_data[[tissue]]), sep = "/")
}), use.names = FALSE)
stopifnot(setequal(external_keys, filtered_cohorts))
stopifnot(length(intersect(external_keys, reference_keys)) == 0L)
stopifnot(setequal(
  unique(ablation_metadata$bank_cohort_key[
    ablation_metadata$analysis_set == "external_query"
  ]),
  filtered_model_cohorts
))
stopifnot(setequal(
  ablation_metadata$sample_id[
    ablation_metadata$cohort_key %in% filtered_cohorts
  ],
  ablation_metadata$sample_id[
    ablation_metadata$bank_cohort_key %in% filtered_model_cohorts
  ]
))

message("test-ablation-data-contract: all tests passed")
