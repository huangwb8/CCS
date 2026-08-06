#!/usr/bin/env Rscript

# Purpose: Verify that the real-data loader is self-contained and audit-ready.
# Input: The PADv20240810 expression RDS, filtered resCCS, and BatchInfo workbook.
# Parameters: Defaults declared by ablation-test-data.R; environment overrides are allowed.
# Output: Assertions covering sample identity, provenance, labels, and external cohorts.

source(file.path("test", "ablation-02", "ablation-test-data.R"))

required_objects <- c(
  "data_all",
  "resCCS",
  "batch_info",
  "cohort_summary",
  "ablation_metadata",
  "external_data",
  "reference_data",
  "filtered_cohorts"
)
stopifnot(all(vapply(required_objects, function(name) {
  exists(name, envir = .GlobalEnv, inherits = FALSE)
}, logical(1))))

required_columns <- c(
  "sample_id",
  "cohort",
  "tissue",
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
stopifnot(length(filtered_cohorts) == 43L)
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

message("test-ablation-data-contract: all tests passed")
