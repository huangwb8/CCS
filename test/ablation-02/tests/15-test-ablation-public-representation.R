#!/usr/bin/env Rscript

# Purpose: Verify the public representation-only ablation orchestration.
# Input: A bounded ACC/STAD slice containing real reference and filtered cohorts.
# Parameters: Three frozen modules, exact retrieval, and disabled readout for speed.
# Output: A CCSAblation object plus the seven planned representation result files.

source(file.path("test", "ablation-02", "01-ablation-test-data.R"))
source(file.path("R", "ablation.R"))

mini_data <- data_all[intersect(c("ACC", "STAD"), names(data_all))]
mini_metadata <- ablation_metadata[
  ablation_metadata$tissue %in% names(mini_data),
  ,
  drop = FALSE
]
mini_modules <- .ablation_module_manifest(resCCS)$modules$module_id[seq_len(3L)]
output_dir <- file.path(
  "test",
  "ablation-02",
  "tmp",
  "public-representation-contract"
)

result <- ablation(
  object = resCCS,
  data = mini_data,
  metadata = mini_metadata,
  experiment = "representation",
  output.dir = output_dir,
  params = list(
    comparison = list(module_ids = mini_modules),
    provenance = list(
      max_reference_samples = 80L,
      max_query_samples = 40L
    ),
    anchors = list(
      primary = "cancer_type",
      primary_role = "independent",
      technical = c("assay_type", "platform_id", "source_system")
    ),
    geometry = list(
      k = c(1L, 3L),
      search = "exact",
      geometry_samples = 80L,
      distance_pairs = 500L
    ),
    validation = list(enabled = FALSE),
    controls = list(null_rp = FALSE, null_perm = TRUE),
    output = list(cover = TRUE, cache_external_d1 = TRUE)
  ),
  seed = 501L,
  verbose = FALSE
)

stopifnot(inherits(result, "CCSAblation"))
stopifnot(identical(result$experiment, "representation"))
stopifnot(identical(result$evidence_level, "confirmatory"))
stopifnot(result$manifest$direct_feature_count == 529L)
stopifnot(result$manifest$tsp_feature_count == 496L)
stopifnot(result$manifest$query_sample_count == 40L)
stopifnot(result$manifest$reference_sample_count == 80L)
stopifnot(all(c("Direct-GSClassifier", "Cohort-d1") %in%
  unique(result$retrieval$per_sample$representation)))
stopifnot(all(result$retrieval$neighbors$query_cohort !=
  result$retrieval$neighbors$reference_cohort))
stopifnot(identical(result$controls$null_perm$status, "not_eligible"))
stopifnot(identical(result$readout$status, "not_run"))
stopifnot(identical(result$learning_curve$status, "not_run"))

expected_files <- c(
  "manifest.rds",
  "native_geometry.rds",
  "retrieval.rds",
  "readout.rds",
  "learning_curve.rds",
  "tradeoffs.rds",
  "audit.csv"
)
stopifnot(all(file.exists(file.path(output_dir, expected_files))))

message("15-test-ablation-public-representation: all tests passed")
