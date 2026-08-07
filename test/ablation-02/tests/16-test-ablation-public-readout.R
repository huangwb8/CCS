#!/usr/bin/env Rscript

# Purpose: Verify public integration of external readout, learning curve, and Null-RP.
# Input: A bounded real ACC/STAD slice and three frozen cohort modules.
# Parameters: One lambda, one full-data repeat, and one random-projection seed.
# Output: Completed validation/control objects and persisted RDS artifacts.

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
  "public-readout-contract"
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
      max_reference_samples = 96L,
      max_query_samples = 40L
    ),
    anchors = list(
      primary = "cancer_type",
      primary_role = "independent",
      technical = c("assay_type", "source_system")
    ),
    geometry = list(
      k = 1L,
      search = "exact",
      geometry_samples = 80L,
      distance_pairs = 300L
    ),
    validation = list(
      enabled = TRUE,
      learning_fractions = 1,
      repeats = 1L,
      inner_folds = 2L,
      lambda = 1,
      nrounds = 5L,
      min_class_n = 5L,
      numCores = 1L
    ),
    controls = list(
      null_rp = TRUE,
      null_rp_rank = 2L,
      null_rp_seeds = 701L,
      null_perm = TRUE
    ),
    tradeoffs = list(
      decoder = TRUE,
      decoder_rank = 2L,
      decoder_lambda = 1
    ),
    output = list(cover = TRUE, cache_external_d1 = TRUE)
  ),
  seed = 702L,
  verbose = FALSE
)

stopifnot(result$readout$status == "complete")
stopifnot(result$learning_curve$status == "complete")
stopifnot(nrow(result$readout$overall) == 2L)
stopifnot(nrow(result$readout$paired_by_cohort) > 0L)
stopifnot(nrow(result$learning_curve$metrics) == 2L)
stopifnot(nrow(result$learning_curve$paired) == 1L)
stopifnot(result$controls$null_rp$status == "complete")
stopifnot(length(result$controls$null_rp$seeds) == 1L)
stopifnot(result$tradeoffs$decoder$status == "complete")
stopifnot(nrow(result$tradeoffs$decoder$summary) == 3L)
stopifnot(file.exists(file.path(output_dir, "readout.rds")))
stopifnot(file.exists(file.path(output_dir, "learning_curve.rds")))

saved_readout <- readRDS(file.path(output_dir, "readout.rds"))
saved_curve <- readRDS(file.path(output_dir, "learning_curve.rds"))
stopifnot(saved_readout$status == "complete")
stopifnot(saved_curve$status == "complete")

message("16-test-ablation-public-readout: all tests passed")
