#!/usr/bin/env Rscript

# Purpose: Run the complete d1 vs Direct-GSClassifier ablation on real CCS data.
# Input: Objects prepared by 01-ablation-test-data.R and functions from R/ablation.R.
# Parameters: External cohorts, paired learning/scaling curves, nulls, and decoder settings.
# Output: Full unfiltered result objects under the shared tmp/ablation-experiment/ cache.

source(file.path("test", "ablation-02", "00.Environment.R"))
source(file.path("test", "ablation-02", "01-ablation-test-data.R"))
source(file.path("test", "ablation-02", "02-ablation-experiment_functions.R"))
source(file.path("R", "ablation.R"))

output_dir <- file.path(
  "test", "ablation-02", "tmp", "ablation-experiment"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ablation_params <- .ae_ablation_params(filtered_cohorts, n_cores)

ablation_result <- ablation(
  object = resCCS,
  data = data_all,
  metadata = ablation_metadata,
  experiment = "representation",
  output.dir = output_dir,
  params = ablation_params,
  seed = 20260805L,
  verbose = TRUE
)

saveRDS(
  ablation_result,
  file.path(output_dir, "ablation-result.rds")
)

luckyBase::LuckyVerbose(
  "02-ablation-experiment: complete; evidence level = ",
  ablation_result$evidence_level,
  "; output = ",
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
)
