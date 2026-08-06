#!/usr/bin/env Rscript

# Purpose: Run the complete d1 vs Direct-GSClassifier ablation on real CCS data.
# Input: Objects prepared by ablation-test-data.R and functions sourced from R/ablation.R.
# Parameters: External filtered cohorts, paired learning curves, nulls, and decoder settings.
# Output: Full unfiltered result objects under tmp/ablation-experiment/.

source(file.path("test", "ablation-02", "00.Environment.R"))
source(file.path("test", "ablation-02", "ablation-test-data.R"))
source(file.path("R", "ablation.R"))

output_dir <- file.path("test", "ablation-02", "tmp", "ablation-experiment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Preserve the complete data audit before any display threshold is applied in Rmd.
saveRDS(
  ablation_data_profile,
  file.path(output_dir, "data-profile.rds")
)

ablation_result <- ablation(
  object = resCCS,
  data = data_all,
  metadata = ablation_metadata,
  experiment = "representation",
  output.dir = output_dir,
  params = list(
    comparison = list(
      module_ids = NULL,
      direct_group = "Direct-GSClassifier",
      cohort_group = "Cohort-d1"
    ),
    provenance = list(
      external_cohorts = filtered_cohorts,
      max_reference_samples = Inf,
      max_query_samples = Inf,
      require_external = TRUE
    ),
    anchors = list(
      primary = "cancer_type",
      primary_role = "independent",
      bank_aligned = "tissue",
      technical = c("assay_type", "platform_id", "source_system"),
      min_reference_cohorts = 2L
    ),
    geometry = list(
      k = c(5L, 15L, 30L),
      search = "annoy",
      n_trees = 75L,
      search_k = 10000L,
      exact_validation_queries = 30L,
      min_annoy_recall = 0.90,
      geometry_samples = 5000L,
      distance_pairs = 100000L
    ),
    validation = list(
      enabled = TRUE,
      learning_fractions = c(0.10, 0.25, 0.50, 1.00),
      repeats = 3L,
      inner_folds = 3L,
      lambda = c(0.1, 1, 10),
      nrounds = 30L,
      min_class_n = 20L,
      numCores = n_cores
    ),
    controls = list(
      null_rp = TRUE,
      null_rp_rank = 100L,
      null_rp_seeds = 20260805L + seq_len(3L),
      null_perm = TRUE
    ),
    tradeoffs = list(
      decoder = TRUE,
      decoder_rank = 50L,
      decoder_lambda = 1,
      decoder_max_reference_samples = 10000L,
      decoder_max_query_samples = 5000L
    ),
    output = list(
      cover = TRUE,
      cache_external_d1 = TRUE
    )
  ),
  seed = 20260805L,
  verbose = TRUE
)

saveRDS(
  ablation_result,
  file.path(output_dir, "ablation-result.rds")
)

luckyBase::LuckyVerbose(
  "ablation-experiment: complete; evidence level = ",
  ablation_result$evidence_level,
  "; output = ",
  normalizePath(output_dir, winslash = "/", mustWork = TRUE)
)
