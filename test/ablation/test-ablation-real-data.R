library(luckyBase)
Plus.library(c("CCS", "digest"))

paramMD5 <- "5ff3a2de76e6cf902e765e8224f9cb66"
source(file.path("test", "ablation", "ablation-test-data.R"))
source(file.path("R", "ablation.R"))

smoke_params <- list(
  general = list(
    rank = 5,
    k = 5,
    n_folds = 2,
    bootstrap = 20,
    max_samples = 400,
    probe = FALSE,
    cover = TRUE
  ),
  cohort = list(
    rank_sensitivity = 5,
    geometry_samples = 200,
    distance_pairs = 1000,
    mechanism_samples = 100,
    rp_seeds = 1701,
    permutation_seeds = 1801
  )
)
output_dir <- file.path("test", "ablation", "tmp", "real-data-smoke")

result_a <- ablation(
  object = resCCS,
  data = data_all,
  experiment = "cohort",
  output.dir = output_dir,
  params = smoke_params,
  seed = 1901,
  verbose = FALSE
)
result_b <- ablation(
  object = resCCS,
  data = data_all,
  experiment = "cohort",
  output.dir = output_dir,
  params = smoke_params,
  seed = 1901,
  verbose = FALSE
)

stopifnot(
  result_a$manifest$module_count == 150,
  result_a$manifest$tsp_feature_count == 496,
  result_a$manifest$direct_feature_count == 529,
  identical(
    as.integer(result_a$manifest$direct_feature_type_count[
      c("gene_pair", "set_pair", "single_bin")
    ]),
    c(496L, 1L, 32L)
  ),
  result_a$manifest$sample_count == 400,
  identical(result_a$experiments$cohort$metrics, result_b$experiments$cohort$metrics),
  identical(
    result_a$experiments$cohort$dimension_free_geometry,
    result_b$experiments$cohort$dimension_free_geometry
  )
)

message("Real CCS ablation smoke test passed.")
