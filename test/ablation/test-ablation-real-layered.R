library(luckyBase)
Plus.library(c("CCS", "digest"))

paramMD5 <- "5ff3a2de76e6cf902e765e8224f9cb66"
source(file.path("test", "ablation", "ablation-test-data.R"))
source(file.path("test", "ablation", "ablation.R"))

result <- ablation(
  object = resCCS,
  data = data_all,
  experiment = c("scaling", "tissue_first"),
  output.dir = file.path("test", "ablation", "tmp", "real-layered-smoke"),
  params = list(
    rank = 5,
    rank_sensitivity = 5,
    k = 5,
    n_folds = 2,
    bootstrap = 5,
    max_samples = 120,
    geometry_samples = 80,
    distance_pairs = 500,
    mechanism_samples = 60,
    rp_seeds = 2701,
    permutation_seeds = 2801,
    probe = FALSE,
    scaling_counts = 25,
    scaling_sequences = 1,
    scaling_embedding_counts = 25,
    scaling_embedding_sequences = 1,
    scaling_embedding_seeds = 2901,
    scaling_subsample_fraction = 1,
    tissue_seeds = 3001,
    tissue_subsample_fraction = 1,
    fidelity_samples = 120,
    gate1 = list(enforce = FALSE),
    cover = TRUE
  ),
  seed = 3101,
  verbose = FALSE
)

stopifnot(
  result$manifest$module_count == 150,
  result$manifest$tsp_feature_count == 496,
  nrow(result$experiments$scaling$embedding$metrics) > 0,
  setequal(
    unique(result$experiments$tissue_first$metrics$group_id),
    c("Two-stage", "One-stage")
  )
)

message("Real CCS layered ablation smoke test passed.")
