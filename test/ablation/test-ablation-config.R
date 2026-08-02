source(file.path("R", "ablation.R"))

# Canonical nested parameters are kept in human-reviewable domains.
config <- .ablation_resolve_config(
  seed = 101,
  params = list(
    general = list(rank = 7, dr = list(min_dist = 0.2)),
    cohort = list(rp_seeds = 201),
    scaling = list(counts = c(2, 4), gate = list(enforce = TRUE)),
    tissue_first = list(seeds = 301),
    metaccs = list(resample_seeds = 401, umap_seeds = 501)
  )
)
stopifnot(
  identical(names(config), c(
    "general", "cohort", "scaling", "tissue_first", "metaccs"
  )),
  identical(config$general$rank, 7),
  identical(config$general$dr$min_dist, 0.2),
  identical(config$general$dr$method, "UWOT"),
  identical(config$cohort$rp_seeds, 201),
  identical(config$scaling$counts, c(2, 4)),
  isTRUE(config$scaling$gate$enforce),
  identical(config$tissue_first$seeds, 301)
)

# Existing flat calls normalize to the same canonical schema.
legacy <- .ablation_resolve_config(
  seed = 101,
  params = list(
    rank = 7,
    dr = list(min_dist = 0.2),
    rp_seeds = 201,
    scaling_counts = c(2, 4),
    gate1 = list(enforce = TRUE),
    tissue_seeds = 301,
    metaccs = list(resample_seeds = 401, umap_seeds = 501)
  )
)
stopifnot(identical(config, legacy))

# Ambiguous mixed schemas and unknown fields fail before an experiment starts.
conflict <- tryCatch(
  .ablation_resolve_config(101, list(rank = 3, general = list(rank = 7))),
  error = identity
)
unknown <- tryCatch(
  .ablation_resolve_config(101, list(cohort = list(typo_seed = 1))),
  error = identity
)
stopifnot(
  inherits(conflict, "error"),
  grepl("provided more than once", conditionMessage(conflict), fixed = TRUE),
  inherits(unknown, "error"),
  grepl("unknown params field", conditionMessage(unknown), fixed = TRUE)
)

message("Ablation config schema tests passed.")
