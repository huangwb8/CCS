local({
  source(file.path("test", "ablation-03", "targets", "functions.R"), local = TRUE)
  cache_root <- tempfile("ablation03-target-isolation-")
  dir.create(file.path(cache_root, "01-representations"), recursive = TRUE)
  on.exit(unlink(cache_root, recursive = TRUE), add = TRUE)
  computed_path <- file.path(cache_root, "01-representations", "structural-inputs.rds")
  computed <- list(marker = "representation_inputs")
  saveRDS(computed, computed_path)

  data_target <- list(input = list(
    biology_inputs = list(structural_anchor_cache = list(marker = "biology")),
    structural_inputs = list(structural_anchor_cache = list(marker = "structural"))
  ))
  .ablation03_materialize_optional_inputs(data_target, list(cache_root = cache_root))

  stopifnot(identical(readRDS(computed_path), computed))
  stopifnot(file.exists(file.path(cache_root, "01-biology", "structural-anchor-cache.rds")))
})

message("targets representation input isolation passed.")
