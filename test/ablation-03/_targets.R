# ablation-03 targets entry point
#
# The project is intentionally serial in this first migration step. Dynamic
# branches and crew/future are introduced only after this dependency graph has
# passed equivalence and resource review.

if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Install the targets package before running ablation-03.", call. = FALSE)
}

targets::tar_source("targets")

cache_root <- .ablation03_cache_root()
store <- .ablation03_store(cache_root)
targets::tar_config_set(store = store)
targets::tar_option_set(
  packages = c("CCS", "digest", "rmarkdown", "yaml"),
  error = "stop",
  memory = "transient"
)

list(
  targets::tar_target(runtime_config, .ablation03_runtime_config()),
  targets::tar_target(data_inputs, .ablation03_read_inputs(runtime_config)),
  targets::tar_target(
    representation_inputs,
    CCS::ablation_prepare_representation_inputs(
      object = data_inputs$object,
      data = data_inputs$data,
      metadata = data_inputs$metadata,
      params = .ablation03_representation_params(),
      seed = runtime_config$seed,
      cache_dir = NULL,
      verbose = FALSE
    )
  ),
  # One serial calculation target is the migration-safe bridge. The returned
  # list is split into named downstream targets so report changes do not force
  # re-preparation of raw inputs once the store is populated.
  targets::tar_target(
    representation_nodes,
    CCS::ablation_run_representation_nodes(
      representation_inputs,
      seed = runtime_config$seed,
      verbose = FALSE
    )
  ),
  targets::tar_target(native_geometry, representation_nodes$native_geometry),
  targets::tar_target(retrieval, representation_nodes$retrieval),
  targets::tar_target(readout, representation_nodes$readout),
  targets::tar_target(learning_curve, representation_nodes$learning_curve),
  targets::tar_target(cohort_scaling, representation_nodes$scaling),
  targets::tar_target(decoder, representation_nodes$decoder),
  targets::tar_target(
    learning_jobs,
    CCS::ablation_make_learning_curve_jobs(representation_inputs)
  ),
  targets::tar_target(
    learning_job_results,
    .ablation03_split_jobs(learning_jobs, learning_curve$metrics)
  ),
  targets::tar_target(
    learning_curve_combined,
    CCS::ablation_combine_learning_curve_jobs(learning_jobs, learning_job_results)
  ),
  targets::tar_target(scaling_jobs, CCS::ablation_make_scaling_jobs(representation_inputs)),
  targets::tar_target(
    biology_inputs,
    .ablation03_require_optional_input(data_inputs, "biology_inputs")
  ),
  targets::tar_target(
    biology_result,
    .ablation03_passthrough_result(biology_inputs, "biology_inputs")
  ),
  targets::tar_target(
    structural_inputs,
    .ablation03_require_optional_input(data_inputs, "structural_inputs")
  ),
  targets::tar_target(
    structural_result,
    .ablation03_passthrough_result(structural_inputs, "structural_inputs")
  ),
  targets::tar_target(
    report_payload,
    list(
      runtime = runtime_config,
      inputs = representation_inputs,
      native_geometry = native_geometry,
      retrieval = retrieval,
      readout = readout,
      learning_curve = learning_curve_combined,
      scaling = cohort_scaling,
      decoder = decoder,
      biology = biology_result,
      structural = structural_result
    )
  )
)

