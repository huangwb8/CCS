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
    context,
    CCS::ablation(
      object = data_inputs$object,
      data = data_inputs$data,
      metadata = data_inputs$metadata,
      params = .ablation03_representation_params(),
      seed = runtime_config$seed,
      step = "context",
      output.dir = file.path(runtime_config$cache_root, "output"),
      cache.root = runtime_config$cache_root,
      verbose = FALSE
    )
  ),
  targets::tar_target(plan, CCS::ablation(step = "plan", input = context, verbose = FALSE)),
  targets::tar_target(run, CCS::ablation(step = "run", input = plan, verbose = FALSE)),
  targets::tar_target(
    result,
    CCS::ablation(
      step = "result",
      input = run,
      output.dir = file.path(runtime_config$cache_root, "output"),
      params = .ablation03_representation_params(),
      verbose = FALSE
    )
  ),
  targets::tar_target(native_geometry, run$value$native_geometry),
  targets::tar_target(retrieval, run$value$retrieval),
  targets::tar_target(readout, run$value$readout),
  targets::tar_target(learning_curve, run$value$learning_curve),
  targets::tar_target(cohort_scaling, run$value$cohort_scaling),
  targets::tar_target(decoder, run$value$tradeoffs$decoder),
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
      inputs = plan$context$analysis,
      result = result,
      native_geometry = native_geometry,
      retrieval = retrieval,
      readout = readout,
      learning_curve = learning_curve,
      scaling = cohort_scaling,
      decoder = decoder,
      biology = biology_result,
      structural = structural_result
    )
  )
)
