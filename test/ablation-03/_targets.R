# ablation-03 targets entry point
#
# targets owns pipeline provenance and crew owns target-level parallel workers.
# The scientific API remains in the installed CCS package; observability is
# written outside the targets store under the fixed cache root.

renv_activation <- file.path(getwd(), "renv", "activate.R")
if (!file.exists(renv_activation)) {
  stop(
    "ablation-03 requires project renv; run scripts/renv-ablation03.R --command init.",
    call. = FALSE
  )
}
source(renv_activation, local = .GlobalEnv)
if (!requireNamespace("renv", quietly = TRUE)) {
  stop("ablation-03 renv activation did not provide renv.", call. = FALSE)
}

if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Install the targets package before running ablation-03.", call. = FALSE)
}
if (!requireNamespace("crew", quietly = TRUE) ||
    !requireNamespace("autometric", quietly = TRUE)) {
  stop("Install crew and autometric in the ablation-03 renv project.", call. = FALSE)
}

targets::tar_source("targets")

cache_root <- .ablation03_cache_root()
store <- .ablation03_store(cache_root)
observability <- .ablation03_observability_config(cache_root)
controller <- .ablation03_crew_controller(observability)
targets::tar_config_set(store = store)
targets::tar_option_set(
  packages = c("CCS", "autometric", "crew", "digest", "rmarkdown", "yaml"),
  error = "stop",
  memory = "transient",
  controller = controller
)

# The main targets process gets its own autometric log. Worker-specific logs
# and CPU/RAM samples are configured above through crew_options_metrics().
if (targets::tar_active()) {
  autometric::log_start(
    path = observability$main_log,
    seconds = observability$metrics_interval
  )
}

list(
  targets::tar_target(runtime_config, .ablation03_runtime_config()),
  targets::tar_target(
    observability_config,
    .ablation03_observability_config(runtime_config$cache_root)
  ),
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
    resource_metrics,
    {
      run
      .ablation03_read_resource_metrics(
        list(observability = observability_config)
      )
    }
  ),
  targets::tar_target(worker_health, .ablation03_worker_health(resource_metrics)),
  targets::tar_target(
    report_payload,
    list(
      runtime = runtime_config,
      observability = observability_config,
      inputs = plan$context$analysis,
      result = result,
      native_geometry = native_geometry,
      retrieval = retrieval,
      readout = readout,
      learning_curve = learning_curve,
      scaling = cohort_scaling,
      decoder = decoder,
      biology = biology_result,
      structural = structural_result,
      resource_metrics = resource_metrics,
      worker_health = worker_health
    )
  )
)
