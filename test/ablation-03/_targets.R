# ablation-03: one targets graph for formal analysis and input-subset tests.
#
# Runs use independent input RDS files and cache roots. All scientific stages
# share the same targets graph, package API and parameter contract.

renv_activation <- file.path(getwd(), "renv", "activate.R")
if (!file.exists(renv_activation)) {
  stop("ablation-03 requires the project renv environment.", call. = FALSE)
}
source(renv_activation, local = .GlobalEnv)
if (!requireNamespace("renv", quietly = TRUE)) {
  stop("ablation-03 renv activation did not provide renv.", call. = FALSE)
}
if (!requireNamespace("targets", quietly = TRUE)) {
  stop("Install targets in the ablation-03 renv project.", call. = FALSE)
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

if (targets::tar_active()) {
  autometric::log_start(
    path = observability$main_log,
    seconds = observability$metrics_interval
  )
}

list(
  targets::tar_target(
    input_rds,
    Sys.getenv('CCS_ABLATION_INPUT_RDS'),
    format = 'file'
  ),
  targets::tar_target(runtime_config, {
    input_rds
    .ablation03_runtime_config()
  }),
  targets::tar_target(
    observability_config,
    .ablation03_observability_config(runtime_config$cache_root)
  ),
  targets::tar_target(
    data_preparation,
    {
      input_rds
      .ablation03_prepare_data_target(runtime_config)
    }
  ),
  targets::tar_target(
    representation_inputs,
    .ablation03_target_stage(
      "01.02.00. 表示输入准备.R",
      dependency = data_preparation,
      cache_root = runtime_config$cache_root
    )
  ),
  targets::tar_target(
    biology_inputs,
    .ablation03_biology_target(
      data_target = data_preparation,
      representation_target = representation_inputs,
      runtime_config = runtime_config
    )
  ),
  targets::tar_target(
    representation_analysis,
    .ablation03_target_stage(
      "02.01.00. 表示分析.R",
      dependency = representation_inputs,
      cache_root = runtime_config$cache_root
    )
  ),
  targets::tar_target(
    biology_sources,
    c("02.02.00. 生物锚点分析.R",
      "02.02.00. 生物锚点分析_functions.R"),
    format = "file"
  ),
  targets::tar_target(
    biology_analysis,
    {
      biology_sources
      result <- .ablation03_target_stage(
        "02.02.00. 生物锚点分析.R",
        dependency = list(representation_analysis, biology_inputs),
        cache_root = runtime_config$cache_root
      )
      report_inputs <- file.path(result$directory, c(
        "anchor_coverage.csv", "anchor_utility.csv",
        "anchor_contrasts.csv", "anchor_inference.csv",
        "anchor_cohort_deltas.csv", "anchor_missing_pairs.csv",
        "anchor_per_query_utility.rds",
        "ablation03-biology.rds"
      ))
      if (!all(file.exists(report_inputs))) {
        stop("Biology target did not produce all report inputs.", call. = FALSE)
      }
      result$report_input_md5 <- unname(tools::md5sum(report_inputs))
      result
    }
  ),
  targets::tar_target(
    structural_analysis,
    .ablation03_target_stage(
      "02.03.00. 结构复现分析.R",
      dependency = list(representation_analysis, biology_inputs),
      cache_root = runtime_config$cache_root
    )
  ),
  targets::tar_target(
    statistical_inference,
    .ablation03_statistical_inference(
      representation_inputs = representation_inputs,
      representation_analysis = representation_analysis,
      biology_analysis = biology_analysis,
      structural_analysis = structural_analysis,
      cache_root = runtime_config$cache_root
    ),
    format = "file"
  ),
  targets::tar_target(
    learning_query_inference,
    .ablation03_learning_query_inference(
      representation_inputs = representation_inputs,
      representation_analysis = representation_analysis,
      cache_root = runtime_config$cache_root
    ),
    format = "file"
  ),
  targets::tar_target(
    geometry_sensitivity,
    .ablation03_geometry_sensitivity(
      representation_inputs = representation_inputs,
      representation_analysis = representation_analysis,
      cache_root = runtime_config$cache_root
    ),
    format = "file"
  ),
  targets::tar_target(
    resource_metrics,
    {
      representation_analysis
      biology_analysis
      structural_analysis
      statistical_inference
      learning_query_inference
      geometry_sensitivity
      .ablation03_read_resource_metrics(list(observability = observability_config))
    }
  ),
  targets::tar_target(worker_health, .ablation03_worker_health(resource_metrics)),
  targets::tar_target(
    report_payload,
    list(
      runtime = runtime_config,
      observability = observability_config,
      data = data_preparation,
      representations = representation_inputs,
      biology_inputs = biology_inputs,
      representation = representation_analysis,
      biology = biology_analysis,
      structural = structural_analysis,
      statistical_inference = statistical_inference,
      learning_query_inference = learning_query_inference,
      geometry_sensitivity = geometry_sensitivity,
      resource_metrics = resource_metrics,
      worker_health = worker_health
    )
  ),
  targets::tar_target(
    data_overview_report_sources,
    c("01.04.00. 数据概览.Rmd",
      "01.01.00. 数据准备_functions.R",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R"),
    format = "file"
  ),
  targets::tar_target(
    data_overview_report,
    {
      data_overview_report_sources
      .ablation03_render_report(
        "01.04.00. 数据概览.Rmd",
        dependency = biology_inputs,
        cache_root = runtime_config$cache_root
      )
    },
    format = "file"
  ),
  targets::tar_target(
    representation_report_sources,
    c("02.01.00. 表示分析.Rmd",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R"),
    format = "file"
  ),
  targets::tar_target(
    representation_report,
    {
      representation_report_sources
      .ablation03_render_report(
        "02.01.00. 表示分析.Rmd",
        dependency = list(representation_analysis, statistical_inference,
          learning_query_inference, geometry_sensitivity),
        cache_root = runtime_config$cache_root
      )
    },
    format = "file"
  ),
  targets::tar_target(
    biology_report_sources,
    c("02.02.00. 生物锚点分析.Rmd",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R"),
    format = "file"
  ),
  targets::tar_target(
    biology_report,
    {
      biology_report_sources
      .ablation03_render_report(
        "02.02.00. 生物锚点分析.Rmd",
        dependency = list(biology_analysis, statistical_inference),
        cache_root = runtime_config$cache_root
      )
    },
    format = "file"
  ),
  targets::tar_target(
    structural_report_sources,
    c("02.03.00. 结构复现分析.Rmd",
      "02.01.00. 表示分析_functions.R",
      "02.03.00. 结构复现分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R"),
    format = "file"
  ),
  targets::tar_target(
    structural_report,
    {
      structural_report_sources
      .ablation03_render_report(
        "02.03.00. 结构复现分析.Rmd",
        dependency = list(structural_analysis, statistical_inference),
        cache_root = runtime_config$cache_root
      )
    },
    format = "file"
  )
)
