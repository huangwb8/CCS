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
targets::tar_source("R")
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
  targets::tar_target(ccs_code_files, .ablation03_ccs_code_files(), format = "file"),
  targets::tar_target(ccs_code_identity, .ablation03_assert_ccs_code(ccs_code_files)),
  targets::tar_target(
    representation_sources,
    c(
      "01.02.00. 表示输入准备.R",
      "02.01.00. 表示分析.R",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/workflow_helpers.R"
    ),
    format = "file"
  ),
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
      dependency = list(data_preparation, ccs_code_identity,
        representation_sources),
      cache_root = runtime_config$cache_root
    )
  ),
  targets::tar_target(
    biology_input_sources,
    c("config/biological-anchors.yml", "01.03.00. 生物输入准备.R",
      "02.02.00. 生物锚点分析_functions.R"), format = "file"
  ),
  targets::tar_target(
    biology_inputs,
    {
      biology_input_sources
      .ablation03_biology_target(
        data_target = data_preparation,
        representation_target = representation_inputs,
        runtime_config = runtime_config
      )
    }
  ),
  targets::tar_target(
    representation_analysis,
    .ablation03_target_stage(
      "02.01.00. 表示分析.R",
      dependency = list(representation_inputs, representation_sources),
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
    {
      result <- .ablation03_target_stage(
        "02.03.00. 结构复现分析.R",
        dependency = list(representation_analysis, biology_inputs),
        cache_root = runtime_config$cache_root
      )
      result$artifact_md5 <- unname(tools::md5sum(result$files))
      result
    }
  ),
  # Root-cause diagnostics preserve the original biology_analysis products.
  targets::tar_target(biology_diagnostic_config_file,
    "config/biology-diagnostics.yml", format = "file"),
  targets::tar_target(biology_diagnostic_config,
    yaml::read_yaml(biology_diagnostic_config_file)),
  targets::tar_target(biology_diagnostic_code_file,
    file.path("..", "..", "R", "ablation_biology.R"), format = "file"),
  targets::tar_target(biology_diagnostic_code_identity,
    .ablation03_bio_code_identity(biology_diagnostic_code_file)),
  targets::tar_target(biology_diagnostic_sources,
    .ablation03_bio_sources(biology_inputs), format = "file"),
  targets::tar_target(biology_diagnostic_representation_file,
    file.path(representation_inputs$directory, "representation-inputs.rds"), format = "file"),
  targets::tar_target(biology_diagnostic_baseline_files,
    c(file.path(biology_inputs$directory, "expression-anchor-cache.rds"),
      file.path(biology_analysis$directory, c("anchor_cohort_deltas.csv",
        "anchor_per_query_utility.rds", "ablation03-biology.rds")),
      file.path(runtime_config$cache_root, "ablation-experiment", "anchor_retrieval.rds")),
    format = "file"),
  targets::tar_target(biology_diagnostic_audit, {
    biology_diagnostic_code_identity
    ccs_code_identity
    biology_diagnostic_baseline_files
    .ablation03_bio_audit(biology_diagnostic_representation_file, biology_inputs, biology_analysis,
      data_preparation, biology_diagnostic_config, biology_diagnostic_sources, runtime_config)
  }, format = "file"),
  targets::tar_target(biology_diagnostic_gene_lists,
    .ablation03_bio_gene_lists(biology_diagnostic_audit, runtime_config), format = "file"),
  targets::tar_target(biology_diagnostic_scores,
    {
      biology_diagnostic_code_identity
      ccs_code_identity
      .ablation03_bio_scores(biology_diagnostic_audit, runtime_config)
    }, format = "file"),
  targets::tar_target(biology_diagnostic_retrieval_specs,
    .ablation03_bio_specs(biology_diagnostic_config), iteration = "list"),
  targets::tar_target(biology_diagnostic_retrieval,
    {
      biology_diagnostic_code_identity
      ccs_code_identity
      .ablation03_bio_retrieval(biology_diagnostic_representation_file, biology_diagnostic_audit,
        biology_diagnostic_config, biology_diagnostic_retrieval_specs, runtime_config)
    },
    pattern = map(biology_diagnostic_retrieval_specs), format = "file"),
  targets::tar_target(biology_diagnostic_readout,
    {
      biology_diagnostic_code_identity
      ccs_code_identity
      .ablation03_bio_readout(biology_diagnostic_representation_file, biology_diagnostic_audit,
        biology_diagnostic_scores, biology_diagnostic_config, runtime_config)
    }, format = "file"),
  targets::tar_target(biology_diagnostic_inference,
    {
    biology_diagnostic_code_identity
    ccs_code_identity
    biology_diagnostic_gene_lists
    .ablation03_bio_inference(biology_diagnostic_representation_file, biology_diagnostic_audit,
      biology_diagnostic_scores, biology_diagnostic_retrieval,
      biology_diagnostic_readout, biology_diagnostic_config, runtime_config)
    }, format = "file"),
  # Additional high-coverage validation leaves all broad/low-coverage targets intact.
  targets::tar_target(biology_high_coverage_config_file,
    "config/biology-high-coverage.yml", format = "file"),
  targets::tar_target(biology_high_coverage_config,
    yaml::read_yaml(biology_high_coverage_config_file)),
  targets::tar_target(biology_high_coverage_sources,
    c("R/biology_high_coverage.R", "targets/biology_high_coverage.R"), format = "file"),
  targets::tar_target(biology_high_coverage_inventory, {
    biology_high_coverage_sources
    biology_diagnostic_code_identity
    ccs_code_identity
    .ablation03_hc_inventory(biology_diagnostic_representation_file,
      biology_diagnostic_sources, biology_inputs, runtime_config)
  }, format = "file"),
  targets::tar_target(biology_high_coverage_measurement_frontier,
    .ablation03_hc_frontier(biology_high_coverage_inventory,
      biology_high_coverage_config, runtime_config), format = "file"),
  targets::tar_target(biology_high_coverage_contracts,
    .ablation03_hc_contracts(biology_high_coverage_measurement_frontier, runtime_config), format = "file"),
  targets::tar_target(biology_high_coverage_scores, {
    biology_high_coverage_sources
    biology_diagnostic_code_identity
    ccs_code_identity
    .ablation03_hc_scores(biology_high_coverage_contracts,
      biology_diagnostic_sources, runtime_config)
  }, format = "file"),
  targets::tar_target(biology_high_coverage_retrieval_specs,
    .ablation03_hc_specs(biology_high_coverage_contracts), iteration = "list"),
  targets::tar_target(biology_high_coverage_retrieval, {
    biology_high_coverage_sources
    biology_diagnostic_code_identity
    ccs_code_identity
    .ablation03_hc_retrieval(biology_diagnostic_representation_file,
      biology_high_coverage_contracts, biology_high_coverage_scores,
      biology_high_coverage_retrieval_specs, runtime_config)
  }, pattern = map(biology_high_coverage_retrieval_specs), format = "file"),
  targets::tar_target(biology_high_coverage_readout_specs,
    .ablation03_hc_readout_specs(biology_high_coverage_contracts), iteration = "list"),
  targets::tar_target(biology_high_coverage_readout, {
    biology_high_coverage_sources
    biology_diagnostic_code_identity
    ccs_code_identity
    .ablation03_hc_readout(biology_diagnostic_representation_file,
      biology_high_coverage_contracts, biology_high_coverage_scores,
      biology_high_coverage_readout_specs, runtime_config)
  }, pattern = map(biology_high_coverage_readout_specs), format = "file"),
  targets::tar_target(biology_high_coverage_sensitivity, {
    biology_high_coverage_sources
    biology_diagnostic_code_identity
    ccs_code_identity
    .ablation03_hc_sensitivity(biology_diagnostic_representation_file,
      biology_high_coverage_contracts, biology_high_coverage_scores,
      biology_high_coverage_retrieval, biology_high_coverage_readout_specs, runtime_config)
  }, pattern = map(biology_high_coverage_readout_specs), format = "file"),
  targets::tar_target(biology_high_coverage_inference, {
    biology_high_coverage_sources
    biology_diagnostic_code_identity
    ccs_code_identity
    .ablation03_hc_inference(biology_high_coverage_contracts,
      biology_high_coverage_scores, biology_high_coverage_retrieval,
      biology_high_coverage_readout, biology_high_coverage_sensitivity,
      biology_diagnostic_baseline_files, runtime_config)
  }, format = "file"),
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
      biology_diagnostic_inference
      biology_high_coverage_inference
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
      biology_diagnostics = biology_diagnostic_inference,
      biology_high_coverage = biology_high_coverage_inference,
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
      "00.Environment.R",
      "01.01.00. 数据准备_functions.R",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/anchor_plot_labels.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R",
      "templates/liquid_glass_theme.css",
      "templates/liquid_glass_lightbox.html"),
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
      "references/figure-s3-cohort-gain.csv",
      "references/figure-s3-provenance.csv",
      "references/figure-s3-importance.md",
      "00.Environment.R",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R",
      "templates/liquid_glass_theme.css",
      "templates/liquid_glass_lightbox.html"),
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
      "00.Environment.R",
      "02.01.00. 表示分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/anchor_plot_labels.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R",
      "templates/liquid_glass_theme.css",
      "templates/liquid_glass_lightbox.html",
      file.path(biology_analysis$directory, "ablation03-biology.rds"),
      file.path(runtime_config$cache_root, "ablation-experiment",
        c("sample-contract.rds", "anchor_retrieval.rds")),
      file.path(runtime_config$cache_root, "01-biology", "expression-anchor-cache.rds")),
    format = "file"
  ),
  targets::tar_target(
    biology_report,
    {
      biology_report_sources
      .ablation03_render_report(
        "02.02.00. 生物锚点分析.Rmd",
        dependency = list(biology_analysis, statistical_inference,
          biology_diagnostic_inference, biology_high_coverage_inference),
        cache_root = runtime_config$cache_root
      )
    },
    format = "file"
  ),
  targets::tar_target(
    structural_report_sources,
    c("02.03.00. 结构复现分析.Rmd",
      "00.Environment.R",
      "02.01.00. 表示分析_functions.R",
      "02.03.00. 结构复现分析_functions.R",
      "scripts/helpers/nature_colors.R",
      "scripts/helpers/nature_theme.R",
      "scripts/helpers/plot_delivery_helpers.R",
      "scripts/helpers/datatables_helper.R",
      "templates/liquid_glass_theme.css",
      "templates/liquid_glass_lightbox.html"),
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
