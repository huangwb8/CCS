#' Run layered ablation experiments for a CCS model
#'
#' @description
#' Evaluate the frozen CCS representation without retraining cohort submodels.
#' The default `"representation"` experiment reconstructs the complete native
#' GSClassifier input, re-encodes filtered cohorts with the frozen model bank,
#' and compares Direct-GSClassifier with Cohort-d1 on independent
#' query-to-reference retrieval, grouped linear readout, paired learning curves,
#' null controls, and feature-type reconstruction. Explicit requests for the
#' legacy scaling, tissue-first, and metaCCS experiments remain available during
#' the API transition. All stochastic operations are paired across comparison
#' groups and recorded in an audit table.
#'
#' @param object A `CCS` object containing the frozen d1 representation.
#' @param data Raw RNA expression data. A tissue/cohort nested list is preferred;
#'   each leaf can be an expression matrix or a list containing `expr`.
#' @param metadata Optional sample annotation with sample, cohort, tissue and
#'   biological-label columns. Common CCS column names are recognized.
#' @param experiment `"representation"` by default. The deprecated single value
#'   `"cohort"` maps to `"representation"` with a warning. Explicit legacy
#'   requests may contain one or more of `"cohort"`, `"scaling"`,
#'   `"tissue_first"`, and `"metaccs"`; `"representation"` cannot be mixed with
#'   legacy experiments.
#' @param output.dir Independent output directory. Existing CCS products are
#'   never overwritten.
#' @param params Named nested list. For `"representation"`, values are merged
#'   onto `.ablation_representation_default_params(seed)` under `comparison`,
#'   `provenance`, `anchors`, `geometry`, `validation`, `controls`, `tradeoffs`,
#'   and `output`. For explicit legacy experiments, values are merged onto
#'   `.ablation_default_params(seed)` under `general`, `cohort`, `scaling`,
#'   `tissue_first`, and `metaccs`; existing flat legacy names remain accepted.
#'   Unknown fields are rejected before computation.
#'
#'   Representation-specific cohort-bank scaling settings:
#'
#'   \describe{
#'     \item{`scaling$enabled`}{Logical; when `TRUE`, evaluate tissue breadth,
#'       within-tissue cohort depth, and matched-size frozen cohort-model banks.
#'       Default: `FALSE`.}
#'     \item{`scaling$module_counts`}{Positive module counts used for matched-size
#'       breadth-heavy versus depth-heavy contrasts. Infeasible sizes are
#'       retained in the exclusion audit. Default:
#'       `c(10L, 25L, 50L, 75L, 100L, 125L, 150L)`.}
#'     \item{`scaling$sequences`}{Number of independently randomized,
#'       reproducible bank-design repeats. Each repeat, rather than its samples
#'       or grid cells, is the uncertainty unit. Default: `10L`.}
#'     \item{`scaling$direct_feature_type`}{Main Direct-GSClassifier diagnostic
#'       contract: `"all"` for the complete native input or `"gene_pair"` for
#'       ordinary TSPs only. Default: `"all"`.}
#'     \item{`scaling$sensitivity_feature_type`}{Optional secondary Direct
#'       contract. `"gene_pair"` reports the explicitly named
#'       Direct-GSClassifier-TSP sensitivity analysis; `"none"` disables it.
#'       Default: `"gene_pair"`.}
#'     \item{`scaling$biology_anchors`}{Independent sample-level biological or
#'       clinical metadata columns used for external-utility neighborhood
#'       consistency. Missing or empty anchors are returned as `not_evaluated`;
#'       cancer type is not substituted. Default: `character()`.}
#'     \item{`scaling$score_reference_samples`}{Fixed stratified reference
#'       sample cap shared by every bank score. This controls repeated neighbor
#'       and readout cost without changing bank-repeat uncertainty. `Inf` keeps
#'       all reference samples. Default: `5000L`.}
#'     \item{`scaling$score_query_samples`}{Fixed stratified external-query cap
#'       shared by every bank score. `Inf` keeps all eligible queries. Default:
#'       `2000L`.}
#'     \item{`scaling$lambda`}{Single pre-specified L2 penalty shared by Direct
#'       and matched-bank lineage diagnostics. Fixing it avoids re-tuning at
#'       every bank composition. Default: `1`.}
#'     \item{`scaling$bootstrap`}{Positive number of repeat-level bootstrap
#'       draws used for breadth/depth slopes, interactions, and matched-size
#'       summaries. Default:
#'       `1000L`.}
#'   }
#'
#'   Shared representation, validation, and metric settings:
#'
#'   \describe{
#'     \item{`general$rank`}{Positive target PCA rank used by Experiment 1 and the d1
#'       scaling curves. It is reduced automatically when a fold has fewer
#'       samples or features. Default: `50L`.}
#'     \item{`cohort$rank_sensitivity`}{Positive integer vector of additional PCA ranks
#'       evaluated in Experiment 1; infeasible values are reduced to the common
#'       maximum rank. Default: `c(25L, 50L, 100L)`.}
#'     \item{`general$k`}{Positive neighborhood size used by kNN similarity, mixing,
#'       purity, fidelity, and stability metrics. It is capped at the available
#'       sample count minus one. Default: `30L`.}
#'     \item{`general$distance`}{Distance label stored in metric and audit outputs.
#'       Current computations remain Euclidean regardless of this value; this is
#'       not yet a distance-method switch. Default: `"euclidean"`.}
#'     \item{`general$n_folds`}{Number of whole-cohort validation folds shared by
#'       Experiments 1 and 2. `Inf` creates one fold per cohort; a finite value is
#'       capped at the number of cohorts and must yield at least two folds.
#'       Default: `Inf`.}
#'     \item{`general$bootstrap`}{Positive number of bootstrap draws used for metric
#'       summaries and paired 95 percent confidence intervals. Default: `1000L`.}
#'     \item{`general$max_samples`}{Maximum number of aligned input samples retained
#'       before all experiments. `Inf` keeps every sample; finite limits use
#'       tissue-by-cohort stratified sampling. Default: `Inf`.}
#'     \item{`cohort$geometry_samples`}{Maximum stratified sample count used for the
#'       dimension-free Direct-GSClassifier-versus-d1 geometry comparison in Experiment 1.
#'       Default: `5000L`.}
#'     \item{`cohort$distance_pairs`}{Maximum number of sampled row pairs used to
#'       estimate the Spearman correlation between distance rankings. Default:
#'       `100000L`.}
#'     \item{`cohort$mechanism_samples`}{Maximum stratified sample count used by the
#'       selective-reconstruction mechanism metrics in Experiment 1. Default:
#'       `1000L`.}
#'     \item{`general$probe`}{Logical; whether Experiments 1 and 2 fit the cross-cohort
#'       linear XGBoost probe and report macro AUROC and balanced accuracy.
#'       Default: `TRUE`.}
#'     \item{`general$probe_label`}{Metadata column decoded by the probe. Labels not
#'       represented in at least two cohorts are excluded. Default: `"tissue"`.}
#'     \item{`general$probe_nrounds`}{Positive XGBoost boosting-round count for the
#'       linear probe. Default: `50L`.}
#'     \item{`general$numCores`}{Positive thread count passed to the XGBoost probe.
#'       Default: `1L`.}
#'   }
#'
#'   Experiment 1 null-control settings:
#'
#'   \describe{
#'     \item{`cohort$rp_density`}{Probability that an entry of the sparse Achlioptas
#'       Null-RP projection matrix is nonzero; it should lie in `(0, 1]`.
#'       Default: `1 / 3`.}
#'     \item{`cohort$rp_seeds`}{Non-empty vector of seeds for independent Null-RP
#'       repeats. Default: `seed + seq_len(20)`.}
#'     \item{`cohort$permutation_seeds`}{Non-empty vector of seeds for independent
#'       within-cohort, whole-module Null-Perm repeats. Default:
#'       `seed + 1000L + seq_len(20)`.}
#'   }
#'
#'   Experiment 2 cohort-axis scaling settings:
#'
#'   \describe{
#'     \item{`scaling$counts`}{Positive module counts evaluated on each nested
#'       sequence; values above the available module bank are capped. Default:
#'       `c(10L, 25L, 50L, 75L, 100L, 125L, 150L)`.}
#'     \item{`scaling$sequences`}{Positive number of tissue-balanced nested
#'       module sequences used for the d1 scaling curves. Default: `100L`.}
#'     \item{`scaling$embedding_counts`}{Module counts at which the more
#'       expensive downstream two-stage embedding and DBSCAN analysis is run;
#'       unavailable counts are capped. Default: `c(25L, 50L, 100L, 150L)`.}
#'     \item{`scaling$embedding_sequences`}{Number of nested sequences, starting
#'       from the generated sequence list, retained for downstream embedding.
#'       Default: `10L`.}
#'     \item{`scaling$embedding_seeds`}{Seeds for paired sample subsampling,
#'       two-stage reduction, clustering, and stability repeats at each selected
#'       module count. Default: `seed + 2000L + seq_len(10)`.}
#'     \item{`scaling$subsample_fraction`}{Fraction in `(0, 1]` sampled within
#'       tissue-by-cohort strata for each downstream scaling repeat. Default:
#'       `0.8`.}
#'   }
#'
#'   Experiment 3 tissue-first settings:
#'
#'   \describe{
#'     \item{`tissue_first$seeds`}{Seeds defining paired stratified subsamples and
#'       reduction repeats shared by the Two-stage and One-stage arms. Default:
#'       `seed + 3000L + seq_len(20)`.}
#'     \item{`tissue_first$subsample_fraction`}{Fraction in `(0, 1]` sampled within
#'       tissue-by-cohort strata for each Experiment 3 repeat. Default: `0.8`.}
#'     \item{`general$fidelity_samples`}{Maximum stratified sample count used when
#'       computing embedding trustworthiness and continuity. Default: `2000L`.}
#'   }
#'
#'   Shared dimensional-reduction settings are supplied as
#'   `general = list(dr = list(...))`:
#'
#'   \describe{
#'     \item{`general$dr$method`}{Reduction method forwarded to the CCS reduction
#'       helpers. Default: `"UWOT"`.}
#'     \item{`general$dr$dimension`}{Two positive target dimensions: the first is the
#'       per-tissue d2 dimension and the second is the final global d3 dimension.
#'       Low-rank blocks use the largest feasible smaller value. Default:
#'       `c(5L, 2L)`.}
#'     \item{`general$dr$n_neighbors`}{Neighborhood size forwarded to UWOT and capped
#'       for small blocks; it must be at least two when supplied. Default: `30L`.}
#'     \item{`general$dr$min_dist`}{UWOT minimum-distance parameter. Default: `0.01`.}
#'     \item{`general$dr$spread`}{UWOT spread parameter. Default: `0.75`.}
#'     \item{`general$dr$set_op_mix_ratio`}{UWOT fuzzy-set intersection/union mixing
#'       parameter. Default: `1`.}
#'     \item{`general$dr$metric`}{Distance metric forwarded to the reduction backend.
#'       Default: `"euclidean"`.}
#'     \item{`general$dr$n_threads`}{Positive reduction thread count, or `NULL` to use
#'       the CCS/backend default. Default: `NULL`.}
#'   }
#'
#'   Shared DBSCAN settings are supplied as `general = list(cluster = list(...))`:
#'
#'   \describe{
#'     \item{`general$cluster$eps`}{Positive DBSCAN neighborhood radius applied after
#'       column-standardizing d3. Default: `0.02`.}
#'     \item{`general$cluster$minPts`}{Positive DBSCAN core-point neighborhood threshold.
#'       Default: `20L`.}
#'   }
#'
#'   Experiment 4 end-to-end metaCCS settings are supplied as
#'   `metaccs = list(...)`:
#'
#'   \describe{
#'     \item{`metaccs$resample_seeds`}{Non-empty, finite, unique seeds defining
#'       tissue-by-cohort stratified sample-composition repeats. Default:
#'       `seed + 4000L + seq_len(10)`.}
#'     \item{`metaccs$umap_seeds`}{Non-empty, finite, unique algorithm seeds run
#'       within each metaCCS resample. Default:
#'       `seed + 5000L + seq_len(5)`.}
#'     \item{`metaccs$subsample_fraction`}{Fraction in `(0, 1]` retained within
#'       each tissue-by-cohort stratum. A value of `1` measures algorithmic rather
#'       than sample-composition variation. Default: `0.8`.}
#'     \item{`metaccs$parameter_mode`}{`"fixed"` uses the shared `general$dr` and
#'       `general$cluster` settings once; `"grid"` evaluates the Cartesian product of
#'       `dr_grid` and `cluster_grid` for both arms. Default: `"fixed"`.}
#'     \item{`metaccs$dr_grid`}{For grid mode, a non-empty list of named partial
#'       `general$dr` lists, or a data frame whose rows are partial configurations.
#'       Each entry is merged onto `general$dr`. `NULL` uses only the base configuration.
#'       Default: `NULL`.}
#'     \item{`metaccs$cluster_grid`}{For grid mode, a non-empty list of named
#'       partial `general$cluster` lists, or a data frame whose rows are partial
#'       configurations. Each entry is merged onto `general$cluster`. `NULL` uses only
#'       the base configuration. Default: `NULL`.}
#'     \item{`metaccs$direct_feature_mode`}{Feature-support rule for the Direct
#'       arm. The only currently supported value, `"tissue_model_union"`, uses
#'       the union of frozen-model GSClassifier input features within each tissue. Default:
#'       `"tissue_model_union"`.}
#'     \item{`metaccs$retain_assignments`}{Logical; retain per-sample raw DBSCAN
#'       cluster and noise assignments for every metaCCS run. Default: `TRUE`.}
#'   }
#'
#'   Gate 1 settings are supplied as `scaling = list(gate = list(...))`:
#'
#'   \describe{
#'     \item{`scaling$gate$enforce`}{Logical; Gate 1 is always calculated before
#'       scaling, but scaling is stopped on failure only when this is `TRUE`.
#'       Keep `FALSE` for exploratory analysis. Use `TRUE` for preregistered
#'       confirmatory or compute-gated analysis. Default: `FALSE`.}
#'     \item{`scaling$gate$primary_metric`}{Experiment 1 metric used for the three
#'       Cohort-minus-baseline decisions. If `"balanced_accuracy"` is not
#'       estimable, the current implementation falls back to `"biology_purity"`.
#'       Default: `"balanced_accuracy"`.}
#'     \item{`scaling$gate$min_gain`}{Minimum required lower bound of the paired 95 percent
#'       confidence interval for each of Cohort-Direct, Cohort-Null-RP, and
#'       Cohort-Null-Perm on the primary metric. `0` requires non-negative
#'       evidence; `0.02` requires at least a 0.02 lower-bound gain. Default: `0`.}
#'     \item{`scaling$gate$purity_tolerance`}{Largest accepted decrease in
#'       biology purity: the Cohort-Direct paired CI lower bound must be at least
#'       the negative tolerance. For example, `0.01` permits a lower bound of
#'       `-0.01`. Default: `0`.}
#'     \item{`scaling$gate$mixing_tolerance`}{Largest accepted decrease in cohort
#'       mixing, interpreted in the same way as `purity_tolerance`. Default: `0`.}
#'   }
#'
#'   Output handling:
#'
#'   \describe{
#'     \item{`general$cover`}{Logical; allow writing into a non-empty `output.dir` and
#'       replacing same-named ablation result files. Default: `FALSE`.}
#'   }
#' @param seed Master random seed.
#' @param verbose Whether to report progress.
#'
#' @return An object of class `CCSAblation`.
#' @author Weibin Huang <hwb2012@@qq.com>
#' @md
#' @export
ablation <- function(
    object,
    data,
    metadata = NULL,
    experiment = "representation",
    output.dir = file.path(getwd(), "ccs-ablation"),
    params = list(),
    seed = 20260727,
    verbose = TRUE
) {
  experiment <- unique(as.character(experiment))
  legacy_only <- c("scaling", "tissue_first", "metaccs")
  unknown <- setdiff(experiment, c("representation", "cohort", legacy_only))
  if (length(unknown) > 0) {
    stop(
      "ablation: unknown experiment(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (identical(experiment, "cohort")) {
    warning(
      "ablation: experiment = 'cohort' is deprecated; using 'representation'.",
      call. = FALSE
    )
    experiment <- "representation"
  }
  if ("representation" %in% experiment && length(experiment) > 1) {
    stop(
      "ablation: representation cannot be mixed with legacy downstream experiments.",
      call. = FALSE
    )
  }
  if (identical(experiment, "representation")) {
    return(.ablation_run_representation(
      object = object,
      data = data,
      metadata = metadata,
      output.dir = output.dir,
      params = params,
      seed = seed,
      verbose = verbose
    ))
  }

  # Preserve explicitly requested downstream legacy experiments during the
  # transition. A mixed legacy request that includes cohort keeps its original
  # Gate-1 dependency rather than silently changing historical behavior.
  .ablation_legacy(
    object = object,
    data = data,
    metadata = metadata,
    experiment = experiment,
    output.dir = output.dir,
    params = params,
    seed = seed,
    verbose = verbose
  )
}


.ablation_legacy <- function(
    object,
    data,
    metadata = NULL,
    experiment = c("cohort", "scaling", "tissue_first", "metaccs"),
    output.dir = file.path(getwd(), "ccs-ablation"),
    params = list(),
    seed = 20260727,
    verbose = TRUE
) {
  # Step 1: Validate the public inputs and merge user overrides into reproducible defaults.
  if (!methods::is(object, "CCS")) {
    stop("ablation: object must be a CCS object.", call. = FALSE)
  }

  choices <- c("cohort", "scaling", "tissue_first", "metaccs")
  experiment <- unique(match.arg(experiment, choices, several.ok = TRUE))
  config <- .ablation_resolve_config(seed, params)

  # Keep each run isolated and never overwrite existing results unless cover is explicit.
  if (dir.exists(output.dir) &&
      length(list.files(output.dir)) > 0 &&
      !config$general$cover) {
    stop(
      paste0(
        "ablation: output.dir is not empty. Use a new directory or set ",
        "params$general$cover = TRUE."
      ),
      call. = FALSE
    )
  }
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # Step 2: Align frozen d1, GSClassifier inputs, and sample annotations.
  if (verbose) {
    luckyBase::LuckyVerbose("ablation: Prepare frozen CCS inputs...")
  }
  prepared <- .ablation_prepare_input(
    object = object,
    data = data,
    metadata = metadata,
    max_samples = config$general$max_samples,
    seed = seed
  )
  manifest <- .ablation_build_manifest(object, prepared, config, seed)
  saveRDS(manifest, file.path(output.dir, "manifest.rds"))
  saveRDS(config, file.path(output.dir, "config.rds"))

  experiments <- list()
  audit_parts <- list()

  # Step 3: Run the requested experiments. Gate 1 is evaluated before scaling,
  # so a scaling-only request still runs the cohort experiment first.
  if ("cohort" %in% experiment || "scaling" %in% experiment) {
    if (verbose) {
      luckyBase::LuckyVerbose("ablation: Experiment 1 - cohort representation...")
    }
    experiments$cohort <- .ablation_experiment_cohort(
      prepared = prepared,
      config = config,
      manifest = manifest,
      seed = seed,
      verbose = verbose
    )
    audit_parts$cohort <- experiments$cohort$audit
    saveRDS(experiments$cohort, file.path(output.dir, "experiment-01-cohort.rds"))
  }

  if ("scaling" %in% experiment) {
    # When enforced, Gate 1 stops scaling unless the cohort representation beats its baselines.
    gate <- .ablation_gate_one(experiments$cohort, config$scaling$gate)
    if (gate$pass || !config$scaling$gate$enforce) {
      if (verbose) {
        luckyBase::LuckyVerbose("ablation: Experiment 2 - cohort-axis scaling...")
      }
      experiments$scaling <- .ablation_experiment_scaling(
        prepared = prepared,
        config = config,
        manifest = manifest,
        seed = seed + 10000,
        verbose = verbose
      )
      experiments$scaling$gate_one <- gate
      audit_parts$scaling <- experiments$scaling$audit
    } else {
      experiments$scaling <- list(
        status = "stopped_by_gate_one",
        gate_one = gate,
        metrics = data.frame(),
        summary = data.frame(),
        audit = data.frame()
      )
      if (verbose) {
        luckyBase::LuckyVerbose(
          "ablation: Experiment 2 stopped because experiment 1 did not pass Gate 1."
        )
      }
    }
    saveRDS(experiments$scaling, file.path(output.dir, "experiment-02-scaling.rds"))
  }

  if ("tissue_first" %in% experiment) {
    if (verbose) {
      luckyBase::LuckyVerbose("ablation: Experiment 3 - tissue-first reduction...")
    }
    experiments$tissue_first <- .ablation_experiment_tissue_first(
      prepared = prepared,
      config = config,
      manifest = manifest,
      seed = seed + 20000,
      verbose = verbose
    )
    audit_parts$tissue_first <- experiments$tissue_first$audit
    saveRDS(
      experiments$tissue_first,
      file.path(output.dir, "experiment-03-tissue-first.rds")
    )
  }

  if ("metaccs" %in% experiment) {
    if (verbose) {
      luckyBase::LuckyVerbose("ablation: Experiment 4 - end-to-end metaCCS...")
    }
    experiments$metaccs <- .ablation_experiment_metaccs(
      prepared = prepared,
      config = config,
      manifest = manifest,
      seed = seed + 30000,
      verbose = verbose
    )
    audit_parts$metaccs <- experiments$metaccs$audit
    saveRDS(
      experiments$metaccs,
      file.path(output.dir, "experiment-04-metaccs.rds")
    )
  }

  audit <- .ablation_rbind(audit_parts)
  # Step 4: Combine experiments and audit fields, then save machine- and reviewer-friendly outputs.
  result <- structure(
    list(
      call = match.call(),
      manifest = manifest,
      config = config,
      experiments = experiments,
      audit = audit,
      output.dir = normalizePath(output.dir, winslash = "/", mustWork = TRUE)
    ),
    class = "CCSAblation"
  )

  utils::write.csv(audit, file.path(output.dir, "audit.csv"), row.names = FALSE)
  saveRDS(result, file.path(output.dir, "ablation-result.rds"))
  if (verbose) {
    luckyBase::LuckyVerbose("ablation: All requested experiments are complete.")
  }
  result
}


# Centralize defaults shared by all experiments. Independent seed ranges let
# each random control or reduction repeat be reproduced in isolation.
.ablation_default_params <- function(seed = 20260727) {
  list(
    general = list(
      rank = 50L,
      k = 30L,
      distance = "euclidean",
      n_folds = Inf,
      bootstrap = 1000L,
      max_samples = Inf,
      probe = TRUE,
      probe_label = "tissue",
      probe_nrounds = 50L,
      numCores = 1L,
      fidelity_samples = 2000L,
      dr = list(
        method = "UWOT",
        dimension = c(5L, 2L),
        n_neighbors = 30L,
        min_dist = 0.01,
        spread = 0.75,
        set_op_mix_ratio = 1,
        metric = "euclidean",
        n_threads = NULL
      ),
      cluster = list(eps = 0.02, minPts = 20L),
      cover = FALSE
    ),
    cohort = list(
      rank_sensitivity = c(25L, 50L, 100L),
      geometry_samples = 5000L,
      distance_pairs = 100000L,
      mechanism_samples = 1000L,
      rp_density = 1 / 3,
      rp_seeds = seed + seq_len(20),
      permutation_seeds = seed + 1000L + seq_len(20)
    ),
    scaling = list(
      counts = c(10L, 25L, 50L, 75L, 100L, 125L, 150L),
      sequences = 100L,
      embedding_counts = c(25L, 50L, 100L, 150L),
      embedding_sequences = 10L,
      embedding_seeds = seed + 2000L + seq_len(10),
      subsample_fraction = 0.8,
      gate = list(
        enforce = FALSE,
        primary_metric = "balanced_accuracy",
        min_gain = 0,
        purity_tolerance = 0,
        mixing_tolerance = 0
      )
    ),
    tissue_first = list(
      seeds = seed + 3000L + seq_len(20),
      subsample_fraction = 0.8
    ),
    metaccs = list(
      resample_seeds = seed + 4000L + seq_len(10),
      umap_seeds = seed + 5000L + seq_len(5),
      subsample_fraction = 0.8,
      parameter_mode = "fixed",
      dr_grid = NULL,
      cluster_grid = NULL,
      direct_feature_mode = "tissue_model_union",
      retain_assignments = TRUE
    )
  )
}


# Resolve the public nested schema while keeping existing flat calls compatible.
.ablation_resolve_config <- function(seed, params = list()) {
  default <- .ablation_default_params(seed)
  override <- .ablation_normalize_params(params, default)
  config <- .ablation_merge_lists(default, override)
  .ablation_validate_config(config)
  config
}


.ablation_normalize_params <- function(params, default = .ablation_default_params()) {
  if (!is.list(params)) {
    stop("ablation: params must be a named list.", call. = FALSE)
  }
  if (length(params) == 0) {
    return(list())
  }
  if (is.null(names(params)) || any(!nzchar(names(params))) || anyDuplicated(names(params))) {
    stop("ablation: params must have unique, non-empty names.", call. = FALSE)
  }

  groups <- names(default)
  legacy_paths <- list(
    rank = c("general", "rank"),
    k = c("general", "k"),
    distance = c("general", "distance"),
    n_folds = c("general", "n_folds"),
    bootstrap = c("general", "bootstrap"),
    max_samples = c("general", "max_samples"),
    probe = c("general", "probe"),
    probe_label = c("general", "probe_label"),
    probe_nrounds = c("general", "probe_nrounds"),
    numCores = c("general", "numCores"),
    fidelity_samples = c("general", "fidelity_samples"),
    dr = c("general", "dr"),
    cluster = c("general", "cluster"),
    cover = c("general", "cover"),
    rank_sensitivity = c("cohort", "rank_sensitivity"),
    geometry_samples = c("cohort", "geometry_samples"),
    distance_pairs = c("cohort", "distance_pairs"),
    mechanism_samples = c("cohort", "mechanism_samples"),
    rp_density = c("cohort", "rp_density"),
    rp_seeds = c("cohort", "rp_seeds"),
    permutation_seeds = c("cohort", "permutation_seeds"),
    scaling_counts = c("scaling", "counts"),
    scaling_sequences = c("scaling", "sequences"),
    scaling_embedding_counts = c("scaling", "embedding_counts"),
    scaling_embedding_sequences = c("scaling", "embedding_sequences"),
    scaling_embedding_seeds = c("scaling", "embedding_seeds"),
    scaling_subsample_fraction = c("scaling", "subsample_fraction"),
    gate1 = c("scaling", "gate"),
    tissue_seeds = c("tissue_first", "seeds"),
    tissue_subsample_fraction = c("tissue_first", "subsample_fraction")
  )
  known <- c(groups, names(legacy_paths))
  unknown <- setdiff(names(params), known)
  if (length(unknown) > 0) {
    stop(
      "ablation: unknown params field(s): ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }

  override <- params[intersect(groups, names(params))]
  for (legacy_name in intersect(names(legacy_paths), names(params))) {
    path <- legacy_paths[[legacy_name]]
    group <- path[1]
    field <- path[2]
    if (!is.null(override[[group]]) && field %in% names(override[[group]])) {
      stop(
        "ablation: params$", legacy_name, " and params$", group, "$", field,
        " were provided more than once.",
        call. = FALSE
      )
    }
    if (is.null(override[[group]])) {
      override[[group]] <- list()
    }
    override[[group]][[field]] <- params[[legacy_name]]
  }
  override <- override[intersect(groups, names(override))]
  .ablation_validate_override(default, override)
  override
}


.ablation_validate_override <- function(default, override, path = "params") {
  if (!is.list(override)) {
    stop("ablation: ", path, " must be a named list.", call. = FALSE)
  }
  if (length(override) == 0) {
    return(invisible(TRUE))
  }
  if (is.null(names(override)) || any(!nzchar(names(override))) || anyDuplicated(names(override))) {
    stop("ablation: ", path, " must have unique, non-empty names.", call. = FALSE)
  }
  unknown <- setdiff(names(override), names(default))
  if (length(unknown) > 0) {
    stop(
      "ablation: unknown params field(s): ",
      paste(paste0(path, "$", unknown), collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (name in intersect(names(default), names(override))) {
    if (is.list(default[[name]])) {
      .ablation_validate_override(
        default[[name]],
        override[[name]],
        paste0(path, "$", name)
      )
    }
  }
  invisible(TRUE)
}


# Recursively merge nested settings so callers only specify leaves that differ.
.ablation_merge_lists <- function(default, override) {
  if (length(override) == 0) {
    return(default)
  }
  for (name in names(override)) {
    if (is.list(default[[name]]) && is.list(override[[name]])) {
      default[[name]] <- .ablation_merge_lists(default[[name]], override[[name]])
    } else {
      default[[name]] <- override[[name]]
    }
  }
  default
}


# Validate the minimum settings that could invalidate an experiment before costly work starts.
.ablation_validate_config <- function(config) {
  if (length(config$general$rank) != 1 || config$general$rank < 1) {
    stop("ablation: params$general$rank must be a positive integer.", call. = FALSE)
  }
  if (length(config$general$k) != 1 || config$general$k < 1) {
    stop("ablation: params$general$k must be a positive integer.", call. = FALSE)
  }
  if (length(config$cohort$rp_seeds) < 1 ||
      length(config$cohort$permutation_seeds) < 1) {
    stop("ablation: both null controls require at least one seed.", call. = FALSE)
  }
  if (length(config$general$bootstrap) != 1 || config$general$bootstrap < 1) {
    stop("ablation: params$general$bootstrap must be a positive integer.", call. = FALSE)
  }
  if (length(config$general$dr$dimension) != 2 ||
      any(config$general$dr$dimension < 1)) {
    stop(
      "ablation: params$general$dr$dimension must contain two positive integers.",
      call. = FALSE
    )
  }
  if (length(config$metaccs$resample_seeds) < 1 ||
      length(config$metaccs$umap_seeds) < 1) {
    stop("ablation: metaccs requires resample_seeds and umap_seeds.", call. = FALSE)
  }
  metaccs_seeds <- c(config$metaccs$resample_seeds, config$metaccs$umap_seeds)
  if (any(!is.finite(metaccs_seeds)) ||
      anyDuplicated(config$metaccs$resample_seeds) ||
      anyDuplicated(config$metaccs$umap_seeds)) {
    stop("ablation: metaccs seeds must be finite and unique within each seed type.", call. = FALSE)
  }
  if (config$metaccs$subsample_fraction <= 0 ||
      config$metaccs$subsample_fraction > 1) {
    stop("ablation: params$metaccs$subsample_fraction must be in (0, 1].", call. = FALSE)
  }
  if (!config$metaccs$parameter_mode %in% c("fixed", "grid")) {
    stop("ablation: params$metaccs$parameter_mode must be fixed or grid.", call. = FALSE)
  }
  if (!identical(config$metaccs$direct_feature_mode, "tissue_model_union")) {
    stop("ablation: the supported Direct feature mode is tissue_model_union.", call. = FALSE)
  }
  .ablation_metaccs_parameter_manifest(config)
  invisible(TRUE)
}


# Recover tissue|cohort module boundaries from d1 column names. Each block remains
# intact during permutation and scaling so one frozen cohort model is never split.
.ablation_module_manifest <- function(object) {
  d1 <- object@Data$Probability$d1
  if (!is.matrix(d1) && !is.data.frame(d1)) {
    stop("ablation: object@Data$Probability$d1 must be a matrix.", call. = FALSE)
  }
  parts <- strsplit(colnames(d1), "|", fixed = TRUE)
  if (any(lengths(parts) < 3)) {
    stop(
      "ablation: d1 columns must follow tissue|cohort|feature naming.",
      call. = FALSE
    )
  }
  module_id <- vapply(parts, function(x) paste(x[1:2], collapse = "|"), character(1))
  module_levels <- unique(module_id)
  blocks <- lapply(module_levels, function(x) which(module_id == x))
  names(blocks) <- module_levels
  module_parts <- strsplit(module_levels, "|", fixed = TRUE)
  modules <- data.frame(
    module_id = module_levels,
    tissue = vapply(module_parts, `[`, character(1), 1),
    cohort = vapply(module_parts, `[`, character(1), 2),
    block_width = lengths(blocks),
    first_column = vapply(blocks, min, integer(1)),
    last_column = vapply(blocks, max, integer(1)),
    stringsAsFactors = FALSE
  )
  list(modules = modules, blocks = blocks)
}


# Read only the ordinary gene-pair features used by frozen models.
.ablation_extract_tsp_features <- function(object, module_manifest) {
  .ablation_frozen_feature_manifest(object, module_manifest)$tsp_features
}


# Read the complete GSClassifier input feature support used by frozen models.
.ablation_extract_direct_features <- function(object, module_manifest) {
  .ablation_frozen_feature_manifest(object, module_manifest)$features
}


# Recover module- and tissue-level frozen GSClassifier feature support once.
.ablation_frozen_feature_manifest <- function(object, module_manifest) {
  if (!identical(object@Repeat$method, "GSClassifier")) {
    stop(
      "ablation: frozen cohort models must use method = 'GSClassifier'.",
      call. = FALSE
    )
  }
  models <- object@Model
  use_embedded <- length(models) > 0 && !identical(models, list(NA))
  path_map <- if (use_embedded) NULL else .ablation_model_path_map(object)

  module_records <- lapply(seq_len(nrow(module_manifest$modules)), function(i) {
    tissue <- module_manifest$modules$tissue[i]
    cohort <- module_manifest$modules$cohort[i]
    module_id <- module_manifest$modules$module_id[i]
    model <- if (use_embedded) {
      models[[tissue]][[cohort]]
    } else {
      path <- path_map[[module_id]]
      if (is.null(path)) {
        stop("ablation: frozen model is missing for module ", module_id, ".", call. = FALSE)
      }
      readRDS(path)
    }
    features <- .ablation_model_features(model)
    if (length(features) == 0) {
      stop(
        "ablation: frozen model has no GSClassifier features for module ",
        module_id,
        ".",
        call. = FALSE
      )
    }
    list(
      features = data.frame(
        module_id = module_id,
        tissue = tissue,
        cohort = cohort,
        feature = features,
        feature_type = .ablation_feature_type(features),
        stringsAsFactors = FALSE
      ),
      break_vectors = .ablation_model_break_vectors(model)
    )
  })
  module_features <- do.call(rbind, lapply(module_records, `[[`, "features"))
  features <- unique(module_features$feature)
  if (length(features) == 0) {
    stop("ablation: no GSClassifier features were found in frozen models.", call. = FALSE)
  }
  feature_manifest <- unique(module_features[, c("feature", "feature_type")])
  feature_manifest <- feature_manifest[match(features, feature_manifest$feature), , drop = FALSE]
  tissue_features <- lapply(
    split(module_features$feature, module_features$tissue),
    unique
  )
  break_vectors <- unlist(
    lapply(module_records, `[[`, "break_vectors"),
    recursive = FALSE,
    use.names = FALSE
  )
  break_hashes <- vapply(break_vectors, digest::digest, character(1), algo = "md5")
  unique_breaks <- break_vectors[!duplicated(break_hashes)]
  if (length(unique_breaks) != 1) {
    stop(
      "ablation: frozen GSClassifier models use heterogeneous breakVec values.",
      call. = FALSE
    )
  }
  tsp_features <- feature_manifest$feature[feature_manifest$feature_type == "gene_pair"]
  list(
    features = features,
    tsp_features = tsp_features,
    feature_manifest = feature_manifest,
    module_features = module_features,
    tissue_features = tissue_features,
    break_vec = unique_breaks[[1]],
    break_vec_hash = digest::digest(unique_breaks[[1]], algo = "md5"),
    direct_feature_hash = digest::digest(
      list(tissue_features = tissue_features, break_vec = unique_breaks[[1]]),
      algo = "md5"
    )
  )
}


# Collect the exact input features used across all repeats and frozen class models.
.ablation_model_features <- function(model) {
  if (is.null(model) || is.null(model$Model)) {
    stop("ablation: malformed frozen cohort model.", call. = FALSE)
  }
  unique(unlist(lapply(model$Model, function(repeat_model) {
    unlist(lapply(repeat_model, function(class_model) {
      if (is.null(class_model)) {
        return(character())
      }
      features <- class_model$bst$feature_names
      if (length(features) == 0) {
        features <- class_model$genes
      }
      features
    }), use.names = FALSE)
  }), use.names = FALSE))
}


.ablation_model_break_vectors <- function(model) {
  vectors <- unlist(lapply(model$Model, function(repeat_model) {
    lapply(repeat_model, function(class_model) {
      if (is.null(class_model)) NULL else class_model$breakVec
    })
  }), recursive = FALSE, use.names = FALSE)
  vectors <- Filter(function(x) length(x) > 0, vectors)
  if (length(vectors) == 0 || any(!vapply(vectors, function(x) {
    is.numeric(x) &&
      length(x) >= 2 &&
      all(is.finite(x)) &&
      all(x >= 0 & x <= 1) &&
      all(diff(x) > 0)
  }, logical(1)))) {
    stop("ablation: malformed GSClassifier breakVec in frozen model.", call. = FALSE)
  }
  lapply(vectors, as.numeric)
}


.ablation_feature_type <- function(features) {
  ifelse(
    grepl(":", features, fixed = TRUE),
    "gene_pair",
    ifelse(grepl("^s[0-9]+s[0-9]+$", features), "set_pair", "single_bin")
  )
}


# Map on-disk modelFit.rds files to tissue|cohort keys for objects without embedded models.
.ablation_model_path_map <- function(object) {
  paths <- list.files(
    object@Repeat$model.dir,
    pattern = "modelFit.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(paths) == 0) {
    stop("ablation: no frozen modelFit.rds files were found.", call. = FALSE)
  }
  keys <- vapply(paths, function(path) {
    pieces <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
    paste(tail(pieces, 3)[1:2], collapse = "|")
  }, character(1))
  stats::setNames(paths, keys)
}


# Align CCS, RNA expression, and metadata, then build GSClassifier features.
.ablation_prepare_input <- function(
    object,
    data,
    metadata = NULL,
    max_samples = Inf,
    seed = 20260727
) {
  # Derive the feature universe from frozen models before organizing expression and metadata.
  module_manifest <- .ablation_module_manifest(object)
  feature_manifest <- .ablation_frozen_feature_manifest(object, module_manifest)
  direct_features <- feature_manifest$features
  tsp_features <- feature_manifest$tsp_features
  flattened <- .ablation_flatten_expression(data)
  metadata <- .ablation_prepare_metadata(metadata, flattened$metadata, object)

  d1 <- as.matrix(object@Data$Probability$d1)
  # Duplicate sample IDs cannot be aligned uniquely; exclude the entire duplicate set and audit it.
  excluded_duplicate_samples <- intersect(
    flattened$excluded_duplicate_samples,
    rownames(d1)
  )
  sample_ids <- Reduce(
    intersect,
    list(rownames(d1), colnames(flattened$expr), metadata$sample_id)
  )
  if (length(sample_ids) < 3) {
    stop("ablation: fewer than three aligned samples are available.", call. = FALSE)
  }
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]

  # For capped runs, sample round-robin across tissue x cohort strata to limit dominance.
  if (is.finite(max_samples) && length(sample_ids) > max_samples) {
    keep <- .ablation_stratified_sample(
      metadata,
      size = as.integer(max_samples),
      seed = seed
    )
    metadata <- metadata[keep, , drop = FALSE]
    sample_ids <- metadata$sample_id
  }

  expr <- flattened$expr[, sample_ids, drop = FALSE]
  direct <- .ablation_gsclassifier_matrix(object, expr, feature_manifest)
  tsp_features <- feature_manifest$tsp_features
  tsp <- direct[, tsp_features, drop = FALSE]
  d1 <- d1[sample_ids, , drop = FALSE]
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
  rownames(metadata) <- metadata$sample_id

  list(
    direct = direct,
    tsp = tsp,
    d1 = d1,
    metadata = metadata,
    module_manifest = module_manifest,
    feature_manifest = feature_manifest,
    direct_features = direct_features,
    tsp_features = tsp_features,
    excluded_duplicate_samples = excluded_duplicate_samples
  )
}


# Normalize either one expression matrix or a tissue/cohort nested list to a common structure.
.ablation_flatten_expression <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) {
    expr <- as.matrix(data)
    duplicate_ids <- unique(colnames(expr)[duplicated(colnames(expr))])
    keep <- !colnames(expr) %in% duplicate_ids
    return(list(
      expr = expr[, keep, drop = FALSE],
      metadata = NULL,
      excluded_duplicate_samples = duplicate_ids
    ))
  }
  if (!is.list(data) || is.null(names(data))) {
    stop("ablation: data must be an expression matrix or named nested list.", call. = FALSE)
  }

  matrices <- list()
  annotations <- list()
  index <- 1L
  for (tissue in names(data)) {
    tissue_data <- data[[tissue]]
    for (cohort in names(tissue_data)) {
      leaf <- tissue_data[[cohort]]
      expr <- if (is.matrix(leaf) || is.data.frame(leaf)) leaf else leaf$expr
      expr <- as.matrix(expr)
      sample_ids <- colnames(expr)
      biology <- if (is.list(leaf) && !is.null(leaf$subtype)) {
        as.character(leaf$subtype)
      } else {
        rep(tissue, length(sample_ids))
      }
      if (length(biology) != length(sample_ids)) {
        stop("ablation: subtype length does not match expression samples.", call. = FALSE)
      }
      matrices[[index]] <- expr
      annotations[[index]] <- data.frame(
        sample_id = sample_ids,
        cohort = cohort,
        tissue = tissue,
        biology = biology,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  metadata <- do.call(rbind, annotations)
  duplicate_ids <- unique(metadata$sample_id[duplicated(metadata$sample_id)])
  # Match CCS::getResData(): join cohorts by gene union and retain missing cells as NA.
  genes <- unique(unlist(lapply(matrices, rownames), use.names = FALSE))
  expr <- do.call(cbind, lapply(matrices, function(x) {
    aligned <- matrix(
      NA_real_,
      nrow = length(genes),
      ncol = ncol(x),
      dimnames = list(genes, colnames(x))
    )
    aligned[rownames(x), ] <- x
    aligned
  }))
  keep <- !metadata$sample_id %in% duplicate_ids
  list(
    expr = expr[, keep, drop = FALSE],
    metadata = metadata[keep, , drop = FALSE],
    excluded_duplicate_samples = duplicate_ids
  )
}


# Normalize common metadata aliases and derive tissue/biology from CCS when possible.
.ablation_prepare_metadata <- function(metadata, derived, object) {
  if (is.null(metadata)) {
    if (is.null(derived)) {
      stop("ablation: metadata is required for matrix input.", call. = FALSE)
    }
    metadata <- derived
  }
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  names_lower <- tolower(colnames(metadata))
  # Keep the closure's name cache synchronized while selecting the first alias for each field.
  rename_one <- function(target, candidates) {
    hit <- match(tolower(candidates), names_lower, nomatch = 0)
    hit <- hit[hit > 0]
    if (length(hit) > 0 && !target %in% colnames(metadata)) {
      colnames(metadata)[hit[1]] <<- target
      names_lower <<- tolower(colnames(metadata))
    }
  }
  rename_one("sample_id", c("sample_id", "sampleids", "id"))
  rename_one("cohort", c("cohort", "dataset", "study"))
  rename_one("tissue", c("tissue", "cancertype", "cancer_type", "tumor_type"))
  rename_one("biology", c("biology", "subtype", "label", "anchor"))

  required <- c("sample_id", "cohort")
  if (!all(required %in% colnames(metadata))) {
    stop("ablation: metadata must contain sample_id and cohort.", call. = FALSE)
  }
  if (!"tissue" %in% colnames(metadata)) {
    metadata$tissue <- unname(object@Data$CancerType[metadata$sample_id])
  }
  if (!"biology" %in% colnames(metadata)) {
    metadata$biology <- metadata$tissue
  }
  core_columns <- c("sample_id", "cohort", "tissue", "biology")
  metadata[, c(core_columns, setdiff(colnames(metadata), core_columns)), drop = FALSE]
}


# Sample round-robin across tissue x cohort strata; a fixed seed reproduces the selected rows.
.ablation_stratified_sample <- function(metadata, size, seed) {
  set.seed(seed)
  groups <- split(seq_len(nrow(metadata)), interaction(
    metadata$tissue,
    metadata$cohort,
    drop = TRUE
  ))
  shuffled <- lapply(groups, sample)
  selected <- integer()
  while (length(selected) < size && any(lengths(shuffled) > 0)) {
    for (name in names(shuffled)) {
      if (length(shuffled[[name]]) > 0 && length(selected) < size) {
        selected <- c(selected, shuffled[[name]][1])
        shuffled[[name]] <- shuffled[[name]][-1]
      }
    }
  }
  sort(selected)
}


# Reconstruct the frozen GSClassifier input space without refitting any model.
.ablation_gsclassifier_matrix <- function(object, expr, feature_manifest) {
  if (!identical(object@Repeat$method, "GSClassifier")) {
    stop(
      "ablation: Direct features require method = 'GSClassifier'.",
      call. = FALSE
    )
  }
  if (length(object@Repeat$geneSet) == 0) {
    stop("ablation: GSClassifier Direct features require a non-empty geneSet.", call. = FALSE)
  }
  matched <- GSClassifier::geneMatch(
    X = expr,
    geneAnnotation = object@Repeat$geneAnnotation,
    geneid = object@Repeat$geneid,
    matchmode = "fix"
  )
  build_features <- getFromNamespace("trainDataProc_X", "GSClassifier")
  transformed <- build_features(
    Xmat = as.matrix(matched$Subset),
    geneSet = object@Repeat$geneSet,
    breakVec = feature_manifest$break_vec
  )$dat$Xbin
  transformed <- as.matrix(transformed)
  missing <- setdiff(feature_manifest$features, colnames(transformed))
  if (length(missing) > 0) {
    stop(
      "ablation: GSClassifier did not reconstruct frozen features: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  direct <- transformed[, feature_manifest$features, drop = FALSE]
  if (any(!is.finite(direct))) {
    stop(
      "ablation: GSClassifier Direct features contain non-finite values.",
      call. = FALSE
    )
  }
  direct
}


# Record input sizes, versions, settings, and sample/feature hashes for audit and comparison.
.ablation_build_manifest <- function(object, prepared, config, seed) {
  sample_hash <- digest::digest(sort(prepared$metadata$sample_id), algo = "md5")
  feature_hash <- digest::digest(
    list(
      tsp = prepared$tsp_features,
      direct = prepared$direct_features,
      modules = prepared$module_manifest$modules,
      d1_columns = colnames(prepared$d1)
    ),
    algo = "md5"
  )
  list(
    version = 2L,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = seed,
    object_class = class(object)[1],
    object_version = tryCatch(as.character(utils::packageVersion("CCS")), error = function(e) NA_character_),
    gsclassifier_version = tryCatch(
      as.character(utils::packageVersion("GSClassifier")),
      error = function(e) NA_character_
    ),
    gsclassifier_feature_builder = "trainDataProc_X",
    model_dir = object@Repeat$model.dir,
    sample_count = nrow(prepared$d1),
    tsp_feature_count = ncol(prepared$tsp),
    direct_feature_count = ncol(prepared$direct),
    direct_feature_type_count = table(factor(
      prepared$feature_manifest$feature_manifest$feature_type,
      levels = c("single_bin", "gene_pair", "set_pair")
    )),
    d1_dimension = dim(prepared$d1),
    module_count = nrow(prepared$module_manifest$modules),
    excluded_duplicate_sample_count = length(prepared$excluded_duplicate_samples),
    sample_manifest_hash = sample_hash,
    feature_manifest_hash = feature_hash,
    direct_feature_hash = prepared$feature_manifest$direct_feature_hash,
    direct_break_vec = prepared$feature_manifest$break_vec,
    direct_break_vec_hash = prepared$feature_manifest$break_vec_hash,
    tissue_feature_count = vapply(
      prepared$feature_manifest$tissue_features,
      length,
      integer(1)
    ),
    tissue_feature_hash = vapply(
      prepared$feature_manifest$tissue_features,
      digest::digest,
      character(1),
      algo = "md5"
    ),
    module_manifest = prepared$module_manifest$modules,
    config_hash = digest::digest(config, algo = "md5")
  )
}


# Use whole cohorts as validation units and greedily balance sample counts across folds.
# A cohort never appears in both training and test data, preventing cohort leakage.
.ablation_grouped_folds <- function(
    cohort,
    n_folds = Inf,
    seed = 20260727,
    label = NULL
) {
  sizes <- sort(table(cohort), decreasing = TRUE)
  n_groups <- length(sizes)
  n_folds <- if (is.finite(n_folds)) min(as.integer(n_folds), n_groups) else n_groups
  if (n_folds < 2) {
    stop("ablation: grouped validation requires at least two cohorts.", call. = FALSE)
  }
  set.seed(seed)
  ordered <- names(sizes)
  tied <- split(ordered, sizes[ordered])
  ordered <- unlist(lapply(tied, sample), use.names = FALSE)
  fold_load <- numeric(n_folds)
  group_fold <- integer(n_groups)
  names(group_fold) <- ordered

  label_table <- NULL
  label_load <- NULL
  if (!is.null(label)) {
    if (length(label) != length(cohort)) {
      stop("ablation: label and cohort must have the same length.", call. = FALSE)
    }
    valid <- !is.na(label) & nzchar(as.character(label))
    label_table <- table(
      factor(as.character(cohort[valid]), levels = ordered),
      as.character(label[valid])
    )
    label_load <- matrix(
      0,
      nrow = n_folds,
      ncol = ncol(label_table),
      dimnames = list(seq_len(n_folds), colnames(label_table))
    )
  }

  for (group in ordered) {
    if (is.null(label_table) || ncol(label_table) == 0) {
      fold <- which.min(fold_load)
    } else {
      group_labels <- label_table[group, , drop = TRUE]
      label_total <- pmax(colSums(label_table), 1)
      label_cost <- vapply(seq_len(n_folds), function(i) {
        sum(((label_load[i, ] + group_labels)^2) / label_total)
      }, numeric(1))
      load_cost <- (fold_load + as.numeric(sizes[group])) / sum(sizes)
      fold <- which.min(label_cost + load_cost)
      label_load[fold, ] <- label_load[fold, ] + group_labels
    }
    group_fold[group] <- fold
    fold_load[fold] <- fold_load[fold] + sizes[group]
  }
  unname(group_fold[as.character(cohort)])
}


# Null-Perm independently shuffles each module within each cohort while moving all
# block columns together. This preserves block and cohort marginals but breaks cross-block pairing.
.ablation_permute_blocks <- function(d1, blocks, cohort, seed) {
  set.seed(seed)
  result <- d1
  cohort_rows <- split(seq_len(nrow(d1)), cohort)
  for (block in blocks) {
    for (rows in cohort_rows) {
      permutation <- rows[sample.int(length(rows))]
      result[rows, block] <- d1[permutation, block, drop = FALSE]
    }
  }
  result
}


# Estimate scaling parameters on training data and apply them unchanged to prevent test leakage.
.ablation_scale_train_apply <- function(train, test) {
  center <- colMeans(train)
  scale <- apply(train, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(train, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(test, 2, center, "-"), 2, scale, "/"),
    center = center,
    scale = scale
  )
}


# Fit fixed-rank PCA on training data and project test data; dimensions cap the achievable rank.
.ablation_fit_pca <- function(train, test, rank_q) {
  scaled <- .ablation_scale_train_apply(train, test)
  rank_q <- min(as.integer(rank_q), ncol(train), nrow(train) - 1L)
  fit <- stats::prcomp(
    scaled$train,
    center = FALSE,
    scale. = FALSE,
    rank. = rank_q
  )
  rank_q <- min(rank_q, ncol(fit$rotation))
  rotation <- fit$rotation[, seq_len(rank_q), drop = FALSE]
  list(
    train = scaled$train %*% rotation,
    test = scaled$test %*% rotation,
    rank = rank_q,
    rotation = rotation,
    center = scaled$center,
    scale = scaled$scale
  )
}


# Generate a sparse Achlioptas matrix for an equal-rank Null-RP without cohort learning.
.ablation_projection_matrix <- function(n_features, rank_q, seed, density = 1 / 3) {
  set.seed(seed)
  values <- sample(
    c(-sqrt(1 / density), 0, sqrt(1 / density)),
    n_features * rank_q,
    replace = TRUE,
    prob = c(density / 2, 1 - density, density / 2)
  )
  matrix(values / sqrt(rank_q), nrow = n_features, ncol = rank_q)
}


# Apply a fixed-seed random projection after optional scaling for auxiliary comparisons.
.ablation_random_projection <- function(
    data,
    rank_q,
    seed,
    density = 1 / 3,
    standardize = TRUE
) {
  x <- as.matrix(data)
  if (standardize) {
    x <- .ablation_scale_train_apply(x, x)$train
  }
  projection <- .ablation_projection_matrix(ncol(x), rank_q, seed, density)
  result <- x %*% projection
  rownames(result) <- rownames(data)
  result
}


# Linear CKA measures global geometric agreement; values closer to 1 indicate greater similarity.
.ablation_linear_cka <- function(x, y) {
  x <- scale(as.matrix(x), center = TRUE, scale = FALSE)
  y <- scale(as.matrix(y), center = TRUE, scale = FALSE)
  numerator <- sum(crossprod(x, y)^2)
  denominator <- sqrt(sum(crossprod(x)^2) * sum(crossprod(y)^2))
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}


# Return row indices of each sample's k nearest neighbors, capping k at n - 1.
.ablation_knn <- function(data, k) {
  data <- as.matrix(data)
  k <- min(as.integer(k), nrow(data) - 1L)
  if (k < 1) {
    stop("ablation: kNN requires at least two samples.", call. = FALSE)
  }
  as.matrix(dbscan::kNN(data, k = k)$id)
}


# Compare neighborhood sets per sample; higher mean Jaccard indicates better local agreement.
.ablation_knn_jaccard <- function(x, y, k) {
  x_nn <- .ablation_knn(x, k)
  y_nn <- .ablation_knn(y, k)
  mean(vapply(seq_len(nrow(x_nn)), function(i) {
    length(intersect(x_nn[i, ], y_nn[i, ])) /
      length(union(x_nn[i, ], y_nn[i, ]))
  }, numeric(1)))
}


# Define effective rank as the exponential entropy of singular-value energy.
.ablation_effective_rank <- function(data) {
  x <- scale(as.matrix(data), center = TRUE, scale = FALSE)
  singular <- svd(x, nu = 0, nv = 0)$d
  probability <- singular^2 / sum(singular^2)
  probability <- probability[probability > 0]
  exp(-sum(probability * log(probability)))
}


# Sample row pairs and correlate Euclidean-distance ranks to avoid a full distance matrix.
.ablation_distance_spearman <- function(x, y, n_pairs, seed) {
  n <- nrow(x)
  total <- n * (n - 1) / 2
  n_pairs <- min(as.integer(n_pairs), total)
  set.seed(seed)
  first <- sample.int(n, n_pairs, replace = TRUE)
  second <- sample.int(n, n_pairs, replace = TRUE)
  same <- first == second
  while (any(same)) {
    second[same] <- sample.int(n, sum(same), replace = TRUE)
    same <- first == second
  }
  distance_x <- sqrt(rowSums((x[first, , drop = FALSE] - x[second, , drop = FALSE])^2))
  distance_y <- sqrt(rowSums((y[first, , drop = FALSE] - y[second, , drop = FALSE])^2))
  suppressWarnings(stats::cor(distance_x, distance_y, method = "spearman"))
}


# Measure cross-cohort mixing within biology and global biological neighborhood purity.
# Higher values indicate stronger cohort removal and biological preservation, respectively.
.ablation_mixing_purity <- function(data, metadata, k) {
  nn <- .ablation_knn(data, k)
  biology <- as.character(metadata$biology)
  cohort <- as.character(metadata$cohort)
  purity <- vapply(seq_len(nrow(nn)), function(i) {
    mean(biology[nn[i, ]] == biology[i], na.rm = TRUE)
  }, numeric(1))

  mixing <- rep(NA_real_, nrow(nn))
  # Condition mixing on biology so tissue differences are not mistaken for cohort separation.
  for (label in unique(biology)) {
    rows <- which(biology == label)
    if (length(rows) < 2) {
      next
    }
    local_k <- min(k, length(rows) - 1L)
    local_nn <- .ablation_knn(data[rows, , drop = FALSE], local_k)
    mixing[rows] <- vapply(seq_along(rows), function(i) {
      mean(cohort[rows[local_nn[i, ]]] != cohort[rows[i]])
    }, numeric(1))
  }
  list(
    cohort_mixing = mean(mixing, na.rm = TRUE),
    biology_purity = mean(purity, na.rm = TRUE),
    per_sample = data.frame(
      sample_id = metadata$sample_id,
      cohort_mixing = mixing,
      biology_purity = purity,
      stringsAsFactors = FALSE
    )
  )
}


# Use a lightweight linear XGBoost probe to test decodable biological information.
# Fit only on training rows and evaluate only test classes observed during training.
.ablation_probe <- function(train, test, train_label, test_label, seed, config) {
  train_label <- as.character(train_label)
  test_label <- as.character(test_label)
  keep_train <- !is.na(train_label)
  train <- train[keep_train, , drop = FALSE]
  train_label <- train_label[keep_train]
  classes <- sort(unique(train_label))
  keep_test <- !is.na(test_label) & test_label %in% classes
  if (length(classes) < 2 || sum(keep_test) < 2) {
    return(c(macro_auroc = NA_real_, balanced_accuracy = NA_real_))
  }
  encoded <- match(train_label, classes) - 1L
  set.seed(seed)
  fit <- xgboost::xgboost(
    data = as.matrix(train),
    label = encoded,
    booster = "gblinear",
    updater = "coord_descent",
    feature_selector = "cyclic",
    objective = "multi:softprob",
    num_class = length(classes),
    nrounds = as.integer(config$general$probe_nrounds),
    eta = 0.1,
    lambda = 1,
    alpha = 0,
    nthread = as.integer(config$general$numCores),
    verbose = 0
  )
  probability <- predict(fit, as.matrix(test[keep_test, , drop = FALSE]))
  probability <- matrix(probability, ncol = length(classes), byrow = TRUE)
  truth <- test_label[keep_test]
  prediction <- classes[max.col(probability, ties.method = "first")]
  recalls <- vapply(classes, function(class) {
    rows <- truth == class
    if (!any(rows)) NA_real_ else mean(prediction[rows] == class)
  }, numeric(1))
  auc <- vapply(seq_along(classes), function(i) {
    .ablation_binary_auc(as.integer(truth == classes[i]), probability[, i])
  }, numeric(1))
  c(
    macro_auroc = mean(auc, na.rm = TRUE),
    balanced_accuracy = mean(recalls, na.rm = TRUE)
  )
}


# Keep labels spanning at least two cohorts so the probe cannot infer them from cohort identity.
.ablation_eligible_probe_labels <- function(metadata, label_column) {
  label <- as.character(metadata[[label_column]])
  cohort <- as.character(metadata$cohort)
  valid <- !is.na(label) & nzchar(label)
  cohort_count <- vapply(
    split(cohort[valid], label[valid]),
    function(x) length(unique(x)),
    integer(1)
  )
  eligible <- names(cohort_count)[cohort_count >= 2]
  label[!label %in% eligible] <- NA_character_
  label
}


# Compute binary AUC from rank sums; return NA when either class is absent.
.ablation_binary_auc <- function(label, score) {
  positive <- sum(label == 1)
  negative <- sum(label == 0)
  if (positive == 0 || negative == 0) {
    return(NA_real_)
  }
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[label == 1]) - positive * (positive + 1) / 2) / (positive * negative)
}


# Compare Direct and candidate distance-rank changes along original feature edges,
# separating biology-discordant edges from biology-matched, cross-cohort edges.
.ablation_selective_reconstruction <- function(
    original,
    direct,
    candidate,
    metadata,
    k,
    max_samples,
    seed
) {
  if (nrow(original) > max_samples) {
    rows <- .ablation_stratified_sample(metadata, max_samples, seed)
    original <- original[rows, , drop = FALSE]
    direct <- direct[rows, , drop = FALSE]
    candidate <- candidate[rows, , drop = FALSE]
    metadata <- metadata[rows, , drop = FALSE]
  }
  k <- min(k, nrow(original) - 1L)
  edges <- .ablation_knn(original, k)
  direct_distance <- as.matrix(stats::dist(direct))
  candidate_distance <- as.matrix(stats::dist(candidate))
  direct_rank <- t(apply(direct_distance, 1, rank, ties.method = "average")) - 1
  candidate_rank <- t(apply(candidate_distance, 1, rank, ties.method = "average")) - 1
  # delta > 0 means an original neighbor ranks farther away in the candidate than in Direct.
  anchor <- rep(seq_len(nrow(edges)), each = ncol(edges))
  neighbor <- as.vector(t(edges))
  delta <- candidate_rank[cbind(anchor, neighbor)] - direct_rank[cbind(anchor, neighbor)]
  biology_same <- metadata$biology[anchor] == metadata$biology[neighbor]
  cohort_same <- metadata$cohort[anchor] == metadata$cohort[neighbor]
  c(
    discordant_rank_change = mean(delta[!biology_same], na.rm = TRUE),
    concordant_cross_cohort_rank_change = mean(
      delta[biology_same & !cohort_same],
      na.rm = TRUE
    )
  )
}


# Collect all metrics for one Experiment 1 group x fold x seed and attach random settings.
.ablation_metric_rows <- function(
    direct_scores,
    candidate_scores,
    original_test,
    metadata_test,
    metadata_train,
    train_scores,
    group_id,
    fold_id,
    seed_type,
    stochastic_seed,
    rank_q,
    config,
    seed
) {
  local_k <- min(config$general$k, nrow(candidate_scores) - 1L)
  # Compute geometry, mixing, probe, and mechanism metrics on the same test sample set.
  mixing <- .ablation_mixing_purity(candidate_scores, metadata_test, local_k)
  mechanism <- .ablation_selective_reconstruction(
    original = original_test,
    direct = direct_scores,
    candidate = candidate_scores,
    metadata = metadata_test,
    k = local_k,
    max_samples = config$cohort$mechanism_samples,
    seed = seed
  )
  probe <- if (config$general$probe) {
    probe_metadata <- rbind(metadata_train, metadata_test)
    probe_label <- .ablation_eligible_probe_labels(
      probe_metadata,
      config$general$probe_label
    )
    train_label <- probe_label[seq_len(nrow(metadata_train))]
    test_label <- probe_label[nrow(metadata_train) + seq_len(nrow(metadata_test))]
    .ablation_probe(
      train = train_scores,
      test = candidate_scores,
      train_label = train_label,
      test_label = test_label,
      seed = seed,
      config = config
    )
  } else {
    c(macro_auroc = NA_real_, balanced_accuracy = NA_real_)
  }
  values <- c(
    linear_cka = .ablation_linear_cka(direct_scores, candidate_scores),
    distance_spearman = .ablation_distance_spearman(
      direct_scores,
      candidate_scores,
      config$cohort$distance_pairs,
      seed
    ),
    knn_jaccard = .ablation_knn_jaccard(direct_scores, candidate_scores, local_k),
    effective_rank = .ablation_effective_rank(candidate_scores),
    cohort_mixing = mixing$cohort_mixing,
    biology_purity = mixing$biology_purity,
    probe,
    mechanism
  )
  data.frame(
    experiment_id = "cohort",
    group_id = group_id,
    fold_id = as.character(fold_id),
    held_out_cohort = paste(unique(metadata_test$cohort), collapse = ";"),
    seed_type = seed_type,
    stochastic_seed = stochastic_seed,
    rank_q = rank_q,
    k = local_k,
    distance = config$general$distance,
    metric_name = names(values),
    metric_value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}


# Experiment 1 compares Direct, frozen Cohort, and two null representations to test
# whether the cohort layer adds biological information beyond compression or randomness.
.ablation_experiment_cohort <- function(prepared, config, manifest, seed, verbose) {
  folds <- .ablation_grouped_folds(
    prepared$metadata$cohort,
    n_folds = config$general$n_folds,
    seed = seed
  )
  # Generate each permutation once and reuse it across ranks and folds for strict pairing.
  permuted <- lapply(config$cohort$permutation_seeds, function(seed_i) {
    .ablation_permute_blocks(
      prepared$d1,
      prepared$module_manifest$blocks,
      prepared$metadata$cohort,
      seed_i
    )
  })
  names(permuted) <- as.character(config$cohort$permutation_seeds)

  minimum_train_size <- min(vapply(
    sort(unique(folds)),
    function(fold) sum(folds != fold),
    integer(1)
  ))
  maximum_rank <- min(
    ncol(prepared$direct),
    ncol(prepared$d1),
    minimum_train_size - 1L
  )
  # Truncate all requested ranks to the largest common rank estimable in every fold.
  rank_targets <- unique(pmin(
    as.integer(c(config$general$rank, config$cohort$rank_sensitivity)),
    maximum_rank
  ))
  main_rank <- min(as.integer(config$general$rank), maximum_rank)
  metrics_by_rank <- lapply(rank_targets, function(rank_target) {
    .ablation_cohort_rank_metrics(
      prepared = prepared,
      config = config,
      folds = folds,
      permuted = permuted,
      rank_target = rank_target,
      seed = seed
    )
  })
  names(metrics_by_rank) <- as.character(rank_targets)
  metrics <- metrics_by_rank[[as.character(main_rank)]]
  sensitivity_parts <- metrics_by_rank[names(metrics_by_rank) != as.character(main_rank)]
  sensitivity_metrics <- .ablation_rbind(sensitivity_parts)

  # Use the main rank for Gate 1 and the remaining ranks only for sensitivity analysis.
  geometry <- .ablation_dimension_free_geometry(prepared, config, seed)
  summary <- .ablation_metric_summary(metrics, config$general$bootstrap, seed)
  contrasts <- .ablation_paired_contrasts(metrics, config$general$bootstrap, seed)
  sensitivity <- list(
    metrics = sensitivity_metrics,
    summary = if (nrow(sensitivity_metrics) > 0) {
      .ablation_metric_summary(
        sensitivity_metrics,
        config$general$bootstrap,
        seed + 500L
      )
    } else {
      data.frame()
    }
  )
  audit <- .ablation_add_audit_hashes(
    .ablation_rbind(c(list(main = metrics), sensitivity_parts)),
    manifest
  )
  list(
    status = "complete",
    folds = folds,
    main_rank = main_rank,
    dimension_free_geometry = geometry,
    metrics = metrics,
    summary = summary,
    contrasts = contrasts,
    rank_sensitivity = sensitivity,
    audit = audit
  )
}


# Run the full leave-cohort-out comparison at rank_q. All groups share folds and test
# samples, while Null-RP/Null-Perm seeds remain within-fold stochastic repeats.
.ablation_cohort_rank_metrics <- function(
    prepared,
    config,
    folds,
    permuted,
    rank_target,
    seed
) {
  metrics <- list()
  index <- 1L
  for (fold in sort(unique(folds))) {
    test_rows <- folds == fold
    train_rows <- !test_rows
    rank_q <- min(
      as.integer(rank_target),
      ncol(prepared$direct),
      ncol(prepared$d1),
      sum(train_rows) - 1L
    )
    direct <- .ablation_fit_pca(
      prepared$direct[train_rows, , drop = FALSE],
      prepared$direct[test_rows, , drop = FALSE],
      rank_q
    )
    cohort <- .ablation_fit_pca(
      prepared$d1[train_rows, , drop = FALSE],
      prepared$d1[test_rows, , drop = FALSE],
      rank_q
    )
    metadata_train <- prepared$metadata[train_rows, , drop = FALSE]
    metadata_test <- prepared$metadata[test_rows, , drop = FALSE]
    direct_test <- prepared$direct[test_rows, , drop = FALSE]

    # Direct is equal-rank PCA of frozen GSClassifier inputs; Cohort uses d1.
    metrics[[index]] <- .ablation_metric_rows(
      direct$test, direct$test, direct_test, metadata_test, metadata_train,
      direct$train, "Direct", fold, NA_character_, NA_integer_, direct$rank,
      config, seed + fold
    )
    index <- index + 1L
    metrics[[index]] <- .ablation_metric_rows(
      direct$test, cohort$test, direct_test, metadata_test, metadata_train,
      cohort$train, "Cohort", fold, NA_character_, NA_integer_, cohort$rank,
      config, seed + 100L + fold
    )
    index <- index + 1L

    scaled_direct <- .ablation_scale_train_apply(
      prepared$direct[train_rows, , drop = FALSE],
      prepared$direct[test_rows, , drop = FALSE]
    )
    # Null-RP tests whether equal-rank random compression alone improves the metrics.
    for (seed_i in config$cohort$rp_seeds) {
      projection <- .ablation_projection_matrix(
        ncol(prepared$direct),
        rank_q,
        seed_i,
        config$cohort$rp_density
      )
      rp_train <- scaled_direct$train %*% projection
      rp_test <- scaled_direct$test %*% projection
      metrics[[index]] <- .ablation_metric_rows(
        direct$test, rp_test, direct_test, metadata_test, metadata_train,
        rp_train, "Null-RP", fold, "projection_seed", seed_i, rank_q,
        config, seed_i + fold
      )
      index <- index + 1L
    }

    # Null-Perm preserves within-cohort block structure but breaks cross-block sample pairing.
    for (seed_i in config$cohort$permutation_seeds) {
      perm_fit <- .ablation_fit_pca(
        permuted[[as.character(seed_i)]][train_rows, , drop = FALSE],
        permuted[[as.character(seed_i)]][test_rows, , drop = FALSE],
        rank_q
      )
      metrics[[index]] <- .ablation_metric_rows(
        direct$test, perm_fit$test, direct_test, metadata_test, metadata_train,
        perm_fit$train, "Null-Perm", fold, "permutation_seed", seed_i,
        perm_fit$rank, config, seed_i + fold
      )
      index <- index + 1L
    }
  }
  do.call(rbind, metrics)
}


# Compare raw GSClassifier inputs and d1 with dimension-free geometry metrics.
.ablation_dimension_free_geometry <- function(prepared, config, seed) {
  rows <- seq_len(nrow(prepared$d1))
  if (length(rows) > config$cohort$geometry_samples) {
    rows <- .ablation_stratified_sample(
      prepared$metadata,
      config$cohort$geometry_samples,
      seed
    )
  }
  direct <- prepared$direct[rows, , drop = FALSE]
  d1 <- prepared$d1[rows, , drop = FALSE]
  k <- min(config$general$k, length(rows) - 1L)
  data.frame(
    metric_name = c("linear_cka", "distance_spearman", "knn_jaccard"),
    metric_value = c(
      .ablation_linear_cka(direct, d1),
      .ablation_distance_spearman(
        direct,
        d1,
        config$cohort$distance_pairs,
        seed
      ),
      .ablation_knn_jaccard(direct, d1, k)
    ),
    sample_count = length(rows),
    stringsAsFactors = FALSE
  )
}


# Average stochastic repeats within each fold, then bootstrap folds for a 95% CI.
# This keeps cohorts as statistical units instead of treating seeds as independent samples.
.ablation_metric_summary <- function(metrics, bootstrap, seed) {
  keys <- interaction(metrics$group_id, metrics$metric_name, drop = TRUE)
  parts <- split(metrics, keys)
  result <- lapply(seq_along(parts), function(i) {
    data <- parts[[i]]
    group_id <- data$group_id[1]
    metric_name <- data$metric_name[1]
    data <- data[is.finite(data$metric_value), , drop = FALSE]
    fold_values <- if (nrow(data) == 0) {
      numeric()
    } else {
      stats::aggregate(
        metric_value ~ fold_id,
        data = data,
        FUN = mean
      )$metric_value
    }
    if (length(fold_values) == 0) {
      interval <- c(NA_real_, NA_real_)
      estimate <- NA_real_
    } else {
      set.seed(seed + i)
      boot <- replicate(
        as.integer(bootstrap),
        mean(sample(fold_values, replace = TRUE))
      )
      interval <- stats::quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      estimate <- mean(fold_values)
    }
    data.frame(
      group_id = group_id,
      metric_name = metric_name,
      estimate = estimate,
      ci_low = interval[1],
      ci_high = interval[2],
      n_folds = length(fold_values),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}


# Bootstrap Cohort-minus-baseline differences on matched folds to preserve pairing.
# Positive values therefore indicate a higher metric for Cohort.
.ablation_paired_contrasts <- function(metrics, bootstrap, seed) {
  data <- metrics[is.finite(metrics$metric_value), , drop = FALSE]
  data <- stats::aggregate(
    metric_value ~ group_id + fold_id + metric_name,
    data = data,
    FUN = mean
  )
  contrast_groups <- list(
    `Cohort-Direct` = c("Cohort", "Direct"),
    `Cohort-Null-RP` = c("Cohort", "Null-RP"),
    `Cohort-Null-Perm` = c("Cohort", "Null-Perm")
  )
  result <- list()
  index <- 1L
  for (metric in unique(data$metric_name)) {
    metric_data <- data[data$metric_name == metric, , drop = FALSE]
    wide <- reshape(
      metric_data[, c("fold_id", "group_id", "metric_value")],
      idvar = "fold_id",
      timevar = "group_id",
      direction = "wide"
    )
    for (contrast in names(contrast_groups)) {
      groups <- contrast_groups[[contrast]]
      first <- paste0("metric_value.", groups[1])
      second <- paste0("metric_value.", groups[2])
      if (!all(c(first, second) %in% colnames(wide))) {
        next
      }
      delta <- wide[[first]] - wide[[second]]
      delta <- delta[is.finite(delta)]
      if (length(delta) == 0) {
        next
      }
      set.seed(seed + index)
      boot <- replicate(
        as.integer(bootstrap),
        mean(sample(delta, replace = TRUE))
      )
      interval <- stats::quantile(
        boot,
        c(0.025, 0.975),
        names = FALSE,
        na.rm = TRUE
      )
      result[[index]] <- data.frame(
        contrast = contrast,
        metric_name = metric,
        estimate = mean(delta),
        ci_low = interval[1],
        ci_high = interval[2],
        n_pairs = length(delta),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  if (length(result) == 0) data.frame() else do.call(rbind, result)
}


# General paired-fold contrast reused for Two-stage minus One-stage in Experiment 3.
.ablation_two_group_contrasts <- function(
    metrics,
    first_group,
    second_group,
    bootstrap,
    seed
) {
  data <- metrics[
    metrics$group_id %in% c(first_group, second_group) &
      is.finite(metrics$metric_value),
    ,
    drop = FALSE
  ]
  data <- stats::aggregate(
    metric_value ~ group_id + fold_id + metric_name,
    data = data,
    FUN = mean
  )
  result <- lapply(seq_along(unique(data$metric_name)), function(i) {
    metric <- unique(data$metric_name)[i]
    metric_data <- data[data$metric_name == metric, , drop = FALSE]
    wide <- reshape(
      metric_data[, c("fold_id", "group_id", "metric_value")],
      idvar = "fold_id",
      timevar = "group_id",
      direction = "wide"
    )
    first <- paste0("metric_value.", first_group)
    second <- paste0("metric_value.", second_group)
    if (!all(c(first, second) %in% colnames(wide))) {
      return(NULL)
    }
    delta <- wide[[first]] - wide[[second]]
    delta <- delta[is.finite(delta)]
    if (length(delta) == 0) {
      return(NULL)
    }
    set.seed(seed + i)
    boot <- replicate(
      as.integer(bootstrap),
      mean(sample(delta, replace = TRUE))
    )
    interval <- stats::quantile(
      boot,
      c(0.025, 0.975),
      names = FALSE,
      na.rm = TRUE
    )
    data.frame(
      contrast = paste(first_group, second_group, sep = "-"),
      metric_name = metric,
      estimate = mean(delta),
      ci_low = interval[1],
      ci_high = interval[2],
      n_pairs = length(delta),
      stringsAsFactors = FALSE
    )
  })
  result <- result[!vapply(result, is.null, logical(1))]
  if (length(result) == 0) data.frame() else do.call(rbind, result)
}


# Attach input hashes and random seeds to each long-format metric row for a uniform audit table.
.ablation_add_audit_hashes <- function(metrics, manifest) {
  cbind(
    metrics[, c(
      "experiment_id", "group_id", "fold_id", "held_out_cohort",
      "seed_type", "stochastic_seed", "rank_q", "k", "distance",
      "metric_name", "metric_value"
    ), drop = FALSE],
    sample_manifest_hash = manifest$sample_manifest_hash,
    feature_manifest_hash = manifest$feature_manifest_hash,
    module_sequence_id = NA_character_,
    module_count = manifest$module_count,
    projection_seed = ifelse(metrics$seed_type == "projection_seed", metrics$stochastic_seed, NA),
    permutation_seed = ifelse(metrics$seed_type == "permutation_seed", metrics$stochastic_seed, NA),
    umap_seed = NA_real_,
    stringsAsFactors = FALSE
  )
}


# Gate 1 requires Cohort to beat Direct and both nulls on the primary metric without
# materially degrading biological purity or cohort mixing.
.ablation_gate_one <- function(cohort_result, gate) {
  summary <- cohort_result$summary
  contrasts <- cohort_result$contrasts
  get_value <- function(group, metric) {
    index <- which(summary$group_id == group & summary$metric_name == metric)
    value <- summary$estimate[index]
    if (length(value) == 0) NA_real_ else value[1]
  }
  max_finite <- function(...) {
    value <- c(...)
    value <- value[is.finite(value)]
    if (length(value) == 0) NA_real_ else max(value)
  }
  get_contrast <- function(contrast, metric, column = "ci_low") {
    index <- which(
      contrasts$contrast == contrast & contrasts$metric_name == metric
    )
    value <- contrasts[[column]][index]
    if (length(value) == 0) NA_real_ else value[1]
  }
  primary <- gate$primary_metric
  cohort_primary <- get_value("Cohort", primary)
  direct_primary <- get_value("Direct", primary)
  null_primary <- max_finite(
    get_value("Null-RP", primary),
    get_value("Null-Perm", primary)
  )
  primary_gains <- c(
    get_contrast("Cohort-Direct", primary),
    get_contrast("Cohort-Null-RP", primary),
    get_contrast("Cohort-Null-Perm", primary)
  )
  # If the probe is undefined in small or sparse classes, fall back to biological purity.
  if (!all(is.finite(primary_gains)) && primary == "balanced_accuracy") {
    primary <- "biology_purity"
    cohort_primary <- get_value("Cohort", primary)
    direct_primary <- get_value("Direct", primary)
    null_primary <- max_finite(
      get_value("Null-RP", primary),
      get_value("Null-Perm", primary)
    )
    primary_gains <- c(
      get_contrast("Cohort-Direct", primary),
      get_contrast("Cohort-Null-RP", primary),
      get_contrast("Cohort-Null-Perm", primary)
    )
  }
  purity_gain <- get_contrast("Cohort-Direct", "biology_purity")
  mixing_gain <- get_contrast("Cohort-Direct", "cohort_mixing")
  purity_ok <- is.finite(purity_gain) &&
    purity_gain >= -gate$purity_tolerance
  mixing_ok <- is.finite(mixing_gain) &&
    mixing_gain >= -gate$mixing_tolerance
  gain_ok <- all(is.finite(primary_gains)) &&
    all(primary_gains >= gate$min_gain)
  list(
    pass = isTRUE(gain_ok && purity_ok && mixing_ok),
    primary_metric = primary,
    cohort = cohort_primary,
    direct = direct_primary,
    strongest_null = null_primary,
    paired_ci_lower = stats::setNames(
      primary_gains,
      c("Cohort-Direct", "Cohort-Null-RP", "Cohort-Null-Perm")
    ),
    gain_ok = gain_ok,
    purity_ok = purity_ok,
    mixing_ok = mixing_ok
  )
}


# Bind audit tables by their union of columns and fill unavailable fields explicitly with NA.
.ablation_rbind <- function(parts) {
  parts <- parts[vapply(parts, nrow, integer(1)) > 0]
  if (length(parts) == 0) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(data) {
    missing <- setdiff(all_names, names(data))
    for (name in missing) data[[name]] <- NA
    data[, all_names, drop = FALSE]
  })
  do.call(rbind, parts)
}


# Experiment 2 follows tissue-stratified nested sequences while adding cohort modules.
# Keeping it independent makes future module-count ablations easy to extend.
.ablation_experiment_scaling <- function(prepared, config, manifest, seed, verbose) {
  sequences <- .ablation_nested_module_sequences(
    prepared$module_manifest$modules,
    config$scaling$sequences,
    seed
  )
  counts <- sort(unique(pmin(config$scaling$counts, manifest$module_count)))
  folds <- .ablation_grouped_folds(
    prepared$metadata$cohort,
    config$general$n_folds,
    seed
  )
  metrics <- list()
  index <- 1L

  # Larger counts strictly contain smaller counts within a sequence, enabling paired increments.
  for (sequence_id in seq_along(sequences)) {
    order <- sequences[[sequence_id]]
    for (module_count in counts) {
      module_ids <- order[seq_len(module_count)]
      columns <- unlist(prepared$module_manifest$blocks[module_ids], use.names = FALSE)
      candidate <- prepared$d1[, columns, drop = FALSE]
      for (fold in sort(unique(folds))) {
        test_rows <- folds == fold
        train_rows <- !test_rows
        rank_q <- min(config$general$rank, ncol(candidate), sum(train_rows) - 1L)
        full_fit <- .ablation_fit_pca(
          prepared$d1[train_rows, , drop = FALSE],
          prepared$d1[test_rows, , drop = FALSE],
          rank_q
        )
        subset_fit <- .ablation_fit_pca(
          candidate[train_rows, , drop = FALSE],
          candidate[test_rows, , drop = FALSE],
          rank_q
        )
        # Use the same rank_q for full and subset fits so dimension does not drive geometry.
        metadata_test <- prepared$metadata[test_rows, , drop = FALSE]
        metadata_train <- prepared$metadata[train_rows, , drop = FALSE]
        local_k <- min(config$general$k, sum(test_rows) - 1L)
        mixing <- .ablation_mixing_purity(subset_fit$test, metadata_test, local_k)
        probe <- if (config$general$probe) {
          probe_metadata <- rbind(metadata_train, metadata_test)
          probe_label <- .ablation_eligible_probe_labels(
            probe_metadata,
            config$general$probe_label
          )
          train_label <- probe_label[seq_len(nrow(metadata_train))]
          test_label <- probe_label[
            nrow(metadata_train) + seq_len(nrow(metadata_test))
          ]
          .ablation_probe(
            train = subset_fit$train,
            test = subset_fit$test,
            train_label = train_label,
            test_label = test_label,
            seed = seed + sequence_id + fold + module_count,
            config = config
          )
        } else {
          c(macro_auroc = NA_real_, balanced_accuracy = NA_real_)
        }
        values <- c(
          cka_to_full = .ablation_linear_cka(full_fit$test, subset_fit$test),
          knn_to_full = .ablation_knn_jaccard(full_fit$test, subset_fit$test, local_k),
          effective_rank = .ablation_effective_rank(subset_fit$test),
          cohort_mixing = mixing$cohort_mixing,
          biology_purity = mixing$biology_purity,
          probe
        )
        metrics[[index]] <- data.frame(
          experiment_id = "scaling",
          group_id = paste0(module_count, "-modules"),
          fold_id = as.character(fold),
          held_out_cohort = paste(unique(metadata_test$cohort), collapse = ";"),
          seed_type = "module_sequence",
          stochastic_seed = seed + sequence_id,
          rank_q = subset_fit$rank,
          k = local_k,
          distance = config$general$distance,
          metric_name = names(values),
          metric_value = as.numeric(values),
          module_sequence_id = sprintf("S%03d", sequence_id),
          module_count = module_count,
          cumulative_dimension = ncol(candidate),
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
  }
  metrics <- do.call(rbind, metrics)
  summary <- .ablation_scaling_summary(metrics, config$general$bootstrap, seed)
  d1_audit <- cbind(
    metrics,
    sample_manifest_hash = manifest$sample_manifest_hash,
    feature_manifest_hash = manifest$feature_manifest_hash,
    projection_seed = NA_real_,
    permutation_seed = NA_real_,
    umap_seed = NA_real_,
    stringsAsFactors = FALSE
  )
  # Beyond the d1 curve, test whether downstream d3 and clustering stabilize at representative counts.
  embedding <- .ablation_scaling_embedding(
    prepared = prepared,
    sequences = sequences,
    config = config,
    manifest = manifest,
    seed = seed + 5000L,
    verbose = verbose
  )
  audit <- .ablation_rbind(list(d1 = d1_audit, embedding = embedding$audit))
  list(
    status = "complete",
    sequences = sequences,
    metrics = metrics,
    summary = summary,
    embedding = embedding,
    audit = audit
  )
}


# Generate a tissue-stratified interleaved order for each repeat. Prefixes cover tissues
# broadly while randomizing the entry order of cohorts within each tissue.
.ablation_nested_module_sequences <- function(modules, n_sequences, seed) {
  lapply(seq_len(as.integer(n_sequences)), function(i) {
    set.seed(seed + i)
    strata <- split(modules$module_id, modules$tissue)
    strata <- lapply(strata, sample)
    module_ids <- unlist(strata, use.names = FALSE)
    keys <- unlist(lapply(strata, function(x) {
      (seq_along(x) - stats::runif(length(x))) / length(x)
    }), use.names = FALSE)
    names(keys) <- module_ids
    names(sort(keys))
  })
}


# Summarize the mean scaling curve and adjacent count increments within matched sequence x fold.
.ablation_scaling_summary <- function(metrics, bootstrap, seed) {
  base <- stats::aggregate(
    metric_value ~ group_id + module_count + cumulative_dimension + metric_name,
    data = metrics,
    FUN = mean,
    na.rm = TRUE
  )
  base <- base[order(base$metric_name, base$module_count), , drop = FALSE]
  increments <- do.call(rbind, lapply(split(metrics, metrics$metric_name), function(data) {
    wide <- reshape(
      data[, c("module_sequence_id", "fold_id", "module_count", "metric_value")],
      idvar = c("module_sequence_id", "fold_id"),
      timevar = "module_count",
      direction = "wide"
    )
    counts <- sort(unique(data$module_count))
    if (length(counts) < 2) return(NULL)
    # Adjacent differences preserve nested pairing, so the CI reflects marginal module gain.
    do.call(rbind, lapply(seq_len(length(counts) - 1L), function(i) {
      from <- counts[i]
      to <- counts[i + 1L]
      delta <- wide[[paste0("metric_value.", to)]] - wide[[paste0("metric_value.", from)]]
      delta <- delta[is.finite(delta)]
      if (length(delta) == 0) {
        return(data.frame(
          metric_name = data$metric_name[1],
          from = from,
          to = to,
          estimate = NA_real_,
          ci_low = NA_real_,
          ci_high = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      set.seed(seed + from * 1000L + to + sum(utf8ToInt(data$metric_name[1])))
      boot <- replicate(
        as.integer(bootstrap),
        mean(sample(delta, replace = TRUE))
      )
      interval <- stats::quantile(
        boot,
        c(0.025, 0.975),
        names = FALSE,
        na.rm = TRUE
      )
      data.frame(
        metric_name = data$metric_name[1],
        from = from,
        to = to,
        estimate = mean(delta, na.rm = TRUE),
        ci_low = interval[1],
        ci_high = interval[2],
        stringsAsFactors = FALSE
      )
    }))
  }))
  saturation <- .ablation_saturation_fit(base)
  list(
    curve = base,
    paired_increments = increments,
    saturation = saturation
  )
}


# Fit S(m) = S_inf - a * exp(-m / tau), where smaller tau indicates earlier saturation.
# Grid-search tau and solve S_inf/a linearly to avoid unstable nonlinear optimization.
.ablation_saturation_fit <- function(curve) {
  result <- lapply(split(curve, curve$metric_name), function(data) {
    data <- data[is.finite(data$metric_value), , drop = FALSE]
    data <- data[order(data$module_count), , drop = FALSE]
    if (nrow(data) < 3 || length(unique(data$module_count)) < 3) {
      return(NULL)
    }
    module_count <- data$module_count
    value <- data$metric_value
    tau_grid <- exp(seq(
      log(max(1, min(module_count) / 5)),
      log(max(module_count) * 10),
      length.out = 200
    ))
    fits <- lapply(tau_grid, function(tau) {
      design <- cbind(1, exp(-module_count / tau))
      coefficients <- qr.solve(design, value)
      fitted <- as.vector(design %*% coefficients)
      c(
        tau = tau,
        s_inf = coefficients[1],
        a = -coefficients[2],
        rmse = sqrt(mean((value - fitted)^2))
      )
    })
    fits <- do.call(rbind, fits)
    best <- fits[which.min(fits[, "rmse"]), ]
    data.frame(
      metric_name = data$metric_name[1],
      s_inf = best["s_inf"],
      a = best["a"],
      tau = best["tau"],
      rmse = best["rmse"],
      stringsAsFactors = FALSE
    )
  })
  result <- result[!vapply(result, is.null, logical(1))]
  if (length(result) == 0) data.frame() else do.call(rbind, result)
}


# Repeat full Two-stage reduction and clustering at representative module counts to test
# whether d1 scaling extends downstream. The same seed reuses one stratified sample subset.
.ablation_scaling_embedding <- function(
    prepared,
    sequences,
    config,
    manifest,
    seed,
    verbose
) {
  sequence_count <- min(
    as.integer(config$scaling$embedding_sequences),
    length(sequences)
  )
  selected_sequences <- sequences[seq_len(sequence_count)]
  counts <- sort(unique(pmin(
    config$scaling$embedding_counts,
    manifest$module_count
  )))
  metrics <- list()
  embeddings <- list()
  index <- 1L

  # sequence x module_count x seed is the repeat unit; retain embeddings for paired stability.
  for (sequence_id in seq_along(selected_sequences)) {
    order <- selected_sequences[[sequence_id]]
    for (module_count in counts) {
      module_ids <- order[seq_len(module_count)]
      columns <- unlist(
        prepared$module_manifest$blocks[module_ids],
        use.names = FALSE
      )
      for (seed_i in config$scaling$embedding_seeds) {
        rows <- .ablation_tissue_subsample(
          prepared$metadata,
          config$scaling$subsample_fraction,
          seed_i
        )
        d1 <- prepared$d1[rows, columns, drop = FALSE]
        metadata <- prepared$metadata[rows, , drop = FALSE]
        d3 <- .ablation_two_stage_embedding(d1, config$general$dr, seed_i)
        clusters <- .ablation_dbscan(d3, config$general$cluster)
        values <- .ablation_embedding_metrics(
          high = d1,
          low = d3,
          metadata = metadata,
          clusters = clusters,
          config = config,
          seed = seed_i
        )
        group_id <- paste0(module_count, "-modules")
        sequence_name <- sprintf("S%03d", sequence_id)
        metrics[[index]] <- data.frame(
          experiment_id = "scaling_embedding",
          group_id = group_id,
          fold_id = as.character(seed_i),
          held_out_cohort = NA_character_,
          seed_type = "umap_seed",
          stochastic_seed = seed_i,
          rank_q = ncol(d3),
          k = min(config$general$k, nrow(d3) - 1L),
          distance = config$general$distance,
          metric_name = names(values),
          metric_value = as.numeric(values),
          module_sequence_id = sequence_name,
          module_count = module_count,
          cumulative_dimension = ncol(d1),
          stringsAsFactors = FALSE
        )
        key <- paste(sequence_name, group_id, seed_i, sep = "|")
        embeddings[[key]] <- list(
          sample_id = metadata$sample_id,
          d3 = d3,
          cluster = clusters
        )
        index <- index + 1L
      }
    }
  }

  metrics <- do.call(rbind, metrics)
  stability <- .ablation_scaling_embedding_stability(embeddings, config$general$k)
  summary <- .ablation_metric_summary(metrics, config$general$bootstrap, seed)
  audit <- cbind(
    metrics,
    sample_manifest_hash = manifest$sample_manifest_hash,
    feature_manifest_hash = manifest$feature_manifest_hash,
    projection_seed = NA_real_,
    permutation_seed = NA_real_,
    umap_seed = metrics$stochastic_seed,
    stringsAsFactors = FALSE
  )
  list(
    metrics = metrics,
    stability = stability,
    summary = summary,
    audit = audit
  )
}


# Compare d3 neighborhoods and clusters across seed pairs at the same sequence and module count.
.ablation_scaling_embedding_stability <- function(embeddings, k) {
  keys <- names(embeddings)
  groups <- split(keys, sub("\\|[^|]+$", "", keys))
  do.call(rbind, lapply(names(groups), function(group) {
    group_keys <- groups[[group]]
    if (length(group_keys) < 2) return(NULL)
    pairs <- utils::combn(group_keys, 2, simplify = FALSE)
    do.call(rbind, lapply(pairs, function(pair) {
      first <- embeddings[[pair[1]]]
      second <- embeddings[[pair[2]]]
      common <- intersect(first$sample_id, second$sample_id)
      first_rows <- match(common, first$sample_id)
      second_rows <- match(common, second$sample_id)
      local_k <- min(k, length(common) - 1L)
      group_parts <- strsplit(group, "|", fixed = TRUE)[[1]]
      data.frame(
        module_sequence_id = group_parts[1],
        group_id = group_parts[2],
        seed_a = sub("^.*\\|", "", pair[1]),
        seed_b = sub("^.*\\|", "", pair[2]),
        neighbor_jaccard = .ablation_knn_jaccard(
          first$d3[first_rows, , drop = FALSE],
          second$d3[second_rows, , drop = FALSE],
          local_k
        ),
        ari = .ablation_ari(
          first$cluster[first_rows],
          second$cluster[second_rows]
        ),
        cluster_jaccard = .ablation_cluster_jaccard(
          first$cluster[first_rows],
          second$cluster[second_rows]
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))
}


# Experiment 3 reuses CCS reduction and ablates only the tissue-first layer: Two-stage
# reduces each tissue to d2 before d3, whereas One-stage maps full d1 directly to d3.
.ablation_experiment_tissue_first <- function(prepared, config, manifest, seed, verbose) {
  seeds <- config$tissue_first$seeds
  results <- list()
  metrics <- list()
  stratified <- list()
  index <- 1L
  # Both arms share the same sample subset and seed within each repeat for strict pairing.
  for (i in seq_along(seeds)) {
    seed_i <- seeds[i]
    rows <- .ablation_tissue_subsample(
      prepared$metadata,
      config$tissue_first$subsample_fraction,
      seed_i
    )
    d1 <- prepared$d1[rows, , drop = FALSE]
    metadata <- prepared$metadata[rows, , drop = FALSE]
    embeddings <- .ablation_tissue_embeddings(d1, config$general$dr, seed_i)
    for (group in names(embeddings)) {
      d3 <- embeddings[[group]]
      clusters <- .ablation_dbscan(d3, config$general$cluster)
      values <- .ablation_embedding_metrics(
        high = d1,
        low = d3,
        metadata = metadata,
        clusters = clusters,
        config = config,
        seed = seed_i
      )
      metrics[[index]] <- data.frame(
        experiment_id = "tissue_first",
        group_id = group,
        fold_id = as.character(seed_i),
        held_out_cohort = NA_character_,
        seed_type = "umap_seed",
        stochastic_seed = seed_i,
        rank_q = ncol(d3),
        k = min(config$general$k, nrow(d3) - 1L),
        distance = config$general$distance,
        metric_name = names(values),
        metric_value = as.numeric(values),
        stringsAsFactors = FALSE
      )
      stratified[[index]] <- .ablation_tissue_stratified_metrics(
        high = d1,
        low = d3,
        metadata = metadata,
        group_id = group,
        seed = seed_i,
        k = config$general$k
      )
      results[[paste(group, seed_i, sep = "|")]] <- list(
        sample_id = metadata$sample_id,
        d3 = d3,
        cluster = clusters
      )
      index <- index + 1L
    }
  }
  metrics <- do.call(rbind, metrics)
  stratified <- do.call(rbind, stratified)
  quartile <- stats::quantile(
    stratified$sample_count,
    0.25,
    names = FALSE
  )
  # Flag bottom-quartile tissues to test whether tissue-first benefits only large tissues.
  stratified$small_tissue <- stratified$sample_count <= quartile
  stability <- .ablation_embedding_stability(results, config$general$k)
  summary <- .ablation_metric_summary(metrics, config$general$bootstrap, seed)
  contrasts <- .ablation_two_group_contrasts(
    metrics,
    first_group = "Two-stage",
    second_group = "One-stage",
    bootstrap = config$general$bootstrap,
    seed = seed
  )
  audit <- .ablation_add_audit_hashes(metrics, manifest)
  audit$umap_seed <- audit$stochastic_seed
  list(
    status = "complete",
    metrics = metrics,
    stratified = stratified,
    stability = stability,
    summary = summary,
    contrasts = contrasts,
    audit = audit
  )
}


# Draw the requested fraction across tissue x cohort strata; retain all rows when fraction >= 1.
.ablation_tissue_subsample <- function(metadata, fraction, seed) {
  if (fraction >= 1) {
    return(seq_len(nrow(metadata)))
  }
  size <- max(3L, floor(nrow(metadata) * fraction))
  .ablation_stratified_sample(metadata, size, seed)
}


# Generate paired Two-stage and One-stage embeddings from the same d1 and seed.
.ablation_tissue_embeddings <- function(d1, dr_config, seed) {
  list(
    `Two-stage` = .ablation_two_stage_embedding(d1, dr_config, seed),
    `One-stage` = .ablation_one_stage_embedding(d1, dr_config, seed)
  )
}


# Two-stage first reduces each tissue block independently, then maps the combined d2
# to the final dimensions; both stages reuse the existing CCS reduction logic.
.ablation_two_stage_embedding <- function(d1, dr_config, seed) {
  dr_fun <- getFromNamespace("drCCSProbability", "CCS")
  reference <- vapply(
    strsplit(colnames(d1), "|", fixed = TRUE),
    `[`,
    character(1),
    1
  )
  dr_args <- dr_config[setdiff(names(dr_config), c("method", "dimension"))]
  two_d2 <- .ablation_reduce_by_reference(
    data = d1,
    method = dr_config$method,
    reference = reference,
    dims = dr_config$dimension[1],
    seed = seed,
    dr_args = dr_args
  )
  two_d3 <- do.call(dr_fun, c(list(
    data = two_d2,
    method = dr_config$method,
    reference = NULL,
    dims = dr_config$dimension[2],
    seed = seed,
    verbose = FALSE
  ), dr_args))
  two_d3
}


# Apply CCS preprocessing and CORE_DR per reference block (tissue here), then repair
# deduplicated samples. Each block gets a derived seed and returns in original row order.
.ablation_reduce_by_reference <- function(
    data,
    method,
    reference,
    dims,
    seed,
    dr_args
) {
  data_for_dr <- getFromNamespace("data_for_dr", "CCS")
  core_dr <- getFromNamespace("CORE_DR", "CCS")
  repair <- getFromNamespace("repairCCS", "CCS")
  reference_levels <- unique(reference)
  set.seed(seed)
  reference_seeds <- sample(
    1:10000,
    length(reference_levels),
    replace = FALSE
  )

  reduced <- lapply(seq_along(reference_levels), function(i) {
    reference_i <- reference_levels[i]
    block <- data[, reference == reference_i, drop = FALSE]
    prepared <- data_for_dr(block, rm.dup.col = FALSE, verbose = FALSE)
    cleaned <- prepared$cleaned$data
    # Reduce target dimensions for small/low-rank blocks and cap UMAP neighbors accordingly.
    available_dims <- min(ncol(cleaned), nrow(cleaned) - 2L)
    dims_i <- min(as.integer(dims), available_dims)
    if (dims_i < 1) {
      stop(
        "ablation: tissue block ", reference_i,
        " has insufficient rank for dimension reduction.",
        call. = FALSE
      )
    }
    args_i <- dr_args
    if (!is.null(args_i$n_neighbors)) {
      args_i$n_neighbors <- min(args_i$n_neighbors, nrow(cleaned) - 1L)
    }
    embedding <- do.call(core_dr, c(list(
      method = method,
      data = cleaned,
      dims = dims_i,
      seed = reference_seeds[i]
    ), args_i))
    rownames(embedding) <- names(prepared$cleaned$md5)
    embedding <- repair(embedding, prepared)
    colnames(embedding) <- paste0(reference_i, "|D", seq_len(dims_i))
    embedding[rownames(data), , drop = FALSE]
  })
  do.call(cbind, reduced)
}


# One-stage control: skip tissue blocks and reduce the complete d1 directly to final dimensions.
.ablation_one_stage_embedding <- function(d1, dr_config, seed) {
  dr_fun <- getFromNamespace("drCCSProbability", "CCS")
  dr_args <- dr_config[setdiff(names(dr_config), c("method", "dimension"))]
  do.call(dr_fun, c(list(
    data = d1,
    method = dr_config$method,
    reference = NULL,
    dims = dr_config$dimension[2],
    seed = seed,
    verbose = FALSE
  ), dr_args))
}


# Run DBSCAN on column-standardized d3; label 0 denotes noise by dbscan convention.
.ablation_dbscan <- function(d3, cluster_config) {
  fit <- do.call(
    dbscan::dbscan,
    c(list(x = scale(d3)), cluster_config)
  )
  stats::setNames(fit$cluster, rownames(d3))
}


# Evaluate high/low-dimensional fidelity, within-tissue retention, cohort/biology, and clusters.
.ablation_embedding_metrics <- function(high, low, metadata, clusters, config, seed) {
  rows <- seq_len(nrow(high))
  if (length(rows) > config$general$fidelity_samples) {
    rows <- .ablation_stratified_sample(
      metadata,
      config$general$fidelity_samples,
      seed
    )
  }
  fidelity <- .ablation_trust_continuity(
    high[rows, , drop = FALSE],
    low[rows, , drop = FALSE],
    min(config$general$k, length(rows) - 1L)
  )
  local_k <- min(config$general$k, nrow(low) - 1L)
  mixing <- .ablation_mixing_purity(low, metadata, local_k)
  tissue_retention <- .ablation_tissue_knn_retention(high, low, metadata$tissue, local_k)
  c(
    trustworthiness = fidelity["trustworthiness"],
    continuity = fidelity["continuity"],
    tissue_knn_retention = tissue_retention,
    cohort_mixing = mixing$cohort_mixing,
    biology_purity = mixing$biology_purity,
    cluster_count = length(setdiff(unique(clusters), 0)),
    cluster_size_entropy = .ablation_cluster_size_entropy(clusters),
    noise_fraction = mean(clusters == 0)
  )
}


# Compute Shannon entropy of non-noise cluster sizes; at fixed cluster count, higher is more balanced.
.ablation_cluster_size_entropy <- function(clusters) {
  sizes <- table(clusters[clusters != 0])
  if (length(sizes) == 0) {
    return(NA_real_)
  }
  probability <- as.numeric(sizes) / sum(sizes)
  -sum(probability * log(probability))
}


# Compute retention, mixing, and purity per tissue so global means do not mask small-tissue loss.
.ablation_tissue_stratified_metrics <- function(
    high,
    low,
    metadata,
    group_id,
    seed,
    k
) {
  parts <- lapply(split(seq_len(nrow(metadata)), metadata$tissue), function(rows) {
    if (length(rows) < 2) {
      return(NULL)
    }
    local_k <- min(k, length(rows) - 1L)
    local_metadata <- metadata[rows, , drop = FALSE]
    mixing <- .ablation_mixing_purity(
      low[rows, , drop = FALSE],
      local_metadata,
      local_k
    )
    data.frame(
      group_id = group_id,
      umap_seed = seed,
      tissue = local_metadata$tissue[1],
      sample_count = length(rows),
      cohort_count = length(unique(local_metadata$cohort)),
      knn_retention = .ablation_knn_jaccard(
        high[rows, , drop = FALSE],
        low[rows, , drop = FALSE],
        local_k
      ),
      cohort_mixing = mixing$cohort_mixing,
      biology_purity = mixing$biology_purity,
      stringsAsFactors = FALSE
    )
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) data.frame() else do.call(rbind, parts)
}


# Compute classical trustworthiness and continuity. The former penalizes false low-dimensional
# neighbors, the latter missing high-dimensional neighbors; values closer to 1 are better.
.ablation_trust_continuity <- function(high, low, k) {
  high_distance <- as.matrix(stats::dist(scale(high)))
  low_distance <- as.matrix(stats::dist(scale(low)))
  diag(high_distance) <- Inf
  diag(low_distance) <- Inf
  high_rank <- t(apply(high_distance, 1, rank, ties.method = "average"))
  low_rank <- t(apply(low_distance, 1, rank, ties.method = "average"))
  high_nn <- t(apply(high_distance, 1, function(x) order(x)[seq_len(k)]))
  low_nn <- t(apply(low_distance, 1, function(x) order(x)[seq_len(k)]))
  n <- nrow(high)
  # Penalize intruders by high-dimensional ranks and missing neighbors by low-dimensional ranks.
  penalty_t <- sum(vapply(seq_len(n), function(i) {
    intruders <- setdiff(low_nn[i, ], high_nn[i, ])
    sum(high_rank[i, intruders] - k)
  }, numeric(1)))
  penalty_c <- sum(vapply(seq_len(n), function(i) {
    missing <- setdiff(high_nn[i, ], low_nn[i, ])
    sum(low_rank[i, missing] - k)
  }, numeric(1)))
  normalizer <- 2 / (n * k * (2 * n - 3 * k - 1))
  c(
    trustworthiness = 1 - normalizer * penalty_t,
    continuity = 1 - normalizer * penalty_c
  )
}


# Compare high- and low-dimensional kNN overlap within each tissue, then average over samples.
.ablation_tissue_knn_retention <- function(high, low, tissue, k) {
  values <- unlist(lapply(split(seq_len(nrow(high)), tissue), function(rows) {
    if (length(rows) < 2) return(NA_real_)
    local_k <- min(k, length(rows) - 1L)
    high_nn <- .ablation_knn(high[rows, , drop = FALSE], local_k)
    low_nn <- .ablation_knn(low[rows, , drop = FALSE], local_k)
    vapply(seq_along(rows), function(i) {
      length(intersect(high_nn[i, ], low_nn[i, ])) / local_k
    }, numeric(1))
  }))
  mean(values, na.rm = TRUE)
}


# Compare seed pairs within each arm after aligning common samples, using neighborhood
# Jaccard, ARI, and cluster-set Jaccard.
.ablation_embedding_stability <- function(results, k) {
  groups <- split(names(results), sub("\\|.*$", "", names(results)))
  do.call(rbind, lapply(names(groups), function(group) {
    keys <- groups[[group]]
    if (length(keys) < 2) return(NULL)
    pairs <- utils::combn(keys, 2, simplify = FALSE)
    values <- lapply(pairs, function(pair) {
      first <- results[[pair[1]]]
      second <- results[[pair[2]]]
      common <- intersect(first$sample_id, second$sample_id)
      first_rows <- match(common, first$sample_id)
      second_rows <- match(common, second$sample_id)
      local_k <- min(k, length(common) - 1L)
      data.frame(
        group_id = group,
        seed_a = sub("^.*\\|", "", pair[1]),
        seed_b = sub("^.*\\|", "", pair[2]),
        neighbor_jaccard = .ablation_knn_jaccard(
          first$d3[first_rows, , drop = FALSE],
          second$d3[second_rows, , drop = FALSE],
          local_k
        ),
        ari = .ablation_ari(first$cluster[first_rows], second$cluster[second_rows]),
        cluster_jaccard = .ablation_cluster_jaccard(
          first$cluster[first_rows],
          second$cluster[second_rows]
        ),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, values)
  }))
}


# Compute adjusted Rand index directly from a contingency table, correcting chance agreement.
.ablation_ari <- function(x, y) {
  table_xy <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  if (sum(table_xy) < 2) {
    return(NA_real_)
  }
  index <- sum(choose2(table_xy))
  row_index <- sum(choose2(rowSums(table_xy)))
  col_index <- sum(choose2(colSums(table_xy)))
  total <- choose2(sum(table_xy))
  expected <- row_index * col_index / total
  maximum <- (row_index + col_index) / 2
  if (!is.finite(expected) || !is.finite(maximum)) {
    return(NA_real_)
  }
  if (maximum == expected) 1 else (index - expected) / (maximum - expected)
}


# Match every non-noise cluster to its best Jaccard counterpart and average both directions,
# reducing sensitivity to arbitrary labels and asymmetric containment.
.ablation_cluster_jaccard <- function(x, y) {
  one_way <- function(a, b) {
    clusters <- setdiff(unique(a), 0)
    if (length(clusters) == 0) return(NA_real_)
    mean(vapply(clusters, function(cluster) {
      a_rows <- which(a == cluster)
      candidates <- setdiff(unique(b[a_rows]), 0)
      if (length(candidates) == 0) return(0)
      max(vapply(candidates, function(other) {
        b_rows <- which(b == other)
        length(intersect(a_rows, b_rows)) / length(union(a_rows, b_rows))
      }, numeric(1)))
    }, numeric(1)))
  }
  values <- c(one_way(x, y), one_way(y, x))
  if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
}


# Expand the fixed configuration or a shared DR x DBSCAN grid into stable IDs.
.ablation_metaccs_parameter_manifest <- function(config) {
  mode <- config$metaccs$parameter_mode
  dr_sets <- .ablation_parameter_sets(
    base = config$general$dr,
    grid = if (identical(mode, "grid")) config$metaccs$dr_grid else NULL,
    label = "dr_grid"
  )
  cluster_sets <- .ablation_parameter_sets(
    base = config$general$cluster,
    grid = if (identical(mode, "grid")) config$metaccs$cluster_grid else NULL,
    label = "cluster_grid"
  )
  pairs <- expand.grid(
    dr_index = seq_along(dr_sets),
    cluster_index = seq_along(cluster_sets),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dr_ids <- vapply(dr_sets, function(x) {
    paste0("dr-", substr(digest::digest(x, algo = "md5"), 1, 8))
  }, character(1))
  cluster_ids <- vapply(cluster_sets, function(x) {
    paste0("cluster-", substr(digest::digest(x, algo = "md5"), 1, 8))
  }, character(1))
  .ablation_validate_parameter_sets(dr_sets, cluster_sets)
  manifest <- data.frame(
    parameter_mode = mode,
    direct_feature_mode = config$metaccs$direct_feature_mode,
    dr_param_set_id = dr_ids[pairs$dr_index],
    cluster_param_set_id = cluster_ids[pairs$cluster_index],
    dr = I(dr_sets[pairs$dr_index]),
    cluster = I(cluster_sets[pairs$cluster_index]),
    stringsAsFactors = FALSE
  )
  manifest[!duplicated(manifest[c("dr_param_set_id", "cluster_param_set_id")]), , drop = FALSE]
}


# Merge each grid row onto the shared base so both arms receive identical budgets.
.ablation_parameter_sets <- function(base, grid, label) {
  if (is.null(grid)) {
    return(list(base))
  }
  if (is.data.frame(grid)) {
    grid <- lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, , drop = FALSE]))
  }
  if (!is.list(grid) || length(grid) == 0 ||
      !all(vapply(grid, is.list, logical(1)))) {
    stop("ablation: params$metaccs$", label, " must be a non-empty list of lists.", call. = FALSE)
  }
  lapply(grid, function(x) {
    unknown <- setdiff(names(x), names(base))
    if (length(unknown) > 0) {
      stop(
        "ablation: unknown ", label, " fields: ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    .ablation_merge_lists(base, x)
  })
}


.ablation_validate_parameter_sets <- function(dr_sets, cluster_sets) {
  for (dr in dr_sets) {
    if (length(dr$dimension) != 2 || any(!is.finite(dr$dimension)) ||
        any(dr$dimension < 1)) {
      stop("ablation: every metaccs DR set requires two positive dimensions.", call. = FALSE)
    }
    if (!is.null(dr$n_neighbors) &&
        (length(dr$n_neighbors) != 1 || !is.finite(dr$n_neighbors) ||
          dr$n_neighbors < 2)) {
      stop("ablation: every metaccs n_neighbors must be at least 2.", call. = FALSE)
    }
    if (!is.null(dr$n_threads) &&
        (length(dr$n_threads) != 1 || !is.finite(dr$n_threads) ||
          dr$n_threads < 1)) {
      stop("ablation: every metaccs n_threads must be positive or NULL.", call. = FALSE)
    }
  }
  for (cluster in cluster_sets) {
    if (length(cluster$eps) != 1 || !is.finite(cluster$eps) || cluster$eps <= 0 ||
        length(cluster$minPts) != 1 || !is.finite(cluster$minPts) ||
        cluster$minPts < 1) {
      stop("ablation: every metaccs DBSCAN set requires eps > 0 and minPts >= 1.", call. = FALSE)
    }
  }
  invisible(TRUE)
}


# Build Direct-GSClassifier blocks from each tissue's frozen feature union.
.ablation_direct_tissue_blocks <- function(direct, tissue_features) {
  missing <- setdiff(unique(unlist(tissue_features, use.names = FALSE)), colnames(direct))
  if (length(missing) > 0) {
    stop(
      "ablation: Direct-GSClassifier lacks frozen tissue features: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  blocks <- lapply(tissue_features, function(features) {
    direct[, features, drop = FALSE]
  })
  if (length(blocks) == 0 || any(lengths(tissue_features) == 0)) {
    stop("ablation: every tissue requires at least one Direct-GSClassifier feature.", call. = FALSE)
  }
  blocks
}


# Pair Direct-GSClassifier and Cohort-d1 through tissue and global reduction. Both arms
# use the same tissue order, feasible dimensions, neighbors, and derived seeds.
.ablation_paired_two_stage_embeddings <- function(
    direct_blocks,
    cohort_d1,
    dr_config,
    seed
) {
  cohort_reference <- vapply(
    strsplit(colnames(cohort_d1), "|", fixed = TRUE),
    `[`,
    character(1),
    1
  )
  tissues <- names(direct_blocks)
  if (!setequal(tissues, unique(cohort_reference))) {
    stop("ablation: Direct-GSClassifier and Cohort-d1 tissue blocks do not match.", call. = FALSE)
  }
  tissues <- tissues[tissues %in% unique(cohort_reference)]
  cohort_blocks <- lapply(tissues, function(tissue) {
    cohort_d1[, cohort_reference == tissue, drop = FALSE]
  })
  names(cohort_blocks) <- tissues
  if (!all(vapply(direct_blocks, rownames, character(nrow(cohort_d1))) ==
      rownames(cohort_d1))) {
    stop("ablation: Direct-GSClassifier and Cohort-d1 sample order must match.", call. = FALSE)
  }

  set.seed(seed)
  stage_seeds <- sample.int(10000000L, length(tissues) + 1L, replace = FALSE)
  dr_args <- dr_config[setdiff(names(dr_config), c("method", "dimension"))]
  direct_d2 <- list()
  cohort_d2 <- list()
  audit <- list()

  for (i in seq_along(tissues)) {
    tissue <- tissues[i]
    direct_prepared <- .ablation_prepare_dr_block(direct_blocks[[tissue]])
    cohort_prepared <- .ablation_prepare_dr_block(cohort_blocks[[tissue]])
    dimensions <- .ablation_common_dr_budget(
      direct_prepared,
      cohort_prepared,
      dr_config$dimension[1],
      dr_args$n_neighbors
    )
    args_i <- dr_args
    if (!is.null(args_i$n_neighbors)) {
      args_i$n_neighbors <- dimensions$n_neighbors
    }
    direct_d2[[tissue]] <- .ablation_run_prepared_dr(
      direct_prepared,
      method = dr_config$method,
      dims = dimensions$dimension,
      seed = stage_seeds[i],
      dr_args = args_i,
      column_prefix = tissue
    )
    cohort_d2[[tissue]] <- .ablation_run_prepared_dr(
      cohort_prepared,
      method = dr_config$method,
      dims = dimensions$dimension,
      seed = stage_seeds[i],
      dr_args = args_i,
      column_prefix = tissue
    )
    audit[[i]] <- data.frame(
      stage = "tissue",
      tissue = tissue,
      stage_seed = stage_seeds[i],
      direct_available_dimension = dimensions$direct_available_dimension,
      cohort_available_dimension = dimensions$cohort_available_dimension,
      common_dimension = dimensions$dimension,
      direct_dimension = dimensions$dimension,
      cohort_dimension = dimensions$dimension,
      n_neighbors = dimensions$n_neighbors,
      stringsAsFactors = FALSE
    )
  }

  direct_d2 <- do.call(cbind, direct_d2)
  cohort_d2 <- do.call(cbind, cohort_d2)
  direct_global <- .ablation_prepare_dr_block(direct_d2)
  cohort_global <- .ablation_prepare_dr_block(cohort_d2)
  global_budget <- .ablation_common_dr_budget(
    direct_global,
    cohort_global,
    dr_config$dimension[2],
    dr_args$n_neighbors
  )
  global_args <- dr_args
  if (!is.null(global_args$n_neighbors)) {
    global_args$n_neighbors <- global_budget$n_neighbors
  }
  global_seed <- stage_seeds[length(stage_seeds)]
  direct_d3 <- .ablation_run_prepared_dr(
    direct_global,
    method = dr_config$method,
    dims = global_budget$dimension,
    seed = global_seed,
    dr_args = global_args,
    column_prefix = "all"
  )
  cohort_d3 <- .ablation_run_prepared_dr(
    cohort_global,
    method = dr_config$method,
    dims = global_budget$dimension,
    seed = global_seed,
    dr_args = global_args,
    column_prefix = "all"
  )
  audit[[length(audit) + 1L]] <- data.frame(
    stage = "global",
    tissue = NA_character_,
    stage_seed = global_seed,
    direct_available_dimension = global_budget$direct_available_dimension,
    cohort_available_dimension = global_budget$cohort_available_dimension,
    common_dimension = global_budget$dimension,
    direct_dimension = global_budget$dimension,
    cohort_dimension = global_budget$dimension,
    n_neighbors = global_budget$n_neighbors,
    stringsAsFactors = FALSE
  )

  direct_high <- do.call(cbind, lapply(tissues, function(tissue) {
    block <- direct_blocks[[tissue]]
    colnames(block) <- paste(tissue, colnames(block), sep = "|")
    block
  }))
  list(
    embeddings = list(`Direct-GSClassifier` = direct_d3, `Cohort-d1` = cohort_d3),
    high = list(`Direct-GSClassifier` = direct_high, `Cohort-d1` = cohort_d1),
    parameter_audit = do.call(rbind, audit)
  )
}


# Apply CCS preprocessing before determining a shared, feasible DR budget.
.ablation_prepare_dr_block <- function(data) {
  prepare <- getFromNamespace("data_for_dr", "CCS")
  list(
    original_rows = rownames(data),
    prepared = prepare(data, rm.dup.col = FALSE, verbose = FALSE)
  )
}


.ablation_common_dr_budget <- function(direct, cohort, dims, n_neighbors = NULL) {
  direct_clean <- direct$prepared$cleaned$data
  cohort_clean <- cohort$prepared$cleaned$data
  direct_available <- min(ncol(direct_clean), nrow(direct_clean) - 2L)
  cohort_available <- min(ncol(cohort_clean), nrow(cohort_clean) - 2L)
  dimension <- min(as.integer(dims), direct_available, cohort_available)
  if (dimension < 1) {
    stop("ablation: paired reduction has insufficient common rank.", call. = FALSE)
  }
  neighbors <- if (is.null(n_neighbors)) {
    NA_integer_
  } else {
    min(as.integer(n_neighbors), nrow(direct_clean) - 1L, nrow(cohort_clean) - 1L)
  }
  list(
    dimension = dimension,
    n_neighbors = neighbors,
    direct_available_dimension = direct_available,
    cohort_available_dimension = cohort_available
  )
}


# Run CORE_DR on an already prepared block, then restore duplicate rows exactly
# as the existing CCS reduction path does.
.ablation_run_prepared_dr <- function(
    prepared,
    method,
    dims,
    seed,
    dr_args,
    column_prefix
) {
  core_dr <- getFromNamespace("CORE_DR", "CCS")
  repair <- getFromNamespace("repairCCS", "CCS")
  embedding <- do.call(core_dr, c(list(
    method = method,
    data = prepared$prepared$cleaned$data,
    dims = dims,
    seed = seed
  ), dr_args))
  embedding <- as.matrix(embedding)
  rownames(embedding) <- names(prepared$prepared$cleaned$md5)
  embedding <- repair(embedding, prepared$prepared)
  embedding <- as.matrix(embedding[prepared$original_rows, , drop = FALSE])
  colnames(embedding) <- paste0(column_prefix, "|D", seq_len(ncol(embedding)))
  embedding
}


# Experiment 4 compares raw end-to-end DBSCAN solutions under shared samples,
# reductions, seeds, and parameter sets. It does not normalize clusters to four labels.
.ablation_experiment_metaccs <- function(prepared, config, manifest, seed, verbose) {
  parameter_manifest <- .ablation_metaccs_parameter_manifest(config)
  direct_blocks <- .ablation_direct_tissue_blocks(
    prepared$direct,
    prepared$feature_manifest$tissue_features
  )
  input_hash <- digest::digest(
    list(direct = prepared$direct, d1 = prepared$d1),
    algo = "md5"
  )
  metrics <- list()
  stratified <- list()
  assignments <- list()
  runs <- list()
  cross_arm <- list()
  parameter_audit <- list()
  metric_index <- 1L
  assignment_index <- 1L
  run_index <- 1L
  parameter_audit_index <- 1L

  for (resample_index in seq_along(config$metaccs$resample_seeds)) {
    resample_seed <- config$metaccs$resample_seeds[resample_index]
    resample_id <- sprintf("R%03d", resample_index)
    rows <- .ablation_tissue_subsample(
      prepared$metadata,
      config$metaccs$subsample_fraction,
      resample_seed
    )
    metadata <- prepared$metadata[rows, , drop = FALSE]
    resample_sample_hash <- digest::digest(sort(metadata$sample_id), algo = "md5")
    cohort_d1 <- prepared$d1[rows, , drop = FALSE]
    direct_i <- lapply(direct_blocks, function(x) x[rows, , drop = FALSE])

    for (umap_seed in config$metaccs$umap_seeds) {
      for (dr_param_set_id in unique(parameter_manifest$dr_param_set_id)) {
        dr_row <- which(parameter_manifest$dr_param_set_id == dr_param_set_id)[1]
        paired <- .ablation_paired_two_stage_embeddings(
          direct_blocks = direct_i,
          cohort_d1 = cohort_d1,
          dr_config = parameter_manifest$dr[[dr_row]],
          seed = umap_seed
        )
        audit_i <- paired$parameter_audit
        audit_i$resample_id <- resample_id
        audit_i$resample_seed <- resample_seed
        audit_i$resample_sample_hash <- resample_sample_hash
        audit_i$umap_seed <- umap_seed
        audit_i$dr_param_set_id <- dr_param_set_id
        audit_i$dr_cache_key <- paste(
          resample_id,
          umap_seed,
          dr_param_set_id,
          sep = "|"
        )
        parameter_audit[[parameter_audit_index]] <- audit_i
        parameter_audit_index <- parameter_audit_index + 1L

        cluster_rows <- which(parameter_manifest$dr_param_set_id == dr_param_set_id)
        for (parameter_row in cluster_rows) {
          cluster_param_set_id <- parameter_manifest$cluster_param_set_id[parameter_row]
          clusters <- lapply(paired$embeddings, function(d3) {
            .ablation_dbscan(d3, parameter_manifest$cluster[[parameter_row]])
          })
          run_signature <- paste(
            resample_id,
            umap_seed,
            dr_param_set_id,
            cluster_param_set_id,
            sep = "|"
          )

          for (group_id in names(paired$embeddings)) {
            d3 <- paired$embeddings[[group_id]]
            high <- paired$high[[group_id]]
            cluster <- clusters[[group_id]]
            values <- c(
              .ablation_embedding_metrics(
                high = high,
                low = d3,
                metadata = metadata,
                clusters = cluster,
                config = config,
                seed = umap_seed
              ),
              .ablation_cluster_biology(cluster, metadata$biology)
            )
            metrics[[metric_index]] <- data.frame(
              experiment_id = "metaccs",
              group_id = group_id,
              fold_id = resample_id,
              held_out_cohort = NA_character_,
              seed_type = "umap_seed",
              stochastic_seed = umap_seed,
              rank_q = ncol(d3),
              k = min(config$general$k, nrow(d3) - 1L),
              distance = config$general$distance,
              metric_name = names(values),
              metric_value = as.numeric(values),
              resample_id = resample_id,
              resample_seed = resample_seed,
              resample_sample_hash = resample_sample_hash,
              umap_seed = umap_seed,
              dr_param_set_id = dr_param_set_id,
              cluster_param_set_id = cluster_param_set_id,
              dr_cache_key = paste(resample_id, umap_seed, dr_param_set_id, sep = "|"),
              stringsAsFactors = FALSE
            )
            stratified_i <- .ablation_tissue_stratified_metrics(
              high = high,
              low = d3,
              metadata = metadata,
              group_id = group_id,
              seed = umap_seed,
              k = config$general$k
            )
            stratified_i$resample_id <- resample_id
            stratified_i$resample_seed <- resample_seed
            stratified_i$resample_sample_hash <- resample_sample_hash
            stratified_i$dr_param_set_id <- dr_param_set_id
            stratified_i$cluster_param_set_id <- cluster_param_set_id
            stratified[[metric_index]] <- stratified_i
            runs[[run_index]] <- list(
              run_id = paste(group_id, run_signature, sep = "|"),
              run_signature = run_signature,
              group_id = group_id,
              resample_id = resample_id,
              resample_seed = resample_seed,
              resample_sample_hash = resample_sample_hash,
              umap_seed = umap_seed,
              dr_param_set_id = dr_param_set_id,
              cluster_param_set_id = cluster_param_set_id,
              sample_id = metadata$sample_id,
              d3 = d3,
              cluster = cluster
            )
            run_index <- run_index + 1L
            if (isTRUE(config$metaccs$retain_assignments)) {
              assignments[[assignment_index]] <- data.frame(
                sample_id = metadata$sample_id,
                group_id = group_id,
                resample_id = resample_id,
                resample_seed = resample_seed,
                resample_sample_hash = resample_sample_hash,
                umap_seed = umap_seed,
                dr_param_set_id = dr_param_set_id,
                cluster_param_set_id = cluster_param_set_id,
                raw_cluster = unname(cluster),
                noise = unname(cluster) == 0,
                stringsAsFactors = FALSE
              )
              assignment_index <- assignment_index + 1L
            }
            metric_index <- metric_index + 1L
          }
          cross_arm[[length(cross_arm) + 1L]] <- data.frame(
            resample_id = resample_id,
            resample_seed = resample_seed,
            resample_sample_hash = resample_sample_hash,
            umap_seed = umap_seed,
            dr_param_set_id = dr_param_set_id,
            cluster_param_set_id = cluster_param_set_id,
            neighbor_jaccard = .ablation_knn_jaccard(
              paired$embeddings$`Direct-GSClassifier`,
              paired$embeddings$`Cohort-d1`,
              min(config$general$k, nrow(metadata) - 1L)
            ),
            ari = .ablation_ari(clusters$`Direct-GSClassifier`, clusters$`Cohort-d1`),
            cluster_jaccard = .ablation_cluster_jaccard(
              clusters$`Direct-GSClassifier`,
              clusters$`Cohort-d1`
            ),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (verbose) {
      luckyBase::LuckyVerbose("ablation: Experiment 4 completed ", resample_id, ".")
    }
  }

  metrics <- do.call(rbind, metrics)
  stratified <- .ablation_rbind(stratified)
  if (nrow(stratified) > 0) {
    tissue_sizes <- stats::aggregate(
      sample_count ~ tissue,
      data = stratified,
      FUN = max
    )
    small_tissue_cutoff <- stats::quantile(
      tissue_sizes$sample_count,
      0.25,
      names = FALSE
    )
    stratified$small_tissue <- stratified$sample_count <= small_tissue_cutoff
  }
  assignments <- .ablation_rbind(assignments)
  cross_arm <- .ablation_rbind(cross_arm)
  parameter_audit <- .ablation_rbind(parameter_audit)
  stability <- .ablation_metaccs_stability(runs, config$general$k)
  stability_contrasts <- .ablation_metaccs_stability_contrasts(stability)
  summary <- .ablation_metaccs_summary(
    metrics,
    bootstrap = config$general$bootstrap,
    seed = seed,
    subsample_fraction = config$metaccs$subsample_fraction
  )
  contrasts <- .ablation_metaccs_contrasts(
    metrics,
    bootstrap = config$general$bootstrap,
    seed = seed,
    subsample_fraction = config$metaccs$subsample_fraction
  )
  audit <- cbind(
    metrics,
    sample_manifest_hash = manifest$sample_manifest_hash,
    feature_manifest_hash = manifest$feature_manifest_hash,
    direct_feature_hash = prepared$feature_manifest$direct_feature_hash,
    input_hash = input_hash,
    module_sequence_id = NA_character_,
    module_count = manifest$module_count,
    projection_seed = NA_real_,
    permutation_seed = NA_real_,
    stringsAsFactors = FALSE
  )
  list(
    status = "complete",
    metrics = metrics,
    summary = summary,
    contrasts = contrasts,
    stratified = stratified,
    stability = stability,
    stability_contrasts = stability_contrasts,
    cross_arm_agreement = cross_arm,
    assignments = assignments,
    parameter_manifest = parameter_manifest,
    parameter_audit = parameter_audit,
    audit = audit
  )
}


# Quantify label-invariant agreement between raw clusters and independent biology.
.ablation_cluster_biology <- function(clusters, biology) {
  biology <- as.character(biology)
  non_noise <- !is.na(clusters) & clusters != 0
  valid <- non_noise & !is.na(biology) & nzchar(biology)
  coverage <- mean(non_noise)
  if (sum(valid) < 2 || length(unique(clusters[valid])) < 1 ||
      length(unique(biology[valid])) < 1) {
    return(c(
      cluster_biology_ari = NA_real_,
      cluster_biology_nmi = NA_real_,
      weighted_cluster_purity = NA_real_,
      non_noise_coverage = coverage
    ))
  }
  contingency <- table(clusters[valid], biology[valid])
  purity <- sum(apply(contingency, 1, max)) / sum(contingency)
  c(
    cluster_biology_ari = .ablation_ari(clusters[valid], biology[valid]),
    cluster_biology_nmi = .ablation_nmi(contingency),
    weighted_cluster_purity = purity,
    non_noise_coverage = coverage
  )
}


.ablation_nmi <- function(contingency) {
  probability <- contingency / sum(contingency)
  row_probability <- rowSums(probability)
  column_probability <- colSums(probability)
  occupied <- which(probability > 0, arr.ind = TRUE)
  information <- sum(vapply(seq_len(nrow(occupied)), function(i) {
    row_i <- occupied[i, 1]
    column_i <- occupied[i, 2]
    value <- probability[row_i, column_i]
    value * log(value / (row_probability[row_i] * column_probability[column_i]))
  }, numeric(1)))
  entropy_row <- -sum(row_probability[row_probability > 0] *
    log(row_probability[row_probability > 0]))
  entropy_column <- -sum(column_probability[column_probability > 0] *
    log(column_probability[column_probability > 0]))
  denominator <- sqrt(entropy_row * entropy_column)
  if (denominator == 0) {
    return(if (identical(dim(contingency), c(1L, 1L))) 1 else NA_real_)
  }
  information / denominator
}


# Compare algorithm/resample runs only within the same arm and parameter pair.
.ablation_metaccs_stability <- function(runs, k) {
  if (length(runs) < 2) {
    return(data.frame())
  }
  group_key <- vapply(runs, function(x) {
    paste(x$group_id, x$dr_param_set_id, x$cluster_param_set_id, sep = "|")
  }, character(1))
  groups <- split(seq_along(runs), group_key)
  result <- list()
  index <- 1L
  for (indices in groups) {
    if (length(indices) < 2) next
    pairs <- utils::combn(indices, 2, simplify = FALSE)
    for (pair in pairs) {
      first <- runs[[pair[1]]]
      second <- runs[[pair[2]]]
      common <- intersect(first$sample_id, second$sample_id)
      if (length(common) < 2) next
      first_rows <- match(common, first$sample_id)
      second_rows <- match(common, second$sample_id)
      local_k <- min(k, length(common) - 1L)
      result[[index]] <- data.frame(
        group_id = first$group_id,
        dr_param_set_id = first$dr_param_set_id,
        cluster_param_set_id = first$cluster_param_set_id,
        run_a = first$run_signature,
        run_b = second$run_signature,
        resample_id_a = first$resample_id,
        resample_id_b = second$resample_id,
        resample_seed_a = first$resample_seed,
        resample_seed_b = second$resample_seed,
        umap_seed_a = first$umap_seed,
        umap_seed_b = second$umap_seed,
        sample_hash_a = first$resample_sample_hash,
        sample_hash_b = second$resample_sample_hash,
        neighbor_jaccard = .ablation_knn_jaccard(
          first$d3[first_rows, , drop = FALSE],
          second$d3[second_rows, , drop = FALSE],
          local_k
        ),
        ari = .ablation_ari(
          first$cluster[first_rows],
          second$cluster[second_rows]
        ),
        cluster_jaccard = .ablation_cluster_jaccard(
          first$cluster[first_rows],
          second$cluster[second_rows]
        ),
        common_sample_count = length(common),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  .ablation_rbind(result)
}


.ablation_metaccs_stability_contrasts <- function(stability) {
  if (nrow(stability) == 0) {
    return(data.frame())
  }
  direct <- stability[stability$group_id == "Direct-GSClassifier", , drop = FALSE]
  cohort <- stability[stability$group_id == "Cohort-d1", , drop = FALSE]
  by <- c("dr_param_set_id", "cluster_param_set_id", "run_a", "run_b")
  paired <- merge(cohort, direct, by = by, suffixes = c(".cohort", ".direct"))
  metrics <- c("neighbor_jaccard", "ari", "cluster_jaccard")
  .ablation_rbind(lapply(metrics, function(metric) {
    data.frame(
      paired[, by, drop = FALSE],
      contrast = "Cohort-d1-Direct-GSClassifier",
      metric_name = metric,
      metric_value = paired[[paste0(metric, ".cohort")]] -
        paired[[paste0(metric, ".direct")]],
      stringsAsFactors = FALSE
    )
  }))
}


# Average UMAP seeds within resample before bootstrapping the actual resample unit.
.ablation_metaccs_summary <- function(
    metrics,
    bootstrap,
    seed,
    subsample_fraction
) {
  data <- metrics[is.finite(metrics$metric_value), , drop = FALSE]
  data <- stats::aggregate(
    metric_value ~ group_id + resample_id + dr_param_set_id +
      cluster_param_set_id + metric_name,
    data = data,
    FUN = mean
  )
  keys <- interaction(
    data$group_id,
    data$dr_param_set_id,
    data$cluster_param_set_id,
    data$metric_name,
    drop = TRUE
  )
  parts <- split(data, keys)
  .ablation_rbind(lapply(seq_along(parts), function(i) {
    part <- parts[[i]]
    interval <- .ablation_bootstrap_mean(part$metric_value, bootstrap, seed + i)
    data.frame(
      group_id = part$group_id[1],
      dr_param_set_id = part$dr_param_set_id[1],
      cluster_param_set_id = part$cluster_param_set_id[1],
      metric_name = part$metric_name[1],
      estimate = mean(part$metric_value),
      ci_low = interval[1],
      ci_high = interval[2],
      n_resamples = nrow(part),
      inference_scope = if (subsample_fraction < 1 && nrow(part) > 1) {
        "resample_variation"
      } else {
        "algorithm_variation_only"
      },
      stringsAsFactors = FALSE
    )
  }))
}


.ablation_metaccs_contrasts <- function(
    metrics,
    bootstrap,
    seed,
    subsample_fraction
) {
  data <- metrics[is.finite(metrics$metric_value), , drop = FALSE]
  data <- stats::aggregate(
    metric_value ~ group_id + resample_id + dr_param_set_id +
      cluster_param_set_id + metric_name,
    data = data,
    FUN = mean
  )
  direct <- data[data$group_id == "Direct-GSClassifier", , drop = FALSE]
  cohort <- data[data$group_id == "Cohort-d1", , drop = FALSE]
  by <- c("resample_id", "dr_param_set_id", "cluster_param_set_id", "metric_name")
  paired <- merge(cohort, direct, by = by, suffixes = c(".cohort", ".direct"))
  paired$delta <- paired$metric_value.cohort - paired$metric_value.direct
  keys <- interaction(
    paired$dr_param_set_id,
    paired$cluster_param_set_id,
    paired$metric_name,
    drop = TRUE
  )
  parts <- split(paired, keys)
  .ablation_rbind(lapply(seq_along(parts), function(i) {
    part <- parts[[i]]
    interval <- .ablation_bootstrap_mean(part$delta, bootstrap, seed + 10000L + i)
    data.frame(
      contrast = "Cohort-d1-Direct-GSClassifier",
      dr_param_set_id = part$dr_param_set_id[1],
      cluster_param_set_id = part$cluster_param_set_id[1],
      metric_name = part$metric_name[1],
      estimate = mean(part$delta),
      ci_low = interval[1],
      ci_high = interval[2],
      n_pairs = nrow(part),
      inference_scope = if (subsample_fraction < 1 && nrow(part) > 1) {
        "resample_variation"
      } else {
        "algorithm_variation_only"
      },
      stringsAsFactors = FALSE
    )
  }))
}


.ablation_bootstrap_mean <- function(values, bootstrap, seed) {
  if (length(values) == 0) {
    return(c(NA_real_, NA_real_))
  }
  set.seed(seed)
  draws <- replicate(
    as.integer(bootstrap),
    mean(sample(values, replace = TRUE))
  )
  stats::quantile(draws, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}


# Predict one frozen cohort module from a precomputed Direct-GSClassifier matrix.
# Reusing the shared 529-feature matrix avoids repeating geneMatch/trainDataProc_X
# for every model while remaining numerically identical to callEnsemble().
.ablation_predict_module_from_direct <- function(direct, model, module_id) {
  direct <- as.matrix(direct)
  repeats <- model$Model
  if (length(repeats) == 0) {
    stop("ablation: frozen module has no ensemble repeats.", call. = FALSE)
  }
  classes <- sort(unique(unlist(lapply(repeats, names), use.names = FALSE)))
  if (length(classes) == 0) {
    stop("ablation: frozen module has no class models.", call. = FALSE)
  }

  probability <- vapply(classes, function(class_id) {
    repeat_probability <- vapply(repeats, function(repeat_model) {
      class_model <- repeat_model[[class_id]]
      if (is.null(class_model) || length(class_model) <= 1) {
        return(rep(0, nrow(direct)))
      }
      features <- class_model$bst$feature_names
      if (length(features) == 0) {
        features <- class_model$genes
      }
      missing <- setdiff(features, colnames(direct))
      if (length(missing) > 0) {
        stop(
          "ablation: Direct matrix is missing frozen features for ",
          module_id,
          ": ",
          paste(missing, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      as.numeric(stats::predict(
        class_model$bst,
        direct[, features, drop = FALSE]
      ))
    }, numeric(nrow(direct)))
    if (is.null(dim(repeat_probability))) {
      repeat_probability
    } else {
      apply(repeat_probability, 1, stats::median)
    }
  }, numeric(nrow(direct)))

  probability <- matrix(
    probability,
    nrow = nrow(direct),
    dimnames = list(
      rownames(direct),
      paste(module_id, classes, sep = "|")
    )
  )
  if (any(!is.finite(probability))) {
    stop("ablation: frozen module prediction contains non-finite values.", call. = FALSE)
  }
  probability
}


# Encode one or more frozen cohort modules from the shared Direct matrix.
# PSOCK is used for cross-platform parallelism; workers only read frozen models
# and return isolated module blocks, while the main process owns column ordering.
.ablation_encode_d1_from_direct <- function(
    object,
    direct,
    module_manifest = .ablation_module_manifest(object),
    module_ids = module_manifest$modules$module_id,
    numCores = 1L,
    verbose = TRUE
) {
  direct <- as.matrix(direct)
  available <- module_manifest$modules$module_id
  module_ids <- as.character(module_ids)
  unknown <- setdiff(module_ids, available)
  if (length(unknown) > 0) {
    stop(
      "ablation: unknown frozen module(s): ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  module_ids <- available[available %in% module_ids]
  numCores <- min(as.integer(numCores), length(module_ids))
  if (!is.finite(numCores) || numCores < 1) {
    stop("ablation: numCores must be a positive integer.", call. = FALSE)
  }

  path_map <- .ablation_model_path_map(object)
  missing_paths <- module_ids[!module_ids %in% names(path_map)]
  if (length(missing_paths) > 0) {
    stop(
      "ablation: modelFit.rds is missing for module(s): ",
      paste(missing_paths, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  worker <- function(module_id) {
    model <- readRDS(path_map[[module_id]])
    .ablation_predict_module_from_direct(direct, model, module_id)
  }

  if (verbose) {
    luckyBase::LuckyVerbose(
      "ablation: encoding ",
      length(module_ids),
      " frozen modules with ",
      numCores,
      " worker(s)..."
    )
  }
  if (numCores == 1L) {
    blocks <- lapply(module_ids, worker)
  } else {
    cluster <- parallel::makePSOCKcluster(numCores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(
      cluster,
      suppressPackageStartupMessages(library(xgboost))
    )
    parallel::clusterExport(
      cluster,
      c(
        "direct",
        "path_map",
        ".ablation_predict_module_from_direct"
      ),
      envir = environment()
    )
    blocks <- parallel::parLapply(cluster, module_ids, function(module_id) {
      model <- readRDS(path_map[[module_id]])
      .ablation_predict_module_from_direct(direct, model, module_id)
    })
  }
  names(blocks) <- module_ids

  encoded <- do.call(cbind, blocks)
  expected_columns <- colnames(object@Data$Probability$d1)[
    unlist(module_manifest$blocks[module_ids], use.names = FALSE)
  ]
  missing_columns <- setdiff(expected_columns, colnames(encoded))
  if (length(missing_columns) > 0) {
    stop(
      "ablation: encoded d1 is missing expected probability columns.",
      call. = FALSE
    )
  }
  encoded <- encoded[, expected_columns, drop = FALSE]
  rownames(encoded) <- rownames(direct)
  encoded
}


# Standardize d1 on the reference boundary and give every cohort module equal
# total squared-distance weight, independent of whether its block has 3 or 4 columns.
.ablation_module_balanced_transform <- function(reference, query, blocks) {
  reference <- as.matrix(reference)
  query <- as.matrix(query)
  if (!identical(colnames(reference), colnames(query))) {
    stop("ablation: reference and query d1 columns must be identical.", call. = FALSE)
  }
  block_columns <- unlist(blocks, use.names = FALSE)
  if (!setequal(block_columns, seq_len(ncol(reference)))) {
    stop("ablation: d1 blocks must cover every column exactly once.", call. = FALSE)
  }
  if (anyDuplicated(block_columns)) {
    stop("ablation: d1 blocks must not overlap.", call. = FALSE)
  }

  center <- colMeans(reference)
  scale <- apply(reference, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  reference_scaled <- sweep(sweep(reference, 2, center, "-"), 2, scale, "/")
  query_scaled <- sweep(sweep(query, 2, center, "-"), 2, scale, "/")

  weights <- numeric(ncol(reference))
  for (block in blocks) {
    weights[block] <- 1 / sqrt(length(block) * length(blocks))
  }
  list(
    reference = sweep(reference_scaled, 2, weights, "*"),
    query = sweep(query_scaled, 2, weights, "*"),
    center = center,
    scale = scale,
    weights = weights,
    distance = "module_balanced_standardized_euclidean"
  )
}


# Evaluate each query only against a separate reference atlas. This avoids the
# held-out-fold bug where cohort mixing is necessarily zero inside one test cohort.
.ablation_query_reference_retrieval <- function(
    reference,
    query,
    reference_metadata,
    query_metadata,
    label_column,
    technical_columns = character(),
    k = c(5L, 15L, 30L),
    search = c("exact", "annoy"),
    seed = 20260727,
    n_trees = 50L,
    search_k = -1L
) {
  reference <- as.matrix(reference)
  query <- as.matrix(query)
  search <- match.arg(search)
  k <- sort(unique(as.integer(k)))
  if (any(k < 1) || max(k) >= nrow(reference)) {
    stop("ablation: retrieval k must be between 1 and n_reference - 1.", call. = FALSE)
  }
  if (!identical(colnames(reference), colnames(query))) {
    stop("ablation: retrieval matrices must have identical columns.", call. = FALSE)
  }
  if (nrow(reference_metadata) != nrow(reference) ||
      nrow(query_metadata) != nrow(query)) {
    stop("ablation: retrieval metadata does not align with matrices.", call. = FALSE)
  }
  required <- c("sample_id", "cohort", label_column)
  if (!all(required %in% colnames(reference_metadata)) ||
      !all(required %in% colnames(query_metadata))) {
    stop("ablation: retrieval metadata is missing required columns.", call. = FALSE)
  }
  technical_columns <- intersect(
    technical_columns,
    intersect(colnames(reference_metadata), colnames(query_metadata))
  )

  max_k <- max(k)
  candidate_k <- min(nrow(reference), max(max_k * 5L, max_k + 20L))
  if (search == "exact") {
    neighbor_id <- matrix(NA_integer_, nrow(query), max_k)
    neighbor_distance <- matrix(NA_real_, nrow(query), max_k)
    for (i in seq_len(nrow(query))) {
      eligible <- reference_metadata$cohort != query_metadata$cohort[i]
      if (sum(eligible) < max_k) {
        stop(
          "ablation: fewer than max(k) cross-cohort reference samples for query ",
          query_metadata$sample_id[i],
          ".",
          call. = FALSE
        )
      }
      distance <- sqrt(rowSums(
        (sweep(reference[eligible, , drop = FALSE], 2, query[i, ], "-"))^2
      ))
      eligible_rows <- which(eligible)
      selected <- order(distance)[seq_len(max_k)]
      neighbor_id[i, ] <- eligible_rows[selected]
      neighbor_distance[i, ] <- distance[selected]
    }
  } else {
    set.seed(seed)
    index <- new(RcppAnnoy::AnnoyEuclidean, ncol(reference))
    index$setSeed(as.integer(seed))
    for (i in seq_len(nrow(reference))) {
      index$addItem(i - 1L, reference[i, ])
    }
    index$build(as.integer(n_trees))
    neighbor_id <- matrix(NA_integer_, nrow(query), max_k)
    neighbor_distance <- matrix(NA_real_, nrow(query), max_k)
    for (i in seq_len(nrow(query))) {
      found <- index$getNNsByVectorList(
        query[i, ],
        as.integer(candidate_k),
        as.integer(search_k),
        TRUE
      )
      ids <- as.integer(found$item) + 1L
      distances <- as.numeric(found$distance)
      keep <- reference_metadata$cohort[ids] != query_metadata$cohort[i]
      ids <- ids[keep]
      distances <- distances[keep]
      if (length(ids) < max_k) {
        stop(
          "ablation: Annoy returned fewer than max(k) cross-cohort neighbors.",
          call. = FALSE
        )
      }
      neighbor_id[i, ] <- ids[seq_len(max_k)]
      neighbor_distance[i, ] <- distances[seq_len(max_k)]
    }
  }

  neighbor_rows <- lapply(seq_len(nrow(query)), function(i) {
    data.frame(
      query_sample = query_metadata$sample_id[i],
      query_cohort = query_metadata$cohort[i],
      query_label = as.character(query_metadata[[label_column]][i]),
      neighbor_rank = seq_len(max_k),
      reference_sample = reference_metadata$sample_id[neighbor_id[i, ]],
      reference_cohort = reference_metadata$cohort[neighbor_id[i, ]],
      reference_label = as.character(reference_metadata[[label_column]][neighbor_id[i, ]]),
      distance = neighbor_distance[i, ],
      stringsAsFactors = FALSE
    )
  })
  neighbors <- do.call(rbind, neighbor_rows)
  neighbors$label_match <- neighbors$query_label == neighbors$reference_label
  for (column in technical_columns) {
    query_value <- rep(query_metadata[[column]], each = max_k)
    reference_value <- reference_metadata[[column]][as.vector(t(neighbor_id))]
    neighbors[[paste0(column, "_match")]] <- query_value == reference_value
  }

  first_match <- vapply(seq_len(nrow(query)), function(i) {
    hit <- which(neighbor_id[i, ] > 0 &
      as.character(reference_metadata[[label_column]][neighbor_id[i, ]]) ==
        as.character(query_metadata[[label_column]][i]))
    if (length(hit) == 0) NA_integer_ else hit[1]
  }, integer(1))
  per_sample <- do.call(rbind, lapply(k, function(k_i) {
    rows <- neighbors$neighbor_rank <= k_i
    selected <- neighbors[rows, , drop = FALSE]
    label_rate <- stats::aggregate(
      label_match ~ query_sample,
      data = selected,
      FUN = mean
    )
    result <- data.frame(
      sample_id = query_metadata$sample_id,
      cohort = query_metadata$cohort,
      label = as.character(query_metadata[[label_column]]),
      k = k_i,
      top1_label_match = as.numeric(
        neighbors$label_match[neighbors$neighbor_rank == 1]
      ),
      top_k_label_rate = label_rate$label_match[
        match(query_metadata$sample_id, label_rate$query_sample)
      ],
      mrr = ifelse(is.na(first_match), 0, 1 / first_match),
      stringsAsFactors = FALSE
    )
    for (column in technical_columns) {
      match_column <- paste0(column, "_match")
      observed <- stats::aggregate(
        selected[[match_column]],
        by = list(query_sample = selected$query_sample),
        FUN = function(x) {
          x <- x[!is.na(x)]
          if (length(x) == 0) NA_real_ else mean(x)
        }
      )
      expected <- vapply(seq_len(nrow(query_metadata)), function(i) {
        pool <- as.character(reference_metadata[[label_column]]) ==
          as.character(query_metadata[[label_column]][i]) &
          reference_metadata$cohort != query_metadata$cohort[i]
        reference_value <- as.character(reference_metadata[[column]][pool])
        query_value <- as.character(query_metadata[[column]][i])
        valid <- !is.na(reference_value) & nzchar(reference_value) &
          !is.na(query_value) & nzchar(query_value)
        if (!any(valid)) {
          NA_real_
        } else {
          mean(reference_value[valid] == query_value)
        }
      }, numeric(1))
      observed_rate <- observed$x[
        match(query_metadata$sample_id, observed$query_sample)
      ]
      result[[paste0(column, "_match_rate")]] <- observed_rate
      result[[paste0(column, "_expected_rate")]] <- expected
      result[[paste0(column, "_match_excess")]] <- observed_rate - expected
    }
    result
  }))

  summary <- stats::aggregate(
    cbind(top1_label_match, top_k_label_rate, mrr) ~ k,
    data = per_sample,
    FUN = mean
  )
  list(
    neighbors = neighbors,
    per_sample = per_sample,
    summary = summary,
    search = search,
    k = k
  )
}


# Quantify approximate-neighbor fidelity on the exact same query/reference task.
.ablation_validate_neighbor_search <- function(
    reference,
    query,
    reference_metadata,
    query_metadata,
    label_column,
    k,
    query_samples = 30L,
    n_trees = 50L,
    search_k = -1L,
    seed = 20260727
) {
  query_samples <- min(as.integer(query_samples), nrow(query))
  if (query_samples < 1) {
    stop("ablation: neighbor validation requires query samples.", call. = FALSE)
  }
  rows <- if (query_samples == nrow(query)) {
    seq_len(nrow(query))
  } else {
    .ablation_stratified_sample(query_metadata, query_samples, seed)
  }
  exact <- .ablation_query_reference_retrieval(
    reference,
    query[rows, , drop = FALSE],
    reference_metadata,
    query_metadata[rows, , drop = FALSE],
    label_column = label_column,
    k = k,
    search = "exact",
    seed = seed
  )
  approximate <- .ablation_query_reference_retrieval(
    reference,
    query[rows, , drop = FALSE],
    reference_metadata,
    query_metadata[rows, , drop = FALSE],
    label_column = label_column,
    k = k,
    search = "annoy",
    seed = seed,
    n_trees = n_trees,
    search_k = search_k
  )
  exact_sets <- split(
    exact$neighbors$reference_sample,
    exact$neighbors$query_sample
  )
  approximate_sets <- split(
    approximate$neighbors$reference_sample,
    approximate$neighbors$query_sample
  )
  recall <- vapply(names(exact_sets), function(sample_id) {
    length(intersect(exact_sets[[sample_id]], approximate_sets[[sample_id]])) /
      length(exact_sets[[sample_id]])
  }, numeric(1))
  list(
    recall = mean(recall),
    per_sample_recall = recall,
    query_sample_count = length(rows),
    k = max(as.integer(k)),
    n_trees = n_trees,
    search_k = search_k,
    seed = seed
  )
}


# Null-Perm has no label-level power when every query cohort has a constant anchor.
.ablation_null_perm_eligibility <- function(metadata, label_column) {
  if (!all(c("cohort", label_column) %in% colnames(metadata))) {
    stop("ablation: Null-Perm metadata is missing cohort or anchor.", call. = FALSE)
  }
  label_count <- vapply(
    split(as.character(metadata[[label_column]]), metadata$cohort),
    function(x) length(unique(x[!is.na(x) & nzchar(x)])),
    integer(1)
  )
  eligible <- names(label_count)[label_count > 1]
  list(
    status = if (length(eligible) > 0) "eligible" else "not_eligible",
    eligible_cohorts = eligible,
    cohort_label_count = label_count,
    reason = if (length(eligible) > 0) {
      NA_character_
    } else {
      "anchor_is_constant_within_every_query_cohort"
    }
  )
}


# Confirmatory conclusions require both a genuinely independent anchor and d1
# generated without exposing the query sample to its own frozen cohort model.
.ablation_evidence_level <- function(
    query_metadata,
    anchor_role,
    provenance_column = "d1_provenance"
) {
  if (!provenance_column %in% colnames(query_metadata)) {
    stop("ablation: query metadata is missing d1 provenance.", call. = FALSE)
  }
  provenance <- unique(as.character(query_metadata[[provenance_column]]))
  qualified <- all(provenance %in% c("external_frozen", "out_of_fold"))
  independent <- identical(anchor_role, "independent")
  reasons <- c(
    if (!independent) "anchor_is_not_independent",
    if (!qualified) "query_d1_provenance_is_not_external_or_out_of_fold"
  )
  list(
    level = if (qualified && independent) "confirmatory" else "descriptive",
    qualified_provenance = qualified,
    anchor_role = anchor_role,
    provenance = provenance,
    reasons = reasons
  )
}


# Defaults are grouped by the scientific question so configuration changes are
# auditable and downstream legacy parameters cannot leak into representation tests.
.ablation_representation_default_params <- function(seed = 20260727) {
  list(
    comparison = list(
      module_ids = NULL,
      direct_group = "Direct-GSClassifier",
      cohort_group = "Cohort-d1"
    ),
    provenance = list(
      external_cohorts = NULL,
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
      n_trees = 50L,
      search_k = 5000L,
      exact_validation_queries = 30L,
      min_annoy_recall = 0.8,
      geometry_samples = 5000L,
      distance_pairs = 100000L
    ),
    validation = list(
      enabled = TRUE,
      learning_fractions = c(0.1, 0.25, 0.5, 1),
      repeats = 3L,
      inner_folds = 3L,
      lambda = c(0.1, 1, 10),
      nrounds = 50L,
      min_class_n = 20L,
      numCores = 1L
    ),
    scaling = list(
      enabled = FALSE,
      module_counts = c(10L, 25L, 50L, 75L, 100L, 125L, 150L),
      sequences = 10L,
      direct_feature_type = "all",
      sensitivity_feature_type = "gene_pair",
      biology_anchors = character(),
      score_reference_samples = 5000L,
      score_query_samples = 2000L,
      lambda = 1,
      bootstrap = 1000L
    ),
    controls = list(
      null_rp = TRUE,
      null_rp_rank = 100L,
      null_rp_seeds = seed + seq_len(3L),
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
      cover = FALSE,
      cache_external_d1 = TRUE
    )
  )
}


.ablation_resolve_representation_config <- function(seed, params) {
  default <- .ablation_representation_default_params(seed)
  .ablation_validate_override(default, params, path = "params")
  config <- .ablation_merge_lists(default, params)
  if (length(config$geometry$k) == 0 || any(config$geometry$k < 1)) {
    stop("ablation: geometry$k must contain positive integers.", call. = FALSE)
  }
  if (!config$geometry$search %in% c("exact", "annoy")) {
    stop("ablation: geometry$search must be exact or annoy.", call. = FALSE)
  }
  if (config$anchors$min_reference_cohorts < 1) {
    stop("ablation: anchors$min_reference_cohorts must be positive.", call. = FALSE)
  }
  if (length(config$scaling$enabled) != 1 ||
      !is.logical(config$scaling$enabled) ||
      is.na(config$scaling$enabled)) {
    stop("ablation: scaling$enabled must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(config$scaling$module_counts) < 2 ||
      any(!is.finite(config$scaling$module_counts)) ||
      any(config$scaling$module_counts < 1) ||
      any(config$scaling$module_counts != as.integer(config$scaling$module_counts))) {
    stop(
      "ablation: scaling$module_counts must contain at least two positive integers.",
      call. = FALSE
    )
  }
  if (length(config$scaling$sequences) != 1 ||
      !is.finite(config$scaling$sequences) ||
      config$scaling$sequences < 2 ||
      config$scaling$sequences != as.integer(config$scaling$sequences)) {
    stop("ablation: scaling$sequences must be an integer of at least two.", call. = FALSE)
  }
  if (length(config$scaling$direct_feature_type) != 1 ||
      !config$scaling$direct_feature_type %in% c("gene_pair", "all")) {
    stop(
      "ablation: scaling$direct_feature_type must be gene_pair or all.",
      call. = FALSE
    )
  }
  if (length(config$scaling$sensitivity_feature_type) != 1 ||
      !config$scaling$sensitivity_feature_type %in% c("gene_pair", "none")) {
    stop(
      paste0(
        "ablation: scaling$sensitivity_feature_type must be gene_pair or none."
      ),
      call. = FALSE
    )
  }
  if (!is.character(config$scaling$biology_anchors) ||
      anyNA(config$scaling$biology_anchors)) {
    stop("ablation: scaling$biology_anchors must be a character vector.", call. = FALSE)
  }
  for (field in c("score_reference_samples", "score_query_samples")) {
    value <- config$scaling[[field]]
    if (length(value) != 1L || (!is.infinite(value) &&
        (!is.finite(value) || value < 2L || value != as.integer(value)))) {
      stop(
        "ablation: scaling$", field, " must be an integer of at least two or Inf.",
        call. = FALSE
      )
    }
  }
  if (length(config$scaling$lambda) != 1 ||
      !is.finite(config$scaling$lambda) ||
      config$scaling$lambda < 0) {
    stop("ablation: scaling$lambda must be one non-negative value.", call. = FALSE)
  }
  if (length(config$scaling$bootstrap) != 1 ||
      !is.finite(config$scaling$bootstrap) ||
      config$scaling$bootstrap < 1) {
    stop("ablation: scaling$bootstrap must be positive.", call. = FALSE)
  }
  if (config$scaling$enabled && !config$validation$enabled) {
    stop(
      "ablation: representation scaling requires validation$enabled = TRUE.",
      call. = FALSE
    )
  }
  for (field in c("max_reference_samples", "max_query_samples")) {
    value <- config$provenance[[field]]
    if (length(value) != 1 || (!is.infinite(value) && value < 1)) {
      stop("ablation: provenance sample caps must be positive or Inf.", call. = FALSE)
    }
  }
  config
}


.ablation_limit_metadata <- function(metadata, size, seed) {
  if (!is.finite(size) || nrow(metadata) <= size) {
    return(metadata)
  }
  metadata[.ablation_stratified_sample(metadata, as.integer(size), seed), , drop = FALSE]
}


# Prepare disjoint reference/query matrices. Existing d1 is used only for the
# reference atlas; every filtered query is re-encoded through the frozen bank.
.ablation_prepare_representation_input <- function(
    object,
    data,
    metadata,
    config,
    output.dir,
    seed,
    verbose
) {
  module_manifest <- .ablation_module_manifest(object)
  feature_manifest <- .ablation_frozen_feature_manifest(object, module_manifest)
  module_ids <- config$comparison$module_ids
  if (is.null(module_ids)) {
    module_ids <- module_manifest$modules$module_id
  }
  module_ids <- as.character(module_ids)
  unknown_modules <- setdiff(module_ids, module_manifest$modules$module_id)
  if (length(unknown_modules) > 0) {
    stop(
      "ablation: comparison$module_ids contains unknown modules: ",
      paste(unknown_modules, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  module_ids <- module_manifest$modules$module_id[
    module_manifest$modules$module_id %in% module_ids
  ]

  flattened <- .ablation_flatten_expression(data)
  metadata <- .ablation_prepare_metadata(metadata, flattened$metadata, object)
  metadata$cohort_key <- paste(metadata$tissue, metadata$cohort, sep = "/")
  filtered <- config$provenance$external_cohorts
  if (is.null(filtered)) {
    filtered <- as.character(object@Data$filtered.cohort)
  }
  filtered <- intersect(filtered, unique(metadata$cohort_key))
  if (config$provenance$require_external && length(filtered) == 0) {
    stop("ablation: no filtered external cohorts are available.", call. = FALSE)
  }

  d1 <- as.matrix(object@Data$Probability$d1)
  expression_ids <- colnames(flattened$expr)
  reference_ids <- Reduce(intersect, list(
    rownames(d1),
    expression_ids,
    metadata$sample_id[!metadata$cohort_key %in% filtered]
  ))
  query_ids <- intersect(
    expression_ids,
    metadata$sample_id[metadata$cohort_key %in% filtered]
  )
  if (length(reference_ids) < 3 || length(query_ids) < 1) {
    stop("ablation: reference or external query samples are unavailable.", call. = FALSE)
  }

  reference_metadata <- metadata[match(reference_ids, metadata$sample_id), , drop = FALSE]
  query_metadata <- metadata[match(query_ids, metadata$sample_id), , drop = FALSE]
  reference_metadata <- .ablation_limit_metadata(
    reference_metadata,
    config$provenance$max_reference_samples,
    seed
  )
  query_metadata <- .ablation_limit_metadata(
    query_metadata,
    config$provenance$max_query_samples,
    seed + 1L
  )
  reference_ids <- reference_metadata$sample_id
  query_ids <- query_metadata$sample_id
  reference_metadata$d1_provenance <- "in_sample"
  query_metadata$d1_provenance <- "external_frozen"

  if (verbose) {
    luckyBase::LuckyVerbose(
      "ablation: building Direct-GSClassifier for ",
      length(reference_ids),
      " reference and ",
      length(query_ids),
      " external query samples..."
    )
  }
  all_ids <- c(reference_ids, query_ids)
  direct <- .ablation_gsclassifier_matrix(
    object,
    flattened$expr[, all_ids, drop = FALSE],
    feature_manifest
  )
  reference_direct <- direct[reference_ids, , drop = FALSE]
  query_direct <- direct[query_ids, , drop = FALSE]

  expected_columns <- colnames(d1)[
    unlist(module_manifest$blocks[module_ids], use.names = FALSE)
  ]
  reference_d1 <- d1[reference_ids, expected_columns, drop = FALSE]
  model_paths <- .ablation_model_path_map(object)[module_ids]
  model_info <- file.info(unname(model_paths))[, c("size", "mtime"), drop = FALSE]
  cache_key <- digest::digest(
    list(
      query_ids = query_ids,
      direct = query_direct,
      module_ids = module_ids,
      model_paths = model_paths,
      model_info = model_info
    ),
    algo = "md5"
  )
  cache_path <- file.path(output.dir, "external-d1-cache.rds")
  query_d1 <- NULL
  if (config$output$cache_external_d1 && file.exists(cache_path)) {
    cache <- readRDS(cache_path)
    if (identical(cache$key, cache_key)) {
      query_d1 <- cache$d1
      if (verbose) {
        luckyBase::LuckyVerbose("ablation: using verified external d1 cache.")
      }
    }
  }
  if (is.null(query_d1)) {
    query_d1 <- .ablation_encode_d1_from_direct(
      object = object,
      direct = query_direct,
      module_manifest = module_manifest,
      module_ids = module_ids,
      numCores = config$validation$numCores,
      verbose = verbose
    )
    if (config$output$cache_external_d1) {
      saveRDS(list(key = cache_key, d1 = query_d1), cache_path)
    }
  }

  selected_blocks <- lapply(module_ids, function(module_id) {
    module_columns <- colnames(d1)[module_manifest$blocks[[module_id]]]
    match(module_columns, expected_columns)
  })
  names(selected_blocks) <- module_ids
  list(
    reference_direct = reference_direct,
    query_direct = query_direct,
    reference_d1 = reference_d1,
    query_d1 = query_d1,
    reference_metadata = reference_metadata,
    query_metadata = query_metadata,
    module_manifest = module_manifest,
    selected_blocks = selected_blocks,
    selected_module_ids = module_ids,
    feature_manifest = feature_manifest,
    excluded_duplicate_samples = flattened$excluded_duplicate_samples,
    cache_key = cache_key,
    filtered_cohorts = filtered
  )
}


.ablation_native_geometry <- function(prepared, transformed, config, seed) {
  metadata <- prepared$reference_metadata
  rows <- seq_len(nrow(metadata))
  if (length(rows) > config$geometry$geometry_samples) {
    rows <- .ablation_stratified_sample(
      metadata,
      config$geometry$geometry_samples,
      seed
    )
  }
  direct <- transformed$direct$reference[rows, , drop = FALSE]
  d1 <- transformed$d1$reference[rows, , drop = FALSE]
  local_k <- min(max(config$geometry$k), length(rows) - 1L)
  metrics <- data.frame(
    metric_name = c(
      "linear_cka",
      "distance_spearman",
      "knn_jaccard",
      "direct_effective_rank",
      "d1_effective_rank"
    ),
    metric_value = c(
      .ablation_linear_cka(direct, d1),
      .ablation_distance_spearman(
        direct,
        d1,
        config$geometry$distance_pairs,
        seed
      ),
      .ablation_knn_jaccard(direct, d1, local_k),
      .ablation_effective_rank(direct),
      .ablation_effective_rank(d1)
    ),
    sample_count = length(rows),
    stringsAsFactors = FALSE
  )

  raw_d1 <- prepared$reference_d1[rows, , drop = FALSE]
  module_diagnostics <- do.call(rbind, lapply(
    names(prepared$selected_blocks),
    function(module_id) {
      block <- prepared$selected_blocks[[module_id]]
      values <- raw_d1[, block, drop = FALSE]
      data.frame(
        module_id = module_id,
        block_width = length(block),
        mean_probability = mean(values),
        mean_column_variance = mean(apply(values, 2, stats::var)),
        mean_row_sum = mean(rowSums(values)),
        row_sum_close_one = mean(abs(rowSums(values) - 1) < 1e-6),
        stringsAsFactors = FALSE
      )
    }
  ))
  module_diagnostics$variance_share <-
    module_diagnostics$mean_column_variance /
    sum(module_diagnostics$mean_column_variance)
  list(
    metrics = metrics,
    module_diagnostics = module_diagnostics,
    d1_probability_contract = list(
      value_range = range(raw_d1),
      row_sum_range = range(unlist(lapply(
        prepared$selected_blocks,
        function(block) rowSums(raw_d1[, block, drop = FALSE])
      ), use.names = FALSE)),
      all_blocks_simplex = all(module_diagnostics$row_sum_close_one == 1),
      distance = transformed$d1$distance
    )
  )
}


.ablation_bind_retrieval <- function(results) {
  neighbors <- do.call(rbind, lapply(names(results), function(name) {
    data <- results[[name]]$neighbors
    data$representation <- name
    data
  }))
  per_sample <- do.call(rbind, lapply(names(results), function(name) {
    data <- results[[name]]$per_sample
    data$representation <- name
    data
  }))
  summary <- do.call(rbind, lapply(names(results), function(name) {
    data <- results[[name]]$summary
    data$representation <- name
    data
  }))
  by_cohort <- stats::aggregate(
    cbind(top1_label_match, top_k_label_rate, mrr) ~ representation + cohort + k,
    data = per_sample,
    FUN = mean
  )

  direct <- per_sample[
    per_sample$representation == "Direct-GSClassifier",
    ,
    drop = FALSE
  ]
  cohort <- per_sample[
    per_sample$representation == "Cohort-d1",
    ,
    drop = FALSE
  ]
  paired <- merge(
    cohort,
    direct,
    by = c("sample_id", "cohort", "label", "k"),
    suffixes = c("_d1", "_direct")
  )
  for (metric in c("top1_label_match", "top_k_label_rate", "mrr")) {
    paired[[paste0("delta_", metric)]] <-
      paired[[paste0(metric, "_d1")]] - paired[[paste0(metric, "_direct")]]
  }
  list(
    neighbors = neighbors,
    per_sample = per_sample,
    summary = summary,
    by_cohort = by_cohort,
    paired = paired
  )
}


# Apply the shared external-query eligibility boundary once so full and
# scaling-only representation runs cannot drift in sample composition.
.ablation_prepare_representation_analysis <- function(
    object,
    data,
    metadata,
    config,
    output.dir,
    seed,
    verbose
) {
  prepared <- .ablation_prepare_representation_input(
    object,
    data,
    metadata,
    config,
    output.dir,
    seed,
    verbose
  )
  anchor <- config$anchors$primary
  if (!anchor %in% colnames(prepared$reference_metadata) ||
      !anchor %in% colnames(prepared$query_metadata)) {
    stop("ablation: primary anchor is missing from metadata.", call. = FALSE)
  }
  reference_cohort_count <- vapply(
    split(
      prepared$reference_metadata$cohort,
      prepared$reference_metadata[[anchor]]
    ),
    function(x) length(unique(x)),
    integer(1)
  )
  eligible_labels <- names(reference_cohort_count)[
    reference_cohort_count >= config$anchors$min_reference_cohorts
  ]
  keep_query <- as.character(prepared$query_metadata[[anchor]]) %in% eligible_labels
  prepared$query_metadata <- prepared$query_metadata[keep_query, , drop = FALSE]
  prepared$query_direct <- prepared$query_direct[keep_query, , drop = FALSE]
  prepared$query_d1 <- prepared$query_d1[keep_query, , drop = FALSE]
  if (nrow(prepared$query_metadata) == 0) {
    stop("ablation: no external query has an eligible reference anchor.", call. = FALSE)
  }
  list(prepared = prepared, anchor = anchor)
}


.ablation_run_representation <- function(
    object,
    data,
    metadata,
    output.dir,
    params,
    seed,
    verbose
) {
  if (!methods::is(object, "CCS")) {
    stop("ablation: object must be a CCS object.", call. = FALSE)
  }
  config <- .ablation_resolve_representation_config(seed, params)
  if (dir.exists(output.dir) && length(list.files(output.dir)) > 0 &&
      !config$output$cover) {
    stop(
      "ablation: output.dir is not empty; set params$output$cover = TRUE.",
      call. = FALSE
    )
  }
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  analysis <- .ablation_prepare_representation_analysis(
    object = object,
    data = data,
    metadata = metadata,
    config = config,
    output.dir = output.dir,
    seed = seed,
    verbose = verbose
  )
  prepared <- analysis$prepared
  anchor <- analysis$anchor

  direct_scaled <- .ablation_scale_train_apply(
    prepared$reference_direct,
    prepared$query_direct
  )
  d1_scaled <- .ablation_module_balanced_transform(
    prepared$reference_d1,
    prepared$query_d1,
    prepared$selected_blocks
  )
  transformed <- list(
    direct = list(
      reference = direct_scaled$train,
      query = direct_scaled$test,
      center = direct_scaled$center,
      scale = direct_scaled$scale,
      distance = "standardized_euclidean"
    ),
    d1 = d1_scaled
  )
  native_geometry <- .ablation_native_geometry(prepared, transformed, config, seed)

  retrieval_results <- list(
    `Direct-GSClassifier` = .ablation_query_reference_retrieval(
      transformed$direct$reference,
      transformed$direct$query,
      prepared$reference_metadata,
      prepared$query_metadata,
      label_column = anchor,
      technical_columns = config$anchors$technical,
      k = config$geometry$k,
      search = config$geometry$search,
      seed = seed,
      n_trees = config$geometry$n_trees,
      search_k = config$geometry$search_k
    ),
    `Cohort-d1` = .ablation_query_reference_retrieval(
      transformed$d1$reference,
      transformed$d1$query,
      prepared$reference_metadata,
      prepared$query_metadata,
      label_column = anchor,
      technical_columns = config$anchors$technical,
      k = config$geometry$k,
      search = config$geometry$search,
      seed = seed,
      n_trees = config$geometry$n_trees,
      search_k = config$geometry$search_k
    )
  )
  search_validation <- if (config$geometry$search == "annoy") {
    validation <- list(
      `Direct-GSClassifier` = .ablation_validate_neighbor_search(
        transformed$direct$reference,
        transformed$direct$query,
        prepared$reference_metadata,
        prepared$query_metadata,
        label_column = anchor,
        k = max(config$geometry$k),
        query_samples = config$geometry$exact_validation_queries,
        n_trees = config$geometry$n_trees,
        search_k = config$geometry$search_k,
        seed = seed + 1000L
      ),
      `Cohort-d1` = .ablation_validate_neighbor_search(
        transformed$d1$reference,
        transformed$d1$query,
        prepared$reference_metadata,
        prepared$query_metadata,
        label_column = anchor,
        k = max(config$geometry$k),
        query_samples = config$geometry$exact_validation_queries,
        n_trees = config$geometry$n_trees,
        search_k = config$geometry$search_k,
        seed = seed + 1000L
      )
    )
    recall <- vapply(validation, `[[`, numeric(1), "recall")
    if (any(recall < config$geometry$min_annoy_recall)) {
      stop(
        "ablation: Annoy recall is below geometry$min_annoy_recall: ",
        paste(names(recall), round(recall, 3), sep = "=", collapse = ", "),
        ". Increase n_trees/search_k or use exact search.",
        call. = FALSE
      )
    }
    validation
  } else {
    list(status = "not_required", search = "exact")
  }
  retrieval <- .ablation_bind_retrieval(retrieval_results)
  retrieval$search_validation <- search_validation
  evidence <- .ablation_evidence_level(
    prepared$query_metadata,
    config$anchors$primary_role
  )
  null_perm <- if (config$controls$null_perm) {
    .ablation_null_perm_eligibility(prepared$query_metadata, anchor)
  } else {
    list(status = "not_run", reason = "disabled")
  }
  null_rp <- if (config$controls$null_rp) {
    rp_rank <- min(
      as.integer(config$controls$null_rp_rank),
      ncol(transformed$direct$reference)
    )
    rp_results <- lapply(config$controls$null_rp_seeds, function(rp_seed) {
      projection <- .ablation_projection_matrix(
        ncol(transformed$direct$reference),
        rp_rank,
        rp_seed
      )
      .ablation_query_reference_retrieval(
        transformed$direct$reference %*% projection,
        transformed$direct$query %*% projection,
        prepared$reference_metadata,
        prepared$query_metadata,
        label_column = anchor,
        technical_columns = config$anchors$technical,
        k = config$geometry$k,
        search = config$geometry$search,
        seed = rp_seed,
        n_trees = config$geometry$n_trees,
        search_k = config$geometry$search_k
      )
    })
    names(rp_results) <- as.character(config$controls$null_rp_seeds)
    rp_summary <- do.call(rbind, lapply(seq_along(rp_results), function(i) {
      data <- rp_results[[i]]$summary
      data$seed <- config$controls$null_rp_seeds[i]
      data
    }))
    list(
      status = "complete",
      rank = rp_rank,
      seeds = config$controls$null_rp_seeds,
      results = rp_results,
      summary = rp_summary
    )
  } else {
    list(status = "not_run", reason = "disabled")
  }
  controls <- list(
    null_rp = null_rp,
    null_perm = null_perm
  )
  if (config$validation$enabled) {
    readout_results <- list(
      `Direct-GSClassifier` = .ablation_linear_readout(
        train = prepared$reference_direct,
        test = prepared$query_direct,
        train_metadata = prepared$reference_metadata,
        test_metadata = prepared$query_metadata,
        label_column = anchor,
        lambda = config$validation$lambda,
        inner_folds = config$validation$inner_folds,
        nrounds = config$validation$nrounds,
        numCores = config$validation$numCores,
        seed = seed + 20000L,
        blocks = NULL
      ),
      `Cohort-d1` = .ablation_linear_readout(
        train = prepared$reference_d1,
        test = prepared$query_d1,
        train_metadata = prepared$reference_metadata,
        test_metadata = prepared$query_metadata,
        label_column = anchor,
        lambda = config$validation$lambda,
        inner_folds = config$validation$inner_folds,
        nrounds = config$validation$nrounds,
        numCores = config$validation$numCores,
        seed = seed + 20000L,
        blocks = prepared$selected_blocks
      )
    )
    overall <- do.call(rbind, lapply(names(readout_results), function(name) {
      data <- readout_results[[name]]$overall
      data$representation <- name
      data$selected_lambda <- readout_results[[name]]$selected_lambda
      data
    }))
    by_cohort <- do.call(rbind, lapply(names(readout_results), function(name) {
      data <- readout_results[[name]]$by_cohort
      data$representation <- name
      data
    }))
    predictions <- do.call(rbind, lapply(names(readout_results), function(name) {
      data <- readout_results[[name]]$predictions
      data$representation <- name
      data
    }))
    direct_cohort <- by_cohort[
      by_cohort$representation == "Direct-GSClassifier",
      ,
      drop = FALSE
    ]
    d1_cohort <- by_cohort[
      by_cohort$representation == "Cohort-d1",
      ,
      drop = FALSE
    ]
    paired_by_cohort <- merge(
      d1_cohort,
      direct_cohort,
      by = "cohort",
      suffixes = c("_d1", "_direct")
    )
    for (metric in c("accuracy", "balanced_accuracy", "macro_auroc")) {
      paired_by_cohort[[paste0("delta_", metric)]] <-
        paired_by_cohort[[paste0(metric, "_d1")]] -
        paired_by_cohort[[paste0(metric, "_direct")]]
    }
    readout <- list(
      status = "complete",
      results = readout_results,
      overall = overall,
      by_cohort = by_cohort,
      paired_by_cohort = paired_by_cohort,
      predictions = predictions
    )
    learning_curve <- .ablation_learning_curve(
      representations = list(
        `Direct-GSClassifier` = list(
          train = prepared$reference_direct,
          test = prepared$query_direct,
          blocks = NULL
        ),
        `Cohort-d1` = list(
          train = prepared$reference_d1,
          test = prepared$query_d1,
          blocks = prepared$selected_blocks
        )
      ),
      train_metadata = prepared$reference_metadata,
      test_metadata = prepared$query_metadata,
      label_column = anchor,
      fractions = config$validation$learning_fractions,
      repeats = config$validation$repeats,
      lambda = config$validation$lambda,
      inner_folds = config$validation$inner_folds,
      nrounds = config$validation$nrounds,
      numCores = config$validation$numCores,
      seed = seed + 30000L
    )
  } else {
    readout <- list(status = "not_run", reason = "disabled")
    learning_curve <- list(status = "not_run", reason = "disabled")
  }

  cohort_scaling <- if (config$scaling$enabled) {
    .ablation_representation_scaling(
      prepared = prepared,
      config = config,
      label_column = anchor,
      seed = seed + 35000L,
      verbose = verbose,
      cache_path = file.path(output.dir, "cohort-scaling-fit-cache.rds")
    )
  } else {
    list(status = "not_run", reason = "disabled")
  }

  feature_counts <- table(factor(
    prepared$feature_manifest$feature_manifest$feature_type,
    levels = c("single_bin", "gene_pair", "set_pair")
  ))
  decoder <- if (config$tradeoffs$decoder) {
    decoder_reference_metadata <- .ablation_limit_metadata(
      prepared$reference_metadata,
      config$tradeoffs$decoder_max_reference_samples,
      seed + 40000L
    )
    decoder_query_metadata <- .ablation_limit_metadata(
      prepared$query_metadata,
      config$tradeoffs$decoder_max_query_samples,
      seed + 40001L
    )
    reference_rows <- match(
      decoder_reference_metadata$sample_id,
      prepared$reference_metadata$sample_id
    )
    query_rows <- match(
      decoder_query_metadata$sample_id,
      prepared$query_metadata$sample_id
    )
    decoded <- .ablation_decode_direct_features(
      reference_d1 = prepared$reference_d1[reference_rows, , drop = FALSE],
      query_d1 = prepared$query_d1[query_rows, , drop = FALSE],
      reference_direct = prepared$reference_direct[reference_rows, , drop = FALSE],
      query_direct = prepared$query_direct[query_rows, , drop = FALSE],
      feature_manifest = prepared$feature_manifest$feature_manifest,
      blocks = prepared$selected_blocks,
      rank = config$tradeoffs$decoder_rank,
      lambda = config$tradeoffs$decoder_lambda
    )
    decoded$reference_sample_count <- length(reference_rows)
    decoded$query_sample_count <- length(query_rows)
    decoded
  } else {
    list(status = "not_run", reason = "disabled")
  }
  tradeoffs <- list(
    status = "complete",
    feature_type_count = feature_counts,
    d1_probability_contract = native_geometry$d1_probability_contract,
    module_diagnostics = native_geometry$module_diagnostics,
    excluded_duplicate_samples = prepared$excluded_duplicate_samples,
    decoder = decoder
  )
  manifest <- list(
    version = 5L,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = seed,
    experiment = "representation",
    groups = c("Direct-GSClassifier", "Cohort-d1"),
    reference_sample_count = nrow(prepared$reference_metadata),
    query_sample_count = nrow(prepared$query_metadata),
    reference_cohort_count = length(unique(prepared$reference_metadata$cohort)),
    query_cohort_count = length(unique(prepared$query_metadata$cohort)),
    direct_feature_count = ncol(prepared$reference_direct),
    tsp_feature_count = length(prepared$feature_manifest$tsp_features),
    d1_feature_count = ncol(prepared$reference_d1),
    module_count = length(prepared$selected_blocks),
    gene_signature_count = unname(as.integer(feature_counts["single_bin"])),
    scaling_direct_feature_count = if (config$scaling$enabled) {
      cohort_scaling$diagnostics$direct_contracts$feature_count[
        cohort_scaling$diagnostics$direct_contracts$contract_role == "main"
      ][1]
    } else {
      NA_integer_
    },
    scaling_schema_version = if (config$scaling$enabled) {
      cohort_scaling$schema_version
    } else {
      NA_integer_
    },
    external_cohorts = prepared$filtered_cohorts,
    anchor = anchor,
    anchor_role = config$anchors$primary_role,
    evidence_level = evidence$level,
    evidence_reasons = evidence$reasons,
    direct_distance = transformed$direct$distance,
    d1_distance = transformed$d1$distance,
    cache_key = prepared$cache_key,
    config = config,
    config_hash = digest::digest(config, algo = "md5")
  )

  audit <- retrieval$summary
  audit$evidence_level <- evidence$level
  audit$anchor <- anchor
  audit$reference_sample_count <- nrow(prepared$reference_metadata)
  audit$query_sample_count <- nrow(prepared$query_metadata)
  audit$config_hash <- manifest$config_hash

  saveRDS(manifest, file.path(output.dir, "manifest.rds"))
  saveRDS(native_geometry, file.path(output.dir, "native_geometry.rds"))
  saveRDS(retrieval, file.path(output.dir, "retrieval.rds"))
  saveRDS(readout, file.path(output.dir, "readout.rds"))
  saveRDS(learning_curve, file.path(output.dir, "learning_curve.rds"))
  saveRDS(cohort_scaling, file.path(output.dir, "cohort_scaling.rds"))
  saveRDS(tradeoffs, file.path(output.dir, "tradeoffs.rds"))
  utils::write.csv(audit, file.path(output.dir, "audit.csv"), row.names = FALSE)

  structure(
    list(
      call = match.call(),
      experiment = "representation",
      evidence_level = evidence$level,
      evidence = evidence,
      manifest = manifest,
      native_geometry = native_geometry,
      retrieval = retrieval,
      readout = readout,
      learning_curve = learning_curve,
      cohort_scaling = cohort_scaling,
      tradeoffs = tradeoffs,
      controls = controls,
      output.dir = normalizePath(output.dir, winslash = "/", mustWork = TRUE)
    ),
    class = "CCSAblation"
  )
}


.ablation_readout_transform <- function(train, test, blocks = NULL) {
  if (is.null(blocks)) {
    scaled <- .ablation_scale_train_apply(train, test)
    return(list(train = scaled$train, test = scaled$test))
  }
  balanced <- .ablation_module_balanced_transform(train, test, blocks)
  list(train = balanced$reference, test = balanced$query)
}


.ablation_xgb_linear_predict <- function(
    train,
    test,
    train_label,
    classes,
    lambda,
    nrounds,
    numCores,
    seed
) {
  encoded <- match(as.character(train_label), classes) - 1L
  class_n <- table(encoded)
  weight <- as.numeric(1 / class_n[as.character(encoded)])
  weight <- weight / mean(weight)
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(train),
    label = encoded,
    weight = weight
  )
  set.seed(seed)
  if (length(classes) == 2L) {
    fit <- xgboost::xgboost(
      data = dtrain,
      booster = "gblinear",
      updater = "coord_descent",
      feature_selector = "cyclic",
      objective = "binary:logistic",
      nrounds = as.integer(nrounds),
      eta = 0.1,
      lambda = lambda,
      alpha = 0,
      nthread = as.integer(numCores),
      verbose = 0
    )
    positive <- as.numeric(stats::predict(fit, as.matrix(test)))
    probability <- cbind(1 - positive, positive)
  } else {
    fit <- xgboost::xgboost(
      data = dtrain,
      booster = "gblinear",
      updater = "coord_descent",
      feature_selector = "cyclic",
      objective = "multi:softprob",
      num_class = length(classes),
      nrounds = as.integer(nrounds),
      eta = 0.1,
      lambda = lambda,
      alpha = 0,
      nthread = as.integer(numCores),
      verbose = 0
    )
    probability <- matrix(
      stats::predict(fit, as.matrix(test)),
      ncol = length(classes),
      byrow = TRUE
    )
  }
  colnames(probability) <- classes
  rownames(probability) <- rownames(test)
  list(
    fit = fit,
    probability = probability,
    prediction = classes[max.col(probability, ties.method = "first")]
  )
}


.ablation_classification_metrics <- function(truth, prediction, probability, classes) {
  truth <- as.character(truth)
  prediction <- as.character(prediction)
  recalls <- vapply(classes, function(class_id) {
    rows <- truth == class_id
    if (!any(rows)) NA_real_ else mean(prediction[rows] == class_id)
  }, numeric(1))
  auc <- vapply(seq_along(classes), function(i) {
    .ablation_binary_auc(as.integer(truth == classes[i]), probability[, i])
  }, numeric(1))
  data.frame(
    sample_count = length(truth),
    class_count = length(unique(truth)),
    accuracy = mean(prediction == truth),
    balanced_accuracy = mean(recalls, na.rm = TRUE),
    macro_auroc = mean(auc, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


# Tune one deterministic linear probe inside grouped reference folds, then fit
# once on all reference samples and evaluate only the untouched external query.
.ablation_linear_readout <- function(
    train,
    test,
    train_metadata,
    test_metadata,
    label_column,
    lambda = c(0.1, 1, 10),
    inner_folds = 3L,
    nrounds = 50L,
    numCores = 1L,
    seed = 20260727,
    blocks = NULL
) {
  train <- as.matrix(train)
  test <- as.matrix(test)
  if (nrow(train) != nrow(train_metadata) || nrow(test) != nrow(test_metadata)) {
    stop("ablation: readout metadata does not align with matrices.", call. = FALSE)
  }
  if (!label_column %in% colnames(train_metadata) ||
      !label_column %in% colnames(test_metadata)) {
    stop("ablation: readout label is missing from metadata.", call. = FALSE)
  }
  train_label <- as.character(train_metadata[[label_column]])
  test_label <- as.character(test_metadata[[label_column]])
  valid_train <- !is.na(train_label) & nzchar(train_label)
  classes <- sort(unique(train_label[valid_train]))
  valid_test <- !is.na(test_label) & test_label %in% classes
  train <- train[valid_train, , drop = FALSE]
  train_metadata <- train_metadata[valid_train, , drop = FALSE]
  train_label <- train_label[valid_train]
  test <- test[valid_test, , drop = FALSE]
  test_metadata <- test_metadata[valid_test, , drop = FALSE]
  test_label <- test_label[valid_test]
  if (length(classes) < 2 || nrow(test) < 2) {
    stop("ablation: readout requires at least two train/test classes.", call. = FALSE)
  }

  fold <- .ablation_grouped_folds(
    train_metadata$cohort,
    n_folds = inner_folds,
    seed = seed,
    label = train_label
  )
  lambda <- sort(unique(as.numeric(lambda)))
  if (length(lambda) == 1L) {
    selected_lambda <- lambda
    inner_cv <- data.frame(
      lambda = lambda,
      balanced_accuracy = NA_real_,
      estimable_folds = 0L,
      selection_mode = "fixed",
      stringsAsFactors = FALSE
    )
  } else {
    cv_rows <- lapply(seq_along(lambda), function(lambda_index) {
      fold_score <- vapply(sort(unique(fold)), function(fold_id) {
        inner_train <- fold != fold_id
        inner_test <- fold == fold_id
        fold_classes <- sort(unique(train_label[inner_train]))
        eligible_test <- inner_test & train_label %in% fold_classes
        if (length(fold_classes) < 2 || sum(eligible_test) < 2) {
          return(NA_real_)
        }
        transformed <- .ablation_readout_transform(
          train[inner_train, , drop = FALSE],
          train[eligible_test, , drop = FALSE],
          blocks
        )
        prediction <- .ablation_xgb_linear_predict(
          transformed$train,
          transformed$test,
          train_label[inner_train],
          fold_classes,
          lambda[lambda_index],
          nrounds,
          numCores,
          seed + lambda_index * 100L + fold_id
        )
        metrics <- .ablation_classification_metrics(
          train_label[eligible_test],
          prediction$prediction,
          prediction$probability,
          fold_classes
        )
        metrics$balanced_accuracy
      }, numeric(1))
      data.frame(
        lambda = lambda[lambda_index],
        balanced_accuracy = mean(fold_score, na.rm = TRUE),
        estimable_folds = sum(is.finite(fold_score)),
        selection_mode = "tuned",
        stringsAsFactors = FALSE
      )
    })
    inner_cv <- do.call(rbind, cv_rows)
    eligible_lambda <- inner_cv$estimable_folds > 0 &
      is.finite(inner_cv$balanced_accuracy)
    if (!any(eligible_lambda)) {
      stop("ablation: no lambda is estimable in grouped inner CV.", call. = FALSE)
    }
    ranking <- order(
      -inner_cv$balanced_accuracy[eligible_lambda],
      -inner_cv$lambda[eligible_lambda]
    )
    selected_lambda <- inner_cv$lambda[eligible_lambda][ranking[1]]
  }

  transformed <- .ablation_readout_transform(train, test, blocks)
  final <- .ablation_xgb_linear_predict(
    transformed$train,
    transformed$test,
    train_label,
    classes,
    selected_lambda,
    nrounds,
    numCores,
    seed + 10000L
  )
  predictions <- data.frame(
    sample_id = test_metadata$sample_id,
    cohort = test_metadata$cohort,
    true_label = test_label,
    predicted_label = final$prediction,
    max_probability = apply(final$probability, 1, max),
    stringsAsFactors = FALSE
  )
  overall <- .ablation_classification_metrics(
    test_label,
    final$prediction,
    final$probability,
    classes
  )
  by_cohort <- do.call(rbind, lapply(
    split(seq_len(nrow(predictions)), predictions$cohort),
    function(rows) {
      metric <- .ablation_classification_metrics(
        test_label[rows],
        final$prediction[rows],
        final$probability[rows, , drop = FALSE],
        classes
      )
      metric$cohort <- predictions$cohort[rows[1]]
      metric
    }
  ))
  rownames(by_cohort) <- NULL
  list(
    status = "complete",
    selected_lambda = selected_lambda,
    inner_cv = inner_cv,
    predictions = predictions,
    probability = final$probability,
    overall = overall,
    by_cohort = by_cohort,
    classes = classes,
    fold = fold,
    train_sample_hash = digest::digest(sort(train_metadata$sample_id), algo = "md5"),
    test_sample_hash = digest::digest(sort(test_metadata$sample_id), algo = "md5")
  )
}


.ablation_sample_training_cohorts <- function(metadata, label_column, fraction, seed) {
  cohort_label <- stats::aggregate(
    metadata[[label_column]],
    by = list(cohort = metadata$cohort),
    FUN = function(x) names(sort(table(x), decreasing = TRUE))[1]
  )
  colnames(cohort_label)[2] <- "label"
  set.seed(seed)
  selected <- unlist(lapply(split(cohort_label$cohort, cohort_label$label), function(x) {
    target <- min(length(x), max(2L, ceiling(length(x) * fraction)))
    sample(x, target)
  }), use.names = FALSE)
  sort(unique(selected))
}


# Build paired curves over shared cohort subsets. Each representation is tuned
# independently inside the same subset, preserving equal search budgets.
.ablation_learning_curve <- function(
    representations,
    train_metadata,
    test_metadata,
    label_column,
    fractions,
    repeats,
    lambda,
    inner_folds,
    nrounds,
    numCores,
    seed
) {
  rows <- list()
  index <- 1L
  test_hash <- digest::digest(sort(test_metadata$sample_id), algo = "md5")
  for (fraction_index in seq_along(fractions)) {
    fraction <- fractions[fraction_index]
    for (repeat_id in seq_len(as.integer(repeats))) {
      subset_seed <- seed + fraction_index * 1000L + repeat_id
      cohorts <- .ablation_sample_training_cohorts(
        train_metadata,
        label_column,
        fraction,
        subset_seed
      )
      train_rows <- train_metadata$cohort %in% cohorts
      subset_hash <- digest::digest(cohorts, algo = "md5")
      for (representation in names(representations)) {
        input <- representations[[representation]]
        fit <- .ablation_linear_readout(
          train = input$train[train_rows, , drop = FALSE],
          test = input$test,
          train_metadata = train_metadata[train_rows, , drop = FALSE],
          test_metadata = test_metadata,
          label_column = label_column,
          lambda = lambda,
          inner_folds = inner_folds,
          nrounds = nrounds,
          numCores = numCores,
          seed = subset_seed,
          blocks = input$blocks
        )
        rows[[index]] <- data.frame(
          representation = representation,
          requested_fraction = fraction,
          realized_cohort_fraction = length(cohorts) /
            length(unique(train_metadata$cohort)),
          repeat_id = repeat_id,
          train_cohort_count = length(cohorts),
          train_sample_count = sum(train_rows),
          selected_lambda = fit$selected_lambda,
          accuracy = fit$overall$accuracy,
          balanced_accuracy = fit$overall$balanced_accuracy,
          macro_auroc = fit$overall$macro_auroc,
          cohort_subset_hash = subset_hash,
          test_sample_hash = test_hash,
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
    }
  }
  metrics <- do.call(rbind, rows)
  direct <- metrics[metrics$representation == "Direct-GSClassifier", , drop = FALSE]
  d1 <- metrics[metrics$representation == "Cohort-d1", , drop = FALSE]
  paired <- merge(
    direct,
    d1,
    by = c("requested_fraction", "repeat_id"),
    suffixes = c("_direct", "_d1")
  )
  paired$delta_balanced_accuracy <-
    paired$balanced_accuracy_d1 - paired$balanced_accuracy_direct
  paired$delta_macro_auroc <- paired$macro_auroc_d1 - paired$macro_auroc_direct
  list(
    status = "complete",
    metrics = metrics,
    paired = paired,
    test_sample_hash = test_hash
  )
}


# Summarize nested-bank contrasts with the randomized module sequence as the
# uncertainty unit. Slopes are expressed per doubling of the cohort-module bank.
.ablation_representation_scaling_summary_v1 <- function(metrics, bootstrap, seed) {
  metric_names <- c(
    "balanced_accuracy_direct",
    "balanced_accuracy_d1",
    "delta_balanced_accuracy",
    "macro_auroc_direct",
    "macro_auroc_d1",
    "delta_macro_auroc"
  )
  point_parts <- lapply(seq_along(metric_names), function(metric_index) {
    metric_name <- metric_names[metric_index]
    do.call(rbind, lapply(
      split(metrics, metrics$module_count),
      function(data) {
        values <- data[[metric_name]]
        values <- values[is.finite(values)]
        interval <- .ablation_bootstrap_mean(
          values,
          bootstrap,
          seed + metric_index * 1000L + data$module_count[1]
        )
        data.frame(
          module_count = data$module_count[1],
          module_fraction = data$module_fraction[1],
          d1_feature_count = data$d1_feature_count[1],
          direct_feature_count = data$direct_feature_count[1],
          metric_name = metric_name,
          estimate = mean(values),
          ci_low = interval[1],
          ci_high = interval[2],
          n_sequences = length(values),
          stringsAsFactors = FALSE
        )
      }
    ))
  })
  pointwise <- do.call(rbind, point_parts)
  rownames(pointwise) <- NULL

  delta_names <- c("delta_balanced_accuracy", "delta_macro_auroc")
  sequence_trends <- do.call(rbind, lapply(
    split(metrics, metrics$sequence_id),
    function(data) {
      data <- data[order(data$module_count), , drop = FALSE]
      do.call(rbind, lapply(delta_names, function(metric_name) {
        values <- data[[metric_name]]
        keep <- is.finite(values)
        fit <- stats::lm(values[keep] ~ log2(data$module_count[keep]))
        data.frame(
          sequence_id = data$sequence_id[1],
          metric_name = metric_name,
          slope_per_doubling = unname(stats::coef(fit)[2]),
          endpoint_change = utils::tail(values[keep], 1) - values[keep][1],
          spearman_rho = if (length(unique(values[keep])) == 1L) {
            0
          } else {
            stats::cor(
              data$module_count[keep],
              values[keep],
              method = "spearman"
            )
          },
          stringsAsFactors = FALSE
        )
      }))
    }
  ))

  statistic_names <- c(
    "slope_per_doubling",
    "endpoint_change",
    "spearman_rho"
  )
  trend <- do.call(rbind, lapply(seq_along(delta_names), function(metric_index) {
    metric_name <- delta_names[metric_index]
    data <- sequence_trends[
      sequence_trends$metric_name == metric_name,
      ,
      drop = FALSE
    ]
    do.call(rbind, lapply(seq_along(statistic_names), function(statistic_index) {
      statistic <- statistic_names[statistic_index]
      values <- data[[statistic]]
      values <- values[is.finite(values)]
      interval <- .ablation_bootstrap_mean(
        values,
        bootstrap,
        seed + 10000L + metric_index * 1000L + statistic_index
      )
      data.frame(
        metric_name = metric_name,
        statistic = statistic,
        estimate = mean(values),
        ci_low = interval[1],
        ci_high = interval[2],
        n_sequences = length(values),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(trend) <- NULL

  list(
    pointwise = pointwise,
    sequence_trends = sequence_trends,
    trend = trend
  )
}


# Hold the external query set and Direct TSP contract fixed while expanding d1
# through tissue-balanced nested cohort-module banks. Direct and the fixed-lambda
# readout seed remain unchanged across the entire curve, so adjacent d1 changes
# are not contaminated by a moving baseline or repeated model selection.
.ablation_representation_scaling_v1 <- function(
    prepared,
    config,
    label_column,
    seed,
    verbose,
    cache_path = NULL
) {
  modules <- prepared$module_manifest$modules
  modules <- modules[
    match(prepared$selected_module_ids, modules$module_id),
    ,
    drop = FALSE
  ]
  counts <- sort(unique(pmin(
    as.integer(config$scaling$module_counts),
    nrow(modules)
  )))
  if (length(counts) < 2) {
    stop(
      "ablation: representation scaling needs at least two distinct module counts.",
      call. = FALSE
    )
  }
  sequences <- .ablation_nested_module_sequences(
    modules,
    config$scaling$sequences,
    seed
  )

  direct_features <- if (config$scaling$direct_feature_type == "gene_pair") {
    prepared$feature_manifest$tsp_features
  } else {
    prepared$feature_manifest$features
  }
  direct_columns <- match(direct_features, colnames(prepared$reference_direct))
  if (anyNA(direct_columns)) {
    stop(
      "ablation: the fixed scaling baseline is missing Direct features.",
      call. = FALSE
    )
  }
  reference_direct <- prepared$reference_direct[, direct_columns, drop = FALSE]
  query_direct <- prepared$query_direct[, direct_columns, drop = FALSE]
  test_hash <- digest::digest(
    sort(prepared$query_metadata$sample_id),
    algo = "md5"
  )
  readout_seed <- seed + 50000L
  fit_cache_key <- digest::digest(
    list(
      prepared_cache_key = prepared$cache_key,
      scaling = config$scaling,
      validation = config$validation[c("inner_folds", "nrounds", "numCores")],
      label_column = label_column,
      seed = seed,
      test_hash = test_hash
    ),
    algo = "md5"
  )
  fit_cache <- list(key = fit_cache_key, direct = NULL, d1 = list())
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$key, fit_cache_key)) {
      fit_cache <- cached
    }
  }
  save_fit_cache <- function() {
    if (!is.null(cache_path)) {
      saveRDS(fit_cache, cache_path)
    }
  }
  direct_fit <- fit_cache$direct
  if (is.null(direct_fit)) {
    direct_result <- .ablation_linear_readout(
      train = reference_direct,
      test = query_direct,
      train_metadata = prepared$reference_metadata,
      test_metadata = prepared$query_metadata,
      label_column = label_column,
      lambda = config$scaling$lambda,
      inner_folds = config$validation$inner_folds,
      nrounds = config$validation$nrounds,
      numCores = config$validation$numCores,
      seed = readout_seed,
      blocks = NULL
    )
    direct_fit <- list(
      overall = direct_result$overall,
      selected_lambda = direct_result$selected_lambda
    )
    fit_cache$direct <- direct_fit
    save_fit_cache()
  }

  rows <- list()
  index <- 1L
  for (sequence_id in seq_along(sequences)) {
    sequence_seed <- seed + sequence_id
    module_order <- sequences[[sequence_id]]
    if (verbose) {
      luckyBase::LuckyVerbose(
        "ablation: representation scaling sequence ",
        sequence_id,
        "/",
        length(sequences),
        "..."
      )
    }
    for (module_count in counts) {
      prefix_ids <- module_order[seq_len(module_count)]
      module_ids <- modules$module_id[modules$module_id %in% prefix_ids]
      selected_module_hash <- digest::digest(sort(module_ids), algo = "md5")
      full_blocks <- prepared$selected_blocks[module_ids]
      full_columns <- unlist(full_blocks, use.names = FALSE)
      block_ends <- cumsum(lengths(full_blocks))
      block_starts <- c(1L, utils::head(block_ends, -1L) + 1L)
      blocks <- Map(seq.int, block_starts, block_ends)
      names(blocks) <- module_ids
      d1_fit <- fit_cache$d1[[selected_module_hash]]
      if (is.null(d1_fit)) {
        d1_result <- .ablation_linear_readout(
          train = prepared$reference_d1[, full_columns, drop = FALSE],
          test = prepared$query_d1[, full_columns, drop = FALSE],
          train_metadata = prepared$reference_metadata,
          test_metadata = prepared$query_metadata,
          label_column = label_column,
          lambda = config$scaling$lambda,
          inner_folds = config$validation$inner_folds,
          nrounds = config$validation$nrounds,
          numCores = config$validation$numCores,
          seed = readout_seed,
          blocks = blocks
        )
        d1_fit <- list(
          overall = d1_result$overall,
          selected_lambda = d1_result$selected_lambda
        )
        fit_cache$d1[[selected_module_hash]] <- d1_fit
        save_fit_cache()
      }
      rows[[index]] <- data.frame(
        sequence_id = sprintf("S%03d", sequence_id),
        module_sequence_seed = sequence_seed,
        stochastic_seed = readout_seed,
        module_count = module_count,
        module_fraction = module_count / nrow(modules),
        d1_feature_count = length(full_columns),
        direct_feature_count = ncol(reference_direct),
        balanced_accuracy_direct = direct_fit$overall$balanced_accuracy,
        balanced_accuracy_d1 = d1_fit$overall$balanced_accuracy,
        delta_balanced_accuracy = d1_fit$overall$balanced_accuracy -
          direct_fit$overall$balanced_accuracy,
        macro_auroc_direct = direct_fit$overall$macro_auroc,
        macro_auroc_d1 = d1_fit$overall$macro_auroc,
        delta_macro_auroc = d1_fit$overall$macro_auroc -
          direct_fit$overall$macro_auroc,
        selected_lambda_direct = direct_fit$selected_lambda,
        selected_lambda_d1 = d1_fit$selected_lambda,
        module_sequence_hash = digest::digest(module_order, algo = "md5"),
        selected_module_hash = selected_module_hash,
        test_sample_hash = test_hash,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  metrics <- do.call(rbind, rows)
  rownames(metrics) <- NULL
  summaries <- .ablation_representation_scaling_summary_v1(
    metrics,
    config$scaling$bootstrap,
    seed + 100000L
  )
  list(
    status = "complete",
    direct_group = if (config$scaling$direct_feature_type == "gene_pair") {
      "Direct-GSClassifier-TSP"
    } else {
      "Direct-GSClassifier"
    },
    d1_group = "Cohort-d1",
    direct_feature_type = config$scaling$direct_feature_type,
    module_counts = counts,
    sequences = sequences,
    metrics = metrics,
    pointwise = summaries$pointwise,
    sequence_trends = summaries$sequence_trends,
    trend = summaries$trend,
    test_sample_hash = test_hash
  )
}


# Build auditable tissue-breadth, within-tissue-depth, and matched-size banks.
# Every stochastic repeat owns an independent composition sequence; samples are
# never treated as scaling replicates.
.ablation_cohort_bank_design <- function(
    modules,
    module_counts,
    repeats,
    seed
) {
  required <- c("module_id", "tissue", "cohort")
  if (!all(required %in% colnames(modules))) {
    stop("ablation: module manifest is missing bank-design fields.", call. = FALSE)
  }
  modules <- modules[, required, drop = FALSE]
  modules[] <- lapply(modules, as.character)
  if (anyNA(modules) || any(!nzchar(as.matrix(modules)))) {
    stop("ablation: module bank contains missing identifiers.", call. = FALSE)
  }
  if (anyDuplicated(modules$module_id)) {
    stop("ablation: module IDs must be unique within the bank.", call. = FALSE)
  }

  requested_counts <- sort(unique(as.integer(module_counts)))
  requested_counts <- requested_counts[
    requested_counts > 0L & requested_counts <= nrow(modules)
  ]
  rows <- list()
  excluded <- list()
  row_index <- 1L
  excluded_index <- 1L

  make_row <- function(
      design_id,
      family,
      role,
      repeat_id,
      level,
      module_ids,
      pair_id = NA_character_,
      parent_design_id = NA_character_
  ) {
    selected <- modules[match(module_ids, modules$module_id), , drop = FALSE]
    cohort_counts <- table(selected$tissue)
    data.frame(
      design_id = design_id,
      design_family = family,
      design_role = role,
      repeat_id = repeat_id,
      level = as.integer(level),
      pair_id = pair_id,
      parent_design_id = parent_design_id,
      module_count = length(module_ids),
      tissue_count = length(unique(selected$tissue)),
      min_cohort_depth = min(as.integer(cohort_counts)),
      max_cohort_depth = max(as.integer(cohort_counts)),
      mean_cohort_depth = mean(as.integer(cohort_counts)),
      module_hash = digest::digest(sort(module_ids), algo = "md5"),
      module_ids = I(list(as.character(module_ids))),
      tissues = I(list(sort(unique(selected$tissue)))),
      cohorts_per_tissue = I(list(cohort_counts)),
      stringsAsFactors = FALSE
    )
  }

  take_round_robin <- function(tissue_modules, tissue_order, count) {
    selected <- character()
    depth <- 1L
    while (length(selected) < count) {
      added <- unlist(lapply(tissue_order, function(tissue) {
        ids <- tissue_modules[[tissue]]
        if (length(ids) >= depth) ids[depth] else character()
      }), use.names = FALSE)
      selected <- c(selected, added)
      depth <- depth + 1L
    }
    selected[seq_len(count)]
  }

  for (repeat_index in seq_len(as.integer(repeats))) {
    set.seed(seed + repeat_index)
    repeat_id <- sprintf("R%03d", repeat_index)
    tissue_modules <- split(modules$module_id, modules$tissue)
    tissue_modules <- lapply(tissue_modules, sample)
    tissue_order <- sample(names(tissue_modules))
    capacity <- lengths(tissue_modules)

    # Breadth adds one cohort from each new tissue, so within-tissue depth is one.
    parent <- NA_character_
    for (breadth in seq_along(tissue_order)) {
      design_id <- sprintf("%s-B%03d", repeat_id, breadth)
      module_ids <- unlist(
        lapply(tissue_order[seq_len(breadth)], function(x) tissue_modules[[x]][1]),
        use.names = FALSE
      )
      rows[[row_index]] <- make_row(
        design_id,
        "breadth",
        "sequence",
        repeat_id,
        breadth,
        module_ids,
        parent_design_id = parent
      )
      row_index <- row_index + 1L
      parent <- design_id
    }

    # Depth holds the eligible tissue set fixed and adds one cohort per tissue.
    depth_tissues <- tissue_order[capacity[tissue_order] >= 2L]
    if (length(depth_tissues) == 0L) {
      excluded[[excluded_index]] <- data.frame(
        repeat_id = repeat_id,
        design_family = "depth",
        requested_module_count = NA_integer_,
        reason = "no_tissue_has_two_independent_cohorts",
        stringsAsFactors = FALSE
      )
      excluded_index <- excluded_index + 1L
    } else {
      max_depth <- min(capacity[depth_tissues])
      parent <- NA_character_
      for (depth in seq_len(max_depth)) {
        design_id <- sprintf("%s-D%03d", repeat_id, depth)
        module_ids <- unlist(lapply(depth_tissues, function(x) {
          tissue_modules[[x]][seq_len(depth)]
        }), use.names = FALSE)
        rows[[row_index]] <- make_row(
          design_id,
          "depth",
          "sequence",
          repeat_id,
          depth,
          module_ids,
          parent_design_id = parent
        )
        row_index <- row_index + 1L
        parent <- design_id
      }
    }

    # Matched banks have identical module counts but deliberately different
    # tissue diversity. Breadth uses round-robin allocation; depth fills the
    # largest tissue banks first.
    tie_break <- stats::runif(length(capacity))
    depth_order <- names(sort(capacity + tie_break * 1e-6, decreasing = TRUE))
    for (module_count in requested_counts) {
      breadth_ids <- take_round_robin(tissue_modules, tissue_order, module_count)
      depth_ids <- unlist(lapply(depth_order, function(x) tissue_modules[[x]]),
        use.names = FALSE
      )[seq_len(module_count)]
      breadth_tissues <- unique(modules$tissue[match(breadth_ids, modules$module_id)])
      depth_tissues_used <- unique(modules$tissue[match(depth_ids, modules$module_id)])
      if (length(breadth_tissues) <= length(depth_tissues_used)) {
        excluded[[excluded_index]] <- data.frame(
          repeat_id = repeat_id,
          design_family = "matched",
          requested_module_count = module_count,
          reason = "no_tissue_diversity_contrast_at_this_size",
          stringsAsFactors = FALSE
        )
        excluded_index <- excluded_index + 1L
        next
      }
      pair_id <- sprintf("%s-M%03d", repeat_id, module_count)
      rows[[row_index]] <- make_row(
        paste0(pair_id, "-B"),
        "matched",
        "breadth_heavy",
        repeat_id,
        module_count,
        breadth_ids,
        pair_id = pair_id
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- make_row(
        paste0(pair_id, "-D"),
        "matched",
        "depth_heavy",
        repeat_id,
        module_count,
        depth_ids,
        pair_id = pair_id
      )
      row_index <- row_index + 1L
    }
  }

  design <- do.call(rbind, rows)
  rownames(design) <- NULL
  exclusions <- if (length(excluded) == 0L) {
    data.frame(
      repeat_id = character(),
      design_family = character(),
      requested_module_count = integer(),
      reason = character(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, excluded)
  }
  design_hash <- digest::digest(
    lapply(seq_len(nrow(design)), function(i) {
      list(
        design_id = design$design_id[i],
        module_ids = design$module_ids[[i]],
        parent_design_id = design$parent_design_id[i]
      )
    }),
    algo = "md5"
  )
  list(design = design, exclusions = exclusions, design_hash = design_hash)
}


# Resolve one bank into raw matrices plus locally indexed module blocks.
.ablation_cohort_bank_matrices <- function(prepared, module_ids) {
  full_blocks <- prepared$selected_blocks[module_ids]
  if (length(full_blocks) != length(module_ids) || any(lengths(full_blocks) == 0L)) {
    stop("ablation: bank design references unavailable module blocks.", call. = FALSE)
  }
  full_columns <- unlist(full_blocks, use.names = FALSE)
  block_ends <- cumsum(lengths(full_blocks))
  block_starts <- c(1L, utils::head(block_ends, -1L) + 1L)
  blocks <- Map(seq.int, block_starts, block_ends)
  names(blocks) <- module_ids
  list(
    reference = prepared$reference_d1[, full_columns, drop = FALSE],
    query = prepared$query_d1[, full_columns, drop = FALSE],
    blocks = blocks
  )
}


# Represent each probability block by its first centered singular direction.
# Absolute between-module correlations then measure redundant cohort evidence
# without allowing wider probability blocks to receive more weight.
.ablation_module_score_matrix <- function(data, blocks) {
  scores <- lapply(blocks, function(block) {
    candidate <- scale(data[, block, drop = FALSE], center = TRUE, scale = FALSE)
    if (!any(candidate != 0)) {
      return(rep(0, nrow(candidate)))
    }
    fit <- base::svd(candidate, nu = 1L, nv = 0L)
    as.numeric(fit$u[, 1L] * fit$d[1L])
  })
  result <- do.call(cbind, scores)
  colnames(result) <- names(blocks)
  result
}


.ablation_neighbor_index_jaccard <- function(x_neighbors, y_neighbors) {
  if (!identical(dim(x_neighbors), dim(y_neighbors))) {
    stop("ablation: neighbor-index matrices must have identical dimensions.",
      call. = FALSE
    )
  }
  mean(vapply(seq_len(nrow(x_neighbors)), function(i) {
    length(intersect(x_neighbors[i, ], y_neighbors[i, ])) /
      length(union(x_neighbors[i, ], y_neighbors[i, ]))
  }, numeric(1)))
}


.ablation_scaling_metric_row <- function(
    design,
    metric_name,
    metric_role,
    estimate,
    status = "evaluated",
    reason = NA_character_,
    query_coverage = "all"
) {
  estimate <- as.numeric(estimate)[1]
  if (status == "evaluated" && !is.finite(estimate)) {
    status <- "not_evaluated"
    reason <- if (is.na(reason)) "non_finite_metric_estimate" else reason
    estimate <- NA_real_
  }
  data.frame(
    design_id = design$design_id,
    design_family = design$design_family,
    design_role = design$design_role,
    repeat_id = design$repeat_id,
    level = design$level,
    pair_id = design$pair_id,
    module_count = design$module_count,
    tissue_count = design$tissue_count,
    metric_name = metric_name,
    metric_role = metric_role,
    query_coverage = query_coverage,
    estimate = estimate,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}


# Calculate bank-level information and representation-rewrite diagnostics.
.ablation_score_bank_geometry <- function(
    design,
    bank,
    balanced,
    full_balanced,
    full_neighbors,
    config,
    seed
) {
  rank_max <- min(nrow(balanced$query) - 1L, ncol(balanced$query))
  effective_rank <- .ablation_effective_rank(balanced$query)
  module_scores <- .ablation_module_score_matrix(bank$reference, bank$blocks)
  redundancy <- if (ncol(module_scores) < 2L) {
    NA_real_
  } else {
    correlation <- stats::cor(module_scores)
    mean(abs(correlation[upper.tri(correlation)]), na.rm = TRUE)
  }
  module_variance <- vapply(bank$blocks, function(block) {
    mean(apply(bank$reference[, block, drop = FALSE], 2, stats::var))
  }, numeric(1))
  variance_share <- module_variance / sum(module_variance)
  local_k <- min(max(config$geometry$k), nrow(balanced$query) - 1L)
  bank_neighbors <- .ablation_knn(balanced$query, local_k)
  values <- c(
    normalized_effective_rank = effective_rank / rank_max,
    module_covariance_redundancy = redundancy,
    module_variance_concentration = sum(variance_share^2),
    cka_to_full_d1 = .ablation_linear_cka(
      balanced$query,
      full_balanced$query
    ),
    distance_spearman_to_full_d1 = .ablation_distance_spearman(
      balanced$query,
      full_balanced$query,
      min(10000L, config$geometry$distance_pairs),
      seed
    ),
    knn_jaccard_to_full_d1 = .ablation_neighbor_index_jaccard(
      bank_neighbors,
      full_neighbors
    )
  )
  roles <- c(
    rep("primary_nonredundancy", 3L),
    rep("diagnostic_mechanism", 3L)
  )
  metrics <- do.call(rbind, lapply(seq_along(values), function(i) {
    .ablation_scaling_metric_row(
      design,
      names(values)[i],
      roles[i],
      values[i]
    )
  }))
  list(metrics = metrics, neighbors = bank_neighbors)
}


.ablation_scaling_coverage <- function(design, bank_tissues, bank_design, metadata) {
  if (!"tissue" %in% colnames(metadata)) {
    return(rep("all", nrow(metadata)))
  }
  current <- as.character(metadata$tissue) %in% bank_tissues
  parent_tissues <- character()
  if (!is.na(design$parent_design_id) && nzchar(design$parent_design_id)) {
    parent_index <- match(design$parent_design_id, bank_design$design_id)
    if (!is.na(parent_index)) parent_tissues <- bank_design$tissues[[parent_index]]
  }
  if (length(parent_tissues) == 0L) {
    return(ifelse(current, "covered", "uncovered"))
  }
  parent <- as.character(metadata$tissue) %in% parent_tissues
  ifelse(current & !parent, "newly_covered", ifelse(parent, "already_covered", "uncovered"))
}


# Convert retrieval output into technical-robustness and lineage diagnostics,
# stratified by whether the bank covers each query tissue.
.ablation_score_bank_retrieval <- function(
    design,
    retrieval,
    coverage,
    technical_columns
) {
  per_sample <- retrieval$per_sample
  per_sample$query_coverage <- rep(coverage, times = length(unique(per_sample$k)))
  strata <- c("all", sort(unique(per_sample$query_coverage)))
  rows <- list()
  index <- 1L
  for (coverage_name in strata) {
    selected <- if (coverage_name == "all") {
      per_sample
    } else {
      per_sample[per_sample$query_coverage == coverage_name, , drop = FALSE]
    }
    if (nrow(selected) == 0L) next
    for (k_i in sort(unique(selected$k))) {
      data <- selected[selected$k == k_i, , drop = FALSE]
      lineage <- c(
        top1 = mean(data$top1_label_match, na.rm = TRUE),
        top_k = mean(data$top_k_label_rate, na.rm = TRUE),
        mrr = mean(data$mrr, na.rm = TRUE)
      )
      for (name in names(lineage)) {
        rows[[index]] <- .ablation_scaling_metric_row(
          design,
          paste0("lineage_", name, "@", k_i),
          "diagnostic_lineage",
          lineage[name],
          query_coverage = coverage_name
        )
        index <- index + 1L
      }
      for (column in technical_columns) {
        metric <- paste0(column, "_match_excess")
        if (!metric %in% colnames(data)) next
        rows[[index]] <- .ablation_scaling_metric_row(
          design,
          paste0("technical_neighbor_excess:", column, "@", k_i),
          "primary_technical",
          mean(data[[metric]], na.rm = TRUE),
          query_coverage = coverage_name
        )
        index <- index + 1L
      }
    }
  }
  do.call(rbind, rows)
}


.ablation_scaling_direct_contracts <- function(prepared, config) {
  types <- c(main = config$scaling$direct_feature_type)
  if (config$scaling$sensitivity_feature_type != "none") {
    types <- c(types, sensitivity = config$scaling$sensitivity_feature_type)
  }
  lapply(seq_along(types), function(i) {
    type <- unname(types[i])
    features <- if (type == "gene_pair") {
      prepared$feature_manifest$tsp_features
    } else {
      prepared$feature_manifest$features
    }
    columns <- match(features, colnames(prepared$reference_direct))
    if (anyNA(columns)) {
      stop("ablation: a Direct scaling contract is missing features.", call. = FALSE)
    }
    list(
      contract_role = names(types)[i],
      group = if (type == "gene_pair") {
        "Direct-GSClassifier-TSP"
      } else {
        "Direct-GSClassifier"
      },
      feature_type = type,
      features = features,
      feature_count = length(columns),
      feature_hash = digest::digest(features, algo = "md5"),
      reference = prepared$reference_direct[, columns, drop = FALSE],
      query = prepared$query_direct[, columns, drop = FALSE]
    )
  })
}


# Estimate breadth/depth slopes, their interaction, and matched-size paired
# differences. Bootstrap summaries resample repeat-level estimates only.
.ablation_representation_scaling_summary <- function(
    metrics,
    design,
    bootstrap,
    seed
) {
  base <- metrics[
    metrics$status == "evaluated" & metrics$query_coverage == "all",
    ,
    drop = FALSE
  ]
  base <- merge(
    base,
    design[, c(
      "design_id", "parent_design_id", "mean_cohort_depth"
    )],
    by = "design_id",
    all.x = TRUE,
    sort = FALSE
  )
  units <- list()
  unit_index <- 1L

  for (family in c("breadth", "depth")) {
    selected <- base[base$design_family == family, , drop = FALSE]
    groups <- split(selected, interaction(
      selected$repeat_id,
      selected$metric_name,
      drop = TRUE
    ))
    for (data in groups) {
      if (nrow(data) < 2L || length(unique(data$level)) < 2L) next
      fit <- stats::lm(estimate ~ level, data = data)
      units[[unit_index]] <- data.frame(
        contrast_type = paste0(family, "_slope"),
        aggregation = "repeat",
        repeat_id = data$repeat_id[1],
        pair_id = NA_character_,
        metric_name = data$metric_name[1],
        component = if (family == "breadth") {
          "per_added_tissue"
        } else {
          "per_added_cohort_per_tissue"
        },
        estimate = unname(stats::coef(fit)[2]),
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_repeats = 1L,
        stringsAsFactors = FALSE
      )
      unit_index <- unit_index + 1L
    }
  }

  matched <- base[base$design_family == "matched", , drop = FALSE]
  breadth <- matched[matched$design_role == "breadth_heavy", , drop = FALSE]
  depth <- matched[matched$design_role == "depth_heavy", , drop = FALSE]
  paired <- merge(
    breadth,
    depth,
    by = c("repeat_id", "pair_id", "metric_name"),
    suffixes = c("_breadth", "_depth")
  )
  if (nrow(paired) > 0L) {
    for (i in seq_len(nrow(paired))) {
      units[[unit_index]] <- data.frame(
        contrast_type = "matched_size",
        aggregation = "pair",
        repeat_id = paired$repeat_id[i],
        pair_id = paired$pair_id[i],
        metric_name = paired$metric_name[i],
        component = "breadth_minus_depth",
        estimate = paired$estimate_breadth[i] - paired$estimate_depth[i],
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_repeats = 1L,
        stringsAsFactors = FALSE
      )
      unit_index <- unit_index + 1L
    }
  }

  marginal_current <- base[
    !is.na(base$parent_design_id) & nzchar(base$parent_design_id),
    ,
    drop = FALSE
  ]
  marginal_parent <- base[, c("design_id", "metric_name", "estimate"), drop = FALSE]
  names(marginal_parent) <- c("parent_design_id", "metric_name", "parent_estimate")
  marginal <- merge(
    marginal_current,
    marginal_parent,
    by = c("parent_design_id", "metric_name")
  )
  if (nrow(marginal) > 0L) {
    for (i in seq_len(nrow(marginal))) {
      units[[unit_index]] <- data.frame(
        contrast_type = "marginal_gain",
        aggregation = "nested_step",
        repeat_id = marginal$repeat_id[i],
        pair_id = marginal$design_id[i],
        metric_name = marginal$metric_name[i],
        component = paste0(marginal$design_family[i], "_child_minus_parent"),
        estimate = marginal$estimate[i] - marginal$parent_estimate[i],
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_repeats = 1L,
        stringsAsFactors = FALSE
      )
      unit_index <- unit_index + 1L
    }
  }

  interaction_groups <- split(matched, interaction(
    matched$repeat_id,
    matched$metric_name,
    drop = TRUE
  ))
  for (data in interaction_groups) {
    if (nrow(data) < 4L || length(unique(data$tissue_count)) < 2L ||
        length(unique(data$mean_cohort_depth)) < 2L) next
    fit <- stats::lm(
      estimate ~ tissue_count * mean_cohort_depth,
      data = data
    )
    coefficient <- stats::coef(fit)["tissue_count:mean_cohort_depth"]
    if (!is.finite(coefficient)) next
    units[[unit_index]] <- data.frame(
      contrast_type = "breadth_depth_interaction",
      aggregation = "repeat",
      repeat_id = data$repeat_id[1],
      pair_id = NA_character_,
      metric_name = data$metric_name[1],
      component = "tissue_count_x_mean_depth",
      estimate = unname(coefficient),
      ci_low = NA_real_,
      ci_high = NA_real_,
      n_repeats = 1L,
      stringsAsFactors = FALSE
    )
    unit_index <- unit_index + 1L
  }

  unit_table <- if (length(units) == 0L) data.frame() else do.call(rbind, units)
  if (nrow(unit_table) == 0L) return(unit_table)
  summary_source <- stats::aggregate(
    estimate ~ contrast_type + repeat_id + metric_name + component,
    data = unit_table,
    FUN = mean
  )
  summary_groups <- split(summary_source, interaction(
    summary_source$contrast_type,
    summary_source$metric_name,
    summary_source$component,
    drop = TRUE
  ))
  summaries <- lapply(seq_along(summary_groups), function(i) {
    data <- summary_groups[[i]]
    interval <- .ablation_bootstrap_mean(
      data$estimate,
      bootstrap,
      seed + i
    )
    data.frame(
      contrast_type = data$contrast_type[1],
      aggregation = "bootstrap_summary",
      repeat_id = NA_character_,
      pair_id = NA_character_,
      metric_name = data$metric_name[1],
      component = data$component[1],
      estimate = mean(data$estimate),
      ci_low = interval[1],
      ci_high = interval[2],
      n_repeats = nrow(data),
      stringsAsFactors = FALSE
    )
  })
  rbind(unit_table, do.call(rbind, summaries))
}


# Score the auditable two-dimensional cohort bank. Primary evidence concerns
# non-redundancy, technical robustness, external biology (when supplied), and
# repeat stability; cancer-type readouts remain explicitly diagnostic.
.ablation_representation_scaling <- function(
    prepared,
    config,
    label_column,
    seed,
    verbose,
    cache_path = NULL
) {
  modules <- prepared$module_manifest$modules
  modules <- modules[match(prepared$selected_module_ids, modules$module_id), , drop = FALSE]
  bank_design <- .ablation_cohort_bank_design(
    modules,
    config$scaling$module_counts,
    config$scaling$sequences,
    seed
  )
  design <- bank_design$design
  design$d1_feature_count <- vapply(design$module_ids, function(module_ids) {
    sum(lengths(prepared$selected_blocks[module_ids]))
  }, integer(1))
  if ("tissue" %in% colnames(prepared$query_metadata)) {
    query_tissue <- as.character(prepared$query_metadata$tissue)
    design$query_covered_count <- vapply(design$tissues, function(tissues) {
      sum(query_tissue %in% tissues)
    }, integer(1))
    design$query_covered_fraction <- design$query_covered_count / length(query_tissue)
    design$query_coverage_hash <- vapply(design$tissues, function(tissues) {
      digest::digest(
        sort(prepared$query_metadata$sample_id[query_tissue %in% tissues]),
        algo = "md5"
      )
    }, character(1))
  } else {
    design$query_covered_count <- NA_integer_
    design$query_covered_fraction <- NA_real_
    design$query_coverage_hash <- NA_character_
  }
  query_hash <- digest::digest(sort(prepared$query_metadata$sample_id), algo = "md5")
  contracts <- .ablation_scaling_direct_contracts(prepared, config)
  score_reference_metadata <- .ablation_limit_metadata(
    prepared$reference_metadata,
    config$scaling$score_reference_samples,
    seed + 101L
  )
  score_query_metadata <- .ablation_limit_metadata(
    prepared$query_metadata,
    config$scaling$score_query_samples,
    seed + 102L
  )
  score_reference_rows <- match(
    score_reference_metadata$sample_id,
    prepared$reference_metadata$sample_id
  )
  score_query_rows <- match(
    score_query_metadata$sample_id,
    prepared$query_metadata$sample_id
  )
  score_reference_hash <- digest::digest(
    sort(score_reference_metadata$sample_id),
    algo = "md5"
  )
  score_query_hash <- digest::digest(
    sort(score_query_metadata$sample_id),
    algo = "md5"
  )
  cache_key <- digest::digest(list(
    prepared_cache_key = prepared$cache_key,
    bank_design_hash = bank_design$design_hash,
    direct_contracts = lapply(contracts, function(x) x$feature_hash),
    query_hash = query_hash,
    score_reference_hash = score_reference_hash,
    score_query_hash = score_query_hash,
    geometry = config$geometry,
    scaling = config$scaling,
    seed = seed
  ), algo = "md5")
  fit_cache <- list(key = cache_key, direct = list(), d1 = list())
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (identical(cached$key, cache_key)) fit_cache <- cached
  }
  save_fit_cache <- function() {
    if (!is.null(cache_path)) saveRDS(fit_cache, cache_path)
  }

  readout_seed <- seed + 50000L
  direct_diagnostics <- lapply(contracts, function(contract) {
    fit <- fit_cache$direct[[contract$feature_hash]]
    if (is.null(fit)) {
      result <- .ablation_linear_readout(
        train = contract$reference[score_reference_rows, , drop = FALSE],
        test = contract$query[score_query_rows, , drop = FALSE],
        train_metadata = score_reference_metadata,
        test_metadata = score_query_metadata,
        label_column = label_column,
        lambda = config$scaling$lambda,
        inner_folds = config$validation$inner_folds,
        nrounds = config$validation$nrounds,
        numCores = config$validation$numCores,
        seed = readout_seed,
        blocks = NULL
      )
      fit <- list(overall = result$overall, selected_lambda = result$selected_lambda)
      fit_cache$direct[[contract$feature_hash]] <- fit
      save_fit_cache()
    }
    data.frame(
      contract_role = contract$contract_role,
      group = contract$group,
      feature_type = contract$feature_type,
      feature_count = contract$feature_count,
      feature_hash = contract$feature_hash,
      balanced_accuracy = fit$overall$balanced_accuracy,
      macro_auroc = fit$overall$macro_auroc,
      selected_lambda = fit$selected_lambda,
      stringsAsFactors = FALSE
    )
  })
  direct_diagnostics <- do.call(rbind, direct_diagnostics)

  full_bank <- .ablation_cohort_bank_matrices(
    prepared,
    prepared$selected_module_ids
  )
  full_bank$reference <- full_bank$reference[score_reference_rows, , drop = FALSE]
  full_bank$query <- full_bank$query[score_query_rows, , drop = FALSE]
  full_balanced <- .ablation_module_balanced_transform(
    full_bank$reference,
    full_bank$query,
    full_bank$blocks
  )
  stability_k <- min(max(config$geometry$k), nrow(score_query_metadata) - 1L)
  full_neighbors <- .ablation_knn(full_balanced$query, stability_k)
  technical_columns <- intersect(
    config$anchors$technical,
    intersect(
      colnames(score_reference_metadata),
      colnames(score_query_metadata)
    )
  )
  biology_anchors <- intersect(
    config$scaling$biology_anchors,
    intersect(
      colnames(score_reference_metadata),
      colnames(score_query_metadata)
    )
  )
  missing_biology <- setdiff(config$scaling$biology_anchors, biology_anchors)
  reasons <- if (length(config$scaling$biology_anchors) == 0L) {
    data.frame(
      evidence_layer = "external_biology",
      status = "not_evaluated",
      reason = "no_independent_biology_anchor_configured",
      stringsAsFactors = FALSE
    )
  } else if (length(missing_biology) > 0L) {
    data.frame(
      evidence_layer = paste0("external_biology:", missing_biology),
      status = "not_evaluated",
      reason = "anchor_missing_from_reference_or_query_metadata",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      evidence_layer = character(),
      status = character(),
      reason = character(),
      stringsAsFactors = FALSE
    )
  }

  metric_parts <- list()
  neighbor_indices <- list()
  metric_index <- 1L
  local_k <- config$geometry$k[config$geometry$k < nrow(score_reference_metadata)]
  if (length(local_k) == 0L) {
    stop("ablation: no scaling retrieval k is feasible.", call. = FALSE)
  }
  for (i in seq_len(nrow(design))) {
    design_row <- design[i, , drop = FALSE]
    module_ids <- design$module_ids[[i]]
    if (verbose) {
      luckyBase::LuckyVerbose(
        "ablation: cohort bank ", i, "/", nrow(design), " (",
        design_row$design_id, ")..."
      )
    }
    bank <- .ablation_cohort_bank_matrices(prepared, module_ids)
    bank$reference <- bank$reference[score_reference_rows, , drop = FALSE]
    bank$query <- bank$query[score_query_rows, , drop = FALSE]
    balanced <- .ablation_module_balanced_transform(
      bank$reference,
      bank$query,
      bank$blocks
    )
    geometry <- .ablation_score_bank_geometry(
      design_row,
      bank,
      balanced,
      full_balanced,
      full_neighbors,
      config,
      seed + i
    )
    neighbor_indices[[design_row$design_id]] <- geometry$neighbors
    metric_parts[[metric_index]] <- geometry$metrics
    metric_index <- metric_index + 1L

    retrieval <- .ablation_query_reference_retrieval(
      balanced$reference,
      balanced$query,
      score_reference_metadata,
      score_query_metadata,
      label_column = label_column,
      technical_columns = technical_columns,
      k = local_k,
      search = config$geometry$search,
      seed = seed + 1000L + i,
      n_trees = config$geometry$n_trees,
      search_k = config$geometry$search_k
    )
    coverage <- .ablation_scaling_coverage(
      design_row,
      design$tissues[[i]],
      design,
      score_query_metadata
    )
    metric_parts[[metric_index]] <- .ablation_score_bank_retrieval(
      design_row,
      retrieval,
      coverage,
      technical_columns
    )
    metric_index <- metric_index + 1L

    for (anchor in biology_anchors) {
      biology <- .ablation_query_reference_retrieval(
        balanced$reference,
        balanced$query,
        score_reference_metadata,
        score_query_metadata,
        label_column = anchor,
        k = local_k,
        search = config$geometry$search,
        seed = seed + 2000L + i,
        n_trees = config$geometry$n_trees,
        search_k = config$geometry$search_k
      )
      for (k_i in sort(unique(biology$per_sample$k))) {
        values <- biology$per_sample$top_k_label_rate[biology$per_sample$k == k_i]
        metric_parts[[metric_index]] <- .ablation_scaling_metric_row(
          design_row,
          paste0("biology_neighbor_consistency:", anchor, "@", k_i),
          "primary_biology",
          mean(values, na.rm = TRUE)
        )
        metric_index <- metric_index + 1L
      }
    }

    # Supervised cancer-type readout is deliberately restricted to matched banks.
    if (design_row$design_family == "matched") {
      module_hash <- design_row$module_hash
      d1_fit <- fit_cache$d1[[module_hash]]
      if (is.null(d1_fit)) {
        result <- .ablation_linear_readout(
          train = bank$reference,
          test = bank$query,
          train_metadata = score_reference_metadata,
          test_metadata = score_query_metadata,
          label_column = label_column,
          lambda = config$scaling$lambda,
          inner_folds = config$validation$inner_folds,
          nrounds = config$validation$nrounds,
          numCores = config$validation$numCores,
          seed = readout_seed,
          blocks = bank$blocks
        )
        d1_fit <- list(overall = result$overall, selected_lambda = result$selected_lambda)
        fit_cache$d1[[module_hash]] <- d1_fit
        save_fit_cache()
      }
      metric_parts[[metric_index]] <- rbind(
        .ablation_scaling_metric_row(
          design_row,
          "lineage_balanced_accuracy_d1",
          "diagnostic_lineage",
          d1_fit$overall$balanced_accuracy
        ),
        .ablation_scaling_metric_row(
          design_row,
          "lineage_macro_auroc_d1",
          "diagnostic_lineage",
          d1_fit$overall$macro_auroc
        )
      )
      metric_index <- metric_index + 1L
    }
  }

  # Neighbor stability compares independent bank compositions at the same design cell.
  stability_groups <- split(seq_len(nrow(design)), interaction(
    design$design_family,
    design$design_role,
    design$level,
    design$module_count,
    drop = TRUE
  ))
  for (indices in stability_groups) {
    if (length(indices) < 2L) next
    for (i in indices) {
      peers <- setdiff(indices, i)
      values <- vapply(peers, function(j) {
        .ablation_neighbor_index_jaccard(
          neighbor_indices[[design$design_id[i]]],
          neighbor_indices[[design$design_id[j]]]
        )
      }, numeric(1))
      metric_parts[[metric_index]] <- .ablation_scaling_metric_row(
        design[i, , drop = FALSE],
        "neighbor_stability",
        "primary_stability",
        mean(values)
      )
      metric_index <- metric_index + 1L
    }
  }

  metrics <- do.call(rbind, metric_parts)
  rownames(metrics) <- NULL
  contrasts <- .ablation_representation_scaling_summary(
    metrics,
    design,
    config$scaling$bootstrap,
    seed + 100000L
  )
  list(
    schema_version = 2L,
    status = "complete",
    reasons = reasons,
    design_hash = bank_design$design_hash,
    design = design,
    design_exclusions = bank_design$exclusions,
    metrics = metrics,
    contrasts = contrasts,
    diagnostics = list(
      direct_contracts = direct_diagnostics,
      lineage_role = "diagnostic",
      mechanism_role = "diagnostic",
      query_hash = query_hash,
      score_reference_hash = score_reference_hash,
      score_query_hash = score_query_hash,
      score_reference_count = nrow(score_reference_metadata),
      score_query_count = nrow(score_query_metadata),
      cache_key = cache_key,
      biology_anchors = biology_anchors
    ),
    module_counts = sort(unique(design$module_count)),
    test_sample_hash = query_hash
  )
}


# Decode the complete Direct feature contract from d1 and evaluate each feature
# type with its native loss. The decoder is trained only on the reference atlas.
.ablation_decode_direct_features <- function(
    reference_d1,
    query_d1,
    reference_direct,
    query_direct,
    feature_manifest,
    blocks,
    rank = 50L,
    lambda = 1
) {
  reference_d1 <- as.matrix(reference_d1)
  query_d1 <- as.matrix(query_d1)
  reference_direct <- as.matrix(reference_direct)
  query_direct <- as.matrix(query_direct)
  if (!identical(colnames(reference_direct), colnames(query_direct))) {
    stop("ablation: decoder Direct columns must be identical.", call. = FALSE)
  }
  feature_manifest <- feature_manifest[
    match(colnames(reference_direct), feature_manifest$feature),
    ,
    drop = FALSE
  ]
  if (any(is.na(feature_manifest$feature_type))) {
    stop("ablation: decoder feature manifest is incomplete.", call. = FALSE)
  }

  balanced <- .ablation_module_balanced_transform(
    reference_d1,
    query_d1,
    blocks
  )
  rank <- min(
    as.integer(rank),
    ncol(balanced$reference),
    nrow(balanced$reference) - 1L
  )
  if (rank < 1) {
    stop("ablation: decoder rank is not estimable.", call. = FALSE)
  }
  if (rank < min(dim(balanced$reference))) {
    pca <- irlba::prcomp_irlba(
      balanced$reference,
      n = rank,
      center = FALSE,
      scale. = FALSE
    )
  } else {
    pca <- stats::prcomp(
      balanced$reference,
      center = FALSE,
      scale. = FALSE,
      rank. = rank
    )
  }
  rotation <- pca$rotation[, seq_len(rank), drop = FALSE]
  z_reference <- balanced$reference %*% rotation
  z_query <- balanced$query %*% rotation

  outcome_center <- colMeans(reference_direct)
  centered_outcome <- sweep(reference_direct, 2, outcome_center, "-")
  penalty <- diag(lambda, ncol(z_reference))
  coefficient <- solve(
    crossprod(z_reference) + penalty,
    crossprod(z_reference, centered_outcome)
  )
  predicted <- sweep(z_query %*% coefficient, 2, outcome_center, "+")
  colnames(predicted) <- colnames(query_direct)
  rownames(predicted) <- rownames(query_direct)

  mean_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else mean(x)
  }
  per_feature <- do.call(rbind, lapply(seq_len(ncol(query_direct)), function(i) {
    truth <- query_direct[, i]
    estimate <- predicted[, i]
    type <- feature_manifest$feature_type[i]
    balanced_accuracy <- NA_real_
    brier <- NA_real_
    spearman <- NA_real_
    mae <- NA_real_
    rmse <- NA_real_
    if (type == "gene_pair") {
      probability <- pmin(1, pmax(0, estimate))
      call <- as.integer(probability >= 0.5)
      recalls <- vapply(c(0, 1), function(class_id) {
        rows <- truth == class_id
        if (!any(rows)) NA_real_ else mean(call[rows] == class_id)
      }, numeric(1))
      balanced_accuracy <- mean_finite(recalls)
      brier <- mean((probability - truth)^2)
    } else if (type == "single_bin") {
      spearman <- suppressWarnings(stats::cor(truth, estimate, method = "spearman"))
      mae <- mean(abs(estimate - truth))
    } else if (type == "set_pair") {
      spearman <- suppressWarnings(stats::cor(truth, estimate, method = "spearman"))
      rmse <- sqrt(mean((estimate - truth)^2))
    }
    data.frame(
      feature = feature_manifest$feature[i],
      feature_type = type,
      balanced_accuracy = balanced_accuracy,
      brier = brier,
      spearman = spearman,
      mae = mae,
      rmse = rmse,
      stringsAsFactors = FALSE
    )
  }))
  summary <- do.call(rbind, lapply(split(per_feature, per_feature$feature_type), function(data) {
    data.frame(
      feature_type = data$feature_type[1],
      feature_count = nrow(data),
      balanced_accuracy = mean_finite(data$balanced_accuracy),
      brier = mean_finite(data$brier),
      spearman = mean_finite(data$spearman),
      mae = mean_finite(data$mae),
      rmse = mean_finite(data$rmse),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  list(
    status = "complete",
    rank = rank,
    lambda = lambda,
    per_feature = per_feature,
    summary = summary,
    prediction = predicted
  )
}
