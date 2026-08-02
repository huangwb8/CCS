#' Run layered ablation experiments for a CCS model
#'
#' @description
#' Evaluate the frozen CCS representation without retraining cohort submodels.
#' Direct features are reconstructed with GSClassifier's native preprocessing
#' contract and the exact feature support stored in the frozen models.
#' The function implements the cohort-representation, cohort-axis scaling,
#' tissue-first and end-to-end metaCCS experiments defined by the CCS layered
#' ablation plans. All
#' stochastic operations are paired across comparison groups and recorded in
#' an audit table.
#'
#' @param object A `CCS` object containing the frozen d1 representation.
#' @param data Raw RNA expression data. A tissue/cohort nested list is preferred;
#'   each leaf can be an expression matrix or a list containing `expr`.
#' @param metadata Optional sample annotation with sample, cohort, tissue and
#'   biological-label columns. Common CCS column names are recognized.
#' @param experiment One or more of `"cohort"`, `"scaling"`,
#'   `"tissue_first"` and `"metaccs"`.
#' @param output.dir Independent output directory. Existing CCS products are
#'   never overwritten.
#' @param params Named nested list recursively merged onto
#'   `.ablation_default_params(seed)`. Callers may provide only the leaves they
#'   want to change, for example
#'   `list(general = list(rank = 25L), scaling = list(gate = list(enforce = TRUE)))`.
#'   The canonical top-level groups are `general`, `cohort`, `scaling`,
#'   `tissue_first`, and `metaccs`. Existing flat parameter names remain accepted
#'   and are normalized to this schema; unknown or ambiguous duplicate fields are
#'   rejected before computation.
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
    experiment = c("cohort", "scaling", "tissue_first", "metaccs"),
    output.dir = file.path(getwd(), "ccs-ablation"),
    params = list(),
    seed = 20260727,
    verbose = TRUE
) {
  # Test
  if (FALSE) {
    # Purpose: Exercise all four ablation experiments with one aligned synthetic CCS fixture.
    # Input: Generated expression, tissue/cohort models, frozen d1 probabilities, and metadata.
    # Parameters: Small paired repeats cover every stochastic path without a production-scale run.
    # Output: A CCSAblation result containing cohort, scaling, tissue_first, and metaccs results.
    luckyBase::Plus.library(c("CCS", "digest"))
    set.seed(101)

    sample_ids <- sprintf("S%03d", seq_len(48))
    cohort <- rep(paste0("C", 1:4), each = 12)
    tissue <- rep(c("T1", "T1", "T2", "T2"), each = 12)
    biology <- rep(rep(c("B1", "B2"), each = 6), 4)

    expr <- matrix(
      stats::rnorm(4 * length(sample_ids)),
      nrow = 4,
      dimnames = list(paste0("g", 1:4), sample_ids)
    )
    expr["g1", biology == "B1"] <- expr["g1", biology == "B1"] + 1
    expr["g4", biology == "B2"] <- expr["g4", biology == "B2"] + 1

    data <- list(
      T1 = list(
        C1 = list(
          expr = expr[, cohort == "C1", drop = FALSE],
          subtype = biology[cohort == "C1"]
        ),
        C2 = list(
          expr = expr[, cohort == "C2", drop = FALSE],
          subtype = biology[cohort == "C2"]
        )
      ),
      T2 = list(
        C3 = list(
          expr = expr[, cohort == "C3", drop = FALSE],
          subtype = biology[cohort == "C3"]
        ),
        C4 = list(
          expr = expr[, cohort == "C4", drop = FALSE],
          subtype = biology[cohort == "C4"]
        )
      )
    )

    tsp_features <- apply(
      utils::combn(rownames(expr), 2),
      2,
      paste,
      collapse = ":"
    )
    tsp <- vapply(strsplit(tsp_features, ":", fixed = TRUE), function(pair) {
      as.integer(expr[pair[1], ] >= expr[pair[2], ])
    }, integer(length(sample_ids)))
    rownames(tsp) <- sample_ids

    # Give each frozen cohort model explicit TSP support. Tissue unions differ so
    # metaCCS can verify that Direct follows the matching tissue model space.
    module_features <- list(
      `T1|C1` = tsp_features[c(1, 2)],
      `T1|C2` = tsp_features[c(3, 4)],
      `T2|C3` = tsp_features[c(5, 6)],
      `T2|C4` = tsp_features[c(3, 6)]
    )
    module_ids <- names(module_features)
    d1 <- do.call(cbind, lapply(seq_along(module_ids), function(i) {
      feature_index <- match(module_features[[i]][1], tsp_features)
      score <- stats::plogis(
        tsp[, feature_index] + stats::rnorm(nrow(tsp), 0, 0.1)
      )
      block <- cbind(`1` = score, `2` = 1 - score)
      colnames(block) <- paste(module_ids[i], colnames(block), sep = "|")
      block
    }))
    rownames(d1) <- sample_ids

    make_test_model <- function(features) {
      class_model <- list(
        bst = NULL,
        breakVec = c(0, 0.5, 1),
        genes = features
      )
      list(
        Repeat = list(),
        Model = list(list(`1` = class_model, `2` = class_model))
      )
    }
    models <- list(
      T1 = list(
        C1 = make_test_model(module_features$`T1|C1`),
        C2 = make_test_model(module_features$`T1|C2`)
      ),
      T2 = list(
        C3 = make_test_model(module_features$`T2|C3`),
        C4 = make_test_model(module_features$`T2|C4`)
      )
    )

    object <- methods::new(
      "CCS",
      Repeat = list(
        method = "GSClassifier",
        geneSet = list(A = c("g1", "g2"), B = c("g3", "g4")),
        geneAnnotation = data.frame(ENSEMBL = rownames(expr)),
        geneid = "ensembl",
        seed = 101,
        model.dir = ""
      ),
      Model = models,
      Data = list(
        Probability = list(
          d1 = d1,
          d2 = matrix(0, length(sample_ids), 4),
          d3 = matrix(0, length(sample_ids), 2)
        ),
        CCS = stats::setNames(
          rep(1:2, length.out = length(sample_ids)),
          sample_ids
        ),
        CancerType = stats::setNames(tissue, sample_ids)
      )
    )

    metadata <- data.frame(
      sample_id = sample_ids,
      cohort = cohort,
      tissue = tissue,
      biology = biology,
      stringsAsFactors = FALSE
    )
    experiment <- c("cohort", "scaling", "tissue_first", "metaccs")
    output.dir <- file.path(tempdir(), "ccs-ablation-example")
    params <- list(
      general = list(
        rank = 3L,
        k = 3L,
        n_folds = 2L,
        bootstrap = 5L,
        probe = FALSE,
        fidelity_samples = length(sample_ids),
        dr = list(
          method = "UWOT",
          dimension = c(3L, 2L),
          n_neighbors = 5L,
          min_dist = 0.1,
          spread = 1,
          set_op_mix_ratio = 1,
          metric = "euclidean",
          n_threads = 1L
        ),
        cluster = list(eps = 0.5, minPts = 3L),
        cover = TRUE
      ),
      cohort = list(
        rank_sensitivity = c(2L, 3L),
        geometry_samples = length(sample_ids),
        distance_pairs = 100L,
        mechanism_samples = 24L,
        rp_seeds = c(701L, 702L),
        permutation_seeds = c(801L, 802L)
      ),
      scaling = list(
        counts = c(2L, 4L),
        sequences = 2L,
        embedding_counts = c(2L, 4L),
        embedding_sequences = 2L,
        embedding_seeds = c(901L, 902L),
        subsample_fraction = 1,
        gate = list(enforce = FALSE)
      ),
      tissue_first = list(
        seeds = c(1001L, 1002L),
        subsample_fraction = 1
      ),
      metaccs = list(
        resample_seeds = 1101L,
        umap_seeds = c(1201L, 1202L),
        subsample_fraction = 1,
        parameter_mode = "fixed",
        direct_feature_mode = "tissue_model_union",
        retain_assignments = TRUE
      )
    )
    seed <- 1301L
    verbose <- FALSE
  }

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
.ablation_grouped_folds <- function(cohort, n_folds = Inf, seed = 20260727) {
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
  for (group in ordered) {
    fold <- which.min(fold_load)
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
