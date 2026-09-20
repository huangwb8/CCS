#' Run representation ablation experiments for a CCS model
#'
#' @description
#' Evaluate the frozen CCS representation without retraining cohort submodels.
#' The representation workflow reconstructs the complete native
#' GSClassifier input and reuses only matching precomputed d1 rows. Filtered
#' query samples absent from the object are reported and excluded,
#' and compares Direct-GSClassifier with Cohort-d1 on independent
#' query-to-reference retrieval, grouped linear readout, paired learning curves,
#' null controls, and feature-type reconstruction. All stochastic operations
#' are paired across comparison groups and recorded in an audit table.
#'
#' @param object A `CCS` object containing the frozen d1 representation.
#' @param data Raw RNA expression data. A tissue/cohort nested list is preferred;
#'   each leaf can be an expression matrix or a list containing `expr`.
#' @param metadata Optional sample annotation with sample, cohort, tissue and
#'   biological-label columns. Common CCS column names are recognized.
#' @param output.dir Independent output directory. Existing CCS products are
#'   never overwritten.
#' @param params Named nested list. Values are merged onto
#'   `.ablation_representation_default_params(seed)` under `comparison`,
#'   `provenance`, `anchors`, `geometry`, `validation`, `controls`, `tradeoffs`,
#'   `scaling`, and `output`. Unknown fields are rejected before computation.
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
#'     \item{`validation$numCores`}{Positive XGBoost thread count used by each
#'       representation readout worker. Default: `1L`.}
#'     \item{`validation$workers`}{Positive PSOCK worker count for independent
#'       learning-curve jobs. The default `1L` preserves serial execution;
#'       callers should divide their total CPU budget between workers and
#'       `validation$numCores` to avoid oversubscription.}
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
#'     \item{`output$cache_direct`}{Logical; persist the reconstructed
#'       Direct-GSClassifier matrix in `direct-feature-cache.rds` and reuse it
#'       only when the expression, sample, model, and feature-contract hashes
#'       all match. Default: `TRUE`.}
#'     \item{`general$cover`}{Logical; allow writing into a non-empty `output.dir` and
#'       replacing same-named ablation result files. Default: `FALSE`.}
#'   }
#' @param seed Master random seed.
#' @param verbose Whether to report progress.
#' @param step Lifecycle stage: `"all"` (default), `"context"`, `"plan"`,
#'   `"run"`, or `"result"`. Staged calls consume the object supplied through
#'   `input` and never write reviewer-facing products unless `step = "all"`.
#' @param input Optional serializable stage input. It must be a context, plan,
#'   runner result, or a list containing one of those objects as required by
#'   `step`.
#' @param cache.root Optional root directory for reusable intermediate cache.
#'   It is kept separate from `output.dir`; when omitted a transient directory
#'   under `tempdir()` is used and the current working directory is untouched.
#'
#' @return An object of class `CCSAblation`.
#' @author Weibin Huang <hwb2012@@qq.com>
#' @md
#' @export
ablation <- function(
    object = NULL,
    data = NULL,
    metadata = NULL,
    output.dir = file.path(getwd(), "ccs-ablation"),
    params = list(),
    seed = 20260727,
    verbose = TRUE,
    step = "all",
    input = NULL,
    cache.root = NULL
) {
  step <- match.arg(step, c("all", "context", "plan", "run", "result"))
  if (!identical(step, "all")) {
    return(.ablation_dispatch_step(
      step = step,
      object = object,
      data = data,
      metadata = metadata,
      output.dir = output.dir,
      params = params,
      seed = seed,
      verbose = verbose,
      input = input,
      cache.root = cache.root
    ))
  }

  .ablation_run_representation(
    object = object,
    data = data,
    metadata = metadata,
    output.dir = output.dir,
    params = params,
    seed = seed,
    verbose = verbose,
    cache.root = cache.root
  )
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

# Keep frozen module IDs while using audited tissue names for bank designs.
.ablation_resolve_bank_tissues <- function(module_manifest, metadata) {
  lookup <- unique(metadata[, c("cohort", "tissue"), drop = FALSE])
  if (anyDuplicated(lookup$cohort)) {
    stop("ablation: cohort-to-tissue mapping is ambiguous.", call. = FALSE)
  }
  modules <- module_manifest$modules
  mapped <- lookup$tissue[match(modules$cohort, lookup$cohort)]
  modules$bank_tissue <- modules$tissue
  use <- !is.na(mapped) & nzchar(mapped)
  modules$tissue[use] <- as.character(mapped[use])
  if (any(modules$tissue == "Undefined")) {
    stop("ablation: bank design contains unresolved tissue labels.", call. = FALSE)
  }
  module_manifest$modules <- modules
  module_manifest
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


# Build a content-addressed key for the expensive Direct reconstruction.  Hash
# only the expression rows that can affect GSClassifier features; hashing the
# complete RNA matrix would recreate the same memory/time bottleneck we are
# trying to avoid.
.ablation_direct_expression_fingerprint <- function(object, expr) {
  gene_candidates <- unique(as.character(unlist(object@Repeat$geneSet, use.names = FALSE)))
  row_ids <- rownames(expr)
  relevant <- rep(TRUE, length(row_ids))
  if (length(gene_candidates) > 0L && any(nzchar(gene_candidates))) {
    annotation <- object@Repeat$geneAnnotation
    geneid <- as.character(object@Repeat$geneid)[1L]
    annotation_names <- if (is.data.frame(annotation)) {
      tolower(colnames(annotation))
    } else {
      character()
    }
    annotation_column <- match(tolower(geneid), annotation_names)
    annotation_values <- if (is.data.frame(annotation) &&
        nzchar(geneid) && !is.na(annotation_column) &&
        nrow(annotation) == length(row_ids)) {
      as.character(annotation[[annotation_column]])
    } else {
      row_ids
    }
    relevant <- annotation_values %in% gene_candidates | row_ids %in% gene_candidates
    if (!any(relevant)) relevant <- rep(TRUE, length(row_ids))
  }
  digest::digest(
    list(
      dim = dim(expr),
      rownames = row_ids,
      colnames = colnames(expr),
      relevant_rows = which(relevant),
      relevant_values = expr[relevant, , drop = FALSE]
    ),
    algo = "xxhash64"
  )
}


.ablation_direct_feature_cache_key <- function(object, expr, feature_manifest, sample_ids) {
  digest::digest(
    list(
      version = 2L,
      method = object@Repeat$method,
      feature_builder = "GSClassifier::trainDataProc_X",
      feature_builder_version = as.character(utils::packageVersion("GSClassifier")),
      geneSet = object@Repeat$geneSet,
      geneAnnotation = object@Repeat$geneAnnotation,
      geneid = object@Repeat$geneid,
      feature_hash = digest::digest(feature_manifest, algo = "md5"),
      sample_ids = as.character(sample_ids),
      expression_hash = .ablation_direct_expression_fingerprint(object, expr)
    ),
    algo = "md5"
  )
}


# Save generated cache files through a temporary sibling and rename so an
# interrupted R session cannot leave a partially written cache that looks
# valid to the next run.
.ablation_atomic_save_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  saveRDS(value, tmp, compress = FALSE)
  backup <- paste0(path, ".bak-", Sys.getpid())
  if (file.exists(path) && !file.rename(path, backup)) {
    stop("ablation: cannot stage existing cache file for replacement: ", path, call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("ablation: cannot finalize cache file: ", path, call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  invisible(path)
}

# Persist one recoverable computation node under a content key. Each key gets
# its own immutable value file and a small state file, so interrupted work is
# visible as "running" and can never be mistaken for a complete cache hit.
.ablation_node_code_identity <- function(node) {
  roots <- switch(
    node,
    retrieval = c(
      ".ablation_scale_train_apply",
      ".ablation_module_balanced_transform",
      ".ablation_query_reference_retrieval",
      ".ablation_validate_neighbor_search",
      ".ablation_bind_retrieval",
      ".ablation_evidence_level",
      ".ablation_null_perm_eligibility",
      ".ablation_projection_matrix"
    ),
    readout = c(".ablation_linear_readout"),
    `learning-curve` = c(".ablation_learning_curve"),
    decoder = c(
      ".ablation_limit_metadata",
      ".ablation_decode_direct_features"
    ),
    `cohort-scaling` = c(".ablation_representation_scaling"),
    stop("ablation: unknown checkpoint node: ", node, call. = FALSE)
  )
  runtime_packages <- switch(
    node,
    retrieval = c("RcppAnnoy"),
    readout = c("xgboost"),
    `learning-curve` = c("xgboost"),
    decoder = c("irlba"),
    `cohort-scaling` = c("xgboost", "RcppAnnoy"),
    stop("ablation: unknown checkpoint node: ", node, call. = FALSE)
  )
  # Bind each node to its own recursive implementation graph. This preserves
  # transitive invalidation without making an unrelated helper change evict all
  # recoverable nodes.
  implementation_environment <- environment(.ablation_node_code_identity)
  functions <- character()
  pending <- roots
  while (length(pending) > 0L) {
    name <- pending[[1L]]
    pending <- pending[-1L]
    if (name %in% functions) next
    if (!exists(name, envir = implementation_environment, inherits = FALSE) ||
        !is.function(get(name, envir = implementation_environment, inherits = FALSE))) {
      stop("ablation: missing checkpoint implementation: ", name, call. = FALSE)
    }
    functions <- c(functions, name)
    fun <- get(name, envir = implementation_environment, inherits = FALSE)
    referenced <- unique(all.names(body(fun), functions = TRUE))
    referenced <- referenced[grepl("^\\.ablation_", referenced)]
    referenced <- referenced[vapply(referenced, function(candidate) {
      exists(candidate, envir = implementation_environment, inherits = FALSE) &&
        is.function(get(candidate, envir = implementation_environment, inherits = FALSE))
    }, logical(1L))]
    pending <- unique(c(pending, referenced))
  }
  functions <- sort(functions)
  definitions <- lapply(functions, function(name) {
    fun <- get(name, envir = implementation_environment, inherits = FALSE)
    list(name = name, formals = formals(fun), body = body(fun))
  })
  list(
    functions = digest::digest(definitions, algo = "md5"),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    digest_version = as.character(utils::packageVersion("digest")),
    runtime = stats::setNames(
      vapply(runtime_packages, function(pkg) {
        as.character(utils::packageVersion(pkg))
      }, character(1L)),
      runtime_packages
    )
  )
}

.ablation_node_cache_key <- function(
    node,
    prepared,
    parameters,
    seed,
    algorithm_revision,
    upstream = list()) {
  digest::digest(
    list(
      schema_version = 1L,
      node = node,
      input_key = prepared$input_key,
      representation_key = prepared$cache_key,
      direct_feature_key = prepared$direct_cache$key,
      parameters = parameters,
      seed = as.integer(seed),
      algorithm_revision = algorithm_revision,
      code = .ablation_node_code_identity(node),
      upstream = upstream
    ),
    algo = "md5"
  )
}

.ablation_inspect_node_cache <- function(path, key, state_path = NULL, node = NULL) {
  miss <- function(reason) list(value = NULL, reason = reason)
  if (!file.exists(path)) return(miss("missing-cache"))

  if (!is.null(state_path)) {
    if (!file.exists(state_path)) return(miss("missing-state"))
    state <- tryCatch(readRDS(state_path), error = function(error) NULL)
    if (is.null(state)) return(miss("unreadable-state"))
    if (!is.list(state) || length(state$schema_version) != 1L ||
        !state$schema_version %in% c(1L, 2L)) {
      return(miss("state-schema-mismatch"))
    }
    if (!identical(state$status, "complete")) {
      status <- if (is.character(state$status) && length(state$status) == 1L) {
        state$status
      } else "not-complete"
      return(miss(paste0("state-", status)))
    }
    if (!identical(state$key, key)) return(miss("state-key-mismatch"))
    if (!is.null(node) && !identical(state$node, node)) {
      return(miss("state-node-mismatch"))
    }
    cache_md5 <- unname(tools::md5sum(path))
    if (!is.character(state$cache_md5) || length(state$cache_md5) != 1L ||
        is.na(cache_md5) || !identical(state$cache_md5, cache_md5)) {
      return(miss("cache-file-hash-mismatch"))
    }
  }

  cached <- tryCatch(readRDS(path), error = function(error) NULL)
  if (is.null(cached)) return(miss("unreadable-cache"))
  if (!is.list(cached) || !identical(cached$schema_version, 1L)) {
    return(miss("cache-schema-mismatch"))
  }
  if (!identical(cached$status, "complete")) return(miss("cache-not-complete"))
  if (!identical(cached$key, key)) return(miss("cache-key-mismatch"))
  if (!is.null(node) && !identical(cached$node, node)) {
    return(miss("cache-node-mismatch"))
  }
  if (is.null(cached$value) || !is.character(cached$value_hash) ||
      length(cached$value_hash) != 1L) {
    return(miss("cache-value-contract-mismatch"))
  }
  if (!identical(digest::digest(cached$value, algo = "md5"), cached$value_hash)) {
    return(miss("cache-value-hash-mismatch"))
  }
  list(value = cached$value, reason = "valid")
}

.ablation_read_node_cache <- function(path, key, state_path = NULL, node = NULL) {
  .ablation_inspect_node_cache(path, key, state_path, node)$value
}

.ablation_read_fit_cache <- function(path, key) {
  if (!file.exists(path)) return(NULL)
  cached <- tryCatch(readRDS(path), error = function(error) NULL)
  fit_is_valid <- function(fit) {
    is.list(fit) && is.data.frame(fit$overall) && nrow(fit$overall) >= 1L &&
      is.numeric(fit$selected_lambda) && length(fit$selected_lambda) == 1L &&
      is.finite(fit$selected_lambda)
  }
  entries_are_valid <- function(entries, allow_single = FALSE) {
    if (is.null(entries)) return(isTRUE(allow_single))
    if (isTRUE(allow_single) && fit_is_valid(entries)) return(TRUE)
    is.list(entries) && all(vapply(entries, fit_is_valid, logical(1L)))
  }
  if (!is.list(cached) || !identical(cached$schema_version, 1L) ||
      !identical(cached$status, "complete") || !identical(cached$key, key) ||
      !is.list(cached$direct) || !is.list(cached$d1)) {
    return(NULL)
  }
  if (!entries_are_valid(cached$direct, allow_single = TRUE) ||
      !entries_are_valid(cached$d1)) return(NULL)
  cached
}

.ablation_runtime_identity <- function() {
  hostname <- unname(Sys.info()[["nodename"]])
  if (is.null(hostname) || is.na(hostname) || !nzchar(hostname)) {
    hostname <- Sys.getenv("COMPUTERNAME", unset = "unknown")
  }
  run_id <- Sys.getenv("CCS_ABLATION_RUN_ID", unset = "")
  if (!nzchar(run_id)) {
    run_id <- paste0("r-", Sys.getpid(), "-", format(Sys.time(), "%Y%m%dT%H%M%S"))
  }
  list(
    run_id = run_id,
    pid = Sys.getpid(),
    hostname = hostname
  )
}

.ablation_pid_is_alive <- function(pid) {
  if (length(pid) != 1L || is.na(pid) || !is.finite(pid) || pid < 1L) {
    return(FALSE)
  }
  pid <- as.integer(pid)
  if (identical(.Platform$OS.type, "windows")) {
    output <- tryCatch(
      suppressWarnings(system2(
        "tasklist", c("/FI", shQuote(paste("PID eq", pid)), "/NH"),
        stdout = TRUE, stderr = FALSE
      )),
      error = function(error) character()
    )
    return(any(grepl(paste0("\\b", pid, "\\b"), output)))
  }
  isTRUE(tryCatch({
    status <- system2("kill", c("-0", as.character(pid)), stdout = FALSE, stderr = FALSE)
    identical(status, 0L)
  }, error = function(error) FALSE))
}

.ablation_running_state_is_active <- function(state) {
  if (!is.list(state) || !identical(state$status, "running")) return(FALSE)
  if (!is.character(state$hostname) || length(state$hostname) != 1L ||
      is.na(state$hostname) || !nzchar(state$hostname) ||
      length(state$pid) != 1L || is.na(state$pid) || !is.finite(state$pid)) {
    return(FALSE)
  }
  identity <- .ablation_runtime_identity()
  if (!identical(state$hostname, identity$hostname)) {
    # A process on another host cannot be checked safely. It remains active
    # until an operator explicitly marks it stale, avoiding time-based theft of
    # a legitimately long computation.
    return(TRUE)
  }
  .ablation_pid_is_alive(state$pid)
}

.ablation_condition_summary <- function(condition, limit = 1000L) {
  message <- conditionMessage(condition)
  message <- gsub("[\r\n\t]+", " ", message)
  message <- gsub("[[:cntrl:]]", "", message)
  substr(message, 1L, as.integer(limit))
}

.ablation_process_memory <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(list(working_set_bytes = NA_real_, peak_working_set_bytes = NA_real_))
  }
  info <- tryCatch({
    handle <- getExportedValue("ps", "ps_handle")()
    getExportedValue("ps", "ps_memory_info")(handle)
  }, error = function(error) NULL)
  if (is.null(info)) {
    return(list(working_set_bytes = NA_real_, peak_working_set_bytes = NA_real_))
  }
  rss <- if ("rss" %in% names(info)) as.numeric(info[["rss"]]) else NA_real_
  peak <- if ("peak_wset" %in% names(info)) {
    as.numeric(info[["peak_wset"]])
  } else rss
  list(working_set_bytes = rss, peak_working_set_bytes = peak)
}

.ablation_cached_node <- function(
    node,
    output.dir,
    key,
    compute,
    verbose = TRUE,
    job = NULL,
    parameter_digest = NULL) {
  if (!grepl("^[a-z][a-z0-9-]*$", node)) {
    stop("ablation: invalid checkpoint node name: ", node, call. = FALSE)
  }
  node_dir <- file.path(output.dir, "checkpoints", node)
  dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(node_dir, paste0(key, ".rds"))
  state_path <- file.path(node_dir, paste0(key, ".state.rds"))
  inspection <- .ablation_inspect_node_cache(
    cache_path, key, state_path = state_path, node = node
  )
  if (!is.null(inspection$value)) {
    if (verbose) luckyBase::LuckyVerbose("ablation: checkpoint hit: ", node)
    return(list(value = inspection$value, status = "hit", lookup_status = "hit",
      reason = inspection$reason, key = key,
      path = file.path("checkpoints", node, basename(cache_path))))
  }

  previous_state <- if (file.exists(state_path)) {
    tryCatch(readRDS(state_path), error = function(error) NULL)
  } else NULL
  if (is.list(previous_state) && identical(previous_state$status, "running")) {
    if (.ablation_running_state_is_active(previous_state)) {
      stop(
        "ablation: checkpoint is active in run ", previous_state$run_id,
        " on ", previous_state$hostname, " (PID ", previous_state$pid, "): ", node,
        call. = FALSE
      )
    }
    previous_state$status <- "stale"
    previous_state$stale_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    previous_state$stale_reason <- "owner-process-not-active"
    .ablation_atomic_save_rds(previous_state, state_path)
    inspection$reason <- "stale-running-state"
  }

  identity <- .ablation_runtime_identity()
  started_at <- Sys.time()
  started_memory <- .ablation_process_memory()
  running_state <- c(
    list(
      schema_version = 2L,
      status = "running",
      node = node,
      key = key,
      started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
      updated_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
      parameter_digest = parameter_digest,
      job = job,
      working_set_start_bytes = started_memory$working_set_bytes
    ),
    identity
  )
  if (is.list(previous_state) && !identical(previous_state$status, "complete")) {
    running_state$recovered_from <- previous_state[c(
      "status", "run_id", "pid", "hostname", "failed_at", "stale_at",
      "stale_reason", "error_class", "error_summary"
    )]
  }
  .ablation_atomic_save_rds(
    running_state,
    state_path
  )
  if (verbose) {
    luckyBase::LuckyVerbose(
      "ablation: checkpoint miss: ", node, " (", inspection$reason, ")"
    )
  }
  value <- tryCatch(
    compute(),
    error = function(error) {
      failed_at <- Sys.time()
      failed_state <- running_state
      failed_state$status <- "failed"
      failed_state$updated_at <- format(failed_at, "%Y-%m-%dT%H:%M:%S%z")
      failed_state$failed_at <- failed_state$updated_at
      failed_state$elapsed_seconds <- as.numeric(difftime(failed_at, started_at, units = "secs"))
      failed_state$error_class <- class(error)[1L]
      failed_state$error_summary <- .ablation_condition_summary(error)
      .ablation_atomic_save_rds(failed_state, state_path)
      stop(error)
    },
    interrupt = function(interrupt) {
      failed_at <- Sys.time()
      failed_state <- running_state
      failed_state$status <- "failed"
      failed_state$updated_at <- format(failed_at, "%Y-%m-%dT%H:%M:%S%z")
      failed_state$failed_at <- failed_state$updated_at
      failed_state$elapsed_seconds <- as.numeric(difftime(failed_at, started_at, units = "secs"))
      failed_state$error_class <- "interrupt"
      failed_state$error_summary <- .ablation_condition_summary(interrupt)
      .ablation_atomic_save_rds(failed_state, state_path)
      stop(interrupt)
    }
  )
  .ablation_atomic_save_rds(
    list(
      schema_version = 1L,
      status = "complete",
      node = node,
      key = key,
      value_hash = digest::digest(value, algo = "md5"),
      value = value
    ),
    cache_path
  )
  completed_at <- Sys.time()
  completed_memory <- .ablation_process_memory()
  complete_state <- running_state
  complete_state$status <- "complete"
  complete_state$updated_at <- format(completed_at, "%Y-%m-%dT%H:%M:%S%z")
  complete_state$completed_at <- complete_state$updated_at
  complete_state$elapsed_seconds <- as.numeric(difftime(completed_at, started_at, units = "secs"))
  complete_state$result_bytes <- as.numeric(utils::object.size(value))
  complete_state$working_set_end_bytes <- completed_memory$working_set_bytes
  complete_state$peak_working_set_bytes <- completed_memory$peak_working_set_bytes
  complete_state$cache_md5 <- unname(tools::md5sum(cache_path))
  .ablation_atomic_save_rds(complete_state, state_path)
  list(value = value, status = "written", lookup_status = "miss",
    reason = inspection$reason, key = key,
    path = file.path("checkpoints", node, basename(cache_path)))
}

.ablation_atomic_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(value, temporary, row.names = FALSE)
  invisible(utils::read.csv(temporary, nrows = 5L, check.names = FALSE))
  backup <- paste0(path, ".bak-", Sys.getpid())
  if (file.exists(path) && !file.rename(path, backup)) {
    stop("ablation: cannot stage CSV output for replacement: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("ablation: cannot finalize CSV output: ", path, call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup, force = TRUE)
  invisible(path)
}


.ablation_read_direct_feature_cache <- function(path, key, sample_ids, feature_manifest) {
  if (!file.exists(path)) return(NULL)
  cached <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(cached) || !identical(cached$key, key)) return(NULL)
  direct <- cached$direct
  if (!is.matrix(direct) ||
      !identical(rownames(direct), as.character(sample_ids)) ||
      !identical(colnames(direct), feature_manifest$features) ||
      any(dim(direct) != c(length(sample_ids), length(feature_manifest$features))) ||
      any(!is.finite(direct))) {
    return(NULL)
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

  if (nrow(query) == 0L) {
    neighbors <- data.frame(
      query_sample = character(), query_cohort = character(),
      query_label = character(), neighbor_rank = integer(),
      reference_sample = character(), reference_cohort = character(),
      reference_label = character(), distance = numeric(),
      label_match = logical(), stringsAsFactors = FALSE
    )
    per_sample <- do.call(rbind, lapply(k, function(k_i) {
      result <- data.frame(
        sample_id = character(), cohort = character(), label = character(),
        k = integer(), top1_label_match = numeric(),
        top_k_label_rate = numeric(), mrr = numeric(),
        stringsAsFactors = FALSE
      )
      for (column in technical_columns) {
        result[[paste0(column, "_match_rate")]] <- numeric()
        result[[paste0(column, "_expected_rate")]] <- numeric()
        result[[paste0(column, "_match_excess")]] <- numeric()
      }
      result
    }))
    summary <- data.frame(
      k = k, top1_label_match = NA_real_,
      top_k_label_rate = NA_real_, mrr = NA_real_
    )
    return(list(
      neighbors = neighbors,
      per_sample = per_sample,
      summary = summary,
      search = search,
      k = k,
      status = "not_estimable",
      reason = "no_estimable_query_cohorts"
    ))
  }

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
      # MRR@k must only credit a relevant neighbor when its first hit is
      # inside the requested top-k list.  Reusing the max-k reciprocal rank
      # for smaller k would leak information from ranks that are not part of
      # that endpoint (and incorrectly make MRR identical for every k).
      mrr = ifelse(
        is.na(first_match) | first_match > k_i,
        0,
        1 / first_match
      ),
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
# generated without exposing either evaluation side to its own fitted model.
.ablation_evidence_level <- function(
    query_metadata,
    anchor_role,
    reference_metadata = query_metadata,
    provenance_column = "d1_provenance"
) {
  if (!provenance_column %in% colnames(query_metadata) ||
      !provenance_column %in% colnames(reference_metadata)) {
    stop("ablation: reference/query metadata is missing d1 provenance.", call. = FALSE)
  }
  query_provenance <- unique(as.character(query_metadata[[provenance_column]]))
  reference_provenance <- unique(as.character(reference_metadata[[provenance_column]]))
  query_qualified <- nrow(query_metadata) > 0L &&
    all(query_provenance %in% c("external_frozen", "out_of_fold"))
  reference_qualified <- nrow(reference_metadata) > 0L &&
    all(reference_provenance %in% c("external_frozen", "out_of_fold"))
  qualified <- query_qualified && reference_qualified
  independent <- identical(anchor_role, "independent")
  reasons <- c(
    if (!independent) "anchor_is_not_independent",
    if (!query_qualified) "query_d1_provenance_is_not_external_or_out_of_fold",
    if (!reference_qualified) "reference_d1_provenance_is_not_external_or_out_of_fold"
  )
  list(
    level = if (qualified && independent) "confirmatory" else "descriptive",
    qualified_provenance = qualified,
    reference_qualified_provenance = reference_qualified,
    query_qualified_provenance = query_qualified,
    anchor_role = anchor_role,
    provenance = unique(c(reference_provenance, query_provenance)),
    reference_provenance = reference_provenance,
    query_provenance = query_provenance,
    reasons = reasons
  )
}


# Defaults are grouped by the scientific question so configuration changes are
# auditable and layered parameters cannot leak into representation tests.
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
      min_reference_cohorts = 2L,
      endpoint_min_reference_cohorts = list(
        cancer_retrieval = 2L,
        technical_excess = 2L,
        cancer_readout = 2L,
        learning_curve = 2L
      )
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
      numCores = 1L,
      workers = 1L
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
      cache_direct = TRUE
    )
  )
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
  endpoint_thresholds <- unlist(config$anchors$endpoint_min_reference_cohorts)
  if (length(endpoint_thresholds) == 0L ||
      any(!is.finite(endpoint_thresholds)) || any(endpoint_thresholds < 1) ||
      any(endpoint_thresholds != as.integer(endpoint_thresholds))) {
    stop(
      "ablation: anchors$endpoint_min_reference_cohorts must contain positive integers.",
      call. = FALSE
    )
  }
  if (length(config$scaling$enabled) != 1 ||
      !is.logical(config$scaling$enabled) ||
      is.na(config$scaling$enabled)) {
    stop("ablation: scaling$enabled must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(config$output$cache_direct) != 1L ||
      !is.logical(config$output$cache_direct) ||
      is.na(config$output$cache_direct)) {
    stop("ablation: output$cache_direct must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(config$validation$workers) != 1L ||
      !is.finite(config$validation$workers) ||
      config$validation$workers < 1L ||
      config$validation$workers != as.integer(config$validation$workers)) {
    stop("ablation: validation$workers must be one positive integer.", call. = FALSE)
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

.ablation_apply_runtime_config <- function(config, analysis = NULL) {
  existing_workers <- config$validation$workers
  if (is.null(existing_workers)) existing_workers <- 1L
  existing_total <- max(
    1L,
    as.integer(config$validation$numCores) * as.integer(existing_workers)
  )
  total_threads <- suppressWarnings(as.integer(Sys.getenv(
    "CCS_ABLATION_CORES",
    unset = as.character(existing_total)
  )))
  workers <- suppressWarnings(as.integer(Sys.getenv(
    "CCS_ABLATION_WORKERS",
    unset = as.character(existing_workers)
  )))
  if (length(total_threads) != 1L || !is.finite(total_threads) || total_threads < 1L) {
    stop("ablation: CCS_ABLATION_CORES must be one positive integer.", call. = FALSE)
  }
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L) {
    stop("ablation: CCS_ABLATION_WORKERS must be one positive integer.", call. = FALSE)
  }
  workers <- min(as.integer(workers), as.integer(total_threads))
  memory_gb <- suppressWarnings(as.numeric(Sys.getenv(
    "CCS_ABLATION_MEMORY_GB",
    unset = "Inf"
  )))
  if (length(memory_gb) != 1L || is.na(memory_gb) || memory_gb <= 0) {
    stop("ablation: CCS_ABLATION_MEMORY_GB must be one positive number.", call. = FALSE)
  }
  worker_bytes <- NA_real_
  if (!is.null(analysis) && is.list(analysis$prepared)) {
    prepared <- analysis$prepared
    matrices <- prepared[c(
      "reference_direct", "query_direct", "reference_d1", "query_d1"
    )]
    worker_bytes <- 3 * sum(vapply(
      matrices,
      function(value) as.numeric(utils::object.size(value)),
      numeric(1L)
    ))
    if (is.finite(memory_gb) && is.finite(worker_bytes) && worker_bytes > 0) {
      memory_workers <- floor(memory_gb * 1024^3 / worker_bytes)
      if (memory_workers < 1L) {
        stop(
          "ablation: CCS_ABLATION_MEMORY_GB cannot accommodate one estimated worker (",
          format(round(worker_bytes / 1024^3, 2), nsmall = 2), " GiB required).",
          call. = FALSE
        )
      }
      workers <- min(workers, memory_workers)
    }
  }
  config$validation$workers <- workers
  config$validation$numCores <- max(1L, floor(total_threads / workers))
  config$validation$memory_gb <- memory_gb
  config$validation$worker_memory_estimate_bytes <- worker_bytes
  config$validation$total_thread_budget <- total_threads
  config
}


.ablation_limit_metadata <- function(metadata, size, seed) {
  if (!is.finite(size) || nrow(metadata) <= size) {
    return(metadata)
  }
  metadata[.ablation_stratified_sample(metadata, as.integer(size), seed), , drop = FALSE]
}


# Prepare disjoint reference/query matrices. Existing d1 rows are reused for
# both partitions; query samples absent from d1 are reported and excluded.
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

  # The ablation API consumes the d1 already prepared in `object`.  Query
  # samples without a matching d1 row are reported and excluded; this
  # boundary keeps raw-data preparation under the caller's control and avoids
  # silently predicting new d1 values inside a downstream analysis.
  precomputed_query_ids <- intersect(query_ids, rownames(d1))
  excluded_query_d1_ids <- setdiff(query_ids, precomputed_query_ids)
  if (length(excluded_query_d1_ids) > 0L) {
    warning(
      "ablation: excluding ", length(excluded_query_d1_ids),
      " query sample(s) without precomputed d1; examples: ",
      paste(utils::head(excluded_query_d1_ids, 5L), collapse = ", "),
      call. = FALSE
    )
    query_ids <- precomputed_query_ids
    query_metadata <- query_metadata[
      match(query_ids, query_metadata$sample_id),
      ,
      drop = FALSE
    ]
  }
  if (length(query_ids) < 1L) {
    stop(
      "ablation: no query samples have precomputed d1 in object.",
      call. = FALSE
    )
  }

  if (verbose) {
    luckyBase::LuckyVerbose(
      "ablation: preparing Direct-GSClassifier for ",
      length(reference_ids),
      " reference and ",
      length(query_ids),
      " external query samples..."
    )
  }
  all_ids <- c(reference_ids, query_ids)
  direct_expr <- flattened$expr[, all_ids, drop = FALSE]
  # The package API can be used as a pure in-memory calculation layer by
  # passing `output.dir = NULL`.  In that mode targets owns persistence and
  # this function must not create package/workflow state on its own.
  direct_cache_enabled <- isTRUE(config$output$cache_direct) &&
    !is.null(output.dir)
  direct_cache_path <- if (direct_cache_enabled) {
    file.path(output.dir, "direct-feature-cache.rds")
  } else {
    NULL
  }
  direct_cache_key <- .ablation_direct_feature_cache_key(
    object = object,
    expr = direct_expr,
    feature_manifest = feature_manifest,
    sample_ids = all_ids
  )
  direct <- if (direct_cache_enabled) {
    .ablation_read_direct_feature_cache(
      path = direct_cache_path,
      key = direct_cache_key,
      sample_ids = all_ids,
      feature_manifest = feature_manifest
    )
  } else {
    NULL
  }
  direct_cache_status <- if (!is.null(direct)) "hit" else "miss"
  if (is.null(direct)) {
    if (verbose) {
      luckyBase::LuckyVerbose(
        "ablation: rebuilding Direct-GSClassifier features (cache miss)..."
      )
    }
    direct <- .ablation_gsclassifier_matrix(
      object,
      direct_expr,
      feature_manifest
    )
    if (direct_cache_enabled) {
      .ablation_atomic_save_rds(
        list(
          schema_version = 1L,
          key = direct_cache_key,
          sample_ids = all_ids,
          feature_manifest = feature_manifest,
          direct = direct
        ),
        direct_cache_path
      )
      direct_cache_status <- "written"
    } else {
      direct_cache_status <- "disabled"
    }
  } else if (verbose) {
    luckyBase::LuckyVerbose(
      "ablation: reusing cached Direct-GSClassifier features from ",
      normalizePath(direct_cache_path, winslash = "/", mustWork = TRUE),
      "."
    )
  }
  reference_direct <- direct[reference_ids, , drop = FALSE]
  query_direct <- direct[query_ids, , drop = FALSE]

  expected_columns <- colnames(d1)[
    unlist(module_manifest$blocks[module_ids], use.names = FALSE)
  ]
  reference_d1 <- d1[reference_ids, expected_columns, drop = FALSE]
  cache_key <- digest::digest(
    list(
      query_ids = query_ids,
      direct = query_direct,
      module_ids = module_ids,
      query_d1 = d1[query_ids, expected_columns, drop = FALSE]
    ),
    algo = "md5"
  )
  query_d1 <- d1[query_ids, expected_columns, drop = FALSE]
  # Prediction caches additionally bind reference values and labels. Keep the
  # representation key stable for the existing independently keyed exact
  # geometry cache, which already includes reference d1 and Direct content.
  input_key <- digest::digest(list(
    cache_key = cache_key, reference_direct = reference_direct,
    reference_d1 = reference_d1, reference_metadata = reference_metadata,
    query_metadata = query_metadata, feature_manifest = feature_manifest
  ), algo = "md5")

  selected_blocks <- lapply(module_ids, function(module_id) {
    module_columns <- colnames(d1)[
      module_manifest$blocks[[module_id]]
    ]
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
    module_manifest = .ablation_resolve_bank_tissues(module_manifest, metadata),
    selected_blocks = selected_blocks,
    selected_module_ids = module_ids,
    feature_manifest = feature_manifest,
    excluded_duplicate_samples = flattened$excluded_duplicate_samples,
    excluded_query_d1_ids = excluded_query_d1_ids,
    direct_cache = list(
      path = direct_cache_path,
      key = direct_cache_key,
      status = direct_cache_status,
      enabled = direct_cache_enabled
    ),
    cache_key = cache_key,
    input_key = input_key,
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


# Cache the exact native-geometry diagnostics separately from downstream
# endpoints. The key binds both representations and the sampling contract.
.ablation_native_geometry_cache_key <- function(prepared, config, seed) {
  digest::digest(
    list(
      schema_version = 1L,
      representation_key = prepared$cache_key,
      direct_key = prepared$direct_cache$key,
      reference_ids = rownames(prepared$reference_d1),
      reference_d1 = prepared$reference_d1,
      geometry = config$geometry,
      seed = seed
    ),
    algo = "md5"
  )
}


.ablation_read_native_geometry_cache <- function(path, key) {
  if (!file.exists(path)) return(NULL)
  cached <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(cached) || !identical(cached$schema_version, 1L) ||
      !identical(cached$key, key) || !is.list(cached$value)) {
    return(NULL)
  }
  cached$value
}


# Promote a legacy exact result only when its companion manifest proves that
# the representation, geometry contract, seed and dimensions all match.
.ablation_promote_legacy_native_geometry <- function(
    output.dir,
    prepared,
    config,
    seed,
    key,
    cache_path
) {
  # Old manifests did not bind reference content. Recompute instead of
  # promoting an unverifiable result after a cache miss.
  return(NULL)

}


.ablation_bind_retrieval <- function(results) {
  neighbors <- do.call(rbind, lapply(names(results), function(name) {
    data <- results[[name]]$neighbors
    data$representation <- rep(name, nrow(data))
    data
  }))
  per_sample <- do.call(rbind, lapply(names(results), function(name) {
    data <- results[[name]]$per_sample
    data$representation <- rep(name, nrow(data))
    data
  }))
  summary <- do.call(rbind, lapply(names(results), function(name) {
    data <- results[[name]]$summary
    data$representation <- rep(name, nrow(data))
    data
  }))
  by_cohort <- if (nrow(per_sample) == 0L) {
    data.frame(
      representation = character(), cohort = character(), k = integer(),
      top1_label_match = numeric(), top_k_label_rate = numeric(), mrr = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    stats::aggregate(
      cbind(top1_label_match, top_k_label_rate, mrr) ~ representation + cohort + k,
      data = per_sample,
      FUN = mean
    )
  }

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


# Build a shared candidate/endpoint qualification table.  Cancer-labelled
# endpoints may consume an estimable view, while label-free endpoints retain
# the complete candidate target set.
.ablation_endpoint_eligibility <- function(
    reference_metadata,
    query_metadata,
    config,
    direction = "reference_bank_to_external_targets"
) {
  required <- c("sample_id", "cohort", "cohort_key", "cancer_type")
  missing_reference <- setdiff(required, colnames(reference_metadata))
  missing_query <- setdiff(required, colnames(query_metadata))
  missing <- c(
    if (length(missing_reference) > 0L) paste0("reference:", missing_reference),
    if (length(missing_query) > 0L) paste0("query:", missing_query)
  )
  if (length(missing) > 0L) {
    stop(
      "ablation: endpoint eligibility metadata is missing: ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  reference_support <- vapply(
    split(reference_metadata$cohort, as.character(reference_metadata$cancer_type)),
    function(x) length(unique(x)),
    integer(1)
  )
  query_cohorts <- unique(query_metadata[, c(
    "cohort_key", "cohort", "cancer_type"
  ), drop = FALSE])
  query_sample_counts <- table(query_metadata$cohort_key)
  endpoint_thresholds <- config$anchors$endpoint_min_reference_cohorts
  if (is.null(endpoint_thresholds)) endpoint_thresholds <- list()
  endpoints <- c(
    "cancer_retrieval", "technical_excess", "cancer_readout",
    "learning_curve", "geometry", "anchor", "structural_reproducibility"
  )
  rows <- do.call(rbind, lapply(seq_len(nrow(query_cohorts)), function(i) {
    label <- as.character(query_cohorts$cancer_type[i])
    support <- unname(reference_support[label])
    if (is.na(support)) support <- 0L
    cancer_qualified <- vapply(
      endpoints,
      function(endpoint) {
        threshold <- endpoint_thresholds[[endpoint]]
        if (is.null(threshold)) threshold <- config$anchors$min_reference_cohorts
        support >= threshold
      },
      logical(1)
    )
    data.frame(
      direction = direction,
      bank_role = if (identical(direction, "reference_bank_to_external_targets"))
        "reference" else "external",
      target_role = if (identical(direction, "reference_bank_to_external_targets"))
        "external" else "reference",
      cohort_key = as.character(query_cohorts$cohort_key[i]),
      cohort = as.character(query_cohorts$cohort[i]),
      cancer_type = label,
      endpoint = endpoints,
      candidate_status = "candidate",
      qualification_status = ifelse(
        endpoints %in% c(
          "cancer_retrieval", "technical_excess", "cancer_readout", "learning_curve"
        ),
        ifelse(cancer_qualified, "estimable", "not_estimable"),
        "estimable"
      ),
      qualification_reason = ifelse(
        endpoints %in% c(
          "cancer_retrieval", "technical_excess", "cancer_readout", "learning_curve"
        ),
        ifelse(
          cancer_qualified,
          "reference_cancer_support_meets_threshold",
          ifelse(support == 0L, "no_reference_cancer_support",
                 "reference_cancer_support_below_threshold")
        ),
        "endpoint_does_not_require_cancer_label"
      ),
      reference_support_cohort_count = as.integer(support),
      target_sample_count = as.integer(unname(query_sample_counts[query_cohorts$cohort_key[i]])),
      fit_target_overlap = 0L,
      stringsAsFactors = FALSE
    )
  }))
  rows$sample_hash <- vapply(rows$cohort_key, function(key) {
    digest::digest(sort(query_metadata$sample_id[query_metadata$cohort_key == key]), algo = "md5")
  }, character(1))
  rows$config_hash <- digest::digest(config, algo = "md5")
  rows
}


.ablation_endpoint_view <- function(prepared, endpoint) {
  audit <- prepared$endpoint_eligibility
  keep_cohorts <- unique(audit$cohort_key[
    audit$endpoint == endpoint & audit$qualification_status == "estimable"
  ])
  rows <- prepared$query_metadata$cohort_key %in% keep_cohorts
  list(
    metadata = prepared$query_metadata[rows, , drop = FALSE],
    direct = prepared$query_direct[rows, , drop = FALSE],
    d1 = prepared$query_d1[rows, , drop = FALSE],
    audit = audit[audit$endpoint == endpoint, , drop = FALSE]
  )
}


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
    object = object,
    data = data,
    metadata = metadata,
    config = config,
    output.dir = output.dir,
    seed = seed,
    verbose = verbose
  )
  anchor <- config$anchors$primary
  if (!anchor %in% colnames(prepared$reference_metadata) ||
      !anchor %in% colnames(prepared$query_metadata)) {
    stop("ablation: primary anchor is missing from metadata.", call. = FALSE)
  }
  prepared$endpoint_eligibility <- .ablation_endpoint_eligibility(
    reference_metadata = prepared$reference_metadata,
    query_metadata = prepared$query_metadata,
    config = config
  )
  prepared$candidate_query_metadata <- prepared$query_metadata
  prepared$candidate_query_direct <- prepared$query_direct
  prepared$candidate_query_d1 <- prepared$query_d1
  prepared$query_views <- setNames(
    lapply(
      unique(prepared$endpoint_eligibility$endpoint),
      function(endpoint) .ablation_endpoint_view(prepared, endpoint)
    ),
    unique(prepared$endpoint_eligibility$endpoint)
  )
  list(prepared = prepared, anchor = anchor)
}


.ablation_run_representation <- function(
    object,
    data,
    metadata,
    output.dir,
    params,
    seed,
    verbose,
    cache.root = NULL
) {
  if (!methods::is(object, "CCS")) {
    stop("ablation: object must be a CCS object.", call. = FALSE)
  }
  config <- .ablation_resolve_representation_config(seed, params)
  cache <- .ablation_resolve_cache_layout(cache.root, output.dir)
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
    output.dir = cache$preparation_dir,
    seed = seed,
    verbose = verbose
  )
  result <- .ablation_run_prepared_representation(
    analysis = analysis, config = config, output.dir = output.dir,
    cache.dir = cache$nodes_dir,
    seed = seed, verbose = verbose
  )
  result$manifest$cache <- cache
  .ablation_atomic_save_rds(result$manifest, file.path(output.dir, "manifest.rds"))
  result$call <- match.call()
  result
}


# Run the same representation calculations from explicitly prepared inputs.
# This boundary lets analysis scripts consume caches without loading raw data.
.ablation_run_prepared_representation <- function(
    analysis, config, output.dir, cache.dir = output.dir, seed, verbose
) {
  # Prepared bundles created before runtime workers were introduced retain the
  # historical single-worker behavior instead of failing on a missing field.
  if (is.null(config$validation$workers)) config$validation$workers <- 1L
  config <- .ablation_apply_runtime_config(config, analysis)
  prepared <- analysis$prepared
  anchor <- analysis$anchor
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

  # Phase 1: put Direct-GSClassifier and Cohort-d1 on their native comparison
  # scales. All later endpoints consume these transformed matrices.
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
  native_geometry_cache_path <- file.path(
    cache.dir,
    "native-geometry-cache.rds"
  )
  native_geometry_cache_key <- .ablation_native_geometry_cache_key(
    prepared,
    config,
    seed
  )
  native_geometry <- .ablation_read_native_geometry_cache(
    native_geometry_cache_path,
    native_geometry_cache_key
  )
  if (is.null(native_geometry)) {
    native_geometry <- .ablation_promote_legacy_native_geometry(
      output.dir = cache.dir,
      prepared = prepared,
      config = config,
      seed = seed,
      key = native_geometry_cache_key,
      cache_path = native_geometry_cache_path
    )
    if (!is.null(native_geometry) && verbose) {
      luckyBase::LuckyVerbose(
        "ablation: promoted matching exact native geometry into cache."
      )
    }
  }
  if (is.null(native_geometry)) {
    if (verbose) {
      luckyBase::LuckyVerbose(
        "ablation: computing exact native geometry (cache miss)..."
      )
    }
    native_geometry <- .ablation_native_geometry(
      prepared,
      transformed,
      config,
      seed
    )
    .ablation_atomic_save_rds(
      list(
        schema_version = 1L,
        key = native_geometry_cache_key,
        value = native_geometry
      ),
      native_geometry_cache_path
    )
  } else if (verbose) {
    luckyBase::LuckyVerbose(
      "ablation: reusing cached exact native geometry from ",
      normalizePath(
        native_geometry_cache_path,
        winslash = "/",
        mustWork = TRUE
      ),
      "."
    )
  }

  # Phase 2: evaluate query-to-reference retrieval. The retrieval view is
  # restricted to estimable cancer-labelled queries, while the complete
  # transformed matrices remain available to structural diagnostics.
  # Cancer-labelled endpoints use their pre-declared estimable view.  The
  # transformed query matrices themselves remain complete candidates so that
  # geometry, continuous anchors and structural diagnostics are not silently
  # truncated by a readout-specific cancer support rule.
  retrieval_node_key <- .ablation_node_cache_key(
    node = "retrieval",
    prepared = prepared,
    parameters = list(
      geometry = config$geometry,
      controls = config$controls,
      anchors = config$anchors,
      endpoint_eligibility = digest::digest(
        prepared$endpoint_eligibility,
        algo = "md5"
      )
    ),
    seed = seed,
    algorithm_revision = "retrieval-v1"
  )
  retrieval_checkpoint <- .ablation_cached_node(
    node = "retrieval",
    output.dir = cache.dir,
    key = retrieval_node_key,
    verbose = verbose,
    compute = function() {
  retrieval_view <- prepared$query_views[["cancer_retrieval"]]
  retrieval_rows <- match(
    retrieval_view$metadata$sample_id,
    prepared$query_metadata$sample_id
  )
  retrieval_metadata <- retrieval_view$metadata
  retrieval_direct <- transformed$direct$query[retrieval_rows, , drop = FALSE]
  retrieval_d1 <- transformed$d1$query[retrieval_rows, , drop = FALSE]

  retrieval_results <- list(
    `Direct-GSClassifier` = .ablation_query_reference_retrieval(
      transformed$direct$reference,
      retrieval_direct,
      prepared$reference_metadata,
      retrieval_metadata,
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
      retrieval_d1,
      prepared$reference_metadata,
      retrieval_metadata,
      label_column = anchor,
      technical_columns = config$anchors$technical,
      k = config$geometry$k,
      search = config$geometry$search,
      seed = seed,
      n_trees = config$geometry$n_trees,
      search_k = config$geometry$search_k
    )
  )
  search_validation <- if (config$geometry$search == "annoy" &&
      nrow(retrieval_metadata) > 0L) {
    validation <- list(
      `Direct-GSClassifier` = .ablation_validate_neighbor_search(
        transformed$direct$reference,
        retrieval_direct,
        prepared$reference_metadata,
        retrieval_metadata,
        label_column = anchor,
        k = max(config$geometry$k),
        query_samples = config$geometry$exact_validation_queries,
        n_trees = config$geometry$n_trees,
        search_k = config$geometry$search_k,
        seed = seed + 1000L
      ),
      `Cohort-d1` = .ablation_validate_neighbor_search(
        transformed$d1$reference,
        retrieval_d1,
        prepared$reference_metadata,
        retrieval_metadata,
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
  } else if (config$geometry$search == "annoy") {
    list(status = "not_estimable", reason = "no_estimable_query_cohorts")
  } else {
    list(status = "not_required", search = "exact")
  }
  retrieval <- .ablation_bind_retrieval(retrieval_results)
  retrieval$search_validation <- search_validation
  # Continuous anchors have no cancer-support requirement. Their neighbors
  # must cover all external candidates, independently of cancer retrieval.
  anchor_results <- lapply(c("Direct-GSClassifier", "Cohort-d1"), function(name) {
    matrices <- if (name == "Direct-GSClassifier") transformed$direct else transformed$d1
    .ablation_query_reference_retrieval(
      matrices$reference, matrices$query,
      prepared$reference_metadata, prepared$query_metadata,
      label_column = anchor, k = config$geometry$k,
      search = config$geometry$search, seed = seed,
      n_trees = config$geometry$n_trees, search_k = config$geometry$search_k
    )
  })
  names(anchor_results) <- c("Direct-GSClassifier", "Cohort-d1")
  anchor_retrieval <- .ablation_bind_retrieval(anchor_results)
  evidence <- .ablation_evidence_level(
    retrieval_metadata,
    config$anchors$primary_role,
    reference_metadata = prepared$reference_metadata
  )

  # Phase 3: run paired null controls against the same retrieval contract.
  # These controls remain separate from the primary comparison in the result.
  null_perm <- if (config$controls$null_perm) {
    .ablation_null_perm_eligibility(retrieval_metadata, anchor)
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
        retrieval_direct %*% projection,
        prepared$reference_metadata,
        retrieval_metadata,
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
  list(
    retrieval = retrieval,
    anchor_retrieval = anchor_retrieval,
    evidence = evidence,
    controls = controls
  )
    }
  )
  retrieval <- retrieval_checkpoint$value$retrieval
  anchor_retrieval <- retrieval_checkpoint$value$anchor_retrieval
  evidence <- retrieval_checkpoint$value$evidence
  controls <- retrieval_checkpoint$value$controls

  # Phase 4: evaluate the supervised cancer readout and learning curves only
  # on the pre-declared estimable query view.
  readout_view <- prepared$query_views[["cancer_readout"]]
  readout_estimable <- config$validation$enabled && nrow(readout_view$metadata) > 0L
  if (readout_estimable) {
    readout_rows <- match(
      readout_view$metadata$sample_id,
      prepared$query_metadata$sample_id
    )
    readout_metadata <- readout_view$metadata
    readout_direct <- prepared$query_direct[readout_rows, , drop = FALSE]
    readout_d1 <- prepared$query_d1[readout_rows, , drop = FALSE]
  }
  readout_node_key <- .ablation_node_cache_key(
    node = "readout",
    prepared = prepared,
    parameters = list(
      enabled = config$validation$enabled,
      lambda = config$validation$lambda,
      inner_folds = config$validation$inner_folds,
      nrounds = config$validation$nrounds,
      numCores = config$validation$numCores,
      anchor = anchor,
      query_view = digest::digest(readout_view$metadata, algo = "md5")
    ),
    seed = seed + 20000L,
    algorithm_revision = "readout-v1"
  )
  readout_checkpoint <- .ablation_cached_node(
    node = "readout",
    output.dir = cache.dir,
    key = readout_node_key,
    verbose = verbose,
    compute = function() {
      if (!readout_estimable) {
        return(list(
          status = if (config$validation$enabled) "not_estimable" else "not_run",
          reason = if (config$validation$enabled) {
            "no_estimable_query_cohorts"
          } else {
            "disabled"
          }
        ))
      }
    readout_results <- list(
      `Direct-GSClassifier` = .ablation_linear_readout(
        train = prepared$reference_direct,
        test = readout_direct,
        train_metadata = prepared$reference_metadata,
        test_metadata = readout_metadata,
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
        test = readout_d1,
        train_metadata = prepared$reference_metadata,
        test_metadata = readout_metadata,
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
    list(
      status = "complete",
      results = readout_results,
      overall = overall,
      by_cohort = by_cohort,
      paired_by_cohort = paired_by_cohort,
      predictions = predictions
    )
    }
  )
  readout <- readout_checkpoint$value

  learning_node_key <- .ablation_node_cache_key(
    node = "learning-curve",
    prepared = prepared,
    parameters = list(
      enabled = config$validation$enabled,
      fractions = config$validation$learning_fractions,
      repeats = config$validation$repeats,
      lambda = config$validation$lambda,
      inner_folds = config$validation$inner_folds,
      nrounds = config$validation$nrounds,
      numCores = config$validation$numCores,
      anchor = anchor,
      query_view = digest::digest(readout_view$metadata, algo = "md5")
    ),
    seed = seed + 30000L,
    algorithm_revision = "learning-curve-v1"
  )
  learning_checkpoint <- .ablation_cached_node(
    node = "learning-curve",
    output.dir = cache.dir,
    key = learning_node_key,
    verbose = verbose,
    compute = function() {
      if (!readout_estimable) return(readout)
      .ablation_learning_curve(
        representations = list(
          `Direct-GSClassifier` = list(
            train = prepared$reference_direct,
            test = readout_direct,
            blocks = NULL
          ),
          `Cohort-d1` = list(
            train = prepared$reference_d1,
            test = readout_d1,
            blocks = prepared$selected_blocks
          )
        ),
        train_metadata = prepared$reference_metadata,
        test_metadata = readout_metadata,
        label_column = anchor,
        fractions = config$validation$learning_fractions,
        repeats = config$validation$repeats,
        lambda = config$validation$lambda,
        inner_folds = config$validation$inner_folds,
        nrounds = config$validation$nrounds,
        numCores = config$validation$numCores,
        workers = config$validation$workers,
        checkpoint_output_dir = cache.dir,
        checkpoint_key = learning_node_key,
        seed = seed + 30000L
      )
    }
  )
  learning_curve <- learning_checkpoint$value

  # Phase 5: optional scaling and decoder diagnostics. These are auxiliary
  # analyses and do not redefine the primary retrieval result.
  cohort_scaling <- if (config$scaling$enabled) {
    .ablation_representation_scaling(
      prepared = prepared,
      config = config,
      label_column = anchor,
      seed = seed + 35000L,
      verbose = verbose,
      cache_path = file.path(cache.dir, "cohort-scaling-fit-cache.rds")
    )
  } else {
    list(status = "not_run", reason = "disabled")
  }

  feature_counts <- table(factor(
    prepared$feature_manifest$feature_manifest$feature_type,
    levels = c("single_bin", "gene_pair", "set_pair")
  ))
  decoder_node_key <- .ablation_node_cache_key(
    node = "decoder",
    prepared = prepared,
    parameters = list(
      enabled = config$tradeoffs$decoder,
      rank = config$tradeoffs$decoder_rank,
      lambda = config$tradeoffs$decoder_lambda,
      max_reference_samples = config$tradeoffs$decoder_max_reference_samples,
      max_query_samples = config$tradeoffs$decoder_max_query_samples
    ),
    seed = seed + 40000L,
    algorithm_revision = "decoder-v1"
  )
  decoder_checkpoint <- .ablation_cached_node(
    node = "decoder",
    output.dir = cache.dir,
    key = decoder_node_key,
    verbose = verbose,
    compute = function() {
      if (!config$tradeoffs$decoder) {
        return(list(status = "not_run", reason = "disabled"))
      }
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
    }
  )
  decoder <- decoder_checkpoint$value
  tradeoffs <- list(
    status = "complete",
    feature_type_count = feature_counts,
    d1_probability_contract = native_geometry$d1_probability_contract,
    module_diagnostics = native_geometry$module_diagnostics,
    excluded_duplicate_samples = prepared$excluded_duplicate_samples,
    decoder = decoder
  )

  # Phase 6: construct the reproducibility manifest and compact audit table.
  # Candidate and estimable query counts are both retained for review.
  endpoint_summary <- stats::aggregate(
    cbind(
      candidate_count = prepared$endpoint_eligibility$candidate_status == "candidate",
      estimable_count = prepared$endpoint_eligibility$qualification_status == "estimable",
      not_estimable_count = prepared$endpoint_eligibility$qualification_status == "not_estimable"
    ) ~ endpoint,
    data = prepared$endpoint_eligibility,
    FUN = sum
  )
  readout_view_for_manifest <- prepared$query_views[["cancer_readout"]]
  manifest <- list(
    version = 7L,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = seed,
    experiment = "representation",
    groups = c("Direct-GSClassifier", "Cohort-d1"),
    reference_sample_count = nrow(prepared$reference_metadata),
    query_sample_count = nrow(readout_view_for_manifest$metadata),
    reference_cohort_count = length(unique(prepared$reference_metadata$cohort)),
    query_cohort_count = length(unique(readout_view_for_manifest$metadata$cohort)),
    external_candidate_sample_count = nrow(prepared$query_metadata),
    external_candidate_cohort_count = length(unique(prepared$query_metadata$cohort)),
    excluded_query_d1_count = length(prepared$excluded_query_d1_ids),
    excluded_query_d1_sample_hash = digest::digest(
      sort(prepared$excluded_query_d1_ids),
      algo = "md5"
    ),
    estimable_cohort_count = length(unique(readout_view_for_manifest$metadata$cohort)),
    not_estimable_cohort_count = sum(
      prepared$endpoint_eligibility$endpoint == "cancer_readout" &
        prepared$endpoint_eligibility$qualification_status == "not_estimable"
    ),
    endpoint_summary = endpoint_summary,
    direct_feature_count = ncol(prepared$reference_direct),
    tsp_feature_count = length(prepared$feature_manifest$tsp_features),
    d1_feature_count = ncol(prepared$reference_d1),
    module_count = length(prepared$selected_blocks),
    direct_feature_cache = prepared$direct_cache,
    node_cache = list(
      retrieval = retrieval_checkpoint[
        c("status", "lookup_status", "reason", "key", "path")
      ],
      readout = readout_checkpoint[
        c("status", "lookup_status", "reason", "key", "path")
      ],
      learning_curve = learning_checkpoint[
        c("status", "lookup_status", "reason", "key", "path")
      ],
      decoder = decoder_checkpoint[
        c("status", "lookup_status", "reason", "key", "path")
      ]
    ),
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
    endpoint_eligibility = prepared$endpoint_eligibility,
    anchor = anchor,
    anchor_role = config$anchors$primary_role,
    evidence_level = evidence$level,
    evidence_reasons = evidence$reasons,
    direct_distance = transformed$direct$distance,
    d1_distance = transformed$d1$distance,
    cache_key = prepared$cache_key,
    input_key = prepared$input_key,
    config = config,
    config_hash = digest::digest(config, algo = "md5")
  )

  audit <- retrieval$summary
  audit$evidence_level <- evidence$level
  audit$anchor <- anchor
  audit$reference_sample_count <- nrow(prepared$reference_metadata)
  audit$query_sample_count <- nrow(readout_view_for_manifest$metadata)
  audit$external_candidate_sample_count <- nrow(prepared$query_metadata)
  audit$external_candidate_cohort_count <- length(unique(prepared$query_metadata$cohort))
  audit$excluded_query_d1_count <- manifest$excluded_query_d1_count
  audit$excluded_query_d1_sample_hash <- manifest$excluded_query_d1_sample_hash
  audit$estimable_cohort_count <- manifest$estimable_cohort_count
  audit$not_estimable_cohort_count <- manifest$not_estimable_cohort_count
  audit$config_hash <- manifest$config_hash

  # Phase 7: persist reviewer-facing products separately, then return the
  # stable CCSAblation object expected by callers of the public API.
  .ablation_atomic_save_rds(manifest, file.path(output.dir, "manifest.rds"))
  .ablation_atomic_save_rds(native_geometry, file.path(output.dir, "native_geometry.rds"))
  .ablation_atomic_save_rds(retrieval, file.path(output.dir, "retrieval.rds"))
  .ablation_atomic_save_rds(anchor_retrieval, file.path(output.dir, "anchor_retrieval.rds"))
  .ablation_atomic_save_rds(
    list(reference = prepared$reference_metadata, query = prepared$query_metadata),
    file.path(output.dir, "sample-contract.rds")
  )
  .ablation_atomic_save_rds(readout, file.path(output.dir, "readout.rds"))
  .ablation_atomic_save_rds(learning_curve, file.path(output.dir, "learning_curve.rds"))
  .ablation_atomic_save_rds(cohort_scaling, file.path(output.dir, "cohort_scaling.rds"))
  .ablation_atomic_save_rds(tradeoffs, file.path(output.dir, "tradeoffs.rds"))
  .ablation_atomic_save_rds(
    prepared$endpoint_eligibility,
    file.path(output.dir, "endpoint_eligibility.rds")
  )
  .ablation_atomic_write_csv(
    prepared$endpoint_eligibility,
    file.path(output.dir, "endpoint_eligibility.csv")
  )
  .ablation_atomic_write_csv(
    data.frame(
      sample_id = prepared$excluded_query_d1_ids,
      reason = rep("missing_precomputed_d1", length(prepared$excluded_query_d1_ids)),
      stringsAsFactors = FALSE
    ),
    file.path(output.dir, "excluded-query-d1.csv")
  )
  .ablation_atomic_write_csv(audit, file.path(output.dir, "audit.csv"))

  structure(
    list(
      call = match.call(),
      experiment = "representation",
      evidence_level = evidence$level,
      evidence = evidence,
      manifest = manifest,
      endpoint_eligibility = prepared$endpoint_eligibility,
      native_geometry = native_geometry,
      retrieval = retrieval,
      anchor_retrieval = anchor_retrieval,
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
    # Scaling/module balancing depends on the fold, not on lambda. Prepare
    # each fold once and reuse those matrices for all lambda candidates.
    fold_inputs <- lapply(sort(unique(fold)), function(fold_id) {
      inner_train <- fold != fold_id
      inner_test <- fold == fold_id
      fold_classes <- sort(unique(train_label[inner_train]))
      eligible_test <- inner_test & train_label %in% fold_classes
      if (length(fold_classes) < 2 || sum(eligible_test) < 2) {
        return(list(fold_id = fold_id, estimable = FALSE))
      }
      list(
        fold_id = fold_id,
        estimable = TRUE,
        train_label = train_label[inner_train],
        test_label = train_label[eligible_test],
        classes = fold_classes,
        transformed = .ablation_readout_transform(
          train[inner_train, , drop = FALSE],
          train[eligible_test, , drop = FALSE],
          blocks
        )
      )
    })
    cv_rows <- lapply(seq_along(lambda), function(lambda_index) {
      fold_score <- vapply(fold_inputs, function(input) {
        if (!isTRUE(input$estimable)) return(NA_real_)
        prediction <- .ablation_xgb_linear_predict(
          input$transformed$train,
          input$transformed$test,
          input$train_label,
          input$classes,
          lambda[lambda_index],
          nrounds,
          numCores,
          seed + lambda_index * 100L + input$fold_id
        )
        metrics <- .ablation_classification_metrics(
          input$test_label,
          prediction$prediction,
          prediction$probability,
          input$classes
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
.ablation_learning_curve_job <- function(
    job,
    representations,
    train_metadata,
    test_metadata,
    label_column,
    fractions,
    lambda,
    inner_folds,
    nrounds,
    numCores,
    seed,
    test_hash
) {
  fraction <- fractions[[job$fraction_index]]
  subset_seed <- seed + job$fraction_index * 1000L + job$repeat_id
  cohorts <- .ablation_sample_training_cohorts(
    train_metadata, label_column, fraction, subset_seed
  )
  train_rows <- train_metadata$cohort %in% cohorts
  input <- representations[[job$representation]]
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
  data.frame(
    representation = job$representation,
    requested_fraction = fraction,
    realized_cohort_fraction = length(cohorts) /
      length(unique(train_metadata$cohort)),
    repeat_id = job$repeat_id,
    train_cohort_count = length(cohorts),
    train_sample_count = sum(train_rows),
    selected_lambda = fit$selected_lambda,
    accuracy = fit$overall$accuracy,
    balanced_accuracy = fit$overall$balanced_accuracy,
    macro_auroc = fit$overall$macro_auroc,
    cohort_subset_hash = digest::digest(cohorts, algo = "md5"),
    test_sample_hash = test_hash,
    stringsAsFactors = FALSE
  )
}

# PSOCK workers receive the large immutable matrices once through clusterExport.
# This wrapper avoids serializing a closure containing them for every job.
.ablation_learning_curve_worker <- function(job) {
  context <- base::get(
    ".ablation_learning_curve_worker_context",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  do.call(
    .ablation_learning_curve_execute_job,
    c(list(job = job), context)
  )
}

.ablation_learning_curve_execute_job <- function(
    job,
    representations,
    train_metadata,
    test_metadata,
    label_column,
    fractions,
    lambda,
    inner_folds,
    nrounds,
    numCores,
    seed,
    test_hash,
    checkpoint_output_dir = NULL,
    checkpoint_key = NULL
) {
  compute <- function() .ablation_learning_curve_job(
    job, representations, train_metadata, test_metadata, label_column,
    fractions, lambda, inner_folds, nrounds, numCores, seed, test_hash
  )
  if (is.null(checkpoint_output_dir) || is.null(checkpoint_key)) return(compute())
  job_key <- digest::digest(
    list(
      parent = checkpoint_key,
      fraction_index = job$fraction_index,
      repeat_id = job$repeat_id,
      representation = job$representation
    ),
    algo = "md5"
  )
  .ablation_cached_node(
    node = "learning-curve-job",
    output.dir = checkpoint_output_dir,
    key = job_key,
    compute = compute,
    verbose = FALSE,
    job = as.list(job),
    parameter_digest = digest::digest(
      list(
        parent = checkpoint_key,
        label_column = label_column,
        fractions = fractions,
        lambda = lambda,
        inner_folds = inner_folds,
        nrounds = nrounds,
        numCores = numCores,
        seed = seed,
        test_hash = test_hash
      ),
      algo = "md5"
    )
  )$value
}

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
    workers = 1L,
    checkpoint_output_dir = NULL,
    checkpoint_key = NULL,
    seed
) {
  started_at <- Sys.time()
  started_memory <- .ablation_process_memory()
  test_hash <- digest::digest(sort(test_metadata$sample_id), algo = "md5")
  workers <- max(1L, as.integer(workers))
  job_table <- do.call(rbind, lapply(seq_along(fractions), function(fraction_index) {
    do.call(rbind, lapply(seq_len(as.integer(repeats)), function(repeat_id) {
      data.frame(
        fraction_index = fraction_index,
        repeat_id = repeat_id,
        representation = names(representations),
        stringsAsFactors = FALSE
      )
    }))
  }))
  jobs <- lapply(seq_len(nrow(job_table)), function(i) job_table[i, , drop = FALSE])
  run_job <- function(job) .ablation_learning_curve_execute_job(
    job, representations, train_metadata, test_metadata, label_column,
    fractions, lambda, inner_folds, nrounds, numCores, seed, test_hash,
    checkpoint_output_dir, checkpoint_key
  )
  if (workers > 1L && length(jobs) > 1L) {
    workers <- min(workers, length(jobs))
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    worker_context <- list(
      representations = representations,
      train_metadata = train_metadata,
      test_metadata = test_metadata,
      label_column = label_column,
      fractions = fractions,
      lambda = lambda,
      inner_folds = inner_folds,
      nrounds = nrounds,
      numCores = numCores,
      seed = seed,
      test_hash = test_hash,
      checkpoint_output_dir = checkpoint_output_dir,
      checkpoint_key = checkpoint_key
    )
    implementation_env <- environment(.ablation_learning_curve)
    function_names <- ls(implementation_env, pattern = "^\\.ablation_", all.names = TRUE)
    function_names <- function_names[vapply(function_names, function(name) {
      is.function(get(name, envir = implementation_env, inherits = FALSE))
    }, logical(1L))]
    parallel::clusterExport(
      cl,
      varlist = function_names,
      envir = implementation_env
    )
    parallel::clusterCall(cl, function(context) {
      assign(
        ".ablation_learning_curve_worker_context",
        context,
        envir = .GlobalEnv
      )
      NULL
    }, worker_context)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(digest))
      suppressPackageStartupMessages(library(xgboost))
      NULL
    })
    rows <- parallel::parLapplyLB(cl, jobs, .ablation_learning_curve_worker)
  } else {
    rows <- lapply(jobs, run_job)
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
  completed_at <- Sys.time()
  completed_memory <- .ablation_process_memory()
  list(
    status = "complete",
    metrics = metrics,
    paired = paired,
    test_sample_hash = test_hash,
    runtime = list(
      started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
      completed_at = format(completed_at, "%Y-%m-%dT%H:%M:%S%z"),
      elapsed_seconds = as.numeric(difftime(completed_at, started_at, units = "secs")),
      workers = workers,
      threads_per_worker = as.integer(numCores),
      total_thread_budget = workers * as.integer(numCores),
      job_count = length(jobs),
      input_bytes = sum(vapply(representations, function(value) {
        as.numeric(utils::object.size(value$train)) +
          as.numeric(utils::object.size(value$test))
      }, numeric(1L))),
      working_set_start_bytes = started_memory$working_set_bytes,
      working_set_end_bytes = completed_memory$working_set_bytes,
      peak_working_set_bytes = completed_memory$peak_working_set_bytes
    )
  )
}


# Percentile bootstrap interval for a finite numeric mean. Preserve the caller's
# RNG state because this helper is used inside deterministic scaling summaries.
.ablation_bootstrap_mean <- function(values, bootstrap, seed) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  bootstrap <- as.integer(bootstrap)
  if (length(values) == 0L || length(bootstrap) != 1L || is.na(bootstrap) ||
      bootstrap < 1L) {
    return(c(NA_real_, NA_real_))
  }
  if (length(values) == 1L) return(rep(values, 2L))
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  estimates <- replicate(
    bootstrap,
    mean(sample(values, length(values), replace = TRUE))
  )
  unname(stats::quantile(estimates, c(0.025, 0.975), names = FALSE, type = 7L))
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
  fit_cache <- list(
    schema_version = 1L,
    status = "complete",
    key = fit_cache_key,
    direct = list(),
    d1 = list()
  )
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached <- .ablation_read_fit_cache(cache_path, fit_cache_key)
    if (!is.null(cached)) fit_cache <- cached
  }
  save_fit_cache <- function() {
    if (!is.null(cache_path)) {
      .ablation_atomic_save_rds(fit_cache, cache_path)
    }
  }
  direct_fit <- fit_cache$direct$fit
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
    fit_cache$direct$fit <- direct_fit
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


# Estimate breadth/depth changes per doubling of module count, their interaction,
# and matched-size paired differences. Bootstrap summaries resample repeat-level
# estimates only.
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
      fit <- stats::lm(estimate ~ log2(module_count), data = data)
      units[[unit_index]] <- data.frame(
        contrast_type = paste0(family, "_slope"),
        aggregation = "repeat",
        repeat_id = data$repeat_id[1],
        pair_id = NA_character_,
        metric_name = data$metric_name[1],
        component = "per_module_count_doubling",
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
    schema_version = 2L,
    prepared_cache_key = prepared$cache_key,
    input_key = digest::digest(list(prepared$reference_direct, prepared$query_direct,
      prepared$reference_d1, prepared$query_d1, prepared$reference_metadata,
      prepared$query_metadata), algo = "md5"),
    bank_design_hash = bank_design$design_hash,
    direct_contracts = lapply(contracts, function(x) x$feature_hash),
    query_hash = query_hash,
    score_reference_hash = score_reference_hash,
    score_query_hash = score_query_hash,
    geometry = config$geometry,
    scaling = config$scaling,
    validation = config$validation,
    seed = seed,
    code = .ablation_node_code_identity("cohort-scaling")
  ), algo = "md5")
  fit_cache <- list(
    schema_version = 1L,
    status = "complete",
    key = cache_key,
    direct = list(),
    d1 = list()
  )
  cache_file <- NULL
  if (!is.null(cache_path)) {
    cache_dir <- sub("\\.rds$", "", cache_path)
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    cached <- .ablation_read_fit_cache(cache_file, cache_key)
    if (!is.null(cached)) fit_cache <- cached
  }
  save_fit_cache <- function() {
    if (!is.null(cache_file)) .ablation_atomic_save_rds(fit_cache, cache_file)
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


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Lifecycle stages and durable cache boundary
# -------------------------------------------------------------------------

.ablation_resolve_cache_layout <- function(cache.root = NULL, output.dir = NULL) {
  transient <- is.null(cache.root)
  if (transient) {
    cache.root <- tempfile("ccs-ablation-cache-", tmpdir = tempdir())
  }
  if (length(cache.root) != 1L || is.na(cache.root) || !nzchar(cache.root)) {
    stop("ablation: cache.root must be one non-empty path.", call. = FALSE)
  }
  root <- suppressWarnings(normalizePath(as.character(cache.root), winslash = "/", mustWork = FALSE))
  if (!is.null(output.dir) && length(output.dir) == 1L && !is.na(output.dir)) {
    output_path <- suppressWarnings(normalizePath(as.character(output.dir), winslash = "/", mustWork = FALSE))
    if (identical(root, output_path)) {
      stop("ablation: cache.root and output.dir must be different directories.", call. = FALSE)
    }
  }
  if (file.exists(root) && !dir.exists(root)) {
    stop("ablation: cache.root points to a file: ", root, call. = FALSE)
  }
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(root)) stop("ablation: cannot create cache.root: ", root, call. = FALSE)
  probe <- tempfile("write-probe-", tmpdir = root)
  if (!isTRUE(file.create(probe))) stop("ablation: cache.root is not writable: ", root, call. = FALSE)
  unlink(probe, force = TRUE)
  directories <- list(
    context_dir = file.path(root, "context"),
    plan_dir = file.path(root, "plan"),
    preparation_dir = file.path(root, "context", "preparation"),
    nodes_dir = file.path(root, "nodes"),
    jobs_dir = file.path(root, "jobs"),
    runner_dir = file.path(root, "runner"),
    state_dir = file.path(root, "state")
  )
  for (path in directories) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!all(vapply(directories, dir.exists, logical(1)))) {
    stop("ablation: cannot create cache layout below cache.root.", call. = FALSE)
  }
  c(list(schema_version = 1L, root = root, default = transient),
    directories)
}

.ablation_stage_value <- function(input, class_name, label) {
  if (!inherits(input, class_name)) {
    stop("ablation: step input must be a ", label, " created by the preceding step.", call. = FALSE)
  }
  input
}

.ablation_make_representation_context <- function(object, data, metadata, params,
                                                   seed, output.dir, cache.root, verbose) {
  if (!methods::is(object, "CCS")) stop("ablation: context step requires a CCS object.", call. = FALSE)
  cache <- .ablation_resolve_cache_layout(cache.root, output.dir)
  config <- .ablation_resolve_representation_config(seed, params)
  analysis <- .ablation_prepare_representation_analysis(
    object = object, data = data, metadata = metadata, config = config,
    output.dir = cache$preparation_dir, seed = seed, verbose = verbose
  )
  config <- .ablation_apply_runtime_config(config, analysis)
  context_key <- digest::digest(
    list(schema_version = 1L, input_key = analysis$prepared$input_key,
         config = config, seed = as.integer(seed)), algo = "md5"
  )
  context <- structure(
    list(schema_version = 1L, stage = "context", context_key = context_key,
         seed = as.integer(seed), config = config, analysis = analysis,
         cache = cache, runtime = list(R = R.version.string,
         created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
    class = c("CCSAblationContext", "CCSAblationStage")
  )
  .ablation_atomic_save_rds(
    list(schema_version = 1L, status = "complete", key = context_key,
         value_hash = digest::digest(context, algo = "md5"), value = context),
    file.path(cache$context_dir, paste0(context_key, ".rds"))
  )
  context
}

.ablation_make_representation_plan <- function(context) {
  context <- .ablation_stage_value(context, "CCSAblationContext", "context")
  node_names <- c("native_geometry", "retrieval", "readout", "learning_curve", "scaling", "decoder")
  node_seed_offsets <- c(10000L, 15000L, 20000L, 30000L, 35000L, 40000L)
  nodes <- data.frame(
    node_id = node_names, parent_key = context$context_key,
    seed = as.integer(context$seed) + node_seed_offsets,
    enabled = c(TRUE, TRUE, context$config$validation$enabled,
      context$config$validation$enabled, context$config$scaling$enabled,
      context$config$tradeoffs$decoder), stringsAsFactors = FALSE
  )
  jobs <- list(
    learning_curve = .ablation_make_learning_curve_jobs(context$config),
    scaling = .ablation_make_scaling_jobs(context$config)
  )
  plan_key <- digest::digest(list(schema_version = 1L, context_key = context$context_key,
    nodes = nodes, jobs = jobs), algo = "md5")
  plan <- structure(list(schema_version = 1L, stage = "plan", plan_key = plan_key,
    context_key = context$context_key, context = context, nodes = nodes,
    jobs = jobs, cache = context$cache), class = c("CCSAblationPlan", "CCSAblationStage"))
  .ablation_atomic_save_rds(
    list(schema_version = 1L, status = "complete", key = plan_key,
      value_hash = digest::digest(plan, algo = "md5"), value = plan),
    file.path(context$cache$plan_dir, paste0(plan_key, ".rds"))
  )
  plan
}

.ablation_run_representation_plan <- function(plan, verbose = FALSE) {
  plan <- .ablation_stage_value(plan, "CCSAblationPlan", "plan")
  context <- plan$context
  value <- .ablation_run_prepared_representation(
    analysis = context$analysis, config = context$config,
    output.dir = context$cache$runner_dir, cache.dir = context$cache$nodes_dir,
    seed = context$seed, verbose = verbose
  )
  run_key <- digest::digest(list(schema_version = 1L, plan_key = plan$plan_key,
    value_hash = digest::digest(value, algo = "md5")), algo = "md5")
  run <- structure(list(schema_version = 1L, stage = "run", run_key = run_key,
    plan_key = plan$plan_key, plan = plan, value = value, status = "complete",
    cache = context$cache), class = c("CCSAblationRun", "CCSAblationStage"))
  .ablation_atomic_save_rds(
    list(schema_version = 1L, status = "complete", key = run_key,
      value_hash = digest::digest(run, algo = "md5"), value = run),
    file.path(context$cache$state_dir, paste0(run_key, ".rds"))
  )
  run
}

.ablation_finalize_representation_stage <- function(run, output.dir, params) {
  run <- .ablation_stage_value(run, "CCSAblationRun", "run result")
  result <- run$value
  if (!inherits(result, "CCSAblation")) stop("ablation: run stage did not return a CCSAblation result.", call. = FALSE)
  if (is.null(output.dir)) return(result)
  cover <- isTRUE(params$output$cover) || isTRUE(params$general$cover)
  if (dir.exists(output.dir) && length(list.files(output.dir)) > 0L && !cover) {
    existing <- file.path(output.dir, "ablation-result.rds")
    reusable <- if (file.exists(existing)) {
      old <- tryCatch(readRDS(existing), error = function(e) NULL)
      is.list(old) && identical(old$manifest$cache$run_key, run$run_key)
    } else FALSE
    if (!reusable) {
      stop("ablation: output.dir is not empty; set the configured cover flag.", call. = FALSE)
    }
  }
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  runner_dir <- run$cache$runner_dir
  if (dir.exists(runner_dir)) {
    runner_files <- list.files(runner_dir, recursive = TRUE, full.names = TRUE)
    runner_files <- runner_files[!dir.exists(runner_files)]
    if (length(runner_files) > 0L) {
      for (source in runner_files) {
        relative <- substring(source, nchar(runner_dir) + 2L)
        destination <- file.path(output.dir, relative)
        dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
        file.copy(source, destination, overwrite = TRUE)
      }
    }
  }
  result$output.dir <- normalizePath(output.dir, winslash = "/", mustWork = TRUE)
  if (is.null(result$manifest$cache)) result$manifest$cache <- list()
  result$manifest$cache$run_key <- run$run_key
  .ablation_atomic_save_rds(result, file.path(output.dir, "ablation-result.rds"))
  if (is.data.frame(result$audit)) .ablation_atomic_write_csv(result$audit, file.path(output.dir, "audit.csv"))
  result
}

.ablation_dispatch_step <- function(step, object, data, metadata,
                                    output.dir, params, seed, verbose, input, cache.root) {
  if (identical(step, "context")) {
    if (!is.null(input)) {
      if (!is.list(input) || is.null(input$object) || is.null(input$data)) {
        stop("ablation: context input must contain object and data.", call. = FALSE)
      }
      object <- input$object; data <- input$data
      if (is.null(metadata)) metadata <- input$metadata
    }
    return(.ablation_make_representation_context(object, data, metadata, params,
      seed, output.dir, cache.root, verbose))
  }
  if (identical(step, "plan")) {
    return(.ablation_make_representation_plan(.ablation_stage_value(input,
      "CCSAblationContext", "context")))
  }
  if (identical(step, "run")) {
    if (inherits(input, "CCSAblationContext")) input <- .ablation_make_representation_plan(input)
    return(.ablation_run_representation_plan(input, verbose = verbose))
  }
  if (identical(step, "result")) return(.ablation_finalize_representation_stage(input, output.dir, params))
  stop("ablation: unsupported lifecycle step: ", step, call. = FALSE)
}

# Internal targets planning helpers
# -------------------------------------------------------------------------
#
# The functions below only build deterministic job descriptions for the
# private `plan` stage.  They do not form a public calculation API, create
# workflow markers, or decide where a targets store lives.

.ablation_make_learning_curve_jobs <- function(config) {
  if (!is.list(config) || !is.list(config$validation)) {
    stop("ablation: job planning requires a resolved configuration.", call. = FALSE)
  }
  config <- config$validation
  representations <- c("Direct-GSClassifier", "Cohort-d1")
  jobs <- expand.grid(
    fraction_index = seq_along(config$learning_fractions),
    repeat_id = seq_len(as.integer(config$repeats)),
    representation = representations,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  jobs <- jobs[order(jobs$fraction_index, jobs$repeat_id, jobs$representation), , drop = FALSE]
  rownames(jobs) <- NULL
  jobs$job_key <- vapply(seq_len(nrow(jobs)), function(i) {
    digest::digest(jobs[i, , drop = FALSE], algo = "md5")
  }, character(1))
  jobs
}


.ablation_make_scaling_jobs <- function(config) {
  if (!is.list(config) || !is.list(config$scaling)) {
    stop("ablation: job planning requires a resolved configuration.", call. = FALSE)
  }
  settings <- config$scaling
  if (!isTRUE(settings$enabled)) {
    return(data.frame(module_count = integer(), repeat_id = integer(), job_key = character()))
  }
  jobs <- expand.grid(
    module_count = as.integer(settings$module_counts),
    repeat_id = seq_len(as.integer(settings$sequences)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  jobs <- jobs[order(jobs$module_count, jobs$repeat_id), , drop = FALSE]
  rownames(jobs) <- NULL
  jobs$job_key <- vapply(seq_len(nrow(jobs)), function(i) {
    digest::digest(jobs[i, , drop = FALSE], algo = "md5")
  }, character(1))
  jobs
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
