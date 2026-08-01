#' Run layered ablation experiments for a CCS model
#'
#' @description
#' Evaluate the frozen CCS representation without retraining cohort submodels.
#' The function implements the cohort-representation, cohort-axis scaling and
#' tissue-first experiments defined by the CCS layered ablation plan. All
#' stochastic operations are paired across comparison groups and recorded in
#' an audit table.
#'
#' @param object A `CCS` object containing the frozen d1 representation.
#' @param data Raw RNA expression data. A tissue/cohort nested list is preferred;
#'   each leaf can be an expression matrix or a list containing `expr`.
#' @param metadata Optional sample annotation with sample, cohort, tissue and
#'   biological-label columns. Common CCS column names are recognized.
#' @param experiment One or more of `"cohort"`, `"scaling"` and
#'   `"tissue_first"`.
#' @param output.dir Independent output directory. Existing CCS products are
#'   never overwritten.
#' @param params Named list overriding `.ablation_default_params()`.
#' @param seed Master random seed.
#' @param verbose Whether to report progress.
#'
#' @return An object of class `CCSAblation`.
#' @author Weibin Huang <hwb2012@@qq.com>
#' @export
ablation <- function(
    object,
    data,
    metadata = NULL,
    experiment = c("cohort", "scaling", "tissue_first"),
    output.dir = file.path(getwd(), "ccs-ablation"),
    params = list(),
    seed = 20260727,
    verbose = TRUE
) {
  # Test
  if (FALSE) {
    # Purpose: Build a self-contained synthetic CCS fixture and run a fast smoke test.
    # Input: Generated expression, frozen d1 probabilities, models, and metadata.
    # Parameters: Small repeat counts keep the example quick and reproducible.
    # Output: A CCSAblation result under the current R session's temporary directory.
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

    module_ids <- c("T1|M1", "T1|M2", "T1|M3", "T2|M4", "T2|M5", "T2|M6")
    d1 <- do.call(cbind, lapply(seq_along(module_ids), function(i) {
      feature_index <- ((i - 1) %% ncol(tsp)) + 1
      score <- stats::plogis(
        tsp[, feature_index] + stats::rnorm(nrow(tsp), 0, 0.1)
      )
      block <- cbind(`1` = score, `2` = 1 - score)
      colnames(block) <- paste(module_ids[i], colnames(block), sep = "|")
      block
    }))
    rownames(d1) <- sample_ids

    make_test_model <- function() {
      class_model <- list(
        bst = NULL,
        breakVec = c(0, 0.5, 1),
        genes = c(rownames(expr), tsp_features)
      )
      list(
        Repeat = list(),
        Model = list(list(`1` = class_model, `2` = class_model))
      )
    }
    models <- list(
      T1 = stats::setNames(rep(list(make_test_model()), 3), paste0("M", 1:3)),
      T2 = stats::setNames(rep(list(make_test_model()), 3), paste0("M", 4:6))
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
    experiment <- "cohort"
    output.dir <- file.path(tempdir(), "ccs-ablation-example")
    params <- list(
      rank = 3L,
      rank_sensitivity = 3L,
      k = 3L,
      n_folds = 2L,
      bootstrap = 5L,
      geometry_samples = length(sample_ids),
      distance_pairs = 100L,
      mechanism_samples = 24L,
      rp_seeds = c(701L, 702L),
      permutation_seeds = c(801L, 802L),
      probe = FALSE,
      cover = TRUE
    )
    seed <- 909L
    verbose <- FALSE
  }

  # Step 1: Validate the public inputs and merge user overrides into reproducible defaults.
  if (!methods::is(object, "CCS")) {
    stop("ablation: object must be a CCS object.", call. = FALSE)
  }

  choices <- c("cohort", "scaling", "tissue_first")
  experiment <- unique(match.arg(experiment, choices, several.ok = TRUE))
  config <- .ablation_merge_lists(.ablation_default_params(seed), params)
  .ablation_validate_config(config)

  # Keep each run isolated and never overwrite existing results unless cover is explicit.
  if (dir.exists(output.dir) && length(list.files(output.dir)) > 0 && !config$cover) {
    stop(
      "ablation: output.dir is not empty. Use a new directory or set params$cover = TRUE.",
      call. = FALSE
    )
  }
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # Step 2: Align frozen d1, RNA/TSP features, and sample annotations, then freeze the manifest.
  if (verbose) {
    luckyBase::LuckyVerbose("ablation: Prepare frozen CCS inputs...")
  }
  prepared <- .ablation_prepare_input(
    object = object,
    data = data,
    metadata = metadata,
    max_samples = config$max_samples,
    seed = seed
  )
  manifest <- .ablation_build_manifest(object, prepared, config, seed)
  saveRDS(manifest, file.path(output.dir, "manifest.rds"))
  saveRDS(config, file.path(output.dir, "config.rds"))

  experiments <- list()
  audit_parts <- list()

  # Step 3: Run the requested experiments. Scaling depends on Gate 1, so a
  # scaling-only request still runs the cohort experiment first.
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
    # Gate 1 prevents interpreting scaling before the cohort representation beats its baselines.
    gate <- .ablation_gate_one(experiments$cohort, config$gate1)
    if (gate$pass || !config$gate1$enforce) {
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
    # Experiment 1: common rank, neighborhood, validation, and geometry settings.
    rank = 50L,
    rank_sensitivity = c(25L, 50L, 100L),
    k = 30L,
    distance = "euclidean",
    n_folds = Inf,
    bootstrap = 1000L,
    max_samples = Inf,
    geometry_samples = 5000L,
    distance_pairs = 100000L,
    mechanism_samples = 1000L,
    probe = TRUE,
    probe_label = "tissue",
    probe_nrounds = 50L,
    numCores = 1L,
    # Repeat random projection and within-cohort module permutation separately.
    rp_density = 1 / 3,
    rp_seeds = seed + seq_len(20),
    permutation_seeds = seed + 1000L + seq_len(20),
    # Experiment 2: nested module-count sequences and downstream embedding repeats.
    scaling_counts = c(10L, 25L, 50L, 75L, 100L, 125L, 150L),
    scaling_sequences = 100L,
    scaling_embedding_counts = c(25L, 50L, 100L, 150L),
    scaling_embedding_sequences = 10L,
    scaling_embedding_seeds = seed + 2000L + seq_len(10),
    scaling_subsample_fraction = 0.8,
    # Experiment 3: paired sampling and reduction settings shared by both arms.
    tissue_seeds = seed + 3000L + seq_len(20),
    tissue_subsample_fraction = 0.8,
    fidelity_samples = 2000L,
    dr = list(
      method = "UWOT",
      dimension = c(5L, 2L),
      n_neighbors = 30L,
      min_dist = 0.01,
      spread = 0.75,
      set_op_mix_ratio = 1
    ),
    cluster = list(eps = 0.02, minPts = 20L),
    # Gate 1 checks the primary effect, biological purity, and cohort mixing together.
    gate1 = list(
      enforce = TRUE,
      primary_metric = "balanced_accuracy",
      min_gain = 0,
      purity_tolerance = 0,
      mixing_tolerance = 0
    ),
    cover = FALSE
  )
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
  if (length(config$rank) != 1 || config$rank < 1) {
    stop("ablation: params$rank must be a positive integer.", call. = FALSE)
  }
  if (length(config$k) != 1 || config$k < 1) {
    stop("ablation: params$k must be a positive integer.", call. = FALSE)
  }
  if (length(config$rp_seeds) < 1 || length(config$permutation_seeds) < 1) {
    stop("ablation: both null controls require at least one seed.", call. = FALSE)
  }
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


# Read only the TSP features used by frozen models; never retrain or modify CCS.
# Prefer object@Model and fall back to the recorded model directory when needed.
.ablation_extract_tsp_features <- function(object, module_manifest) {
  models <- object@Model
  use_embedded <- length(models) > 0 && !identical(models, list(NA))
  path_map <- if (use_embedded) NULL else .ablation_model_path_map(object)

  features <- unlist(lapply(seq_len(nrow(module_manifest$modules)), function(i) {
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
    .ablation_model_features(model)
  }), use.names = FALSE)
  features <- unique(features[grepl(":", features, fixed = TRUE)])
  if (length(features) == 0) {
    stop("ablation: no TSP pair features were found in frozen models.", call. = FALSE)
  }
  features
}


# Collect gene pairs used across all repeats and class-specific frozen models.
.ablation_model_features <- function(model) {
  if (is.null(model) || is.null(model$Model)) {
    stop("ablation: malformed frozen cohort model.", call. = FALSE)
  }
  unique(unlist(lapply(model$Model, function(repeat_model) {
    unlist(lapply(repeat_model, function(class_model) {
      if (is.null(class_model)) character() else class_model$genes
    }), use.names = FALSE)
  }), use.names = FALSE))
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


# Align CCS, RNA expression, and metadata to one sample set and build the shared TSP matrix.
.ablation_prepare_input <- function(
    object,
    data,
    metadata = NULL,
    max_samples = Inf,
    seed = 20260727
) {
  # Derive the feature universe from frozen models before organizing expression and metadata.
  module_manifest <- .ablation_module_manifest(object)
  tsp_features <- .ablation_extract_tsp_features(object, module_manifest)
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
  pair_genes <- unique(unlist(strsplit(tsp_features, ":", fixed = TRUE)))
  # Reuse the CCS geneMatch rule only when raw row names do not cover model genes.
  if (!all(pair_genes %in% rownames(expr))) {
    matched <- GSClassifier::geneMatch(
      X = expr,
      geneAnnotation = object@Repeat$geneAnnotation,
      geneid = object@Repeat$geneid,
      matchmode = "fix"
    )
    expr <- matched$Subset
  }
  missing_genes <- setdiff(pair_genes, rownames(expr))
  if (length(missing_genes) > 0) {
    stop(
      "ablation: RNA data lack TSP genes: ",
      paste(missing_genes, collapse = ", "),
      call. = FALSE
    )
  }

  tsp <- .ablation_tsp_matrix(expr, tsp_features)
  d1 <- d1[sample_ids, , drop = FALSE]
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
  rownames(metadata) <- metadata$sample_id

  list(
    tsp = tsp,
    d1 = d1,
    metadata = metadata,
    module_manifest = module_manifest,
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
  # Intersect genes across cohorts before binding so every sample shares one feature space.
  genes <- Reduce(intersect, lapply(matrices, rownames))
  expr <- do.call(cbind, lapply(matrices, function(x) x[genes, , drop = FALSE]))
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


# Encode each geneA:geneB pair as a binary TSP, with geneA >= geneB mapped to 1.
.ablation_tsp_matrix <- function(expr, features) {
  pairs <- strsplit(features, ":", fixed = TRUE)
  tsp <- vapply(pairs, function(pair) {
    as.integer(expr[pair[1], ] >= expr[pair[2], ])
  }, integer(ncol(expr)))
  rownames(tsp) <- colnames(expr)
  colnames(tsp) <- features
  tsp
}


# Record input sizes, versions, settings, and sample/feature hashes for audit and comparison.
.ablation_build_manifest <- function(object, prepared, config, seed) {
  sample_hash <- digest::digest(sort(prepared$metadata$sample_id), algo = "md5")
  feature_hash <- digest::digest(
    list(
      tsp = prepared$tsp_features,
      modules = prepared$module_manifest$modules,
      d1_columns = colnames(prepared$d1)
    ),
    algo = "md5"
  )
  list(
    version = 1L,
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = seed,
    object_class = class(object)[1],
    object_version = tryCatch(as.character(utils::packageVersion("CCS")), error = function(e) NA_character_),
    model_dir = object@Repeat$model.dir,
    sample_count = nrow(prepared$d1),
    tsp_feature_count = ncol(prepared$tsp),
    d1_dimension = dim(prepared$d1),
    module_count = nrow(prepared$module_manifest$modules),
    excluded_duplicate_sample_count = length(prepared$excluded_duplicate_samples),
    sample_manifest_hash = sample_hash,
    feature_manifest_hash = feature_hash,
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
    nrounds = as.integer(config$probe_nrounds),
    eta = 0.1,
    lambda = 1,
    alpha = 0,
    nthread = as.integer(config$numCores),
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


# Compare Direct and candidate distance-rank changes along original TSP neighbor edges,
# separating biology-discordant edges from biology-matched, cross-cohort edges.
.ablation_selective_reconstruction <- function(
    tsp,
    direct,
    candidate,
    metadata,
    k,
    max_samples,
    seed
) {
  if (nrow(tsp) > max_samples) {
    rows <- .ablation_stratified_sample(metadata, max_samples, seed)
    tsp <- tsp[rows, , drop = FALSE]
    direct <- direct[rows, , drop = FALSE]
    candidate <- candidate[rows, , drop = FALSE]
    metadata <- metadata[rows, , drop = FALSE]
  }
  k <- min(k, nrow(tsp) - 1L)
  edges <- .ablation_knn(tsp, k)
  direct_distance <- as.matrix(stats::dist(direct))
  candidate_distance <- as.matrix(stats::dist(candidate))
  direct_rank <- t(apply(direct_distance, 1, rank, ties.method = "average")) - 1
  candidate_rank <- t(apply(candidate_distance, 1, rank, ties.method = "average")) - 1
  # delta > 0 means the TSP neighbor ranks farther away in the candidate than in Direct.
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
    tsp_test,
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
  local_k <- min(config$k, nrow(candidate_scores) - 1L)
  # Compute geometry, mixing, probe, and mechanism metrics on the same test sample set.
  mixing <- .ablation_mixing_purity(candidate_scores, metadata_test, local_k)
  mechanism <- .ablation_selective_reconstruction(
    tsp = tsp_test,
    direct = direct_scores,
    candidate = candidate_scores,
    metadata = metadata_test,
    k = local_k,
    max_samples = config$mechanism_samples,
    seed = seed
  )
  probe <- if (config$probe) {
    probe_metadata <- rbind(metadata_train, metadata_test)
    probe_label <- .ablation_eligible_probe_labels(
      probe_metadata,
      config$probe_label
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
      config$distance_pairs,
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
    distance = config$distance,
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
    n_folds = config$n_folds,
    seed = seed
  )
  # Generate each permutation once and reuse it across ranks and folds for strict pairing.
  permuted <- lapply(config$permutation_seeds, function(seed_i) {
    .ablation_permute_blocks(
      prepared$d1,
      prepared$module_manifest$blocks,
      prepared$metadata$cohort,
      seed_i
    )
  })
  names(permuted) <- as.character(config$permutation_seeds)

  minimum_train_size <- min(vapply(
    sort(unique(folds)),
    function(fold) sum(folds != fold),
    integer(1)
  ))
  maximum_rank <- min(
    ncol(prepared$tsp),
    ncol(prepared$d1),
    minimum_train_size - 1L
  )
  # Truncate all requested ranks to the largest common rank estimable in every fold.
  rank_targets <- unique(pmin(
    as.integer(c(config$rank, config$rank_sensitivity)),
    maximum_rank
  ))
  main_rank <- min(as.integer(config$rank), maximum_rank)
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
  summary <- .ablation_metric_summary(metrics, config$bootstrap, seed)
  contrasts <- .ablation_paired_contrasts(metrics, config$bootstrap, seed)
  sensitivity <- list(
    metrics = sensitivity_metrics,
    summary = if (nrow(sensitivity_metrics) > 0) {
      .ablation_metric_summary(sensitivity_metrics, config$bootstrap, seed + 500L)
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
      ncol(prepared$tsp),
      ncol(prepared$d1),
      sum(train_rows) - 1L
    )
    direct <- .ablation_fit_pca(
      prepared$tsp[train_rows, , drop = FALSE],
      prepared$tsp[test_rows, , drop = FALSE],
      rank_q
    )
    cohort <- .ablation_fit_pca(
      prepared$d1[train_rows, , drop = FALSE],
      prepared$d1[test_rows, , drop = FALSE],
      rank_q
    )
    metadata_train <- prepared$metadata[train_rows, , drop = FALSE]
    metadata_test <- prepared$metadata[test_rows, , drop = FALSE]
    tsp_test <- prepared$tsp[test_rows, , drop = FALSE]

    # Direct is equal-rank PCA of raw TSPs; Cohort is equal-rank PCA of frozen d1.
    metrics[[index]] <- .ablation_metric_rows(
      direct$test, direct$test, tsp_test, metadata_test, metadata_train,
      direct$train, "Direct", fold, NA_character_, NA_integer_, direct$rank,
      config, seed + fold
    )
    index <- index + 1L
    metrics[[index]] <- .ablation_metric_rows(
      direct$test, cohort$test, tsp_test, metadata_test, metadata_train,
      cohort$train, "Cohort", fold, NA_character_, NA_integer_, cohort$rank,
      config, seed + 100L + fold
    )
    index <- index + 1L

    scaled_tsp <- .ablation_scale_train_apply(
      prepared$tsp[train_rows, , drop = FALSE],
      prepared$tsp[test_rows, , drop = FALSE]
    )
    # Null-RP tests whether equal-rank random compression alone improves the metrics.
    for (seed_i in config$rp_seeds) {
      projection <- .ablation_projection_matrix(
        ncol(prepared$tsp),
        rank_q,
        seed_i,
        config$rp_density
      )
      rp_train <- scaled_tsp$train %*% projection
      rp_test <- scaled_tsp$test %*% projection
      metrics[[index]] <- .ablation_metric_rows(
        direct$test, rp_test, tsp_test, metadata_test, metadata_train,
        rp_train, "Null-RP", fold, "projection_seed", seed_i, rank_q,
        config, seed_i + fold
      )
      index <- index + 1L
    }

    # Null-Perm preserves within-cohort block structure but breaks cross-block sample pairing.
    for (seed_i in config$permutation_seeds) {
      perm_fit <- .ablation_fit_pca(
        permuted[[as.character(seed_i)]][train_rows, , drop = FALSE],
        permuted[[as.character(seed_i)]][test_rows, , drop = FALSE],
        rank_q
      )
      metrics[[index]] <- .ablation_metric_rows(
        direct$test, perm_fit$test, tsp_test, metadata_test, metadata_train,
        perm_fit$train, "Null-Perm", fold, "permutation_seed", seed_i,
        perm_fit$rank, config, seed_i + fold
      )
      index <- index + 1L
    }
  }
  do.call(rbind, metrics)
}


# Compare raw TSP and d1 with dimension-free geometry metrics as evidence beyond PCA.
.ablation_dimension_free_geometry <- function(prepared, config, seed) {
  rows <- seq_len(nrow(prepared$d1))
  if (length(rows) > config$geometry_samples) {
    rows <- .ablation_stratified_sample(
      prepared$metadata,
      config$geometry_samples,
      seed
    )
  }
  tsp <- prepared$tsp[rows, , drop = FALSE]
  d1 <- prepared$d1[rows, , drop = FALSE]
  k <- min(config$k, length(rows) - 1L)
  data.frame(
    metric_name = c("linear_cka", "distance_spearman", "knn_jaccard"),
    metric_value = c(
      .ablation_linear_cka(tsp, d1),
      .ablation_distance_spearman(tsp, d1, config$distance_pairs, seed),
      .ablation_knn_jaccard(tsp, d1, k)
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
    config$scaling_sequences,
    seed
  )
  counts <- sort(unique(pmin(config$scaling_counts, manifest$module_count)))
  folds <- .ablation_grouped_folds(prepared$metadata$cohort, config$n_folds, seed)
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
        rank_q <- min(config$rank, ncol(candidate), sum(train_rows) - 1L)
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
        local_k <- min(config$k, sum(test_rows) - 1L)
        mixing <- .ablation_mixing_purity(subset_fit$test, metadata_test, local_k)
        probe <- if (config$probe) {
          probe_metadata <- rbind(metadata_train, metadata_test)
          probe_label <- .ablation_eligible_probe_labels(
            probe_metadata,
            config$probe_label
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
          distance = config$distance,
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
  summary <- .ablation_scaling_summary(metrics, config$bootstrap, seed)
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
    as.integer(config$scaling_embedding_sequences),
    length(sequences)
  )
  selected_sequences <- sequences[seq_len(sequence_count)]
  counts <- sort(unique(pmin(
    config$scaling_embedding_counts,
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
      for (seed_i in config$scaling_embedding_seeds) {
        rows <- .ablation_tissue_subsample(
          prepared$metadata,
          config$scaling_subsample_fraction,
          seed_i
        )
        d1 <- prepared$d1[rows, columns, drop = FALSE]
        metadata <- prepared$metadata[rows, , drop = FALSE]
        d3 <- .ablation_two_stage_embedding(d1, config$dr, seed_i)
        clusters <- .ablation_dbscan(d3, config$cluster)
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
          k = min(config$k, nrow(d3) - 1L),
          distance = config$distance,
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
  stability <- .ablation_scaling_embedding_stability(embeddings, config$k)
  summary <- .ablation_metric_summary(metrics, config$bootstrap, seed)
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
  seeds <- config$tissue_seeds
  results <- list()
  metrics <- list()
  stratified <- list()
  index <- 1L
  # Both arms share the same sample subset and seed within each repeat for strict pairing.
  for (i in seq_along(seeds)) {
    seed_i <- seeds[i]
    rows <- .ablation_tissue_subsample(
      prepared$metadata,
      config$tissue_subsample_fraction,
      seed_i
    )
    d1 <- prepared$d1[rows, , drop = FALSE]
    metadata <- prepared$metadata[rows, , drop = FALSE]
    embeddings <- .ablation_tissue_embeddings(d1, config$dr, seed_i)
    for (group in names(embeddings)) {
      d3 <- embeddings[[group]]
      clusters <- .ablation_dbscan(d3, config$cluster)
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
        k = min(config$k, nrow(d3) - 1L),
        distance = config$distance,
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
        k = config$k
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
  stability <- .ablation_embedding_stability(results, config$k)
  summary <- .ablation_metric_summary(metrics, config$bootstrap, seed)
  contrasts <- .ablation_two_group_contrasts(
    metrics,
    first_group = "Two-stage",
    second_group = "One-stage",
    bootstrap = config$bootstrap,
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
  if (length(rows) > config$fidelity_samples) {
    rows <- .ablation_stratified_sample(metadata, config$fidelity_samples, seed)
  }
  fidelity <- .ablation_trust_continuity(
    high[rows, , drop = FALSE],
    low[rows, , drop = FALSE],
    min(config$k, length(rows) - 1L)
  )
  local_k <- min(config$k, nrow(low) - 1L)
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
  index <- sum(choose2(table_xy))
  row_index <- sum(choose2(rowSums(table_xy)))
  col_index <- sum(choose2(colSums(table_xy)))
  total <- choose2(sum(table_xy))
  expected <- row_index * col_index / total
  maximum <- (row_index + col_index) / 2
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
  mean(c(one_way(x, y), one_way(y, x)), na.rm = TRUE)
}
