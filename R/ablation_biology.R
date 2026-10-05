# Continuous biology diagnostics. Pure scientific functions are called by the
# ablation-03 targets graph; targets owns execution, persistence and recovery.

.ablation_bio_seed <- function(config, id) {
  as.integer((config$seed + sum(utf8ToInt(id) * seq_along(utf8ToInt(id)))) %%
    .Machine$integer.max)
}

.ablation_bio_known <- function(x) {
  !is.na(x) & nzchar(trimws(x)) &
    !tolower(trimws(x)) %in% c("unknown", "unknow", "undefined", "na", "n/a", "")
}

.ablation_bio_bind <- function(x) {
  x <- Filter(function(y) is.data.frame(y) && nrow(y) > 0L, x)
  if (!length(x)) return(data.frame())
  rownames_result <- do.call(rbind, x)
  rownames(rownames_result) <- NULL
  rownames_result
}

.ablation_bio_legacy_scores <- function(cache) {
  genes <- sort(unique(unlist(cache$anchors, use.names = FALSE)))
  values <- lapply(genes, function(gene) {
    x <- unlist(lapply(cache$cohorts[cache$reference_cohorts], function(cohort) {
      if (!gene %in% rownames(cohort$expression)) return(numeric())
      cohort$expression[gene, cohort$sample_id %in% cache$reference_sample_ids]
    }), use.names = FALSE)
    x <- x[is.finite(x)]
    c(mean = mean(x), sd = stats::sd(x))
  })
  scale <- do.call(rbind, values)
  rownames(scale) <- genes
  scale <- scale[is.finite(scale[, "sd"]) & scale[, "sd"] > 0, , drop = FALSE]
  .ablation_bio_bind(lapply(cache$cohorts, function(cohort) {
    .ablation_bio_bind(lapply(names(cache$anchors), function(anchor) {
      use <- intersect(cache$anchors[[anchor]],
        intersect(rownames(cohort$expression), rownames(scale)))
      if (length(use) < 2L) return(NULL)
      z <- sweep(sweep(cohort$expression[use, , drop = FALSE], 1L,
        scale[use, "mean"], "-"), 1L, scale[use, "sd"], "/")
      data.frame(sample_id = cohort$sample_id, anchor = anchor,
        score = colMeans(z, na.rm = TRUE))
    }))
  }))
}

.ablation_bio_baseline_check <- function(cache, neighbors, old, tolerance) {
  scores <- .ablation_bio_legacy_scores(cache)
  rows <- .ablation_bio_bind(lapply(names(cache$anchors), function(anchor) {
    sc <- scores[scores$anchor == anchor, ]
    value <- stats::setNames(sc$score, sc$sample_id)
    n <- neighbors[neighbors$neighbor_rank <= 15L, ]
    n$utility <- exp(-abs(value[n$query_sample] - value[n$reference_sample]))
    n <- n[is.finite(n$utility), ]
    key <- interaction(n$representation, n$query_sample, drop = TRUE)
    complete <- names(which(table(key) == 15L))
    n <- n[as.character(key) %in% complete, ]
    per_query <- stats::aggregate(utility ~ representation + query_sample +
      query_cohort, n, mean)
    d <- per_query[per_query$representation == "Direct-GSClassifier", ]
    b <- per_query[per_query$representation == "Cohort-d1", ]
    pairs <- merge(d, b, by = c("query_sample", "query_cohort"),
      suffixes = c("_direct", "_d1"))
    means <- stats::aggregate(cbind(utility_direct, utility_d1) ~ query_cohort,
      pairs, mean)
    counts <- stats::aggregate(query_sample ~ query_cohort, pairs, length)
    means$query_count <- counts$query_sample[match(means$query_cohort,
      counts$query_cohort)]
    means$anchor <- anchor
    means$delta_d1_minus_direct <- means$utility_d1 - means$utility_direct
    means
  }))
  compared <- merge(rows, old, by = c("anchor", "query_cohort"),
    suffixes = c("_new", "_old"))
  error <- max(abs(compared$delta_d1_minus_direct_new -
    compared$delta_d1_minus_direct_old))
  if (nrow(compared) != nrow(old) || nrow(compared) != nrow(rows) ||
      any(compared$query_count_new != compared$query_count_old) ||
      !is.finite(error) || error > tolerance) {
    stop("biology diagnostics: zero-change baseline does not reproduce.", call. = FALSE)
  }
  list(status = "PASS", cells = nrow(compared), max_absolute_error = error,
    tolerance = tolerance)
}

# Includes every gene in frozen geneSet support if set pairs are present. This
# conservative superset avoids claiming disjointness from opaque set names.
.ablation_bio_input_genes <- function(feature_manifest, gene_set) {
  features <- feature_manifest$feature_manifest
  plain <- features$feature[features$feature_type == "single_bin"]
  pairs <- features$feature[features$feature_type == "gene_pair"]
  genes <- unique(c(plain, unlist(strsplit(pairs, ":", fixed = TRUE))))
  sets <- features$feature[features$feature_type == "set_pair"]
  complete <- !length(sets) || (is.list(gene_set) && length(gene_set) > 0L &&
    all(lengths(gene_set) > 0L))
  if (length(sets) && complete) genes <- unique(c(genes, unlist(gene_set)))
  list(genes = sort(as.character(genes)), complete = complete,
    reason = if (complete) "complete_conservative_frozen_geneSet_support" else
      "frozen_geneSet_support_unavailable",
    set_pair_count = length(sets))
}

.ablation_bio_eligibility <- function(prepared, config, restriction = NULL) {
  ref <- prepared$reference_metadata
  query <- prepared$query_metadata
  evidence <- function(m) .ablation_bio_known(m$cancer_type) &
    m$metadata_status == "confirmed" & .ablation_bio_known(m$evidence_basis)
  ref_ok <- evidence(ref)
  query_ok <- evidence(query) & query$d1_provenance == "external_frozen"
  query_ok[is.na(query_ok)] <- FALSE
  ref_ok[is.na(ref_ok)] <- FALSE
  pool <- vector("list", nrow(query))
  details <- lapply(seq_len(nrow(query)), function(i) {
    selected <- ref_ok & ref$cancer_type == query$cancer_type[i] &
      ref$cohort_key != query$cohort_key[i]
    if (!is.null(restriction)) {
      selected <- selected & .ablation_bio_known(ref[[restriction]]) &
        .ablation_bio_known(query[[restriction]][i]) &
        ref[[restriction]] == query[[restriction]][i]
    }
    selected[is.na(selected)] <- FALSE
    ids <- which(selected)
    count <- length(unique(ref$cohort_key[ids]))
    valid <- query_ok[i] && count >= config$min_reference_cohorts &&
      length(ids) >= config$k
    pool[[i]] <<- if (valid) ids else integer()
    data.frame(sample_id = query$sample_id[i], cohort = query$cohort_key[i],
      cancer_type = query$cancer_type[i], reference_count = length(ids),
      reference_cohort_count = count, eligible = valid,
      reason = if (!query_ok[i]) "unconfirmed_cancer_or_external_provenance" else
        if (!is.null(restriction) && !.ablation_bio_known(query[[restriction]][i]))
          "unknown_technical_label" else
        if (count < config$min_reference_cohorts) "fewer_than_two_reference_cohorts" else
        if (length(ids) < config$k) "fewer_than_k_reference_samples" else "supported")
  })
  list(queries = .ablation_bio_bind(details), pools = pool)
}

# A single atlas pass supplies availability-defined lists and all score inputs.
# No gene-list decision uses expression values or observed arm effects.
.ablation_bio_audit <- function(atlas, cache, prepared, gene_set, signatures,
    config, preprocessing = NULL) {
  ref <- prepared$reference_metadata
  query <- prepared$query_metadata
  ids <- c(ref$sample_id, query$sample_id)
  if (anyDuplicated(ids) || length(intersect(ref$sample_id, query$sample_id))) {
    stop("biology diagnostics: duplicate or overlapping sample identities.", call. = FALSE)
  }
  modules <- prepared$module_manifest$modules
  if ("cohort" %in% names(modules) &&
      length(intersect(as.character(modules$cohort), as.character(query$cohort)))) {
    stop("biology diagnostics: frozen training cohort overlaps external query.", call. = FALSE)
  }
  for (arm in c("direct", "d1")) {
    if (!identical(rownames(prepared[[paste0("reference_", arm)]]), ref$sample_id) ||
        !identical(rownames(prepared[[paste0("query_", arm)]]), query$sample_id)) {
      stop("biology diagnostics: matrix and sample contract do not align.", call. = FALSE)
    }
  }
  metadata <- rbind(ref, query)
  cohorts <- list()
  for (tissue in names(atlas)) for (name in names(atlas[[tissue]])) {
    source_key <- paste(tissue, name, sep = "/")
    source_ids <- metadata$sample_id[metadata$cohort_key == source_key]
    if (!length(source_ids)) next
    raw <- atlas[[tissue]][[name]]
    if (is.list(raw) && !is.matrix(raw) && !is.data.frame(raw)) raw <- raw$expr
    keep <- which(colnames(raw) %in% source_ids)
    if (!length(keep)) next
    key <- unique(metadata$cohort_key[match(colnames(raw)[keep], metadata$sample_id)])
    if (length(key) != 1L || key %in% names(cohorts) ||
        is.null(rownames(raw)) || anyDuplicated(rownames(raw)) ||
        anyDuplicated(colnames(raw)[keep])) {
      stop("biology diagnostics: ambiguous expression cohort or gene/sample IDs.", call. = FALSE)
    }
    cohorts[[key]] <- list(raw = raw, columns = keep)
  }
  found_ids <- unlist(lapply(cohorts, function(x) colnames(x$raw)[x$columns]))
  if (!setequal(ids, found_ids) || anyDuplicated(found_ids)) {
    stop("biology diagnostics: expression source does not cover sample contract exactly.", call. = FALSE)
  }
  background <- sort(Reduce(intersect, lapply(cohorts, function(x) rownames(x$raw))))
  anchors <- cache$anchors
  split_names <- c(ifn = "IFN\u03b3 signaling", il6 = "IL6-JAK-STAT3 signaling")
  for (name in names(split_names)) {
    genes <- signatures[["IFN-IL6"]][[split_names[[name]]]]
    if (length(genes)) anchors[[name]] <- sort(unique(as.character(genes)))
  }
  input <- .ablation_bio_input_genes(prepared$feature_manifest, gene_set)
  original <- anchors
  for (name in config$primary_anchors) {
    anchors[[paste0(name, "_disjoint")]] <- if (input$complete)
      setdiff(anchors[[name]], input$genes) else character()
  }
  common <- lapply(anchors, intersect, y = background)
  coverage <- .ablation_bio_bind(lapply(names(common), function(name) {
    base <- sub("_disjoint$", "", name)
    count <- length(common[[name]])
    supported <- count >= config$min_signature_genes &&
      count >= config$min_signature_fraction * length(original[[base]]) &&
      (!grepl("_disjoint$", name) || input$complete)
    data.frame(anchor = name, original_count = length(original[[base]]),
      common_count = count, coverage = count / length(original[[base]]),
      input_overlap_count = length(intersect(original[[base]], input$genes)),
      estimable = supported, reason = if (supported) "supported" else
        if (grepl("_disjoint$", name) && !input$complete) input$reason else
          "below_50_percent_or_eight_genes")
  }))
  rank_enabled <- "common_rank" %in% config$score_definitions
  rank_ok <- rank_enabled && length(background) >= config$min_background_genes
  signature_genes <- sort(unique(unlist(common)))
  raw_signature <- matrix(NA_real_, length(signature_genes), length(ids),
    dimnames = list(signature_genes, ids))
  rank_scores <- if (rank_enabled) matrix(NA_real_, length(ids), length(common),
    dimnames = list(ids, names(common))) else NULL
  explore <- isTRUE(config$low_coverage_exploration)
  exploratory_rank <- if (explore && rank_enabled) matrix(NA_real_, length(ids),
    length(config$primary_anchors), dimnames = list(ids, config$primary_anchors)) else NULL
  audit <- .ablation_bio_bind(lapply(names(cohorts), function(key) {
    cohort <- cohorts[[key]]
    sample_ids <- colnames(cohort$raw)[cohort$columns]
    cols <- match(sample_ids, ids)
    raw_signature[, cols] <<- as.matrix(cohort$raw[signature_genes,
      cohort$columns, drop = FALSE])
    bg <- as.matrix(cohort$raw[background, cohort$columns, drop = FALSE])
    storage.mode(bg) <- "double"
    ties <- numeric(ncol(bg))
    for (j in seq_len(ncol(bg))) {
      finite <- is.finite(bg[, j])
      ties[j] <- if (any(finite)) 1 - length(unique(bg[finite, j])) / sum(finite) else NA_real_
      if ((!rank_ok && !explore) || !all(finite) || !nrow(bg)) next
      ranks <- (rank(bg[, j], ties.method = "average") - 0.5) / nrow(bg)
      for (name in names(common)) {
        if (rank_ok && coverage$estimable[match(name, coverage$anchor)]) {
          rank_scores[cols[j], name] <<- mean(ranks[match(common[[name]], background)])
        }
        if (!is.null(exploratory_rank) && name %in% colnames(exploratory_rank) &&
            length(common[[name]])) exploratory_rank[cols[j], name] <<-
          mean(ranks[match(common[[name]], background)])
      }
    }
    m <- metadata[match(sample_ids, metadata$sample_id), ]
    provenance <- if (!is.null(preprocessing)) preprocessing[
      preprocessing$cohort_key == key, , drop = FALSE] else NULL
    data.frame(cohort = key, role = if (all(sample_ids %in% ref$sample_id))
      "reference" else "query", sample_count = length(cols),
      available_genes = nrow(cohort$raw), background_genes = length(background),
      background_coverage = length(background) / nrow(cohort$raw),
      rank_estimable = rank_ok, nonfinite_fraction = mean(!is.finite(bg)),
      complete_background_samples = sum(colSums(!is.finite(bg)) == 0L),
      zero_fraction = mean(bg == 0, na.rm = TRUE), mean_tie_fraction = mean(ties),
      expression_median = stats::median(bg[is.finite(bg)]),
      expression_q01 = if (any(is.finite(bg))) unname(stats::quantile(bg[is.finite(bg)], .01)) else NA_real_,
      expression_q99 = if (any(is.finite(bg))) unname(stats::quantile(bg[is.finite(bg)], .99)) else NA_real_,
      expression_unit = if (!is.null(provenance) && nrow(provenance)) provenance$unit[1L] else "unknown",
      transformation = if (!is.null(provenance) && nrow(provenance)) provenance$transformation[1L] else "unknown",
      preprocessing_evidence = if (!is.null(provenance) && nrow(provenance)) provenance$evidence[1L] else
        "no_cohort_specific_preprocessing_record_found",
      confirmed_cancer_fraction = mean(.ablation_bio_known(m$cancer_type) &
        m$metadata_status == "confirmed"),
      known_assay_fraction = mean(.ablation_bio_known(m$assay_type)),
      known_source_fraction = mean(.ablation_bio_known(m$source_system)),
      known_platform_fraction = mean(.ablation_bio_known(m$platform_id)))
  }))
  eligibility <- .ablation_bio_eligibility(prepared, config)
  old <- prepared$endpoint_eligibility
  old <- old[old$endpoint == "cancer_retrieval", ]
  eligibility$queries$original_endpoint_status <- old$qualification_status[
    match(eligibility$queries$cohort, old$cohort_key)]
  list(audit = audit, coverage = coverage, background = background,
    common_signatures = common, input_support = input, eligibility = eligibility,
    training_cohort_check = if ("cohort" %in% names(modules)) "PASS" else "unavailable_manifest_cohort",
    raw_signature = raw_signature, rank_scores = rank_scores,
    exploratory_rank = exploratory_rank,
    reference_ids = ref$sample_id, query_ids = query$sample_id,
    rank_status = if (!rank_enabled) "disabled_by_design" else
      if (rank_ok) "estimable" else "background_below_5000_genes",
    config = config)
}

.ablation_bio_raw_scores <- function(raw, signatures, fit_ids) {
  train <- raw[, fit_ids, drop = FALSE]
  center <- rowMeans(train)
  scale <- apply(train, 1L, stats::sd)
  z <- sweep(sweep(raw, 1L, center, "-"), 1L, scale, "/")
  result <- matrix(NA_real_, ncol(raw), length(signatures),
    dimnames = list(colnames(raw), names(signatures)))
  for (name in names(signatures)) {
    genes <- signatures[[name]]
    if (!length(genes) || any(!is.finite(scale[match(genes, rownames(raw))])) ||
        any(scale[match(genes, rownames(raw))] <= 0)) next
    result[, name] <- colMeans(z[genes, , drop = FALSE])
  }
  result
}

.ablation_bio_scores <- function(audit) {
  raw <- .ablation_bio_raw_scores(audit$raw_signature,
    audit$common_signatures, audit$reference_ids)
  raw[, !audit$coverage$estimable] <- NA_real_
  definitions <- list(common_raw = raw, common_rank = audit$rank_scores)[
    audit$config$score_definitions]
  unscaled <- definitions
  scales <- .ablation_bio_bind(lapply(names(definitions), function(name) {
    x <- definitions[[name]]
    fit <- x[audit$reference_ids, , drop = FALSE]
    sd <- apply(fit, 2L, stats::sd)
    center <- colMeans(fit)
    for (j in seq_len(ncol(x))) {
      if (!is.finite(sd[j]) || sd[j] <= 0) x[, j] <- NA_real_ else
        x[, j] <- (x[, j] - center[j]) / sd[j]
    }
    definitions[[name]] <<- x
    data.frame(score_definition = name, anchor = colnames(x),
      reference_mean = center, reference_sd = sd,
      status = ifelse(is.finite(sd) & sd > 0, "estimable", "not_estimable"),
      reason = ifelse(is.finite(sd) & sd > 0, "supported",
        "signature_coverage_nonfinite_or_reference_zero_variance"))
  }))
  list(values = definitions, scales = scales, unscaled = unscaled,
    exploratory = if (isTRUE(audit$config$low_coverage_exploration))
      .ablation_bio_exploratory_scores(audit) else NULL)
}

# Fixed residual signatures, complete-sample scores and reference-only scales.
# Strict targets above are never altered by this exploratory product.
.ablation_bio_exploratory_scores <- function(audit) {
  signatures <- audit$common_signatures[audit$config$primary_anchors]
  raw <- audit$raw_signature
  train <- raw[, audit$reference_ids, drop = FALSE]
  finite <- is.finite(train)
  train[!finite] <- NA_real_
  count <- rowSums(finite)
  center <- rowMeans(train, na.rm = TRUE)
  scale <- apply(train, 1L, stats::sd, na.rm = TRUE)
  valid <- count >= 2L & is.finite(scale) & scale > 0
  scale[!valid] <- NA_real_
  z <- sweep(sweep(raw, 1L, center, "-"), 1L, scale, "/")
  z[!is.finite(z)] <- NA_real_
  raw_scores <- matrix(NA_real_, ncol(raw), length(signatures),
    dimnames = list(colnames(raw), names(signatures)))
  for (name in names(signatures)) if (length(signatures[[name]]))
    raw_scores[, name] <- colMeans(z[signatures[[name]], , drop = FALSE])
  definitions <- list(common_raw = raw_scores, common_rank = audit$exploratory_rank)[
    audit$config$score_definitions]
  unscaled <- definitions
  scales <- .ablation_bio_bind(lapply(names(definitions), function(name) {
    x <- definitions[[name]]
    fit <- x[audit$reference_ids, , drop = FALSE]
    fit[!is.finite(fit)] <- NA_real_
    n <- colSums(is.finite(fit))
    mu <- colMeans(fit, na.rm = TRUE)
    sd <- apply(fit, 2L, stats::sd, na.rm = TRUE)
    good <- n >= 2L & is.finite(sd) & sd > 0
    for (j in seq_len(ncol(x))) {
      if (!good[j]) x[, j] <- NA_real_ else x[, j] <- (x[, j] - mu[j]) / sd[j]
    }
    definitions[[name]] <<- x
    data.frame(score_definition = name, anchor = colnames(x),
      reference_valid_count = n, reference_mean = mu, reference_sd = sd,
      status = ifelse(good, "residual_proxy", "not_estimable"),
      reason = ifelse(good, "low_coverage_exploratory_only",
        "missing_signature_or_invalid_reference_scale"))
  }))
  list(values = definitions, unscaled = unscaled, scales = scales,
    gene_scales = data.frame(gene = rownames(raw), reference_valid_count = count,
      reference_mean = center, reference_sd = scale, estimable = valid))
}

.ablation_bio_transform <- function(prepared, distance) {
  d <- .ablation_scale_train_apply(prepared$reference_direct, prepared$query_direct)
  b <- .ablation_module_balanced_transform(prepared$reference_d1,
    prepared$query_d1, prepared$selected_blocks)
  if (distance == "module_unscaled") {
    b$reference <- sweep(sweep(prepared$reference_d1, 2L, b$center, "-"),
      2L, b$weights, "*")
    b$query <- sweep(sweep(prepared$query_d1, 2L, b$center, "-"),
      2L, b$weights, "*")
  }
  if (distance == "tissue_balanced") {
    modules <- prepared$module_manifest$modules
    tissues <- modules$tissue[match(names(prepared$selected_blocks), modules$module_id)]
    if (anyNA(tissues) || any(!.ablation_bio_known(tissues))) return(NULL)
    tissue_n <- table(tissues)
    weights <- numeric(ncol(prepared$reference_d1))
    for (i in seq_along(prepared$selected_blocks)) {
      block <- prepared$selected_blocks[[i]]
      weights[block] <- 1 / sqrt(length(tissue_n) * tissue_n[tissues[i]] * length(block))
    }
    b$reference <- sweep(b$reference, 2L, weights / b$weights, "*")
    b$query <- sweep(b$query, 2L, weights / b$weights, "*")
    b$weights <- weights
  }
  list(Direct = list(reference = d$train, query = d$test),
    d1 = list(reference = b$reference, query = b$query), d1_transform = b)
}

# Exact squared Euclidean distances in bounded batches, with ID tie-breaking.
.ablation_bio_exact <- function(reference, query, k, batch_size = 64L) {
  ids <- matrix(NA_integer_, nrow(query), k)
  distances <- matrix(NA_real_, nrow(query), k)
  ties <- numeric(nrow(query))
  ref_norm <- rowSums(reference^2)
  for (start in seq.int(1L, nrow(query), by = batch_size)) {
    rows <- start:min(nrow(query), start + batch_size - 1L)
    squared <- sweep(sweep(-2 * tcrossprod(query[rows, , drop = FALSE], reference),
      2L, ref_norm, "+"), 1L, rowSums(query[rows, , drop = FALSE]^2), "+")
    squared <- pmax(squared, 0)
    for (j in seq_along(rows)) {
      selected <- order(squared[j, ], rownames(reference))[seq_len(k)]
      ids[rows[j], ] <- selected
      distances[rows[j], ] <- sqrt(squared[j, selected])
      ties[rows[j]] <- sum(abs(squared[j, ] - squared[j, selected[k]]) <= 1e-10)
    }
  }
  list(ids = ids, distances = distances, kth_ties = ties)
}

.ablation_bio_retrieve <- function(prepared, config, pool = "all",
    distance = "baseline", arm = "d1", restriction = NULL) {
  spec <- list(pool = pool, distance = distance, arm = arm, restriction = restriction)
  transformed <- .ablation_bio_transform(prepared, distance)
  if (is.null(transformed)) return(list(status = "not_estimable",
    reason = "unresolved_module_tissue", neighbors = data.frame(), validation = data.frame(),
    spec = spec, validation_lists = list()))
  x <- transformed[[arm]]
  ref <- prepared$reference_metadata
  query <- prepared$query_metadata
  eligibility <- .ablation_bio_eligibility(prepared, config, restriction)
  pools <- if (pool == "same_cancer") eligibility$pools else
    lapply(query$cohort_key, function(key) which(ref$cohort_key != key))
  groups <- split(seq_len(nrow(query)), vapply(pools, function(ids)
    digest::digest(ids), character(1)))
  group_indices <- list()
  validations <- list()
  neighbors <- list()
  for (group in seq_along(groups)) {
    rows <- groups[[group]]
    candidates <- pools[[rows[1L]]]
    if (length(candidates) < config$k) next
    reference <- x$reference[candidates, , drop = FALSE]
    q <- x$query[rows, , drop = FALSE]
    set.seed(.ablation_bio_seed(config, paste(pool, restriction, group)))
    validation_rows <- unlist(lapply(split(seq_along(rows), query$cohort_key[rows]),
      function(i) sample(i, min(length(i), config$retrieval$validation_per_cohort))),
      use.names = FALSE)
    exact <- .ablation_bio_exact(reference, q[validation_rows, , drop = FALSE],
      config$k, config$retrieval$exact_batch_size)
    index <- methods::new(RcppAnnoy::AnnoyEuclidean, ncol(reference))
    index$setSeed(.ablation_bio_seed(config, paste(pool, distance, arm, group)))
    for (i in seq_len(nrow(reference))) index$addItem(i - 1L, reference[i, ])
    index$build(as.integer(config$retrieval$n_trees))
    query_index <- function(query_rows, budget) {
      found <- lapply(query_rows, function(i) index$getNNsByVectorList(q[i, ],
        as.integer(min(config$k + 1L, nrow(reference))), as.integer(budget), TRUE))
      list(ids = do.call(rbind, lapply(found, function(f) (as.integer(f$item) + 1L)[seq_len(config$k)])),
        distances = do.call(rbind, lapply(found, function(f) as.numeric(f$distance)[seq_len(config$k)])),
        boundary_tie = vapply(found, function(f) {
          distances <- as.numeric(f$distance)
          abs(distances[config$k] - distances[config$k - 1L]) <= config$retrieval$tie_tolerance ||
            (length(distances) > config$k &&
              abs(distances[config$k + 1L] - distances[config$k]) <= config$retrieval$tie_tolerance)
        }, logical(1)))
    }
    gate <- FALSE
    for (budget in config$retrieval$search_budgets) {
      found <- query_index(validation_rows, budget)
      strict <- tie_aware <- numeric(length(validation_rows))
      for (j in seq_along(validation_rows)) {
        strict[j] <- length(intersect(found$ids[j, ], exact$ids[j, ])) / config$k
        true_distance <- sqrt(rowSums(sweep(reference[found$ids[j, ], , drop = FALSE],
          2L, q[validation_rows[j], ], "-")^2))
        tie_aware[j] <- mean(true_distance <= exact$distances[j, config$k] +
          config$retrieval$tie_tolerance)
      }
      v <- data.frame(sample_id = query$sample_id[rows[validation_rows]],
        cohort = query$cohort_key[rows[validation_rows]], arm = arm, pool = pool,
        distance = distance, restriction = if (is.null(restriction)) "none" else restriction,
        search_k = budget, strict_recall = strict, recall = tie_aware,
        kth_ties = exact$kth_ties)
      by_cohort <- tapply(v$strict_recall, v$cohort, mean)
      gate <- mean(v$strict_recall) >= config$retrieval$min_mean_recall &&
        all(by_cohort >= config$retrieval$min_cohort_recall)
      if (gate) break
    }
    selected <- if (gate) query_index(seq_len(nrow(q)), budget) else
      .ablation_bio_exact(reference, q, config$k, config$retrieval$exact_batch_size)
    tie_rows <- if (gate) which(selected$boundary_tie) else integer()
    if (length(tie_rows)) {
      tied <- .ablation_bio_exact(reference, q[tie_rows, , drop = FALSE],
        config$k, config$retrieval$exact_batch_size)
      selected$ids[tie_rows, ] <- tied$ids
      selected$distances[tie_rows, ] <- tied$distances
    }
    # Validation stores both lists for score-specific utility error auditing.
    v$method <- if (gate) "annoy_validated" else "exact_fallback"
    v$gate_pass <- gate
    v$exact_tie_query_count <- length(tie_rows)
    validations[[group]] <- v
    group_indices[[group]] <- list(query_ids = query$sample_id[rows[validation_rows]],
      exact_ids = matrix(ref$sample_id[candidates[as.vector(exact$ids)]],
        nrow(exact$ids), config$k),
      approximate_ids = matrix(ref$sample_id[candidates[as.vector(found$ids)]],
        nrow(found$ids), config$k))
    neighbors[[group]] <- data.frame(query_sample = rep(query$sample_id[rows], each = config$k),
      query_cohort = rep(query$cohort_key[rows], each = config$k),
      reference_sample = ref$sample_id[candidates[as.vector(t(selected$ids))]],
      reference_cohort = ref$cohort_key[candidates[as.vector(t(selected$ids))]],
      neighbor_rank = rep(seq_len(config$k), length(rows)),
      neighbor_distance = as.vector(t(selected$distances)), arm = arm,
      pool = pool, distance = distance,
      restriction = if (is.null(restriction)) "none" else restriction,
      method = rep(if (gate) ifelse(seq_len(nrow(q)) %in% tie_rows,
        "exact_tie_boundary", "annoy_validated") else "exact_fallback", each = config$k))
  }
  list(status = if (length(neighbors)) "estimable" else "not_estimable",
    reason = if (length(neighbors)) "supported" else "no_eligible_queries",
    neighbors = .ablation_bio_bind(neighbors), validation = .ablation_bio_bind(validations),
    validation_lists = group_indices, eligibility = eligibility$queries, spec = spec)
}

.ablation_bio_ridge_design <- function(train, test) {
  transform <- .ablation_scale_train_apply(train, test)
  gram <- crossprod(transform$train) / nrow(transform$train)
  list(transform = transform, eig = eigen(gram, symmetric = TRUE))
}

.ablation_bio_ridge <- function(train, target, test, lambda, design = NULL) {
  if (is.null(design)) design <- .ablation_bio_ridge_design(train, test)
  transform <- design$transform
  center <- mean(target)
  scale <- stats::sd(target)
  if (!is.finite(scale) || scale <= 0) return(NULL)
  x <- transform$train
  y <- (target - center) / scale
  eigenvectors <- design$eig$vectors
  beta <- eigenvectors %*% (crossprod(eigenvectors, crossprod(x, y) / nrow(x)) /
    (pmax(design$eig$values, 0) + lambda))
  list(prediction = as.numeric(transform$test %*% beta) * scale + center,
    fit = list(beta = beta, feature_center = transform$center,
      feature_scale = transform$scale, target_center = center,
      target_scale = scale, lambda = lambda))
}

.ablation_bio_readout <- function(prepared, audit, scores, config) {
  ref <- prepared$reference_metadata
  query <- prepared$query_metadata
  cohorts <- sort(unique(ref$cohort_key))
  if (length(cohorts) < config$readout$folds) {
    return(list(status = "not_estimable", reason = "fewer_than_five_reference_cohorts",
      predictions = data.frame(), cv = data.frame(), fits = list()))
  }
  set.seed(.ablation_bio_seed(config, "reference-folds"))
  folds <- stats::setNames(rep(seq_len(config$readout$folds), length.out = length(cohorts)),
    sample(cohorts))
  assignments <- folds[ref$cohort_key]
  designs <- list()
  for (arm in c("Direct", "d1")) {
    x <- prepared[[if (arm == "Direct") "reference_direct" else "reference_d1"]]
    q <- prepared[[if (arm == "Direct") "query_direct" else "query_d1"]]
    designs[[paste(arm, "full")]] <- .ablation_bio_ridge_design(x, q)
    for (fold in seq_len(config$readout$folds)) {
      designs[[paste(arm, fold)]] <- .ablation_bio_ridge_design(
        x[assignments != fold, , drop = FALSE], x[assignments == fold, , drop = FALSE])
    }
  }
  cv <- list()
  predictions <- list()
  fits <- list()
  i <- 0L
  # Raw gene scaling is refitted inside every CV training fold. Rank scoring
  # has no learned transform; its target scale is still training-fold-only.
  for (definition in config$score_definitions) for (anchor in config$primary_anchors) {
    if (!audit$coverage$estimable[match(anchor, audit$coverage$anchor)]) next
    final_y <- scores$unscaled[[definition]][ref$sample_id, anchor]
    query_y <- scores$unscaled[[definition]][query$sample_id, anchor]
    if (any(!is.finite(final_y)) || !is.finite(stats::sd(final_y)) ||
        stats::sd(final_y) <= 0) next
    tuning <- list()
    for (fold in seq_len(config$readout$folds)) {
      train <- which(assignments != fold)
      test <- which(assignments == fold)
      y <- if (definition == "common_raw") .ablation_bio_raw_scores(
        audit$raw_signature, audit$common_signatures, ref$sample_id[train])[, anchor] else
          scores$unscaled[[definition]][, anchor]
      scale <- stats::sd(y[ref$sample_id[train]])
      if (any(!is.finite(y[ref$sample_id])) || !is.finite(scale) || scale <= 0) next
      for (arm in c("Direct", "d1")) for (lambda in config$readout$lambda) {
        x <- prepared[[if (arm == "Direct") "reference_direct" else "reference_d1"]]
        fit <- .ablation_bio_ridge(x[train, , drop = FALSE], y[ref$sample_id[train]],
          x[test, , drop = FALSE], lambda, designs[[paste(arm, fold)]])
        if (is.null(fit)) next
        errors <- data.frame(cohort = ref$cohort_key[test],
          mae = abs(fit$prediction - y[ref$sample_id[test]]) / scale,
          mean_mae = abs(mean(y[ref$sample_id[train]]) - y[ref$sample_id[test]]) / scale)
        result <- stats::aggregate(cbind(mae, mean_mae) ~ cohort, errors, mean)
        result$fold <- fold
        result$arm <- arm
        result$lambda <- lambda
        result$anchor <- anchor
        result$score_definition <- definition
        tuning[[length(tuning) + 1L]] <- result
      }
    }
    tuning <- .ablation_bio_bind(tuning)
    cv[[length(cv) + 1L]] <- tuning
    if (!nrow(tuning) || length(unique(tuning$fold)) != config$readout$folds) next
    selected <- stats::aggregate(mae ~ arm + lambda, tuning, mean)
    best <- lapply(split(selected, selected$arm), function(x) x$lambda[which.min(x$mae)])
    direct_cv <- min(selected$mae[selected$arm == "Direct"])
    d1_cv <- min(selected$mae[selected$arm == "d1"])
    nonlinear <- d1_cv > direct_cv || d1_cv >= mean(tuning$mean_mae)
    for (reader in c("ridge", if (nonlinear) "xgboost")) for (arm in c("Direct", "d1")) {
      train_x <- prepared[[if (arm == "Direct") "reference_direct" else "reference_d1"]]
      query_x <- prepared[[if (arm == "Direct") "query_direct" else "query_d1"]]
      if (reader == "ridge") {
        fitted <- .ablation_bio_ridge(train_x, final_y, query_x, best[[arm]],
          designs[[paste(arm, "full")]])
        predicted <- fitted$prediction
        fit <- fitted$fit
      } else {
        transform <- .ablation_scale_train_apply(train_x, query_x)
        parameters <- config$readout$xgboost
        rounds <- parameters$nrounds
        parameters$nrounds <- NULL
        parameters$seed <- .ablation_bio_seed(config, paste(definition, anchor, "xgboost"))
        set.seed(parameters$seed)
        fit <- xgboost::xgb.train(params = parameters,
          data = xgboost::xgb.DMatrix(transform$train,
            label = (final_y - mean(final_y)) / stats::sd(final_y)),
          nrounds = rounds, verbose = 0)
        predicted <- as.numeric(stats::predict(fit,
          xgboost::xgb.DMatrix(transform$test))) * stats::sd(final_y) + mean(final_y)
        fit <- list(model = xgboost::xgb.save.raw(fit),
          feature_center = transform$center, feature_scale = transform$scale,
          target_center = mean(final_y), target_scale = stats::sd(final_y),
          params = parameters, nrounds = rounds)
      }
      i <- i + 1L
      fits[[paste(definition, anchor, reader, arm, sep = ":")]] <- fit
      cancer_mean <- vapply(as.character(query$cancer_type), function(label) {
        supported <- ref$cancer_type == label & ref$metadata_status == "confirmed"
        if (sum(supported, na.rm = TRUE) < config$k ||
            length(unique(ref$cohort_key[which(supported)])) < config$min_reference_cohorts)
          return(NA_real_)
        mean(final_y[which(supported)])
      }, numeric(1))
      predictions[[i]] <- data.frame(sample_id = query$sample_id,
        cohort = query$cohort_key, score_definition = definition, anchor = anchor,
        reader = reader, arm = arm, truth = query_y, prediction = predicted,
        reference_sd = stats::sd(final_y), reference_mean = mean(final_y),
        cancer_mean = cancer_mean, nonlinear_trigger = nonlinear)
    }
  }
  list(status = if (length(predictions)) "estimable" else "not_estimable",
    reason = if (length(predictions)) "supported" else "no_valid_targets_or_folds",
    predictions = .ablation_bio_bind(predictions), cv = .ablation_bio_bind(cv),
    folds = data.frame(cohort = names(folds), fold = as.integer(folds)), fits = fits,
    training_ids = ref$sample_id, query_ids = query$sample_id)
}

.ablation_bio_utility <- function(neighbors, scores, prepared, config) {
  ref <- prepared$reference_metadata
  query <- prepared$query_metadata
  empty <- data.frame(experiment = character(), pool = character(),
    distance = character(), restriction = character(), arm = character(),
    query_sample = character(), query_cohort = character(), utility = numeric(),
    abs_delta = numeric(), cancer_match = numeric(), score_definition = character(),
    anchor = character())
  if (!nrow(neighbors)) return(empty)
  pieces <- list()
  index <- 0L
  for (definition in names(scores$values)) {
    x <- scores$values[[definition]]
    for (anchor in colnames(x)) {
      n <- neighbors
      n$abs_delta <- abs(x[match(n$query_sample, rownames(x)), anchor] -
        x[match(n$reference_sample, rownames(x)), anchor])
      n$utility <- exp(-n$abs_delta)
      n$cancer_match <- ref$cancer_type[match(n$reference_sample, ref$sample_id)] ==
        query$cancer_type[match(n$query_sample, query$sample_id)]
      n <- n[is.finite(n$utility), ]
      if (!nrow(n)) next
      n$experiment <- paste(n$pool, n$distance, n$restriction, sep = ":")
      key <- interaction(n$experiment, n$arm, n$query_sample, drop = TRUE)
      n <- n[as.character(key) %in% names(which(table(key) == config$k)), ]
      if (!nrow(n)) next
      part <- stats::aggregate(n[, c("utility", "abs_delta", "cancer_match"), drop = FALSE],
        n[, c("experiment", "pool", "distance", "restriction", "arm", "query_sample",
          "query_cohort"), drop = FALSE],
        function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_)
      part$score_definition <- definition
      part$anchor <- anchor
      index <- index + 1L
      pieces[[index]] <- part
    }
  }
  result <- .ablation_bio_bind(pieces)
  if (nrow(result)) result else empty
}

# One shared cohort index matrix is used for every column of a comparison.
# Missing cells remain explicit; a replicate cannot silently change its cohort
# denominator by omitting unavailable cohorts.
.ablation_bio_infer_matrix <- function(values, config, id, query_counts,
    test = TRUE) {
  n <- nrow(values)
  cohort_ids <- rownames(values)
  if (is.null(cohort_ids)) cohort_ids <- as.character(seq_len(n))
  set.seed(.ablation_bio_seed(config, paste("cohort-bootstrap", paste(cohort_ids, collapse = ";"))))
  draws <- matrix(sample.int(max(1L, n), config$bootstrap * max(1L, n), TRUE),
    nrow = config$bootstrap)
  indices_by_anchor <- list()
  results <- lapply(seq_len(ncol(values)), function(j) {
    x <- values[, j]
    finite <- is.finite(x)
    estimate <- if (any(finite)) mean(x[finite]) else NA_real_
    # Distinct effective populations receive distinct explicit cohort lists;
    # every comparison on the same cohort list reuses the same seed/indices.
    selected <- which(finite)
    set.seed(.ablation_bio_seed(config, paste("cohort-bootstrap",
      paste(cohort_ids[selected], collapse = ";"))))
    indices <- if (length(selected)) matrix(sample.int(length(selected),
      config$bootstrap * length(selected), TRUE), nrow = config$bootstrap) else
        matrix(integer(), config$bootstrap, 0L)
    indices_by_anchor[[colnames(values)[j]]] <<- list(cohort_ids = cohort_ids[selected],
      indices = indices)
    bootstrap <- if (length(selected)) apply(indices, 1L,
      function(i) mean(x[selected[i]])) else numeric()
    valid <- sum(is.finite(bootstrap))
    supported <- sum(finite) >= config$min_inference_cohorts &&
      valid >= config$min_valid_bootstrap
    ci <- if (supported) unname(stats::quantile(bootstrap[is.finite(bootstrap)],
      c(.025, .975))) else c(NA_real_, NA_real_)
    p <- NA_real_
    if (test && supported) {
      x <- x[finite]
      if (length(x) <= 18L) {
        null <- vapply(0:(2^length(x) - 1L), function(mask) {
          signs <- ifelse(as.integer(intToBits(mask))[seq_along(x)] == 1L, 1, -1)
          mean(x * signs)
        }, numeric(1))
        p <- mean(abs(null) >= abs(estimate) - 1e-14)
      } else {
        signs <- matrix(sample(c(-1, 1), length(x) * config$bootstrap, TRUE),
          nrow = length(x))
        p <- (1 + sum(abs(colMeans(signs * x)) >= abs(estimate) - 1e-14)) /
          (config$bootstrap + 1L)
      }
    }
    data.frame(anchor = colnames(values)[j], estimate = estimate,
      ci_low = ci[1L], ci_high = ci[2L], p_value = p,
      p_method = if (!test) "not_tested" else if (!supported) "not_estimable" else
        if (sum(finite) <= 18L) "exact_paired_sign_flip" else "monte_carlo_paired_sign_flip",
      cohort_count = sum(finite), query_count = sum(query_counts[finite, j]),
      n_boot = config$bootstrap, valid_bootstrap = valid,
      seed = .ablation_bio_seed(config, paste("cohort-bootstrap",
        paste(cohort_ids[selected], collapse = ";"))),
      status = if (supported) "estimable" else "not_estimable",
      reason = if (supported) "supported" else if (sum(finite) < config$min_inference_cohorts)
        "fewer_than_three_cohorts" else "fewer_than_1900_valid_bootstrap")
  })
  result <- .ablation_bio_bind(results)
  result$q_value <- stats::p.adjust(result$p_value, "BH")
  list(summary = result, draws = draws, indices_by_anchor = indices_by_anchor)
}

.ablation_bio_summarize <- function(paired, config, id, anchors, effect = "delta",
    test = TRUE, local_intervals = TRUE) {
  if (!nrow(paired)) {
    empty <- matrix(numeric(), 0L, length(anchors), dimnames = list(NULL, anchors))
    return(list(summary = .ablation_bio_infer_matrix(empty, config, id, empty, test)$summary,
      cohort = data.frame(anchor = character(), cohort = character(), value = numeric(),
        n = integer(), eligible_inference = logical(), reason = character()),
      local = data.frame(anchor = character(), cohort = character(), query_count = integer(),
        estimate = numeric(), ci_low = numeric(), ci_high = numeric(), scope = character(),
        p_value = numeric()), draws = matrix(integer(), 0L, 0L), cohort_ids = character()))
  }
  paired$value <- paired[[effect]]
  means <- stats::aggregate(value ~ anchor + cohort, paired, mean)
  counts <- stats::aggregate(sample_id ~ anchor + cohort, paired, length)
  means$n <- counts$sample_id[match(paste(means$anchor, means$cohort),
    paste(counts$anchor, counts$cohort))]
  means$eligible_inference <- means$n >= config$min_query_per_cohort
  means$reason <- ifelse(means$eligible_inference, "supported", "fewer_than_twenty_queries")
  cohorts <- sort(unique(means$cohort[means$eligible_inference]))
  values <- matrix(NA_real_, length(cohorts), length(anchors),
    dimnames = list(cohorts, anchors))
  sizes <- values
  sizes[] <- 0L
  use <- means[means$eligible_inference, ]
  positions <- cbind(match(use$cohort, cohorts), match(use$anchor, anchors))
  good <- !is.na(positions[, 2L])
  values[positions[good, , drop = FALSE]] <- use$value[good]
  sizes[positions[good, , drop = FALSE]] <- use$n[good]
  inferred <- .ablation_bio_infer_matrix(values, config, id, sizes, test)
  local <- if (!local_intervals) data.frame() else .ablation_bio_bind(lapply(split(paired,
    interaction(paired$anchor, paired$cohort, drop = TRUE)), function(x) {
    set.seed(.ablation_bio_seed(config, paste(id, x$anchor[1L], x$cohort[1L], "local")))
    ci <- if (nrow(x) >= config$min_query_per_cohort) {
      boot <- replicate(config$bootstrap, mean(sample(x$value, nrow(x), TRUE)))
      unname(stats::quantile(boot, c(.025, .975)))
    } else c(NA_real_, NA_real_)
    data.frame(anchor = x$anchor[1L], cohort = x$cohort[1L], query_count = nrow(x),
      estimate = mean(x$value), ci_low = ci[1L], ci_high = ci[2L],
      scope = "fixed_cohort_paired_query_pointwise", p_value = NA_real_)
  }))
  list(summary = inferred$summary, cohort = means, local = local,
    draws = inferred$draws, indices_by_anchor = inferred$indices_by_anchor,
    cohort_ids = cohorts, values = values)
}

.ablation_bio_pair_utility <- function(per_query) {
  keys <- c("experiment", "score_definition", "anchor", "query_sample", "query_cohort")
  d <- per_query[per_query$arm == "Direct", c(keys, "utility", "abs_delta", "cancer_match")]
  b <- per_query[per_query$arm == "d1", c(keys, "utility", "abs_delta", "cancer_match")]
  paired <- merge(d, b, by = keys, suffixes = c("_direct", "_d1"))
  paired$delta <- paired$utility_d1 - paired$utility_direct
  names(paired)[names(paired) == "query_sample"] <- "sample_id"
  names(paired)[names(paired) == "query_cohort"] <- "cohort"
  paired
}

.ablation_bio_exploratory <- function(prepared, audit, scores, retrieval, config) {
  anchors <- config$primary_anchors
  neighbors <- .ablation_bio_bind(lapply(retrieval, `[[`, "neighbors"))
  neighbors <- neighbors[neighbors$restriction == "none", ]
  per_query <- .ablation_bio_utility(neighbors, scores, prepared, config)
  for (distance in c("module_unscaled", "tissue_balanced")) {
    direct <- per_query[per_query$arm == "Direct" & per_query$distance == "baseline" &
      per_query$score_definition == "common_rank", ]
    if (!nrow(direct)) next
    direct$distance <- rep(distance, nrow(direct))
    direct$experiment <- paste(direct$pool, distance, "none", sep = ":")
    per_query <- rbind(per_query, direct)
  }
  per_query <- per_query[per_query$distance == "baseline" |
    per_query$score_definition == "common_rank", ]
  paired <- .ablation_bio_pair_utility(per_query)
  rows <- cohorts <- common_queries <- draws <- list()
  add <- function(p, id, family) {
    result <- .ablation_bio_summarize(p, config, id, anchors,
      test = FALSE, local_intervals = FALSE)
    summary <- result$summary
    counts <- lengths(audit$common_signatures[summary$anchor])
    single <- counts == 1L
    summary$ci_low[single] <- summary$ci_high[single] <- NA_real_
    summary$p_value <- summary$q_value <- NA_real_
    summary$comparison <- id
    summary$family <- family
    summary$residual_gene_count <- counts
    summary$measurement <- ifelse(single, "single_gene_description", "residual_multigene_proxy")
    summary$scope <- "conditional_on_fixed_residual_genes_bank_and_reference_scale"
    summary$interval_method <- ifelse(single, "not_applicable_single_gene",
      "cohort_percentile_bootstrap_95_percent")
    summary$status[summary$status == "estimable"] <- ifelse(
      single[summary$status == "estimable"], "descriptive_only", "exploratory")
    summary$reason[summary$status %in% c("descriptive_only", "exploratory")] <-
      "low_coverage_not_full_anchor_evidence"
    result$cohort$comparison <- rep(id, nrow(result$cohort))
    rows[[id]] <<- summary
    cohorts[[id]] <<- result$cohort
    common_queries[[id]] <<- p[, c("sample_id", "cohort", "anchor")]
    draws[[id]] <<- result$indices_by_anchor
  }
  grid_ids <- unlist(lapply(config$score_definitions, function(definition)
    paste(definition, c("all", "same_cancer"), "baseline", "none", sep = ":")))
  grid <- paired[paste(paired$score_definition, paired$experiment, sep = ":") %in% grid_ids, ]
  common <- .ablation_bio_bind(lapply(split(grid, grid$anchor), function(x) {
    groups <- split(x, paste(x$score_definition, x$experiment, sep = ":"))
    if (!setequal(names(groups), grid_ids)) return(NULL)
    ids <- Reduce(intersect, lapply(groups, `[[`, "sample_id"))
    x[x$sample_id %in% ids, ]
  }))
  if (!nrow(common)) common <- grid[FALSE, ]
  for (id in grid_ids) add(common[paste(common$score_definition, common$experiment,
    sep = ":") == id, ], paste0("grid:", id), "matched_score_pool_grid")
  changes <- list()
  for (definition in config$score_definitions) changes[[paste0(sub("common_", "", definition),
    "_pool")]] <- paste(definition, c("all", "same_cancer"), "baseline", "none", sep = ":")
  if (all(c("common_raw", "common_rank") %in% config$score_definitions)) {
    changes$all_score <- paste(c("common_raw", "common_rank"), "all:baseline:none", sep = ":")
    changes$same_score <- paste(c("common_raw", "common_rank"), "same_cancer:baseline:none", sep = ":")
  }
  for (id in names(changes)) {
    keys <- paste(common$score_definition, common$experiment, sep = ":")
    joined <- merge(common[keys == changes[[id]][1L], ], common[keys == changes[[id]][2L], ],
      by = c("sample_id", "cohort", "anchor"), suffixes = c("_a", "_b"))
    joined$delta <- joined$delta_b - joined$delta_a
    add(joined, paste0("change:", id), "matched_difference_in_differences")
  }
  for (definition in config$score_definitions) add(paired[paired$score_definition == definition &
    paired$experiment == "all:baseline:none", ], paste0("coverage:", definition), "broad_population")
  for (pool in c("all", "same_cancer")) for (distance in c("module_unscaled", "tissue_balanced")) {
    part <- paired[paired$score_definition == "common_rank" &
      paired$experiment == paste(pool, distance, "none", sep = ":"), ]
    base <- paired[paired$score_definition == "common_rank" &
      paired$experiment == paste(pool, "baseline", "none", sep = ":"), ]
    joined <- merge(base, part, by = c("sample_id", "cohort", "anchor"), suffixes = c("_a", "_b"))
    joined$delta <- joined$delta_b - joined$delta_a
    add(joined, paste0("distance_change:", pool, ":", distance), "matched_distance_changes")
    ids <- paste(joined$sample_id, joined$anchor)
    add(part[paste(part$sample_id, part$anchor) %in% ids, ],
      paste0("distance:", pool, ":", distance), "matched_distance_levels")
  }
  metadata <- rbind(prepared$reference_metadata, prepared$query_metadata)
  availability <- .ablation_bio_bind(lapply(names(scores$values), function(definition) {
    x <- scores$values[[definition]]
    .ablation_bio_bind(lapply(anchors, function(anchor) {
      .ablation_bio_bind(lapply(split(seq_len(nrow(metadata)), metadata$cohort_key), function(i) {
        valid <- is.finite(x[match(metadata$sample_id[i], rownames(x)), anchor])
        data.frame(cohort = metadata$cohort_key[i[1L]],
          role = if (all(metadata$sample_id[i] %in% audit$reference_ids)) "reference" else "query",
          score_definition = definition, anchor = anchor, sample_count = length(i),
          valid_count = sum(valid), missing_count = sum(!valid))
      }))
    }))
  }))
  list(inference = .ablation_bio_bind(rows), cohort = .ablation_bio_bind(cohorts),
    scales = scores$scales, gene_scales = scores$gene_scales, availability = availability,
    common_queries = common_queries, bootstrap_indices = draws)
}

.ablation_bio_prediction_metrics <- function(readout, config) {
  predictions <- readout$predictions
  if (!nrow(predictions)) return(data.frame())
  .ablation_bio_bind(lapply(split(predictions, interaction(predictions$score_definition,
    predictions$anchor, predictions$reader, predictions$cohort, drop = TRUE)), function(x) {
    ids <- Reduce(intersect, lapply(split(x, x$arm), function(a)
      a$sample_id[is.finite(a$prediction) & is.finite(a$truth)]))
    if (!length(ids)) return(NULL)
    x <- x[x$sample_id %in% ids, ]
    .ablation_bio_bind(lapply(split(x, x$arm), function(a) {
      scale <- a$reference_sd[1L]
      correlation <- if (nrow(a) >= 2L && stats::sd(a$truth) > 0 &&
          stats::sd(a$prediction) > 0) stats::cor(a$truth, a$prediction,
            method = "spearman") else NA_real_
      cancer <- is.finite(a$cancer_mean)
      data.frame(score_definition = a$score_definition[1L], anchor = a$anchor[1L],
        reader = a$reader[1L], cohort = a$cohort[1L], arm = a$arm[1L],
        query_count = nrow(a), mae = mean(abs(a$prediction - a$truth)) / scale,
        correlation = correlation,
        correlation_reason = if (is.finite(correlation)) "supported" else "truth_or_prediction_zero_variance",
        reference_mean_mae = mean(abs(a$reference_mean - a$truth)) / scale,
        cancer_mean_count = sum(cancer),
        cancer_mean_mae = if (any(cancer)) mean(abs(a$cancer_mean[cancer] - a$truth[cancer])) / scale else NA_real_,
        model_mae_cancer_subset = if (any(cancer)) mean(abs(a$prediction[cancer] - a$truth[cancer])) / scale else NA_real_)
    }))
  }))
}

.ablation_bio_module_diagnostics <- function(prepared, retrieval) {
  ref <- prepared$reference_metadata
  query <- prepared$query_metadata
  blocks <- prepared$selected_blocks
  modules <- prepared$module_manifest$modules
  module_means <- list()
  saturation <- list()
  contributions <- list()
  for (role in c("reference", "query")) {
    x <- prepared[[paste0(role, "_d1")]]
    module_means[[role]] <- sapply(blocks, function(block)
      rowMeans(x[, block, drop = FALSE]))
    scale <- apply(prepared$reference_d1, 2L, stats::sd)
    scale[!is.finite(scale) | scale == 0] <- 1
    z <- sweep(sweep(x, 2L, colMeans(prepared$reference_d1), "-"), 2L, scale, "/")
    saturation[[role]] <- .ablation_bio_bind(lapply(names(blocks), function(name) {
      block <- blocks[[name]]
      data.frame(role = role, module = name,
        tissue = modules$tissue[match(name, modules$module_id)],
        sample_count = nrow(x), column_count = length(block),
        near_zero_fraction = mean(x[, block, drop = FALSE] <= 1e-6),
        near_one_fraction = mean(x[, block, drop = FALSE] >= 1 - 1e-6),
        mean_absolute_z = mean(abs(z[, block, drop = FALSE])))
    }))
  }
  for (result in retrieval) {
    n <- result$neighbors
    if (!nrow(n) || n$arm[1L] != "d1" || n$restriction[1L] != "none") next
    transform <- .ablation_bio_transform(prepared, n$distance[1L])$d1
    q <- match(n$query_sample, query$sample_id)
    r <- match(n$reference_sample, ref$sample_id)
    # Aggregate to cohort/module, rather than exporting patient-level terms.
    for (name in names(blocks)) {
      block <- blocks[[name]]
      squared <- rowSums((transform$query[q, block, drop = FALSE] -
        transform$reference[r, block, drop = FALSE])^2)
      value <- stats::aggregate(squared, list(cohort = n$query_cohort), mean)
      names(value)[2L] <- "mean_squared_contribution"
      value$module <- name
      value$distance <- n$distance[1L]
      value$pool <- n$pool[1L]
      contributions[[length(contributions) + 1L]] <- value
    }
  }
  list(saturation = .ablation_bio_bind(saturation),
    module_correlations = lapply(module_means, stats::cor),
    distance_contributions = .ablation_bio_bind(contributions))
}

.ablation_bio_inference <- function(prepared, audit, scores, retrieval, readout, config) {
  primary_score <- config$primary_score_definition
  if (is.null(primary_score)) primary_score <- if ("common_rank" %in%
    config$score_definitions) "common_rank" else "common_raw"
  neighbors <- .ablation_bio_bind(lapply(retrieval, `[[`, "neighbors"))
  per_query <- .ablation_bio_utility(neighbors, scores, prepared, config)
  # Direct remains at its baseline distance in each d1 distance comparison.
  for (distance in c("module_unscaled", "tissue_balanced")) for (pool in c("all", "same_cancer")) {
    direct <- per_query[per_query$arm == "Direct" & per_query$distance == "baseline" &
      per_query$pool == pool & per_query$restriction == "none" &
      per_query$score_definition == primary_score, ]
    if (!nrow(direct)) next
    direct$distance <- distance
    direct$experiment <- paste(pool, distance, "none", sep = ":")
    per_query <- rbind(per_query, direct)
  }
  # Distance/technical branches use only the declared primary score.
  per_query <- per_query[per_query$distance == "baseline" & per_query$restriction == "none" |
    per_query$score_definition == primary_score, ]
  paired <- .ablation_bio_pair_utility(per_query)
  rows <- cohorts <- locals <- lists <- draws <- list()
  add <- function(p, id, definition, experiment, family, anchors, effect = "delta") {
    result <- .ablation_bio_summarize(p, config, id, anchors, effect)
    result$summary$comparison <- id
    result$summary$score_definition <- definition
    result$summary$experiment <- experiment
    result$summary$family <- family
    result$cohort$comparison <- rep(id, nrow(result$cohort))
    result$local$comparison <- rep(id, nrow(result$local))
    rows[[id]] <<- result$summary
    cohorts[[id]] <<- result$cohort
    locals[[id]] <<- result$local
    lists[[id]] <<- p[, c("sample_id", "cohort", "anchor")]
    draws[[id]] <<- list(cohort_ids = result$cohort_ids, indices = result$draws,
      by_anchor = result$indices_by_anchor)
  }
  anchors <- config$primary_anchors
  grid_ids <- unlist(lapply(config$score_definitions, function(definition)
    paste(definition, c("all", "same_cancer"), "baseline", "none", sep = ":")),
    use.names = FALSE)
  grid <- paired[paired$anchor %in% anchors &
    paste(paired$score_definition, paired$experiment, sep = ":") %in% grid_ids, ]
  common_grid <- .ablation_bio_bind(lapply(split(grid, grid$anchor), function(x) {
    groups <- split(x, paste(x$score_definition, x$experiment, sep = ":"))
    if (!setequal(names(groups), grid_ids)) return(NULL)
    ids <- Reduce(intersect, lapply(groups, `[[`, "sample_id"))
    x <- x[x$sample_id %in% ids, ]
    # Retain small cohorts descriptively; summarize() applies the inference gate.
    x
  }))
  if (!nrow(common_grid)) common_grid <- grid[FALSE, , drop = FALSE]
  for (id in grid_ids) {
    part <- common_grid[paste(common_grid$score_definition, common_grid$experiment,
      sep = ":") == id, ]
    add(part, paste0("grid:", id), sub(":.*$", "", id),
      sub("^[^:]+:", "", id), "four_anchor_grid", anchors)
  }
  # Difference-in-differences always joins at the paired patient level first.
  changes <- list()
  for (definition in config$score_definitions) {
    changes[[paste0(sub("common_", "", definition), "_pool")]] <-
      paste(definition, c("all", "same_cancer"), "baseline", "none", sep = ":")
  }
  if (all(c("common_raw", "common_rank") %in% config$score_definitions)) {
    changes$all_score <- paste(c("common_raw", "common_rank"),
      "all:baseline:none", sep = ":")
    changes$same_score <- paste(c("common_raw", "common_rank"),
      "same_cancer:baseline:none", sep = ":")
  }
  for (id in names(changes)) {
    comparison <- changes[[id]]
    a <- common_grid[paste(common_grid$score_definition, common_grid$experiment,
      sep = ":") == comparison[1L], ]
    b <- common_grid[paste(common_grid$score_definition, common_grid$experiment,
      sep = ":") == comparison[2L], ]
    p <- merge(a, b, by = c("sample_id", "cohort", "anchor"), suffixes = c("_a", "_b"))
    p$delta <- p$delta_b - p$delta_a
    add(p, paste0("change:", id), paste(comparison, collapse = " -> "), id,
      "prespecified_grid_changes", anchors)
  }
  for (definition in config$score_definitions) {
    p <- paired[paired$score_definition == definition &
      paired$experiment == "all:baseline:none" & paired$anchor %in% anchors, ]
    add(p, paste0("coverage:", definition), definition, "all:baseline:none",
      "broad_score_sensitivity", anchors)
  }
  baseline <- paired[paired$score_definition == primary_score &
    paired$experiment == "same_cancer:baseline:none" & paired$anchor %in% anchors, ]
  for (distance in c("module_unscaled", "tissue_balanced")) for (pool in c("same_cancer", "all")) {
    p <- paired[paired$score_definition == primary_score &
      paired$experiment == paste(pool, distance, "none", sep = ":") & paired$anchor %in% anchors, ]
    base <- paired[paired$score_definition == primary_score &
      paired$experiment == paste(pool, "baseline", "none", sep = ":") & paired$anchor %in% anchors, ]
    joined <- merge(base, p, by = c("sample_id", "cohort", "anchor"), suffixes = c("_a", "_b"))
    if (!nrow(joined)) next
    joined$delta <- joined$delta_b - joined$delta_a
    add(joined, paste0("distance_change:", pool, ":", distance), primary_score,
      paste(pool, distance, sep = ":"), "prespecified_distance_changes", anchors)
    ids <- paste(joined$sample_id, joined$anchor)
    p <- p[paste(p$sample_id, p$anchor) %in% ids, ]
    add(p, paste0("distance:", pool, ":", distance), primary_score,
      paste(pool, distance, sep = ":"), "distance_levels", anchors)
  }
  for (definition in config$score_definitions) for (family in c("split_ifn_il6", "direct_input_disjoint")) {
    supplementary <- if (family == "split_ifn_il6") c("ifn", "il6") else paste0(anchors, "_disjoint")
    p <- paired[paired$score_definition == definition & paired$experiment == "all:baseline:none" &
      paired$anchor %in% supplementary, ]
    add(p, paste("supplement", definition, family, sep = ":"), definition,
      "all:baseline:none", family, supplementary)
  }
  for (restriction in config$technical_restrictions) {
    p <- paired[paired$score_definition == primary_score &
      paired$experiment == paste("same_cancer", "baseline", restriction, sep = ":") &
      paired$anchor %in% anchors, ]
    if (!nrow(p)) next
    query <- prepared$query_metadata
    p$stratum <- query[[restriction]][match(p$sample_id, query$sample_id)]
    for (label in unique(p$stratum)) {
      x <- p[p$stratum == label, ]
      b <- merge(baseline, x, by = c("sample_id", "cohort", "anchor"), suffixes = c("_a", "_b"))
      if (!nrow(b)) next
      b$delta <- b$delta_b - b$delta_a
      add(b, paste("technical_change", restriction, label, sep = ":"), primary_score,
        paste(restriction, label, sep = ":"), "technical_source_changes", anchors)
    }
  }
  metrics <- .ablation_bio_prediction_metrics(readout, config)
  if (nrow(metrics)) for (definition in unique(metrics$score_definition)) for (reader in unique(metrics$reader)) {
    x <- metrics[metrics$score_definition == definition & metrics$reader == reader, ]
    d <- x[x$arm == "Direct", ]
    b <- x[x$arm == "d1", ]
    p <- merge(d, b, by = c("cohort", "anchor"), suffixes = c("_direct", "_d1"))
    p <- p[p$query_count_direct >= config$min_query_per_cohort, ]
    if (!nrow(p)) next
    for (measure in c("mae", "correlation")) {
      delta <- p[[paste0(measure, "_d1")]] - p[[paste0(measure, "_direct")]]
      cohort_ids <- sort(unique(p$cohort))
      values <- sizes <- matrix(NA_real_, length(cohort_ids), length(anchors),
        dimnames = list(cohort_ids, anchors))
      indices <- cbind(match(p$cohort, cohort_ids), match(p$anchor, anchors))
      values[indices] <- delta
      sizes[indices] <- p$query_count_direct
      result <- .ablation_bio_infer_matrix(values, config,
        paste("readout", definition, reader, measure, sep = ":"), sizes,
        test = measure == "mae")$summary
      result$comparison <- paste("readout", definition, reader, measure, sep = ":")
      result$score_definition <- definition
      result$experiment <- paste(reader, measure, sep = ":")
      result$family <- if (measure == "mae") "readout_mae_four_anchors_by_reader" else "readout_correlation_descriptive"
      rows[[result$comparison[1L]]] <- result
    }
  }
  inference <- .ablation_bio_bind(rows)
  if (nrow(inference)) {
    family_id <- ifelse(inference$family %in% c("four_anchor_grid", "distance_levels",
      "broad_score_sensitivity"), inference$comparison, inference$family)
    families <- interaction(inference$score_definition, family_id, drop = TRUE)
    inference$q_value <- ave(inference$p_value, families, FUN = function(p) stats::p.adjust(p, "BH"))
  }
  cohort <- .ablation_bio_bind(cohorts)
  # Symmetric leave-one-cohort and leave-one-source sensitivities, including
  # GSE21983, retain every omission regardless of observed effect direction.
  leave_one <- .ablation_bio_bind(lapply(split(cohort, interaction(cohort$comparison,
    cohort$anchor, drop = TRUE)), function(x) {
    x <- x[x$eligible_inference, ]
    if (!nrow(x)) return(NULL)
    query <- prepared$query_metadata
    lookup <- unique(query[, c("cohort_key", "source_system")])
    x$source <- lookup$source_system[match(x$cohort, lookup$cohort_key)]
    .ablation_bio_bind(lapply(c("cohort", "source"), function(unit) {
      groups <- sort(unique(x[[unit]]))
      .ablation_bio_bind(lapply(groups, function(label) {
        remaining <- x[x[[unit]] != label, ]
        data.frame(comparison = x$comparison[1L], anchor = x$anchor[1L],
          omitted_unit = unit, omitted = label, remaining_cohorts = nrow(remaining),
          estimate = if (nrow(remaining)) mean(remaining$value) else NA_real_)
      }))
    }))
  }))
  source_inference <- .ablation_bio_bind(lapply(split(cohort,
    interaction(cohort$comparison, cohort$anchor, drop = TRUE)), function(x) {
    x <- x[x$eligible_inference, ]
    query <- prepared$query_metadata
    lookup <- unique(query[, c("cohort_key", "source_system")])
    x$source <- lookup$source_system[match(x$cohort, lookup$cohort_key)]
    if (!nrow(x) || any(!.ablation_bio_known(x$source))) return(NULL)
    groups <- split(seq_len(nrow(x)), x$source)
    set.seed(.ablation_bio_seed(config, paste(x$comparison[1L], x$anchor[1L], "source")))
    boot <- replicate(config$bootstrap, {
      selected <- sample(names(groups), length(groups), TRUE)
      indices <- unlist(groups[selected], use.names = FALSE)
      mean(x$value[indices])
    })
    ci <- if (length(groups) >= 3L) unname(stats::quantile(boot, c(.025, .975))) else c(NA_real_, NA_real_)
    data.frame(comparison = x$comparison[1L], anchor = x$anchor[1L], estimate = mean(x$value),
      source_count = length(groups), cohort_count = nrow(x), ci_low = ci[1L], ci_high = ci[2L],
      status = if (length(groups) >= 3L) "estimable" else "not_estimable",
      reason = if (length(groups) >= 3L) "source_cluster_bootstrap_conditional" else "fewer_than_three_sources")
  }))
  distributions <- .ablation_bio_bind(lapply(names(scores$values), function(definition) {
    x <- scores$values[[definition]]
    metadata <- rbind(prepared$reference_metadata, prepared$query_metadata)
    .ablation_bio_bind(lapply(colnames(x), function(anchor) {
      .ablation_bio_bind(lapply(split(seq_len(nrow(metadata)), metadata$cohort_key), function(ids) {
        values <- x[match(metadata$sample_id[ids], rownames(x)), anchor]
        finite <- values[is.finite(values)]
        data.frame(cohort = metadata$cohort_key[ids[1L]], anchor = anchor,
          score_definition = definition, sample_count = length(values), valid_count = length(finite),
          mean = if (length(finite)) mean(finite) else NA_real_,
          sd = if (length(finite) > 1L) stats::sd(finite) else NA_real_,
          min = if (length(finite)) min(finite) else NA_real_,
          max = if (length(finite)) max(finite) else NA_real_)
      }))
    }))
  }))
  validation_error <- .ablation_bio_bind(lapply(retrieval, function(result) {
    .ablation_bio_bind(lapply(result$validation_lists, function(v) {
      if (is.null(v)) return(NULL)
      .ablation_bio_bind(lapply(names(scores$values), function(definition) {
        x <- scores$values[[definition]]
        .ablation_bio_bind(lapply(config$primary_anchors, function(anchor) {
          q <- x[match(v$query_ids, rownames(x)), anchor]
          error <- function(ids) rowMeans(exp(-abs(matrix(
            x[match(as.vector(ids), rownames(x)), anchor], nrow(ids), ncol(ids)) - q)))
          delta <- error(v$approximate_ids) - error(v$exact_ids)
          data.frame(arm = result$neighbors$arm[1L], pool = result$neighbors$pool[1L],
            distance = result$neighbors$distance[1L], restriction = result$neighbors$restriction[1L],
            cohort = prepared$query_metadata$cohort_key[match(v$query_ids, prepared$query_metadata$sample_id)],
            score_definition = definition, anchor = anchor, utility_error = delta)
        }))
      }))
    }))
  }))
  validation_error <- if (nrow(validation_error)) stats::aggregate(
    list(utility_error = validation_error$utility_error),
    validation_error[, c("arm", "pool", "distance", "restriction", "cohort",
      "score_definition", "anchor"), drop = FALSE],
    function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_) else data.frame()
  neighbor_summary <- .ablation_bio_bind(lapply(split(neighbors,
    interaction(neighbors$pool, neighbors$distance, neighbors$restriction, drop = TRUE)), function(n) {
    d <- split(n$reference_sample[n$arm == "Direct"], n$query_sample[n$arm == "Direct"])
    b <- split(n$reference_sample[n$arm == "d1"], n$query_sample[n$arm == "d1"])
    shared <- intersect(names(d), names(b))
    if (!length(shared)) return(NULL)
    overlap <- vapply(shared, function(id) length(intersect(d[[id]], b[[id]])) /
      length(union(d[[id]], b[[id]])), numeric(1))
    .ablation_bio_bind(lapply(split(shared, prepared$query_metadata$cohort_key[
      match(shared, prepared$query_metadata$sample_id)]), function(ids) {
      part <- n[n$query_sample %in% ids, ]
      data.frame(pool = n$pool[1L], distance = n$distance[1L], restriction = n$restriction[1L],
        cohort = part$query_cohort[1L], query_count = length(ids),
        mean_neighbor_jaccard = mean(overlap[match(ids, shared)]),
        reference_cohort_count = length(unique(part$reference_cohort)))
    }))
  }))
  neighbor_sources <- if (nrow(neighbors)) stats::aggregate(reference_sample ~
    arm + pool + distance + restriction + query_cohort + reference_cohort,
    neighbors, length) else data.frame()
  if (nrow(neighbor_sources)) names(neighbor_sources)[ncol(neighbor_sources)] <- "neighbor_pair_count"
  floor <- if (nrow(per_query)) stats::aggregate(utility ~ score_definition + anchor +
    pool + distance + restriction + arm + query_cohort, per_query,
    function(x) mean(x < 1e-8)) else data.frame()
  if (nrow(floor)) names(floor)[ncol(floor)] <- "utility_floor_fraction"
  arm_summary <- if (nrow(per_query)) stats::aggregate(
    per_query[, c("utility", "abs_delta", "cancer_match"), drop = FALSE],
    per_query[, c("score_definition", "anchor", "pool", "distance", "restriction",
      "arm", "query_cohort"), drop = FALSE],
    function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_) else data.frame()
  branch_status <- .ablation_bio_bind(lapply(seq_along(retrieval), function(i) {
    result <- retrieval[[i]]
    spec <- result$spec
    data.frame(domain = "retrieval", experiment = paste(spec$pool, spec$arm,
      spec$distance, if (is.null(spec$restriction)) "none" else spec$restriction, sep = ":"),
      status = result$status, reason = result$reason)
  }))
  if (!nrow(branch_status)) branch_status <- data.frame(domain = character(), experiment = character(),
    status = character(), reason = character())
  score_status <- scores$scales
  branch_status <- rbind(branch_status, data.frame(domain = "score",
    experiment = paste(score_status$score_definition, score_status$anchor, sep = ":"),
    status = score_status$status, reason = score_status$reason))
  readout_status <- .ablation_bio_bind(lapply(config$score_definitions, function(definition) {
    .ablation_bio_bind(lapply(config$primary_anchors, function(anchor) {
      .ablation_bio_bind(lapply(c("ridge", "xgboost"), function(reader) {
        x <- readout$predictions
        valid <- nrow(x) > 0L && any(x$score_definition == definition & x$anchor == anchor & x$reader == reader)
        data.frame(domain = "readout", experiment = paste(definition, anchor, reader, sep = ":"),
          status = if (valid) "estimable" else if (reader == "xgboost") "not_run" else "not_estimable",
          reason = if (valid) "supported" else if (reader == "xgboost")
            "reference_cv_trigger_not_met_or_target_unavailable" else "target_or_reference_folds_unavailable")
      }))
    }))
  }))
  branch_status <- rbind(branch_status, readout_status)
  # Explicitly retain every external cohort for every declared comparison.
  qualification <- .ablation_bio_bind(lapply(seq_len(nrow(inference)), function(i) {
    row <- inference[i, ]
    part <- lists[[row$comparison]]
    .ablation_bio_bind(lapply(unique(prepared$query_metadata$cohort_key), function(key) {
      total <- sum(prepared$query_metadata$cohort_key == key)
      count <- if (is.null(part)) 0L else length(unique(part$sample_id[
        part$cohort == key & part$anchor == row$anchor]))
      metric_reason <- "supported"
      if (grepl("^readout:", row$comparison) && nrow(metrics)) {
        reader_measure <- strsplit(row$experiment, ":", fixed = TRUE)[[1L]]
        m <- metrics[metrics$cohort == key & metrics$anchor == row$anchor &
          metrics$score_definition == row$score_definition & metrics$reader == reader_measure[1L], ]
        if (nrow(m) == 2L) {
          count <- min(m$query_count)
          if (!all(is.finite(m[[reader_measure[2L]]]))) metric_reason <-
            if (reader_measure[2L] == "correlation") "truth_or_prediction_zero_variance" else "nonfinite_mae"
        }
      }
      effective <- if (metric_reason == "supported") count else 0L
      data.frame(comparison = row$comparison, score_definition = row$score_definition,
        anchor = row$anchor, cohort = key, total_queries = total,
        paired_queries = count, effective_queries = effective,
        eligible_inference = effective >= config$min_query_per_cohort,
        reason = if (metric_reason != "supported") metric_reason else
          if (count >= config$min_query_per_cohort) "supported" else if (count > 0L)
          "fewer_than_twenty_queries_descriptive_only" else "no_complete_paired_queries")
    }))
  }))
  list(inference = inference, cohort = cohort, local = .ablation_bio_bind(locals),
    paired_queries = paired, common_queries = lists, bootstrap_indices = draws,
    per_query = per_query, prediction_metrics = metrics, leave_one = leave_one,
    source_inference = source_inference, score_distributions = distributions,
    retrieval_validation = .ablation_bio_bind(lapply(retrieval, `[[`, "validation")),
    retrieval_utility_error = validation_error, neighbor_summary = neighbor_summary,
    neighbor_sources = neighbor_sources, utility_floor = floor, arm_summary = arm_summary,
    branch_status = branch_status, experiment_qualification = qualification,
    module_diagnostics = .ablation_bio_module_diagnostics(prepared, retrieval),
    exploratory = if (isTRUE(config$low_coverage_exploration))
      .ablation_bio_exploratory(prepared, audit, scores$exploratory, retrieval, config) else NULL)
}
