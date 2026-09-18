# Helpers for independent-cohort structural reproducibility.
# The biological entities are anchor-specific low/high states defined within
# each cohort; representation geometry is compared only after entity alignment.

.asr_resolve_module_table <- function(
    module_manifest,
    resolution_audit,
    external_cohorts) {
  modules <- module_manifest$modules
  modules$bank_tissue <- modules$tissue
  mapped_tissue <- resolution_audit$resolved_tissue[
    match(modules$cohort, resolution_audit$cohort)
  ]
  modules$tissue <- ifelse(
    modules$bank_tissue == "Undefined" & !is.na(mapped_tissue),
    mapped_tissue,
    modules$bank_tissue
  )
  unresolved <- is.na(modules$tissue) | !nzchar(modules$tissue) |
    modules$tissue == "Undefined"
  if (any(unresolved)) {
    stop(
      "structural reproducibility: unresolved module tissue: ",
      paste(modules$module_id[unresolved], collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  modules$cohort_key <- paste(modules$tissue, modules$cohort, sep = "/")
  modules$bank_role <- ifelse(
    modules$cohort_key %in% external_cohorts,
    "external",
    "reference"
  )
  modules
}

.asr_extract_anchor_cache <- function(
    data,
    anchors,
    sample_ids,
    cohort_key_lookup = NULL) {
  required_genes <- sort(unique(unlist(anchors, use.names = FALSE)))
  sample_ids <- unique(as.character(sample_ids))
  cohorts <- list()
  for (tissue in names(data)) {
    for (cohort in names(data[[tissue]])) {
      expression <- data[[tissue]][[cohort]]
      if (is.list(expression) && !is.matrix(expression) &&
          !is.data.frame(expression)) {
        expression <- expression$expr
      }
      ids <- colnames(expression)
      genes <- rownames(expression)
      if (is.null(ids) || is.null(genes)) next
      keep <- ids %in% sample_ids
      gene_index <- match(required_genes, genes)
      if (!any(keep) || all(is.na(gene_index))) next
      expression <- as.matrix(expression[
        gene_index[!is.na(gene_index)],
        keep,
        drop = FALSE
      ])
      rownames(expression) <- required_genes[!is.na(gene_index)]
      cohort_key <- if (!is.null(cohort_key_lookup) &&
          cohort %in% names(cohort_key_lookup)) {
        unname(cohort_key_lookup[[cohort]])
      } else {
        paste(tissue, cohort, sep = "/")
      }
        if (cohort_key %in% names(cohorts)) {
          stop(
            "structural reproducibility: duplicated cohort_key: ",
            cohort_key,
            call. = FALSE
          )
        }
        cohorts[[cohort_key]] <- list(
        tissue = tissue,
        cohort = cohort,
        cohort_key = cohort_key,
        expression = expression,
        sample_id = ids[keep]
      )
    }
  }
  if (!length(cohorts)) {
    stop("structural reproducibility: no anchor samples were extracted.",
      call. = FALSE
    )
  }
  list(
    anchors = anchors,
    required_genes = required_genes,
    sample_ids = sort(unique(unlist(lapply(
      cohorts,
      function(cohort) cohort$sample_id
    ), use.names = FALSE))),
    cohorts = cohorts
  )
}

.asr_subset_bank <- function(fit_d1, target_d1, blocks, module_ids) {
  module_ids <- names(blocks)[names(blocks) %in% module_ids]
  if (!length(module_ids)) {
    stop("structural reproducibility: selected bank has no modules.",
      call. = FALSE
    )
  }
  selected_blocks <- blocks[module_ids]
  columns <- unlist(selected_blocks, use.names = FALSE)
  if (any(is.na(columns)) || max(columns) > ncol(fit_d1) ||
      max(columns) > ncol(target_d1)) {
    stop("structural reproducibility: module block is outside d1 columns.",
      call. = FALSE
    )
  }
  widths <- lengths(selected_blocks)
  block_ends <- cumsum(widths)
  block_starts <- c(1L, utils::head(block_ends, -1L) + 1L)
  local_blocks <- Map(seq.int, block_starts, block_ends)
  names(local_blocks) <- module_ids
  list(
    fit = fit_d1[, columns, drop = FALSE],
    target = target_d1[, columns, drop = FALSE],
    blocks = local_blocks,
    module_ids = module_ids,
    feature_count = length(columns)
  )
}

.asr_matched_bank_design <- function(module_table, n_repeats, seed) {
  reference <- module_table[module_table$bank_role == "reference", , drop = FALSE]
  external <- module_table[module_table$bank_role == "external", , drop = FALSE]
  shared_tissues <- intersect(unique(reference$tissue), unique(external$tissue))
  shared_tissues <- sort(shared_tissues)
  if (!length(shared_tissues)) {
    stop("structural reproducibility: banks have no shared tissue.",
      call. = FALSE
    )
  }
  rows <- lapply(seq_len(as.integer(n_repeats)), function(repeat_id) {
    set.seed(seed + repeat_id - 1L)
    tissue_rows <- lapply(shared_tissues, function(tissue) {
      reference_ids <- reference$module_id[reference$tissue == tissue]
      external_ids <- external$module_id[external$tissue == tissue]
      matched_n <- min(length(reference_ids), length(external_ids))
      if (!matched_n) return(NULL)
      rbind(
        data.frame(
            repeat_id = repeat_id,
            seed = as.integer(seed + repeat_id - 1L),
            n_repeats = as.integer(n_repeats),
          bank_role = "reference",
          tissue = tissue,
          module_id = sample(reference_ids, matched_n),
          stringsAsFactors = FALSE
        ),
        data.frame(
            repeat_id = repeat_id,
            seed = as.integer(seed + repeat_id - 1L),
            n_repeats = as.integer(n_repeats),
          bank_role = "external",
          tissue = tissue,
          module_id = sample(external_ids, matched_n),
          stringsAsFactors = FALSE
        )
      )
    })
    do.call(rbind, tissue_rows)
  })
  design <- do.call(rbind, rows)
  rownames(design) <- NULL
  design
}

.asr_prepare_direction_base <- function(
    fit_direct,
    target_direct,
    fit_metadata,
    target_metadata,
    anchor_cache,
    fit_sample_ids,
    target_sample_ids,
    tail_fraction,
    min_entity_n) {
  fit_sample_ids <- Reduce(intersect, list(
    rownames(fit_direct),
    fit_metadata$sample_id,
    anchor_cache$sample_ids,
    fit_sample_ids
  ))
  target_sample_ids <- Reduce(intersect, list(
    rownames(target_direct),
    target_metadata$sample_id,
    anchor_cache$sample_ids,
    target_sample_ids
  ))
  if (length(fit_sample_ids) < 2L || length(target_sample_ids) < 2L) {
    stop("structural reproducibility: too few common fit or target samples.",
      call. = FALSE
    )
  }
  fit_metadata <- fit_metadata[
    match(fit_sample_ids, fit_metadata$sample_id),
    ,
    drop = FALSE
  ]
  target_metadata <- target_metadata[
    match(target_sample_ids, target_metadata$sample_id),
    ,
    drop = FALSE
  ]
  if (!"cohort_key" %in% colnames(fit_metadata)) {
    fit_metadata$cohort_key <- paste(
      fit_metadata$tissue,
      fit_metadata$cohort,
      sep = "/"
    )
  }
  if (!"cohort_key" %in% colnames(target_metadata)) {
    target_metadata$cohort_key <- paste(
      target_metadata$tissue,
      target_metadata$cohort,
      sep = "/"
    )
  }
  if (!"cancer_type" %in% colnames(target_metadata)) {
    target_metadata$cancer_type <- as.character(target_metadata$biology)
  }
  direct_scaled <- .ablation_scale_train_apply(
    fit_direct[fit_sample_ids, , drop = FALSE],
    target_direct[target_sample_ids, , drop = FALSE]
  )
  scores <- .asr_compute_anchor_scores(
    anchor_cache,
    scaling_cohorts = unique(fit_metadata$cohort_key),
    scoring_cohorts = unique(target_metadata$cohort_key),
    scaling_sample_ids = fit_sample_ids,
    scoring_sample_ids = target_sample_ids
  )
  states <- .asr_assign_anchor_states(
    scores,
    tail_fraction = tail_fraction,
    min_entity_n = min_entity_n
  )
  direct_geometry <- .asr_build_cohort_geometries(
    direct_scaled$test,
    target_metadata,
    states
  )
  list(
    fit_sample_ids = fit_sample_ids,
    target_sample_ids = target_sample_ids,
    fit_metadata = fit_metadata,
    target_metadata = target_metadata,
    scores = scores,
    states = states,
    direct_geometry = direct_geometry
  )
}

.asr_evaluate_direction <- function(
    base,
    fit_d1,
    target_d1,
    blocks,
    direction,
    bank_role,
    target_role,
    min_shared_entities,
    n_boot,
    seed) {
  fit_d1 <- fit_d1[base$fit_sample_ids, , drop = FALSE]
  target_d1 <- target_d1[base$target_sample_ids, , drop = FALSE]
  d1_scaled <- .ablation_module_balanced_transform(
    fit_d1,
    target_d1,
    blocks
  )
  d1_geometry <- .asr_build_cohort_geometries(
    d1_scaled$query,
    base$target_metadata,
    base$states
  )
  cohort_metadata <- unique(base$target_metadata[, c(
    "cohort_key", "cancer_type", "assay_type", "source_system"
  )])
  pair_comparisons <- .asr_compare_cohort_pairs(
    base$direct_geometry$geometries,
    d1_geometry$geometries,
    cohort_metadata,
    min_shared_entities = min_shared_entities
  )
  pair_comparisons$direction <- direction
  pair_comparisons$bank_role <- bank_role
  pair_comparisons$target_role <- target_role
  summary <- .asr_summarize_pairs(pair_comparisons, n_boot, seed)
  summary$direction <- direction
  summary$bank_role <- bank_role
  summary$target_role <- target_role
  summary$bank_module_count <- length(blocks)
  summary$bank_feature_count <- ncol(fit_d1)
  summary$fit_sample_count <- nrow(fit_d1)
  summary$target_sample_count <- nrow(target_d1)
  cohorts <- sort(unique(c(
    pair_comparisons$cohort_a,
    pair_comparisons$cohort_b
  )))
  matrices <- list(
    `Direct-GSClassifier` = .asr_similarity_matrix(
      pair_comparisons,
      "direct_similarity",
      cohorts
    ),
    `Cohort-d1` = .asr_similarity_matrix(
      pair_comparisons,
      "d1_similarity",
      cohorts
    )
  )
  list(
    pair_comparisons = pair_comparisons,
    summary = summary,
    similarity_matrices = matrices,
    direct_prototypes = base$direct_geometry$prototypes,
    d1_prototypes = d1_geometry$prototypes
  )
}

.asr_compute_anchor_scores <- function(
    cache,
    scaling_cohorts,
    scoring_cohorts = names(cache$cohorts),
    scaling_sample_ids = NULL,
    scoring_sample_ids = NULL) {
  scaling_keys <- intersect(names(cache$cohorts), scaling_cohorts)
  scoring_keys <- intersect(names(cache$cohorts), scoring_cohorts)
  if (!length(scaling_keys)) {
    stop("structural reproducibility: no cohorts for anchor scaling.",
      call. = FALSE
    )
  }
  if (!length(scoring_keys)) {
    stop("structural reproducibility: no cohorts for anchor scoring.",
      call. = FALSE
    )
  }

  genes <- sort(unique(unlist(cache$anchors, use.names = FALSE)))
  global_stats <- lapply(genes, function(gene) {
    values <- unlist(lapply(cache$cohorts[scaling_keys], function(cohort) {
      if (!gene %in% rownames(cohort$expression)) return(numeric())
      ids <- cohort$sample_id
      keep <- if (is.null(scaling_sample_ids)) {
        rep(TRUE, length(ids))
      } else {
        ids %in% scaling_sample_ids
      }
      as.numeric(cohort$expression[gene, keep, drop = TRUE])
    }), use.names = FALSE)
    values <- values[is.finite(values)]
    if (length(values) < 2L) {
      return(c(mean = NA_real_, sd = NA_real_, n = length(values)))
    }
    c(mean = mean(values), sd = stats::sd(values), n = length(values))
  })
  global_stats <- as.data.frame(do.call(rbind, global_stats))
  global_stats$gene_id <- genes
  global_stats <- global_stats[
    is.finite(global_stats$mean) & is.finite(global_stats$sd) &
      global_stats$sd > 0,
    ,
    drop = FALSE
  ]

  score_rows <- unlist(lapply(cache$cohorts[scoring_keys], function(cohort) {
    ids <- cohort$sample_id
    keep <- if (is.null(scoring_sample_ids)) {
      rep(TRUE, length(ids))
    } else {
      ids %in% scoring_sample_ids
    }
    if (!any(keep)) return(NULL)
    lapply(names(cache$anchors), function(anchor) {
      anchor_genes <- intersect(
        cache$anchors[[anchor]],
        intersect(rownames(cohort$expression), global_stats$gene_id)
      )
      if (length(anchor_genes) < 2L) return(NULL)
      scaling <- global_stats[
        match(anchor_genes, global_stats$gene_id),
        ,
        drop = FALSE
      ]
      values <- cohort$expression[anchor_genes, keep, drop = FALSE]
      values <- sweep(values, 1L, scaling$mean, FUN = "-")
      values <- sweep(values, 1L, scaling$sd, FUN = "/")
      data.frame(
        sample_id = ids[keep],
        cohort_key = cohort$cohort_key,
        anchor = anchor,
        score = as.numeric(colMeans(values, na.rm = TRUE)),
        gene_count = length(anchor_genes),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
  if (!length(score_rows)) {
    stop("structural reproducibility: no finite anchor score rows.",
      call. = FALSE
    )
  }
  scores <- do.call(rbind, score_rows)
  scores <- scores[is.finite(scores$score), , drop = FALSE]
  rownames(scores) <- NULL
  scores
}

.asr_assign_anchor_states <- function(
    scores,
    tail_fraction = 1 / 3,
    min_entity_n = 8L) {
  if (!is.finite(tail_fraction) || tail_fraction <= 0 ||
      tail_fraction >= 0.5) {
    stop("structural reproducibility: tail_fraction must be in (0, 0.5).",
      call. = FALSE
    )
  }
  groups <- split(scores, interaction(
    scores$cohort_key,
    scores$anchor,
    drop = TRUE,
    lex.order = TRUE
  ))
  state_rows <- lapply(groups, function(part) {
    part <- part[order(part$score, part$sample_id), , drop = FALSE]
    tail_n <- floor(nrow(part) * tail_fraction)
    if (tail_n < min_entity_n) return(NULL)
    # Ties at a tail boundary must not be split by arbitrary sample IDs.
    if (part$score[tail_n] >= part$score[tail_n + 1L] ||
        part$score[nrow(part) - tail_n] >= part$score[nrow(part) - tail_n + 1L]) {
      return(NULL)
    }
    selected <- rbind(
      transform(part[seq_len(tail_n), , drop = FALSE], state = "low"),
      transform(
        part[nrow(part) - rev(seq_len(tail_n)) + 1L, , drop = FALSE],
        state = "high"
      )
    )
    selected$entity <- paste(selected$anchor, selected$state, sep = ":")
    selected$cohort_sample_count <- nrow(part)
    selected$entity_sample_count <- tail_n
    selected
  })
  states <- do.call(rbind, state_rows)
  if (is.null(states) || !nrow(states)) {
    stop("structural reproducibility: no biological state met min_entity_n.",
      call. = FALSE
    )
  }
  rownames(states) <- NULL
  states
}

.asr_build_cohort_geometries <- function(representation, metadata, states) {
  sample_ids <- Reduce(intersect, list(
    rownames(representation),
    metadata$sample_id,
    states$sample_id
  ))
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
  representation <- representation[sample_ids, , drop = FALSE]
  states <- states[states$sample_id %in% sample_ids, , drop = FALSE]
  groups <- split(states, interaction(
    states$cohort_key,
    states$entity,
    drop = TRUE,
    lex.order = TRUE
  ))
  prototypes <- lapply(groups, function(part) {
    ids <- intersect(part$sample_id, rownames(representation))
    data.frame(
      cohort_key = part$cohort_key[1L],
      entity = part$entity[1L],
      sample_count = length(ids),
      centroid = I(list(colMeans(representation[ids, , drop = FALSE]))),
      stringsAsFactors = FALSE
    )
  })
  prototypes <- do.call(rbind, prototypes)
  cohort_prototypes <- split(prototypes, prototypes$cohort_key, drop = TRUE)
  geometries <- lapply(cohort_prototypes, function(part) {
    centroid <- do.call(rbind, part$centroid)
    rownames(centroid) <- part$entity
    as.matrix(stats::dist(centroid, method = "euclidean"))
  })
  list(geometries = geometries, prototypes = prototypes)
}

.asr_compare_cohort_pairs <- function(
    direct_geometries,
    d1_geometries,
    cohort_metadata,
    min_shared_entities = 4L) {
  cohorts <- Reduce(intersect, list(
    names(direct_geometries),
    names(d1_geometries),
    cohort_metadata$cohort_key
  ))
  if (length(cohorts) < 2L) {
    stop("structural reproducibility: fewer than two target cohorts are estimable.",
      call. = FALSE
    )
  }
  pairs <- utils::combn(sort(cohorts), 2L, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    entities <- Reduce(intersect, list(
      rownames(direct_geometries[[pair[1L]]]),
      rownames(direct_geometries[[pair[2L]]]),
      rownames(d1_geometries[[pair[1L]]]),
      rownames(d1_geometries[[pair[2L]]])
    ))
    if (length(entities) < min_shared_entities) return(NULL)
    upper <- upper.tri(matrix(0, length(entities), length(entities)))
    direct_a <- direct_geometries[[pair[1L]]][entities, entities][upper]
    direct_b <- direct_geometries[[pair[2L]]][entities, entities][upper]
    d1_a <- d1_geometries[[pair[1L]]][entities, entities][upper]
    d1_b <- d1_geometries[[pair[2L]]][entities, entities][upper]
    direct_similarity <- stats::cor(
      direct_a,
      direct_b,
      method = "spearman"
    )
    d1_similarity <- stats::cor(d1_a, d1_b, method = "spearman")
    cancer_a <- cohort_metadata$cancer_type[
      match(pair[1L], cohort_metadata$cohort_key)
    ]
    cancer_b <- cohort_metadata$cancer_type[
      match(pair[2L], cohort_metadata$cohort_key)
    ]
    data.frame(
      cohort_a = pair[1L],
      cohort_b = pair[2L],
      cancer_type_a = cancer_a,
      cancer_type_b = cancer_b,
      same_cancer_type = identical(cancer_a, cancer_b),
      shared_entity_count = length(entities),
      compared_distance_count = sum(upper),
      direct_similarity = direct_similarity,
      d1_similarity = d1_similarity,
      delta_d1_minus_direct = d1_similarity - direct_similarity,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) {
    stop("structural reproducibility: no cohort pair shares enough entities.",
      call. = FALSE
    )
  }
  result <- do.call(rbind, rows)
  result <- result[
    is.finite(result$direct_similarity) & is.finite(result$d1_similarity),
    ,
    drop = FALSE
  ]
  rownames(result) <- NULL
  result
}

.asr_similarity_matrix <- function(pair_data, value_column, cohorts) {
  result <- matrix(
    NA_real_,
    nrow = length(cohorts),
    ncol = length(cohorts),
    dimnames = list(cohorts, cohorts)
  )
  diag(result) <- 1
  for (i in seq_len(nrow(pair_data))) {
    result[pair_data$cohort_a[i], pair_data$cohort_b[i]] <-
      pair_data[[value_column]][i]
    result[pair_data$cohort_b[i], pair_data$cohort_a[i]] <-
      pair_data[[value_column]][i]
  }
  result
}

.asr_node_bootstrap <- function(pair_data, n_boot, seed) {
  cohorts <- sort(unique(c(pair_data$cohort_a, pair_data$cohort_b)))
  if (n_boot < 1L) {
    return(c(
      ci_low = NA_real_,
      ci_high = NA_real_,
      valid_bootstrap = 0
    ))
  }
  set.seed(seed)
  estimates <- replicate(n_boot, {
    sampled <- sample(cohorts, length(cohorts), replace = TRUE)
    weights <- table(sampled)
    pair_weights <- as.numeric(weights[pair_data$cohort_a]) *
      as.numeric(weights[pair_data$cohort_b])
    keep <- is.finite(pair_weights) & pair_weights > 0
    if (!any(keep)) return(NA_real_)
    stats::weighted.mean(
      pair_data$delta_d1_minus_direct[keep],
      pair_weights[keep]
    )
  })
  estimates <- estimates[is.finite(estimates)]
  c(
    ci_low = unname(stats::quantile(estimates, 0.025)),
    ci_high = unname(stats::quantile(estimates, 0.975)),
    valid_bootstrap = length(estimates)
  )
}

.asr_summarize_pairs <- function(
    pair_data,
    n_boot,
    seed,
    min_inference_cohorts = 4L,
    min_inference_pairs = 3L) {
  scopes <- list(
    all_cohort_pairs = pair_data,
    same_cancer_type = pair_data[pair_data$same_cancer_type, , drop = FALSE]
  )
  rows <- lapply(seq_along(scopes), function(i) {
    part <- scopes[[i]]
    if (!nrow(part)) return(NULL)
    cohort_count <- length(unique(c(part$cohort_a, part$cohort_b)))
    inference_status <- if (
      cohort_count >= min_inference_cohorts &&
        nrow(part) >= min_inference_pairs && n_boot > 0L
    ) "estimable" else "not_estimable"
    interval <- if (identical(inference_status, "estimable")) {
      .asr_node_bootstrap(part, n_boot, seed + i - 1L)
    } else {
      c(ci_low = NA_real_, ci_high = NA_real_, valid_bootstrap = 0)
    }
    data.frame(
      scope = names(scopes)[i],
      cohort_count = cohort_count,
      cohort_pair_count = nrow(part),
      mean_direct_similarity = mean(part$direct_similarity),
      median_direct_similarity = stats::median(part$direct_similarity),
      mean_d1_similarity = mean(part$d1_similarity),
      median_d1_similarity = stats::median(part$d1_similarity),
      mean_delta = mean(part$delta_d1_minus_direct),
      median_delta = stats::median(part$delta_d1_minus_direct),
      fraction_d1_ge_direct = mean(part$delta_d1_minus_direct >= 0),
      mean_delta_ci_low = interval[["ci_low"]],
      mean_delta_ci_high = interval[["ci_high"]],
      valid_bootstrap = as.integer(interval[["valid_bootstrap"]]),
      inference_method = ifelse(
        inference_status == "estimable",
        "cohort_node_bootstrap",
        "not_estimable"
      ),
      inference_status = inference_status,
      min_inference_cohorts = as.integer(min_inference_cohorts),
      min_inference_pairs = as.integer(min_inference_pairs),
      p_value = NA_real_,
      p_value_reason = paste(
        "Cohort-pair rows share cohort nodes; a paired Wilcoxon test would",
        "treat dependent dyads as independent."
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
