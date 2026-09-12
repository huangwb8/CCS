# Helpers for independent-cohort structural reproducibility.
# The biological entities are anchor-specific low/high states defined within
# each cohort; representation geometry is compared only after entity alignment.

.asr_compute_anchor_scores <- function(cache, external_cohorts) {
  reference_keys <- setdiff(names(cache$cohorts), external_cohorts)
  if (!length(reference_keys)) {
    stop("structural reproducibility: no reference cohorts for anchor scaling.",
      call. = FALSE
    )
  }

  genes <- sort(unique(unlist(cache$anchors, use.names = FALSE)))
  global_stats <- lapply(genes, function(gene) {
    values <- unlist(lapply(cache$cohorts[reference_keys], function(cohort) {
      if (!gene %in% rownames(cohort$expression)) return(numeric())
      as.numeric(cohort$expression[gene, , drop = TRUE])
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

  score_rows <- unlist(lapply(cache$cohorts, function(cohort) {
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
      values <- cohort$expression[anchor_genes, , drop = FALSE]
      values <- sweep(values, 1L, scaling$mean, FUN = "-")
      values <- sweep(values, 1L, scaling$sd, FUN = "/")
      data.frame(
        sample_id = cohort$sample_id,
        cohort_key = cohort$cohort_key,
        anchor = anchor,
        score = as.numeric(colMeans(values, na.rm = TRUE)),
        gene_count = length(anchor_genes),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
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

.asr_summarize_pairs <- function(pair_data, n_boot, seed) {
  scopes <- list(
    all_cohort_pairs = pair_data,
    same_cancer_type = pair_data[pair_data$same_cancer_type, , drop = FALSE]
  )
  rows <- lapply(seq_along(scopes), function(i) {
    part <- scopes[[i]]
    if (!nrow(part)) return(NULL)
    interval <- .asr_node_bootstrap(part, n_boot, seed + i - 1L)
    data.frame(
      scope = names(scopes)[i],
      cohort_count = length(unique(c(part$cohort_a, part$cohort_b))),
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
      inference_method = "cohort_node_bootstrap",
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
