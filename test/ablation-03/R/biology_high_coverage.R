# Analysis-specific high-coverage validation. Selection accepts measurement
# metadata and availability only; the installed CCS supplies established math.

.hc_bio <- function(name) getFromNamespace(paste0(".ablation_bio_", name), "CCS")
.hc_bind <- function(x) .hc_bio("bind")(x)
.hc_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

.hc_inventory <- function(atlas, anchors, reference_metadata, query_metadata) {
  fields <- c("sample_id", "cohort_key", "cancer_type", "metadata_status",
    "evidence_basis", "source_system", "assay_type", "platform_id")
  ref <- reference_metadata[, fields, drop = FALSE]
  query <- query_metadata[, c(fields, "d1_provenance"), drop = FALSE]
  ref$role <- "reference"
  ref$d1_provenance <- "reference"
  query$role <- "query"
  metadata <- rbind(ref, query[, names(ref), drop = FALSE])
  if (anyDuplicated(metadata$sample_id) ||
      length(intersect(ref$cohort_key, query$cohort_key))) {
    stop("High-coverage reference/query identities overlap.", call. = FALSE)
  }
  genes <- sort(unique(unlist(anchors, use.names = FALSE)))
  cohorts <- list()
  for (tissue in names(atlas)) for (name in names(atlas[[tissue]])) {
    key <- paste(tissue, name, sep = "/")
    m <- metadata[metadata$cohort_key == key, , drop = FALSE]
    if (!nrow(m)) next
    cohort_fields <- c("role", "source_system")
    if (any(vapply(m[, cohort_fields, drop = FALSE], function(x) length(unique(x)) != 1L, logical(1))))
      stop("Conflicting measurement labels within cohort: ", key, call. = FALSE)
    x <- atlas[[tissue]][[name]]
    if (is.list(x) && !is.data.frame(x)) x <- x$expr
    if (is.null(rownames(x)) || anyDuplicated(rownames(x)) ||
        anyDuplicated(colnames(x)) || !all(m$sample_id %in% colnames(x))) {
      stop("Ambiguous high-coverage expression IDs: ", key, call. = FALSE)
    }
    present <- intersect(genes, rownames(x))
    raw <- as.matrix(x[present, m$sample_id, drop = FALSE])
    storage.mode(raw) <- "double"
    cohorts[[key]] <- list(genes = sort(rownames(x)), raw = raw,
      sample_ids = m$sample_id, role = m$role[1L])
  }
  if (!setequal(names(cohorts), unique(metadata$cohort_key))) {
    stop("Expression atlas does not cover the measurement contract.", call. = FALSE)
  }
  list(metadata = metadata, cohorts = cohorts, anchors = anchors)
}

.hc_supported_samples <- function(inventory, genes, config) {
  m <- inventory$metadata
  known <- .hc_bio("known")
  ok <- known(m$cancer_type) & m$metadata_status == "confirmed" & known(m$evidence_basis) &
    (m$role == "reference" | m$d1_provenance == "external_frozen")
  ok[is.na(ok)] <- FALSE
  finite_ids <- unlist(lapply(inventory$cohorts, function(x) {
    if (!all(genes %in% rownames(x$raw))) return(character())
    x$sample_ids[colSums(!is.finite(x$raw[genes, , drop = FALSE])) == 0L]
  }), use.names = FALSE)
  m <- m[ok & m$sample_id %in% finite_ids, , drop = FALSE]
  ref <- m[m$role == "reference", , drop = FALSE]
  query <- m[m$role == "query", , drop = FALSE]
  reference_counts <- table(ref$cancer_type)
  reference_cohorts <- vapply(split(ref$cohort_key, ref$cancer_type), function(x)
    length(unique(x)), integer(1))
  supported_labels <- names(reference_counts)[reference_counts >= config$k &
    reference_cohorts[names(reference_counts)] >= config$min_reference_cohorts]
  query <- query[query$cancer_type %in% supported_labels, , drop = FALSE]
  sizes <- table(query$cohort_key)
  query <- query[query$cohort_key %in% names(sizes)[sizes >= config$min_query_per_cohort], ]
  ref <- ref[ref$cancer_type %in% query$cancer_type, , drop = FALSE]
  rbind(ref, query)
}

# Availability seeds are deterministic and never scored by biological effect.
# Each seed defines all compatible cohorts. This bounded search is not a proof
# of globally maximal breadth; all seed rules are frozen and reported.
.hc_contract <- function(inventory, anchor, tier, config) {
  signature <- sort(unique(inventory$anchors[[anchor]]))
  need <- max(config$min_signature_genes, ceiling(tier * length(signature)))
  available <- lapply(inventory$cohorts, function(x) intersect(signature, x$genes))
  eligible <- available[lengths(available) >= need]
  frequency <- table(factor(unlist(eligible), levels = signature))
  prefix <- function(g) g[order(-as.numeric(frequency[g]), g)][seq_len(need)]
  seeds <- if (length(eligible)) c(list(prefix(signature)), lapply(eligible, prefix)) else list()
  if (length(eligible) > 1L) {
    pairs <- utils::combn(seq_along(eligible), 2L)
    for (j in seq_len(ncol(pairs))) {
      g <- intersect(eligible[[pairs[1L, j]]], eligible[[pairs[2L, j]]])
      if (length(g) >= need) seeds[[length(seeds) + 1L]] <- prefix(g)
    }
  }
  seeds <- lapply(seeds, sort)
  seeds <- seeds[!duplicated(vapply(seeds, .hc_hash, character(1)))]
  best <- NULL
  best_key <- rep(-1, 5L)
  for (genes in seeds) {
    selected <- .hc_supported_samples(inventory, genes, config)
    if (!nrow(selected)) next
    keys <- sort(unique(selected$cohort_key))
    common <- sort(Reduce(intersect, available[keys]))
    # Freeze the actual intersection, then enforce complete samples for it.
    repeat {
      restricted <- inventory
      restricted$cohorts <- inventory$cohorts[keys]
      restricted$metadata <- inventory$metadata[inventory$metadata$cohort_key %in% keys, ]
      s <- .hc_supported_samples(restricted, common, config)
      new_keys <- sort(unique(s$cohort_key))
      if (!length(new_keys)) break
      new_common <- sort(Reduce(intersect, available[new_keys]))
      stable <- identical(keys, new_keys) && identical(common, new_common)
      keys <- new_keys
      common <- new_common
      if (stable) break
    }
    if (!length(new_keys)) next
    selected <- s
    q <- selected[selected$role == "query", , drop = FALSE]
    r <- selected[selected$role == "reference", , drop = FALSE]
    sources <- sort(unique(q$source_system[.hc_bio("known")(q$source_system)]))
    key <- c(length(unique(q$cohort_key)), length(sources), nrow(q), length(common),
      length(unique(r$cohort_key)))
    different <- which(key != best_key)
    if (!length(different) || key[different[1L]] < best_key[different[1L]]) next
    background <- sort(Reduce(intersect, lapply(inventory$cohorts[keys], `[[`, "genes")))
    best <- list(anchor = anchor, tier = tier, genes = common, background = background,
      metadata = selected, query_cohorts = sort(unique(q$cohort_key)),
      reference_cohorts = sort(unique(r$cohort_key)), source_groups = sources)
    best_key <- key
  }
  if (is.null(best)) best <- list(anchor = anchor, tier = tier, genes = character(),
    background = character(), metadata = inventory$metadata[FALSE, ],
    query_cohorts = character(), reference_cohorts = character(), source_groups = character())
  best$original_genes <- signature
  best$gene_hash <- .hc_hash(best$genes)
  best$cohort_hash <- .hc_hash(list(best$query_cohorts, best$reference_cohorts))
  best$contract_hash <- .hc_hash(best)
  best
}

.hc_frontier <- function(inventory, config) {
  # Explicit whitelist makes effect removal/permutation/addition byte invariant.
  inventory$metadata <- inventory$metadata[, c("sample_id", "cohort_key", "cancer_type",
    "metadata_status", "evidence_basis", "source_system", "assay_type", "platform_id",
    "role", "d1_provenance"), drop = FALSE]
  inventory$cohorts <- lapply(inventory$cohorts, function(x)
    x[c("genes", "raw", "sample_ids", "role")])
  if (!identical(config$selection_rule, "deterministic_availability_seed_search") ||
      !is.null(config$non_inferiority_margin)) {
    stop("Unrecognized selection rule or unsupported non-inferiority margin.", call. = FALSE)
  }
  contracts <- list()
  rows <- list()
  eligibility <- list()
  for (anchor in config$primary_anchors) for (tier in config$coverage_tiers) {
    contract <- .hc_contract(inventory, anchor, tier, config)
    id <- paste(anchor, as.integer(round(100 * tier)), sep = "-")
    contracts[[id]] <- contract
    m <- contract$metadata
    nq <- length(contract$query_cohorts)
    ns <- length(contract$source_groups)
    rows[[id]] <- data.frame(id = id, anchor = anchor, tier = tier,
      original_genes = length(contract$original_genes), common_genes = length(contract$genes),
      coverage = length(contract$genes) / length(contract$original_genes),
      query_cohorts = nq, query_count = sum(m$role == "query"),
      reference_cohorts = length(contract$reference_cohorts),
      reference_count = sum(m$role == "reference"), source_groups = ns,
      unknown_source_queries = sum(m$role == "query" & !.hc_bio("known")(m$source_system)),
      background_genes = length(contract$background),
      rank_eligible = length(contract$background) >= config$min_background_genes,
      breadth_eligible = nq >= config$primary_min_query_cohorts &&
        ns >= config$primary_min_source_groups,
      gene_hash = contract$gene_hash, cohort_hash = contract$cohort_hash,
      contract_hash = contract$contract_hash)
    eligibility[[id]] <- .hc_bind(lapply(names(inventory$cohorts), function(key) {
      x <- inventory$cohorts[[key]]
      m <- inventory$metadata[inventory$metadata$cohort_key == key, ]
      complete <- if (length(contract$genes) && all(contract$genes %in% rownames(x$raw)))
        sum(colSums(!is.finite(x$raw[contract$genes, , drop = FALSE])) == 0L) else 0L
      selected <- key %in% c(contract$query_cohorts, contract$reference_cohorts)
      known <- .hc_bio("known")
      provenance_ok <- known(m$cancer_type) & m$metadata_status == "confirmed" &
        known(m$evidence_basis) & (m$role == "reference" | m$d1_provenance == "external_frozen")
      provenance_ok[is.na(provenance_ok)] <- FALSE
      present <- length(intersect(contract$original_genes, x$genes))
      data.frame(id = id, anchor = anchor, tier = tier, cohort_key = key, role = x$role,
        cancer_type = paste(sort(unique(m$cancer_type)), collapse = " | "), source_system = m$source_system[1L],
        assay_type = paste(sort(unique(m$assay_type)), collapse = " | "),
        platform_id = paste(sort(unique(m$platform_id)), collapse = " | "),
        total_samples = nrow(m), original_signature_present = present,
        complete_signature_samples = complete,
        frozen_sample_count = sum(contract$metadata$cohort_key == key), selected = selected,
        reason = if (selected) "selected_frozen_measurement_contract" else
          if (!any(provenance_ok)) "unconfirmed_cancer_or_external_provenance" else
          if (present < ceiling(tier * length(contract$original_genes))) "signature_below_tier" else
          if (x$role == "query" && complete < config$min_query_per_cohort)
            "incomplete_signature_or_fewer_than_twenty_queries" else
              "not_in_compatible_supported_measurement_set")
    }))
  }
  frontier <- .hc_bind(rows)
  frontier$selected <- FALSE
  frontier$evidence_level <- "high_coverage_sensitivity"
  for (anchor in config$primary_anchors) {
    idx <- which(frontier$anchor == anchor)
    valid <- idx[frontier$breadth_eligible[idx]]
    level <- "primary_design"
    if (!length(valid)) {
      valid <- idx[frontier$query_cohorts[idx] >= config$limited_min_query_cohorts]
      level <- "limited_external_validation"
    }
    if (!length(valid)) {
      valid <- idx[frontier$query_cohorts[idx] >= config$min_inference_cohorts]
      level <- "high_coverage_sensitivity"
    }
    if (!length(valid)) valid <- idx
    chosen <- valid[which.max(frontier$tier[valid])]
    frontier$selected[chosen] <- TRUE
    frontier$evidence_level[chosen] <- level
  }
  list(frontier = frontier, contracts = contracts, eligibility = .hc_bind(eligibility), config = config,
    source_md5 = inventory$source_md5,
    measurement_hash = .hc_hash(list(frontier, lapply(contracts, function(x)
      x[c("genes", "background", "query_cohorts", "reference_cohorts", "metadata")]))))
}

.hc_subset <- function(prepared, reference_ids, query_ids) {
  p <- prepared
  for (role in c("reference", "query")) {
    ids <- if (role == "reference") reference_ids else query_ids
    m <- prepared[[paste0(role, "_metadata")]]
    idx <- match(ids, m$sample_id)
    if (anyNA(idx)) stop("Frozen measurement samples missing from representation.")
    p[[paste0(role, "_metadata")]] <- m[idx, , drop = FALSE]
    for (arm in c("direct", "d1"))
      p[[paste0(role, "_", arm)]] <- prepared[[paste0(role, "_", arm)]][idx, , drop = FALSE]
  }
  p
}

.hc_score_contract <- function(atlas, contract, config) {
  m <- contract$metadata
  ids <- m$sample_id
  genes <- contract$genes
  bg <- contract$background
  raw <- matrix(NA_real_, length(genes), length(ids), dimnames = list(genes, ids))
  ranked <- matrix(NA_real_, length(ids), 1L, dimnames = list(ids, contract$anchor))
  complete_bg <- stats::setNames(rep(FALSE, length(ids)), ids)
  for (key in unique(m$cohort_key)) {
    parts <- strsplit(key, "/", fixed = TRUE)[[1L]]
    x <- atlas[[parts[1L]]][[paste(parts[-1L], collapse = "/")]]
    if (is.list(x) && !is.data.frame(x)) x <- x$expr
    use <- m$sample_id[m$cohort_key == key]
    raw[, use] <- as.matrix(x[genes, use, drop = FALSE])
    if (length(bg) < config$min_background_genes) next
    for (id in use) {
      values <- as.numeric(x[bg, id])
      if (!all(is.finite(values))) next
      complete_bg[id] <- TRUE
      ranks <- (rank(values, ties.method = "average") - 0.5) / length(bg)
      ranked[id, 1L] <- mean(ranks[match(genes, bg)])
    }
  }
  ref <- ids[m$role == "reference"]
  query <- ids[m$role == "query"]
  # Raw reference scaling excludes every incomplete signature sample.
  if (length(ids) && any(!is.finite(raw))) stop("Frozen signature has incomplete samples.")
  audit <- list(raw_signature = raw, common_signatures = stats::setNames(list(genes), contract$anchor),
    coverage = data.frame(anchor = contract$anchor, estimable = length(genes) >= config$min_signature_genes),
    rank_scores = ranked, reference_ids = ref, query_ids = query, config = config)
  # Rank is standardized only on complete-background reference samples.
  scores <- if (length(ref) >= 2L && length(genes)) {
    raw_scores <- .hc_bio("raw_scores")(raw, audit$common_signatures, ref)
    unscaled <- list(common_raw = raw_scores, common_rank = ranked)
    scales <- list()
    values <- lapply(names(unscaled), function(definition) {
      x <- unscaled[[definition]]
      fit <- x[ref, 1L]
      fit <- fit[is.finite(fit)]
      mu <- mean(fit)
      sd <- stats::sd(fit)
      scales[[definition]] <<- data.frame(score_definition = definition, anchor = contract$anchor,
        reference_count = length(fit), reference_mean = mu, reference_sd = sd,
        status = if (is.finite(sd) && sd > 0) "estimable" else "not_estimable")
      if (is.finite(sd) && sd > 0) (x - mu) / sd else { x[] <- NA_real_; x }
    })
    names(values) <- names(unscaled)
    list(values = values, unscaled = unscaled, scales = .hc_bind(scales))
  } else list(values = list(common_raw = ranked, common_rank = ranked),
    unscaled = list(common_raw = ranked, common_rank = ranked), scales = data.frame())
  list(audit = audit, scores = scores, complete_background = complete_bg)
}

.hc_score_pool <- function(prepared, contract, scored, definition, config) {
  x <- scored$scores$values[[definition]][, 1L]
  ids <- names(x)[is.finite(x)]
  m <- contract$metadata[contract$metadata$sample_id %in% ids, , drop = FALSE]
  ref <- m[m$role == "reference", ]
  query <- m[m$role == "query", ]
  if (!nrow(ref) || !nrow(query)) return(.hc_subset(prepared, ref$sample_id, character()))
  p <- .hc_subset(prepared, ref$sample_id, query$sample_id)
  support <- .hc_bio("eligibility")(p, config)$queries
  query <- query[query$sample_id %in% support$sample_id[support$eligible], ]
  counts <- table(query$cohort_key)
  query <- query[query$cohort_key %in% names(counts)[counts >= config$min_query_per_cohort], ]
  .hc_subset(prepared, ref$sample_id, query$sample_id)
}

.hc_summary <- function(paired, config, id, anchor, test = TRUE) {
  .hc_bio("summarize")(paired, config, id, anchor, test = test, local_intervals = FALSE)
}

.hc_leave_one <- function(cohort, metadata, config, anchor, id) {
  c <- cohort[cohort$eligible_inference, , drop = FALSE]
  if (!nrow(c)) return(data.frame())
  source_map <- unique(metadata[metadata$role == "query", c("cohort_key", "source_system")])
  if (anyDuplicated(source_map$cohort_key)) stop("Conflicting source labels within cohort.")
  c$source <- source_map$source_system[match(c$cohort, source_map$cohort_key)]
  rows <- list()
  for (unit in c("cohort", "source")) {
    labels <- sort(unique(c[[unit]]))
    if (unit == "source") labels <- labels[.hc_bio("known")(labels)]
    for (label in labels) {
      keep <- if (unit == "source") is.na(c$source) | c$source != label else c$cohort != label
      part <- c[keep, , drop = FALSE]
      values <- matrix(part$value, ncol = 1L, dimnames = list(part$cohort, anchor))
      sizes <- matrix(part$n, ncol = 1L, dimnames = list(part$cohort, anchor))
      result <- .hc_bio("infer_matrix")(values, config, id, sizes, test = FALSE)$summary
      result$unit <- unit
      result$omitted <- label
      result$comparison <- id
      rows[[length(rows) + 1L]] <- result
    }
  }
  .hc_bind(rows)
}
