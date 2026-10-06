# I/O and targets wiring. Computations live in analysis-specific R/ helpers.

.ablation03_hc_save <- function(value, filename, runtime) {
  directory <- file.path(runtime$cache_root, "biology-high-coverage")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, filename)
  saveRDS(value, path, version = 3)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.ablation03_hc_inventory <- function(representations, sources, biology, runtime) {
  p <- readRDS(representations)$analysis$prepared
  cache <- .ablation03_read_artifact(biology, "expression-anchor-cache.rds")
  hashes <- unname(tools::md5sum(sources))
  if (!identical(hashes, c(cache$source$md5, cache$signature$md5)))
    stop("High-coverage expression/signature provenance changed.", call. = FALSE)
  modules <- p$module_manifest$modules
  if (!"cohort" %in% names(modules) ||
      length(intersect(as.character(modules$cohort), as.character(p$query_metadata$cohort))))
    stop("High-coverage query overlaps frozen bank or training manifest is unavailable.")
  atlas <- readRDS(sources[1L])
  result <- .hc_inventory(atlas, cache$anchors, p$reference_metadata, p$query_metadata)
  rm(atlas)
  invisible(gc())
  result$source_md5 <- hashes
  .ablation03_hc_save(result, "measurement-inventory.rds", runtime)
}

.ablation03_hc_frontier <- function(inventory, config, runtime) {
  result <- .hc_frontier(readRDS(inventory), config)
  .ablation03_hc_save(result, "measurement-frontier.rds", runtime)
}

.ablation03_hc_contracts <- function(frontier, runtime) {
  result <- readRDS(frontier)
  .ablation03_hc_save(result, "frozen-contracts.rds", runtime)
}

.ablation03_hc_scores <- function(contracts, sources, runtime) {
  frozen <- readRDS(contracts)
  if (!identical(frozen$source_md5, unname(tools::md5sum(sources))))
    stop("High-coverage source changed after freezing the measurement contract.")
  atlas <- readRDS(sources[1L])
  identities <- vapply(frozen$contracts, function(x) .hc_hash(x[c("anchor", "genes",
    "background", "metadata")]), character(1))
  unique_ids <- names(identities)[!duplicated(identities)]
  values <- lapply(frozen$contracts[unique_ids], .hc_score_contract, atlas = atlas, config = frozen$config)
  result <- values[match(identities, identities[unique_ids])]
  names(result) <- names(frozen$contracts)
  rm(atlas)
  invisible(gc())
  .ablation03_hc_save(result, "scores.rds", runtime)
}

.ablation03_hc_specs <- function(contracts) {
  f <- readRDS(contracts)
  specs <- list()
  # Equivalent frozen gene/cohort/sample/background contracts are executed once
  # in the DAG; every prespecified tier remains represented in the frontier.
  identities <- vapply(f$contracts, function(x) .hc_hash(x[c("anchor", "genes",
    "background", "metadata")]), character(1))
  for (id in names(identities)[!duplicated(identities)]) for (definition in f$config$score_definitions) {
    specs[[length(specs) + 1L]] <- list(id = id, definition = definition,
      tier_ids = names(identities)[identities == identities[id]])
  }
  specs
}

.ablation03_hc_retrieval <- function(representations, contracts, scores, spec, runtime) {
  f <- readRDS(contracts)
  c <- f$contracts[[spec$id]]
  config <- f$config
  config$primary_anchors <- c$anchor
  prepared <- readRDS(representations)$analysis$prepared
  scored <- readRDS(scores)[[spec$id]]
  p <- .hc_score_pool(prepared, c, scored, spec$definition, config)
  result <- list(spec = spec, inference = list(), paired = list(), cohort = list(),
    retrieval = list(), query_ids = p$query_metadata$sample_id,
    reference_ids = p$reference_metadata$sample_id)
  score <- scored$scores
  score$values <- score$values[spec$definition]
  for (pool in c("all", "same_cancer")) {
    id <- paste(spec$id, spec$definition, pool, sep = ":")
    retrieved <- if (nrow(p$query_metadata) && nrow(p$reference_metadata) >= config$k)
      lapply(c("Direct", "d1"), function(arm) .hc_bio("retrieve")(p, config,
        pool = pool, arm = arm)) else list()
    utility <- .hc_bio("utility")(.hc_bind(lapply(retrieved, `[[`, "neighbors")), score, p, config)
    paired <- .hc_bio("pair_utility")(utility)
    summary <- .hc_summary(paired, config, id, c$anchor)
    summary$summary$score_definition <- spec$definition
    summary$summary$pool <- pool
    summary$summary$comparison <- id
    if (!length(c$genes) || !nrow(c$metadata)) summary$summary$reason <- "no_high_coverage_measurement_contract" else
      if (spec$definition == "common_rank" && length(c$background) < config$min_background_genes)
        summary$summary$reason <- "background_below_5000_genes" else
          if (!nrow(p$query_metadata)) summary$summary$reason <- "no_complete_score_supported_queries"
    result$inference[[pool]] <- summary$summary
    result$cohort[[pool]] <- summary$cohort
    result$paired[[pool]] <- paired
    result$retrieval[[pool]] <- retrieved
  }
  result$inference <- .hc_bind(result$inference)
  .ablation03_hc_save(result, paste0("retrieval-", spec$id, "-", spec$definition, ".rds"), runtime)
}

.ablation03_hc_readout_specs <- function(contracts) {
  f <- readRDS(contracts)
  as.list(f$frontier$id[f$frontier$selected])
}

.ablation03_hc_readout <- function(representations, contracts, scores, id, runtime) {
  f <- readRDS(contracts)
  c <- f$contracts[[id]]
  scored <- readRDS(scores)[[id]]
  prepared <- readRDS(representations)$analysis$prepared
  outputs <- list()
  for (definition in f$config$score_definitions) {
    config <- f$config
    config$primary_anchors <- c$anchor
    config$score_definitions <- definition
    p <- .hc_score_pool(prepared, c, scored, definition, config)
    a <- scored$audit
    a$reference_ids <- p$reference_metadata$sample_id
    a$query_ids <- p$query_metadata$sample_id
    outputs[[definition]] <- if (nrow(p$query_metadata)) .hc_bio("readout")(p, a, scored$scores, config) else
      list(status = "not_estimable", reason = "no_complete_score_supported_queries",
        predictions = data.frame(), cv = data.frame(), fits = list())
  }
  .ablation03_hc_save(list(id = id, results = outputs), paste0("readout-", id, ".rds"), runtime)
}

.ablation03_hc_sensitivity <- function(representations, contracts, scores, retrieval, id, runtime) {
  f <- readRDS(contracts)
  c <- f$contracts[[id]]
  config <- f$config
  config$primary_anchors <- c$anchor
  definition <- config$primary_score_definition
  scored <- readRDS(scores)[[id]]
  prepared <- readRDS(representations)$analysis$prepared
  p <- .hc_score_pool(prepared, c, scored, definition, config)
  score <- scored$scores
  score$values <- score$values[definition]
  existing <- lapply(as.character(retrieval), readRDS)
  baseline <- Filter(function(r) id %in% r$spec$tier_ids && r$spec$definition == definition, existing)[[1L]]
  rows <- cohorts <- status <- list()
  for (change in c("module_unscaled", "tissue_balanced", "assay_type", "source_system", "platform_id")) {
    if (!nrow(p$query_metadata)) {
      row <- .hc_summary(data.frame(), config, change, c$anchor, test = FALSE)$summary
      row$change <- change
      rows[[change]] <- row
      next
    }
    distance <- change %in% c("module_unscaled", "tissue_balanced")
    modified <- if (distance) list(.hc_bio("retrieve")(p, config, pool = "same_cancer",
      arm = "d1", distance = change)) else lapply(c("Direct", "d1"), function(arm)
        .hc_bio("retrieve")(p, config, pool = "same_cancer", arm = arm, restriction = change))
    if (distance) {
      direct <- baseline$retrieval$same_cancer[[1L]]
      direct$neighbors$distance <- change
      modified <- c(list(direct), modified)
    }
    utility <- .hc_bio("utility")(.hc_bind(lapply(modified, `[[`, "neighbors")), score, p, config)
    paired <- .hc_bio("pair_utility")(utility)
    base <- baseline$paired$same_cancer
    matched <- merge(base, paired, by = c("sample_id", "cohort", "anchor"), suffixes = c("_base", "_modified"))
    matched$delta <- matched$delta_modified - matched$delta_base
    summary <- .hc_summary(matched, config, change, c$anchor, test = FALSE)
    summary$summary$change <- change
    rows[[change]] <- summary$summary
    cc <- summary$cohort
    cc$change <- rep(change, nrow(cc))
    cohorts[[change]] <- cc
    status[[change]] <- data.frame(change = change, paired_queries = nrow(matched),
      reason = if (nrow(matched)) "matched_patient_difference_in_differences" else "no_eligible_matched_queries")
  }
  .ablation03_hc_save(list(id = id, inference = .hc_bind(rows), cohort = .hc_bind(cohorts),
    status = .hc_bind(status)), paste0("sensitivity-", id, ".rds"), runtime)
}

.ablation03_hc_inference <- function(contracts, scores, retrieval, readout, sensitivity, baseline_files, runtime) {
  frozen <- readRDS(contracts)
  score_list <- readRDS(scores)
  baseline_hashes <- tools::md5sum(baseline_files)
  retrieved <- lapply(as.character(retrieval), readRDS)
  readers <- lapply(as.character(readout), readRDS)
  sensitivity_results <- lapply(as.character(sensitivity), readRDS)
  rows <- cohorts <- primary_pairs <- list()
  paired_grid <- list()
  for (r in retrieved) for (id in r$spec$tier_ids) {
    row <- r$inference
    row$id <- id
    frontier <- frozen$frontier[match(id, frozen$frontier$id), ]
    row$tier <- frontier$tier
    row$coverage <- frontier$coverage
    row$common_genes <- frontier$common_genes
    row$original_genes <- frontier$original_genes
    row$selected <- frontier$selected
    row$evidence_level <- frontier$evidence_level
    row$source_groups <- frontier$source_groups
    row$background_genes <- frontier$background_genes
    row$primary <- row$selected & row$score_definition == frozen$config$primary_score_definition &
      row$pool == frozen$config$primary_pool
    row$comparison <- paste(id, row$score_definition, row$pool, sep = ":")
    rows[[length(rows) + 1L]] <- row
    for (pool in names(r$cohort)) {
      paired_grid[[paste(id, r$spec$definition, pool, sep = ":")]] <- r$paired[[pool]]
      c <- r$cohort[[pool]]
      c$id <- rep(id, nrow(c))
      c$score_definition <- rep(r$spec$definition, nrow(c))
      c$pool <- rep(pool, nrow(c))
      cohorts[[length(cohorts) + 1L]] <- c
      if (frontier$selected && r$spec$definition == frozen$config$primary_score_definition &&
          pool == frozen$config$primary_pool) primary_pairs[[id]] <- r$paired[[pool]]
    }
  }
  grid <- .hc_bind(rows)
  # Four anchors always remain in each declared testing family, including NA.
  for (key in unique(paste(grid$tier, grid$score_definition, grid$pool))) {
    idx <- which(paste(grid$tier, grid$score_definition, grid$pool) == key)
    grid$q_value[idx] <- stats::p.adjust(grid$p_value[idx], method = "BH", n = 4L)
  }
  primary <- grid[grid$primary, ]
  primary$q_value <- stats::p.adjust(primary$p_value, "BH", n = 4L)
  primary$family <- "four_selected_high_coverage_rank_same_cancer_anchors"
  primary$effective_source_groups <- vapply(seq_len(nrow(primary)), function(i) {
    row <- primary[i, ]
    used <- .hc_bind(cohorts)
    used <- used[used$id == row$id & used$score_definition == row$score_definition &
      used$pool == row$pool & used$eligible_inference, ]
    m <- frozen$contracts[[row$id]]$metadata
    source <- m$source_system[m$role == "query" & m$cohort_key %in% used$cohort]
    length(unique(source[.hc_bio("known")(source)]))
  }, integer(1))
  primary$publication_breadth_eligible <- primary$cohort_count >= frozen$config$primary_min_query_cohorts &
    primary$effective_source_groups >= frozen$config$primary_min_source_groups
  grid$q_value[grid$primary] <- primary$q_value
  cohort <- .hc_bind(cohorts)
  changes <- list()
  contrasts <- list(raw_pool = c("common_raw:all", "common_raw:same_cancer"),
    rank_pool = c("common_rank:all", "common_rank:same_cancer"),
    all_score = c("common_raw:all", "common_rank:all"),
    same_score = c("common_raw:same_cancer", "common_rank:same_cancer"))
  for (id in names(frozen$contracts)) for (change in names(contrasts)) {
    keys <- paste(id, contrasts[[change]], sep = ":")
    a <- paired_grid[[keys[1L]]]
    b <- paired_grid[[keys[2L]]]
    paired <- merge(a, b, by = c("sample_id", "cohort", "anchor"), suffixes = c("_a", "_b"))
    paired$delta <- paired$delta_b - paired$delta_a
    row <- .hc_summary(paired, frozen$config, paste(id, change), frozen$contracts[[id]]$anchor,
      test = FALSE)$summary
    row$id <- id
    row$change <- change
    changes[[length(changes) + 1L]] <- row
  }
  leave_one <- .hc_bind(lapply(frozen$frontier$id[frozen$frontier$selected], function(id) {
    c <- frozen$contracts[[id]]
    p <- cohort[cohort$id == id & cohort$score_definition == frozen$config$primary_score_definition &
      cohort$pool == frozen$config$primary_pool, ]
    .hc_leave_one(p, c$metadata, frozen$config, c$anchor, id)
  }))
  source_summary <- .hc_bind(lapply(frozen$frontier$id[frozen$frontier$selected], function(id) {
    c <- frozen$contracts[[id]]
    p <- cohort[cohort$id == id & cohort$score_definition == frozen$config$primary_score_definition &
      cohort$pool == frozen$config$primary_pool & cohort$eligible_inference, ]
    if (!nrow(p)) return(NULL)
    m <- c$metadata[c$metadata$role == "query", ]
    source <- unique(m[, c("cohort_key", "source_system")])
    p$source_system <- source$source_system[match(p$cohort, source$cohort_key)]
    .hc_bind(lapply(split(p, p$source_system), function(x) data.frame(id = id, anchor = c$anchor,
      source_system = x$source_system[1L], cohort_count = nrow(x), query_count = sum(x$n),
      cohort_equal_mean = mean(x$value), negative_cohorts = sum(x$value < 0),
      positive_cohorts = sum(x$value > 0), known_source = .hc_bio("known")(x$source_system[1L]),
      scope = "descriptive_source_direction_not_independent_source_inference")))
  }))
  selected <- frozen$frontier$id[frozen$frontier$selected]
  core_contract <- Reduce(intersect, lapply(frozen$contracts[selected], `[[`, "query_cohorts"))
  # A common validation population must survive every anchor's score/size gates.
  core <- Reduce(intersect, lapply(selected, function(id) sort(unique(cohort$cohort[
    cohort$id == id & cohort$score_definition == frozen$config$primary_score_definition &
      cohort$pool == frozen$config$primary_pool & cohort$eligible_inference]))))
  core_rows <- .hc_bind(lapply(selected, function(id) {
    c <- frozen$contracts[[id]]
    paired <- primary_pairs[[id]]
    if (is.null(paired)) paired <- data.frame()
    if (nrow(paired)) paired <- paired[paired$cohort %in% core, ]
    row <- .hc_summary(paired, frozen$config, paste0("core:", id), c$anchor, test = FALSE)$summary
    row$comparison <- "C_core_fixed_contract"
    row
  }))
  predictions <- .hc_bind(lapply(readers, function(r) .hc_bind(lapply(r$results, function(x) {
    if (!nrow(x$predictions)) return(NULL)
    metrics <- .hc_bio("prediction_metrics")(x, frozen$config)
    metrics$id <- r$id
    metrics
  }))))
  readout_rows <- list()
  if (nrow(predictions)) for (definition in frozen$config$score_definitions) for (reader in c("ridge", "xgboost")) {
    x <- predictions[predictions$score_definition == definition & predictions$reader == reader, ]
    d <- x[x$arm == "Direct", ]
    b <- x[x$arm == "d1", ]
    pairs <- merge(d, b, by = c("anchor", "cohort", "score_definition", "reader"),
      suffixes = c("_direct", "_d1"))
    for (anchor in frozen$config$primary_anchors) {
      part <- pairs[pairs$anchor == anchor, ]
      values <- matrix(part$mae_d1 - part$mae_direct, ncol = 1L,
        dimnames = list(part$cohort, anchor))
      counts <- matrix(part$query_count_d1, ncol = 1L, dimnames = list(part$cohort, anchor))
      eligible <- part$query_count_d1 >= frozen$config$min_query_per_cohort
      row <- .hc_bio("infer_matrix")(values[eligible, , drop = FALSE], frozen$config,
        paste("readout", definition, reader, anchor), counts[eligible, , drop = FALSE])$summary
      row$score_definition <- definition
      row$reader <- reader
      readout_rows[[length(readout_rows) + 1L]] <- row
    }
  }
  readout_inference <- .hc_bind(readout_rows)
  if (nrow(readout_inference)) for (key in unique(paste(readout_inference$score_definition, readout_inference$reader))) {
    idx <- which(paste(readout_inference$score_definition, readout_inference$reader) == key)
    readout_inference$q_value[idx] <- stats::p.adjust(readout_inference$p_value[idx], "BH", n = 4L)
  }
  status <- .hc_bind(lapply(readers, function(r) .hc_bind(lapply(names(r$results), function(definition) {
    x <- r$results[[definition]]
    data.frame(id = r$id, branch = "readout", score_definition = definition,
      status = x$status, reason = x$reason)
  }))))
  eligibility <- frozen$eligibility
  gene_lists <- .hc_bind(lapply(frozen$contracts, function(c) data.frame(anchor = c$anchor,
    tier = c$tier, gene = c$original_genes, retained = c$original_genes %in% c$genes,
    gene_hash = c$gene_hash, cohort_hash = c$cohort_hash)))
  scales <- .hc_bind(lapply(names(score_list), function(id) {
    x <- score_list[[id]]$scores$scales
    x$id <- rep(id, nrow(x))
    x
  }))
  result <- list(config = frozen$config, frontier = frozen$frontier, primary = primary,
    grid = grid, cohort = cohort, leave_one = leave_one, core = core_rows,
    source_summary = source_summary,
    core_cohorts = core, core_contract_cohorts = core_contract,
    readout_metrics = predictions, readout_inference = readout_inference,
    status = status, eligibility = eligibility, gene_lists = gene_lists, scales = scales,
    changes = .hc_bind(changes),
    sensitivity = .hc_bind(lapply(sensitivity_results, `[[`, "inference")),
    sensitivity_status = .hc_bind(lapply(sensitivity_results, `[[`, "status")),
    measurement_hash = frozen$measurement_hash, baseline_hashes = baseline_hashes,
    source_md5 = frozen$source_md5, runtime_identity = runtime[c("R", "package", "renv")],
    inference_scope = frozen$config$inference_scope)
  if (!identical(baseline_hashes, tools::md5sum(baseline_files))) stop("Broad baseline changed.")
  .ablation03_hc_save(result, "inference.rds", runtime)
}
