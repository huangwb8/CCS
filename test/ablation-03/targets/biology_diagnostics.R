# Wiring and ordinary I/O only. The installed CCS namespace owns calculations.

.ablation03_bio_function <- function(name) getFromNamespace(
  paste0(".ablation_bio_", name), "CCS")

.ablation03_bio_code_identity <- function(source_file) {
  definitions <- parse(source_file, keep.source = FALSE)
  for (expr in definitions) {
    if (!is.call(expr) || !identical(expr[[1L]], as.name("<-")) ||
        !is.symbol(expr[[2L]])) next
    name <- as.character(expr[[2L]])
    installed <- getFromNamespace(name, "CCS")
    if (!identical(deparse(expr[[3L]][[3L]]), deparse(body(installed))) ||
        !identical(expr[[3L]][[2L]], formals(installed))) {
      stop("Biology diagnostic source and installed CCS differ: ", name, call. = FALSE)
    }
  }
  unname(tools::md5sum(source_file))
}

.ablation03_bio_save <- function(value, filename, runtime) {
  directory <- file.path(runtime$cache_root, "biology-diagnostics")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, filename)
  saveRDS(value, path, version = 3)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.ablation03_bio_sources <- function(biology_inputs) {
  cache <- .ablation03_read_artifact(biology_inputs, "expression-anchor-cache.rds")
  c(cache$source$path, cache$signature$path)
}

.ablation03_bio_audit <- function(representations, biology, baseline,
    data, config, sources, runtime) {
  prepared <- readRDS(representations)$analysis$prepared
  cache <- .ablation03_read_artifact(biology, "expression-anchor-cache.rds")
  neighbors <- .ablation03_read_artifact(list(directory = file.path(runtime$cache_root,
    "ablation-experiment")), "anchor_retrieval.rds")$neighbors
  old <- .ablation03_read_artifact(baseline, "ablation03-biology.rds")$cohort_deltas
  baseline_paths <- c(file.path(baseline$directory, c("anchor_cohort_deltas.csv",
    "anchor_per_query_utility.rds", "ablation03-biology.rds")),
    file.path(runtime$cache_root, "ablation-experiment", "anchor_retrieval.rds"))
  hashes <- tools::md5sum(baseline_paths)
  check <- .ablation03_bio_function("baseline_check")(cache, neighbors, old,
    config$baseline_tolerance)
  if (!identical(digest::digest(file = sources[1L], algo = "md5"), cache$source$md5) ||
      !identical(digest::digest(file = sources[2L], algo = "md5"), cache$signature$md5)) {
    stop("Diagnostic expression/signature provenance differs from frozen baseline.", call. = FALSE)
  }
  signatures <- readRDS(sources[2L])
  atlas <- readRDS(sources[1L])
  result <- .ablation03_bio_function("audit")(atlas, cache, prepared,
    data$input$object@Repeat$geneSet, signatures, config,
    preprocessing = data$input$expression_preprocessing)
  # The installed package's legacy status label embeds the original threshold.
  if (startsWith(result$rank_status, "background_below_")) {
    result$rank_status <- paste0("background_below_", config$min_background_genes, "_genes")
  }
  rm(atlas)
  invisible(gc())
  result$baseline_check <- check
  result$baseline_hashes <- hashes
  result$source_md5 <- unname(tools::md5sum(sources))
  result$list_identity <- digest::digest(list(config$experiment_id,
    result$background, result$common_signatures), algo = "sha256")
  if (!identical(hashes, tools::md5sum(baseline_paths))) {
    stop("Biology baseline was modified during audit.", call. = FALSE)
  }
  .ablation03_bio_save(result, "audit.rds", runtime)
}

.ablation03_bio_gene_lists <- function(audit, runtime) {
  a <- readRDS(audit)
  .ablation03_bio_save(list(config = a$config, background = a$background,
    signatures = a$common_signatures, input_support = a$input_support,
    list_identity = a$list_identity), "frozen-gene-lists-v1.rds", runtime)
}

.ablation03_bio_scores <- function(audit, runtime) {
  .ablation03_bio_save(.ablation03_bio_function("scores")(readRDS(audit)),
    "scores.rds", runtime)
}

.ablation03_bio_specs <- function(config) {
  specs <- list()
  for (pool in c("all", "same_cancer")) for (arm in c("Direct", "d1")) {
    specs[[length(specs) + 1L]] <- list(pool = pool, arm = arm,
      distance = "baseline", restriction = NULL)
  }
  for (distance in c("module_unscaled", "tissue_balanced")) for (pool in c("same_cancer", "all")) {
    specs[[length(specs) + 1L]] <- list(pool = pool, arm = "d1",
      distance = distance, restriction = NULL)
  }
  for (restriction in config$technical_restrictions) for (arm in c("Direct", "d1")) {
    specs[[length(specs) + 1L]] <- list(pool = "same_cancer", arm = arm,
      distance = "baseline", restriction = restriction)
  }
  specs
}

.ablation03_bio_retrieval <- function(representations, audit, config, spec, runtime) {
  readRDS(audit)$baseline_check
  prepared <- readRDS(representations)$analysis$prepared
  result <- do.call(.ablation03_bio_function("retrieve"),
    c(list(prepared = prepared, config = config), spec))
  result$spec <- spec
  id <- paste(spec$pool, spec$arm, spec$distance,
    if (is.null(spec$restriction)) "none" else spec$restriction, sep = "-")
  .ablation03_bio_save(result, paste0("retrieval-", id, ".rds"), runtime)
}

.ablation03_bio_readout <- function(representations, audit, scores, config, runtime) {
  prepared <- readRDS(representations)$analysis$prepared
  .ablation03_bio_save(.ablation03_bio_function("readout")(prepared,
    readRDS(audit), readRDS(scores), config), "readout.rds", runtime)
}

.ablation03_bio_inference <- function(representations, audit, scores,
    retrieval, readout, config, runtime) {
  prepared <- readRDS(representations)$analysis$prepared
  a <- readRDS(audit)
  result <- .ablation03_bio_function("inference")(prepared, a,
    readRDS(scores), lapply(as.character(retrieval), readRDS), readRDS(readout), config)
  if (!identical(a$baseline_hashes, tools::md5sum(names(a$baseline_hashes)))) {
    stop("Baseline products changed during diagnostic execution.", call. = FALSE)
  }
  # IFN and IL6 are now primary anchors; omit duplicate supplementary tests.
  if (all(c("ifn", "il6") %in% config$primary_anchors)) {
    for (name in names(result)) {
      value <- result[[name]]
      if (is.data.frame(value) && "comparison" %in% names(value)) {
        result[[name]] <- value[
          !grepl("^supplement:.*:split_ifn_il6$", value$comparison), , drop = FALSE]
      } else if (is.list(value) && !is.null(names(value))) {
        result[[name]] <- value[
          !grepl("^supplement:.*:split_ifn_il6$", names(value))]
      }
    }
  }
  result$config <- config
  result$list_identity <- a$list_identity
  .ablation03_bio_save(result, "inference.rds", runtime)
}
