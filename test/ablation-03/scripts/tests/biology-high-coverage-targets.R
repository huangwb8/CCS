#!/usr/bin/env Rscript
# Same formal DAG, synthetic expression/probabilities and isolated cache/output.
main <- function() {
  repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists("ablation-03.Rproj")) repo <- normalizePath("../..", winslash = "/")
  project <- file.path(repo, "test", "ablation-03")
  setwd(project)
  args <- commandArgs(trailingOnly = TRUE)
  run_root <- if (length(args) >= 2L) normalizePath(args[2L], winslash = "/", mustWork = TRUE) else
    file.path(project, "tmp", "tests", paste0("high-coverage-", format(Sys.time(), "%Y%m%d-%H%M%S")))
  test_root <- normalizePath(file.path(project, "tmp", "tests"), winslash = "/", mustWork = TRUE)
  if (!startsWith(tolower(run_root), paste0(tolower(test_root), "/"))) stop("Test run must remain under tmp/tests.")
  dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
  # Keep the project library pinned while avoiding shared sandbox lock waits.
  Sys.setenv(RENV_PATHS_SANDBOX = file.path(run_root, "renv-sandbox"))
  source("renv/activate.R")
  source("R/biology_high_coverage.R")
  source("targets/biology_high_coverage.R")
  source_bundle <- if (length(args)) args[1L] else
    file.path(project, "tmp", "exam-20260923-01", "input-subset.rds")
  input_path <- file.path(run_root, "synthetic-input.rds")
  atlas_path <- file.path(run_root, "synthetic-atlas.rds")
  signature_path <- Sys.getenv("CCS_GENE_SIGNATURE_RDS",
    "E:/RCloud/database/Signature/report/GeneSignature-HWB.rds")
  if (!file.exists(input_path)) {
    inputs <- readRDS(source_bundle)
    cohorts <- inputs$metadata[!duplicated(inputs$metadata$cohort_key), ]
    set.seed(20261006)
    metadata <- do.call(rbind, lapply(seq_len(nrow(cohorts)), function(i) {
      x <- cohorts[rep(i, 24L), ]
      x$sample_id <- paste0("hc-fixture-", i, "-", seq_len(nrow(x)))
      x$cancer_type <- if (i %% 2L == 0L) "BRCA" else "CRC"
      x$metadata_status <- "confirmed"
      x$evidence_basis <- "synthetic_fixture_not_patient_data"
      x$evidence_url <- NA_character_
      x$source_system <- paste0("fixture-source-", i %% 3L)
      x$assay_type <- "RNAseq"
      x$platform_id <- if (i %% 4L == 0L) "unknown" else "fixture-platform"
      x$duplicate_sample_id_global <- FALSE
      x$duplicate_within_cohort <- FALSE
      x
    }))
    metadata$row_order <- seq_len(nrow(metadata))
    ids <- metadata$sample_id
    ref_ids <- ids[metadata$analysis_set == "reference_atlas"]
    synthetic <- function(object) {
      for (name in names(object@Data$Probability)) {
        x <- object@Data$Probability[[name]]
        if (is.null(dim(x))) next
        use <- if (name == "d2") ref_ids else ids
        object@Data$Probability[[name]] <- matrix(runif(length(use) * ncol(x)),
          length(use), dimnames = list(use, colnames(x)))
      }
      object@Data$CCS <- rep("fixture", length(ref_ids))
      object@Data$CancerType <- metadata$cancer_type[metadata$analysis_set == "reference_atlas"]
      object
    }
    inputs$object <- inputs$resCCS_ablation <- synthetic(inputs$object)
    inputs$resCCS_full <- synthetic(inputs$resCCS_full)
    inputs$full_d1 <- inputs$resCCS_full@Data$Probability$d1
    inputs$metadata <- inputs$ablation_metadata <- metadata
    source("02.02.00. 生物锚点分析_functions.R")
    anchors <- .biology_select_anchors(yaml::read_yaml("raw/config/biological-anchors.yml"), readRDS(signature_path))
    input_genes <- unique(as.character(unlist(inputs$object@Repeat$geneSet)))
    genes <- unique(c(input_genes, unlist(anchors), paste0("fixture-background-", 1:5002)))
    data <- atlas <- list()
    for (i in seq_len(nrow(cohorts))) {
      m <- metadata[metadata$cohort_key == cohorts$cohort_key[i], ]
      x <- matrix(rnorm(length(genes) * nrow(m)), length(genes), dimnames = list(genes, m$sample_id))
      # Missing patterns differ between cohorts; no algorithm outcomes used.
      if (i %% 7L == 0L) x <- x[!rownames(x) %in% unlist(lapply(anchors, function(g) g[1:2])), , drop = FALSE]
      if (i %% 9L == 0L) x[intersect(rownames(x), unlist(anchors))[1L], 1L] <- NA_real_
      tissue <- cohorts$tissue[i]
      name <- cohorts$cohort[i]
      raw <- matrix(rnorm(length(input_genes) * nrow(m)), length(input_genes),
        dimnames = list(input_genes, m$sample_id))
      data[[tissue]][[name]] <- list(expr = raw, subtype = rep("fixture", nrow(m)))
      atlas[[tissue]][[name]] <- list(expr = x)
    }
    inputs$data <- inputs$data_all <- data
    inputs$biology_inputs <- inputs$structural_inputs <- inputs$data_profile <- NULL
    saveRDS(inputs, input_path)
    saveRDS(atlas, atlas_path)
    rm(inputs, atlas, data)
    invisible(gc())
  }
  files <- c("_targets.yaml", list.files("raw", recursive = TRUE, full.names = TRUE),
    list.files("reports", recursive = TRUE, full.names = TRUE),
    list.files(".", pattern = "\\.html$", full.names = TRUE))
  before <- tools::md5sum(files)
  local_bin <- file.path(Sys.getenv("USERPROFILE"), ".local", "bin")
  Sys.setenv(PATH = paste(local_bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
    CCS_ABLATION_INPUT_RDS = input_path, CCS_FULL_EXPRESSION_RDS = atlas_path,
    CCS_GENE_SIGNATURE_RDS = signature_path,
    CCS_ABLATION_CACHE_ROOT = file.path(run_root, "cache"),
    CCS_ABLATION_OUTPUT_ROOT = file.path(run_root, "output"),
    TAR_CONFIG = file.path(run_root, "_targets.yaml"),
    CCS_ABLATION_TARGET_WORKERS = "2", CCS_ABLATION_CORES = "1")
  store <- file.path(run_root, "cache", "targets")
  targets::tar_config_set(store = store)
  targets::tar_make(names = c(biology_high_coverage_inference, biology_report), store = store)
  frozen <- readRDS(targets::tar_read(biology_high_coverage_contracts, store = store))
  inventory <- readRDS(targets::tar_read(biology_high_coverage_inventory, store = store))
  result <- readRDS(targets::tar_read(biology_high_coverage_inference, store = store))
  # One cohort may legitimately contain several platforms; keep sample labels.
  mixed <- inventory$metadata
  mixed$platform_id[2L] <- "fixture-second-platform"
  mixed_inventory <- .hc_inventory(readRDS(atlas_path), inventory$anchors,
    mixed[mixed$role == "reference", ], mixed[mixed$role == "query", ])
  stopifnot(nrow(mixed_inventory$metadata) == nrow(inventory$metadata),
    "fixture-second-platform" %in% mixed_inventory$metadata$platform_id)
  stopifnot(nrow(result$grid) == 64L, nrow(result$primary) == 4L,
    any(result$primary$status == "estimable"), any(result$readout_inference$status == "estimable"),
    file.exists(file.path(run_root, "output", "02.02.00. 生物锚点分析.html")),
    identical(before, tools::md5sum(files)))
  changed <- inventory
  changed$metadata$delta <- rnorm(nrow(changed$metadata))
  changed$utility_direct <- rnorm(10)
  changed$utility_d1 <- rnorm(10)
  stopifnot(identical(serialize(.hc_frontier(changed, frozen$config), NULL),
    serialize(.hc_frontier(inventory, frozen$config), NULL)))
  for (c in frozen$contracts) {
    stopifnot(all(c$genes %in% c$original_genes),
      length(c$genes) / length(c$original_genes) >= c$tier,
      !length(intersect(c$query_cohorts, c$reference_cohorts)))
  }
  # Edge gates are shared across anchors; exercise one complete tier while the
  # formal DAG above verifies all four anchors and all prespecified tiers.
  gate_config <- frozen$config
  gate_config$primary_anchors <- frozen$config$primary_anchors[1L]
  gate_config$coverage_tiers <- max(frozen$config$coverage_tiers)
  bad <- inventory
  bad$metadata$cancer_type[bad$metadata$role == "query"] <- "unknown"
  stopifnot(all(.hc_frontier(bad, gate_config)$frontier$query_count == 0))
  no_ref <- inventory
  ref_keys <- unique(no_ref$metadata$cohort_key[no_ref$metadata$role == "reference"])
  no_ref$metadata$metadata_status[no_ref$metadata$role == "reference" &
    no_ref$metadata$cohort_key != ref_keys[1L]] <- "unknown"
  stopifnot(all(.hc_frontier(no_ref, gate_config)$frontier$query_count == 0))
  low_bg <- inventory
  low_bg$cohorts <- lapply(low_bg$cohorts, function(x) { x$genes <- rownames(x$raw); x })
  stopifnot(!any(.hc_frontier(low_bg, gate_config)$frontier$rank_eligible))
  one_source <- inventory
  one_source$metadata$source_system <- "single-source"
  stopifnot(!any(.hc_frontier(one_source, gate_config)$frontier$breadth_eligible))
  # Restore only the downstream target; verify measured upstream reuse.
  times <- targets::tar_meta(store = store, fields = c(name, time))
  targets::tar_invalidate(biology_high_coverage_inference, store = store)
  targets::tar_make(names = c(biology_high_coverage_inference, biology_report), store = store)
  after <- targets::tar_meta(store = store, fields = c(name, time))
  upstream <- c("biology_high_coverage_inventory", "biology_high_coverage_contracts",
    "biology_high_coverage_scores", "representation_inputs", "biology_diagnostic_inference")
  stopifnot(identical(times$time[match(upstream, times$name)], after$time[match(upstream, after$name)]),
    identical(before, tools::md5sum(files)))
  saveRDS(list(status = "PASS", formal_entrypoint = "_targets.R", store = store,
    subject_md5 = tools::md5sum(c("R/biology_high_coverage.R", "targets/biology_high_coverage.R",
      "config/biology-high-coverage.yml", "_targets.R", "02.02.00. 生物锚点分析.Rmd",
      "scripts/tests/biology-high-coverage-targets.R")),
    baseline_unchanged = TRUE, outcome_blind_contract = TRUE, upstream_reused = upstream,
    assertions = "four_anchors_all_tiers_raw_rank_readout_missingness_source_reference_background_recovery"),
    file.path(run_root, "verification.rds"))
  cat("High-coverage same-DAG synthetic verification PASS\n")
}
main()
