#!/usr/bin/env Rscript
# A synthetic atlas and synthetic sample contract retain the frozen public
# bank structure. All scientific stages use the formal _targets.R unchanged.
main <- function() {
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (file.exists(file.path(repo, "ablation-03.Rproj"))) {
  repo <- normalizePath(file.path(repo, "..", ".."), winslash = "/", mustWork = TRUE)
}
project <- file.path(repo, "test", "ablation-03")
setwd(project)
source("renv/activate.R")
local_bin <- file.path(Sys.getenv("USERPROFILE"), ".local", "bin")
if (file.exists(file.path(local_bin, "pandoc.exe")))
  Sys.setenv(PATH = paste(local_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
source("targets/functions.R")
run_root <- if (length(args) >= 2L) normalizePath(args[2L], winslash = "/", mustWork = TRUE) else
  file.path(project, "tmp", "tests", paste0("biology-", format(Sys.time(), "%Y%m%d-%H%M%S")))
test_root <- normalizePath(file.path(project, "tmp", "tests"), winslash = "/", mustWork = TRUE)
if (!startsWith(tolower(run_root), paste0(tolower(test_root), "/"))) stop("Test run must remain under tmp/tests.")
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
source_bundle <- if (length(args)) args[1L] else
  file.path(project, "tmp", "exam-20260923-01", "input-subset.rds")
if (!file.exists(source_bundle)) stop("Provide a local frozen input bundle for the synthetic bank fixture.")
if (length(args) < 2L) {
inputs <- readRDS(source_bundle)
config <- yaml::read_yaml("config/biology-diagnostics.yml")
set.seed(config$seed)
original <- inputs$metadata
cohorts <- original[!duplicated(original$cohort_key), ]
# Retain all bank cohorts, three reference and supported query cohorts with 36 samples,
# and descriptive small cohorts. Only input size differs from the formal DAG.
reference_labels <- table(cohorts$cancer_type[cohorts$analysis_set == "reference_atlas"])
supported_labels <- names(reference_labels)[reference_labels >= 2L]
large_query_cohorts <- which(cohorts$analysis_set == "external_query" &
  cohorts$cancer_type %in% supported_labels)[1:3]
large_reference_cohorts <- which(cohorts$analysis_set == "reference_atlas")[1:3]
metadata <- do.call(rbind, lapply(seq_len(nrow(cohorts)), function(i) {
  count <- if (i %in% c(large_query_cohorts, large_reference_cohorts)) 36L else
    if (cohorts$analysis_set[i] == "reference_atlas") 8L else 2L
  x <- cohorts[rep(i, count), ]
  x$sample_id <- paste0("synthetic-", i, "-", seq_len(nrow(x)))
  x$metadata_status <- "confirmed"
  x$evidence_basis <- "synthetic_fixture_not_patient_data"
  x$evidence_url <- NA_character_
  x$platform_id <- if (i %% 3L == 0L) "fixture-platform" else "unknown"
  x$assay_type <- if (i %% 2L == 0L) "RNAseq" else "Microarray"
  x$source_system <- paste0("fixture-source-", i %% 3L)
  x$duplicate_sample_id_global <- FALSE
  x$duplicate_within_cohort <- FALSE
  x
}))
metadata$row_order <- seq_len(nrow(metadata))
ids <- metadata$sample_id
reference_ids <- ids[metadata$analysis_set == "reference_atlas"]
object <- inputs$object
full <- inputs$resCCS_full
synthetic_object <- function(object) {
  for (name in names(object@Data$Probability)) {
    x <- object@Data$Probability[[name]]
    if (is.null(dim(x))) next
    selected <- if (name == "d2") reference_ids else ids
    object@Data$Probability[[name]] <- matrix(stats::runif(length(selected) * ncol(x)),
      length(selected), dimnames = list(selected, colnames(x)))
  }
  object@Data$CCS <- rep("fixture", length(reference_ids))
  object@Data$CancerType <- metadata$cancer_type[metadata$analysis_set == "reference_atlas"]
  object
}
inputs$object <- inputs$resCCS_ablation <- synthetic_object(object)
inputs$resCCS_full <- synthetic_object(full)
inputs$full_d1 <- inputs$resCCS_full@Data$Probability$d1
inputs$metadata <- inputs$ablation_metadata <- metadata
signature_path <- Sys.getenv("CCS_GENE_SIGNATURE_RDS",
  unset = "E:/RCloud/database/Signature/report/GeneSignature-HWB.rds")
signatures <- readRDS(signature_path)
source("02.02.00. 生物锚点分析_functions.R")
anchor_genes <- unlist(.biology_select_anchors(yaml::read_yaml("config/biological-anchors.yml"), signatures))
input_genes <- unique(as.character(unlist(inputs$object@Repeat$geneSet)))
genes <- unique(c(input_genes, anchor_genes, paste0("fixture-gene-", seq_len(5002L))))
data <- list()
atlas <- list()
for (i in seq_len(nrow(cohorts))) {
  m <- metadata[metadata$cohort_key == cohorts$cohort_key[i], ]
  x <- matrix(round(stats::rnorm(length(genes) * nrow(m)), 1L), length(genes),
    dimnames = list(genes, m$sample_id))
  signature_rows <- intersect(unique(anchor_genes), genes)
  x[signature_rows, ] <- stats::rnorm(length(signature_rows) * nrow(m))
  # A different monotone scale plus ties; no result-driven preprocessing.
  if (i %% 5L == 0L) x <- x * 100 + 300
  x["fixture-gene-5001", ] <- 1
  if (i == 1L) x <- x[rownames(x) != "fixture-gene-5002", , drop = FALSE]
  tissue <- cohorts$tissue[i]
  name <- cohorts$cohort[i]
  data[[tissue]][[name]] <- list(expr = x[input_genes, , drop = FALSE],
    subtype = rep("fixture", nrow(m)))
  atlas[[tissue]][[name]] <- list(expr = x)
}
inputs$data <- inputs$data_all <- data
inputs$biology_inputs <- inputs$structural_inputs <- inputs$data_profile <- NULL
input_path <- file.path(run_root, "synthetic-input.rds")
atlas_path <- file.path(run_root, "synthetic-atlas.rds")
saveRDS(inputs, input_path)
saveRDS(atlas, atlas_path)
rm(inputs, atlas, data, object, full)
invisible(gc())
} else {
  input_path <- file.path(run_root, "synthetic-input.rds")
  atlas_path <- file.path(run_root, "synthetic-atlas.rds")
  signature_path <- Sys.getenv("CCS_GENE_SIGNATURE_RDS",
    unset = "E:/RCloud/database/Signature/report/GeneSignature-HWB.rds")
  stopifnot(file.exists(input_path), file.exists(atlas_path))
}
formal_files <- c(list.files("reports", recursive = TRUE, full.names = TRUE),
  list.files("raw", recursive = TRUE, full.names = TRUE),
  list.files(".", pattern = "\\.html$", full.names = TRUE))
before <- tools::md5sum(formal_files)
old_config <- readBin("_targets.yaml", "raw", n = file.info("_targets.yaml")$size)
on.exit(writeBin(old_config, "_targets.yaml"), add = TRUE)
Sys.setenv(CCS_ABLATION_INPUT_RDS = input_path,
  CCS_FULL_EXPRESSION_RDS = atlas_path,
  CCS_GENE_SIGNATURE_RDS = signature_path,
  CCS_ABLATION_CACHE_ROOT = file.path(run_root, "cache"),
  CCS_ABLATION_OUTPUT_ROOT = file.path(run_root, "output"),
  CCS_ABLATION_TARGET_WORKERS = "2", CCS_ABLATION_CORES = "1")
store <- file.path(run_root, "cache", "targets")
targets::tar_config_set(store = store)
started <- proc.time()
targets::tar_make(names = c(biology_diagnostic_inference, biology_report),
  store = store, script = "_targets.R")
elapsed <- proc.time() - started
audit <- readRDS(file.path(run_root, "cache", "biology-diagnostics", "audit.rds"))
result <- readRDS(file.path(run_root, "cache", "biology-diagnostics", "inference.rds"))
stopifnot(audit$baseline_check$status == "PASS", nrow(result$inference) > 0L,
  !is.null(audit$rank_scores),
  identical(audit$config$score_definitions, c("common_raw", "common_rank")),
  any(grepl("common_rank", result$inference$comparison)),
  is.list(result$exploratory), nrow(result$exploratory$inference) > 0L,
  all(is.na(result$exploratory$inference$p_value)),
  all(is.na(result$exploratory$inference$q_value)),
  identical(before, tools::md5sum(formal_files)),
  file.exists(file.path(run_root, "output", "02.02.00. 生物锚点分析.html")))
read_bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)
make <- function() targets::tar_make(names = c(biology_diagnostic_inference, biology_report), store = store)
meta <- function() targets::tar_meta(store = store, fields = c(name, time))
assert_rebuilt <- function(before, after, score = TRUE) {
  branches <- grep("^biology_diagnostic_retrieval_[[:xdigit:]]+$", before$name, value = TRUE)
  selected <- c(if (score) "biology_diagnostic_scores", "biology_diagnostic_readout",
    "biology_diagnostic_inference", branches)
  stopifnot(length(branches) == 14L,
    all(after$time[match(selected, after$name)] > before$time[match(selected, before$name)]))
}
# A harmless comment changes only scientific source identity, while the
# installed function bodies and the frozen representation remain identical.
code_path <- file.path(repo, "R", "ablation_biology.R")
code_bytes <- read_bytes(code_path)
expected_code <- code_bytes
on.exit({
  if (identical(read_bytes(code_path), expected_code)) writeBin(code_bytes, code_path)
  else warning("Concurrent source edit preserved; identity fixture comment was not removed.")
}, add = TRUE)
representation_path <- file.path(run_root, "cache", "01-representations", "representation-inputs.rds")
representation_bytes <- read_bytes(representation_path)
expected_representation <- representation_bytes
on.exit({
  if (identical(read_bytes(representation_path), expected_representation))
    writeBin(representation_bytes, representation_path)
  else warning("Concurrent representation edit preserved.")
}, add = TRUE)
times <- meta()
representation_hash <- tools::md5sum(representation_path)
expected_code <- c(code_bytes, charToRaw("\n# Synthetic targets identity propagation assertion.\n"))
writeBin(expected_code, code_path)
make()
assert_rebuilt(times, meta())
stopifnot(identical(representation_hash, tools::md5sum(representation_path)))
stopifnot(identical(read_bytes(code_path), expected_code))
writeBin(code_bytes, code_path)
expected_code <- code_bytes
make()
# Changing an unused audit field tests file content propagation without
# changing any scientific values or the frozen sample/bank contract.
times <- meta()
representation <- readRDS(representation_path)
representation$diagnostic_fixture_note <- "file content propagation assertion"
saveRDS(representation, representation_path, version = 3)
expected_representation <- read_bytes(representation_path)
make()
after <- meta()
assert_rebuilt(times, after, score = FALSE)
identity_target <- "biology_diagnostic_code_identity"
stopifnot(identical(times$time[match(identity_target, times$name)],
  after$time[match(identity_target, after$name)]),
  identical(read_bytes(representation_path), expected_representation))
writeBin(representation_bytes, representation_path)
expected_representation <- representation_bytes
make()
times <- targets::tar_meta(store = store, fields = c(name, time))
targets::tar_invalidate(biology_diagnostic_inference, store = store)
targets::tar_make(names = c(biology_diagnostic_inference, biology_report), store = store)
after <- targets::tar_meta(store = store, fields = c(name, time))
upstream <- c("representation_inputs", "biology_diagnostic_audit", "biology_diagnostic_scores",
  "biology_diagnostic_readout")
stopifnot(identical(times$time[match(upstream, times$name)], after$time[match(upstream, after$name)]),
  identical(before, tools::md5sum(formal_files)))
saveRDS(list(status = "PASS", elapsed = elapsed, run_root = run_root,
  baseline_check = audit$baseline_check, upstream_reused = upstream,
  identity_propagation = "PASS", representation_file_propagation = "PASS",
  assertions = "synthetic_same_DAG_report_and_local_recovery_formal_paths_unchanged"),
  file.path(run_root, "verification.rds"))
cat("Synthetic targets diagnostic integration PASS: ", run_root, "\n")
}
main()
