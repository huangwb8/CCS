#!/usr/bin/env Rscript
# Prepare expression anchors for both biological and structural analyses.

options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") {
  utf8_locale <- Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")
  if (!nzchar(utf8_locale)) stop("biology cache: UTF-8 locale is required for signature names.")
}

bootstrap <- c("00-workflow_functions.R",
  "test/ablation-03/00-workflow_functions.R")
bootstrap <- bootstrap[file.exists(bootstrap)][1L]
if (is.na(bootstrap)) stop("Run from ablation-03 or the repository root.", call. = FALSE)
source(bootstrap, local = TRUE)
ablation_dir <- .wf_root
out_dir <- .wf_output("01-biology")
cache_path <- file.path(out_dir, "expression-anchor-cache.rds")
result_dir <- .wf_output("01-representations")
config_path <- file.path(ablation_dir, "config", "biological-anchors.yml")
full_path <- Sys.getenv(
  "CCS_FULL_EXPRESSION_RDS",
  unset = "E:/Sync/@Analysis/PanCan_Data/Level 1/PanCan_CancerSample_DataListForCCS_GEO+cBioPortal+UCXCXenav20240809.rds"
)
sig_path <- Sys.getenv(
  "CCS_GENE_SIGNATURE_RDS",
  unset = "E:/RCloud/database/Signature/report/GeneSignature-HWB.rds"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(full_path)) stop("ablation-03 cache: expression RDS not found: ", full_path, call. = FALSE)
if (!file.exists(sig_path)) stop("ablation-03 cache: signature RDS not found: ", sig_path, call. = FALSE)

source(file.path(ablation_dir, "02-ablation03-representation_functions.R"))
source(file.path(ablation_dir, "03-ablation03-biology_functions.R"))
.wf_validate(result_dir)
contract <- readRDS(file.path(result_dir, "sample-contract.rds"))
target_ids <- sort(unique(c(contract$reference$sample_id, contract$query$sample_id)))
reference_ids <- as.character(contract$reference$sample_id)
query_ids <- as.character(contract$query$sample_id)
if (length(intersect(reference_ids, query_ids))) stop("biology: fit/query overlap.")
metadata <- rbind(contract$reference, contract$query)
cohort_lookup <- unique(metadata[, c("cohort", "cohort_key")])
if (anyDuplicated(cohort_lookup$cohort)) stop("biology: ambiguous cohort identity.")
luckyBase::Plus.library("yaml")
anchor_config <- yaml::read_yaml(config_path)
signatures <- readRDS(sig_path)
anchors <- .biology_select_anchors(anchor_config, signatures)
anchor_names <- names(anchors)
required_genes <- sort(unique(unlist(anchors, use.names = FALSE)))
sample_key_hash <- digest::digest(paste(target_ids, collapse = "\n"), algo = "md5", serialize = FALSE)
source_hash <- digest::digest(file = full_path, algo = "md5")
signature_hash <- digest::digest(file = sig_path, algo = "md5")
config_hash <- digest::digest(file = config_path, algo = "md5")

atlas <- readRDS(full_path)
cohorts <- list(); coverage <- list(); missing <- list(); index <- 0L
for (tissue in names(atlas)) {
  for (cohort in names(atlas[[tissue]])) {
    raw <- atlas[[tissue]][[cohort]]
    if (is.list(raw) && !is.matrix(raw) && !is.data.frame(raw)) raw <- raw$expr
    ids <- colnames(raw)
    genes <- rownames(raw)
    if (is.null(ids) || is.null(genes)) next
    key <- cohort_lookup$cohort_key[match(cohort, cohort_lookup$cohort)]
    if (is.na(key)) next
    allowed <- metadata$sample_id[metadata$cohort_key == key]
    keep <- ids %in% allowed
    gene_idx <- match(required_genes, genes)
    if (!any(keep)) next
    # Subset before conversion to avoid materialising a second full atlas matrix.
    subset_expr <- as.matrix(raw[gene_idx[!is.na(gene_idx)], keep, drop = FALSE])
    rownames(subset_expr) <- required_genes[!is.na(gene_idx)]
    resolved_tissue <- sub("/.*$", "", key)
    index <- index + 1L
    cohorts[[key]] <- list(
      tissue = resolved_tissue, cohort = cohort, cohort_key = key,
      expression = subset_expr, sample_id = ids[keep], gene_id = rownames(subset_expr)
    )
    for (anchor in anchor_names) {
      found <- intersect(anchors[[anchor]], rownames(subset_expr))
      coverage[[length(coverage) + 1L]] <- data.frame(
        tissue = resolved_tissue, cohort = cohort, cohort_key = key, anchor = anchor,
        genes_required = length(anchors[[anchor]]), genes_found = length(found),
        coverage = length(found) / length(anchors[[anchor]]), sample_count = sum(keep),
        id_type = ifelse(any(grepl("^ENSG", found, ignore.case = TRUE)), "ENSEMBL", "SYMBOL"),
        status = ifelse(length(found) >= 2L, "estimable", "not_estimable"),
        reason = ifelse(length(found) >= 2L, NA_character_, "fewer_than_two_genes"),
        stringsAsFactors = FALSE
      )
      missing[[length(missing) + 1L]] <- data.frame(
        cohort_key = key, anchor = anchor,
        missing_genes = paste(setdiff(anchors[[anchor]], found), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
}
if (!length(cohorts)) stop("ablation-03 cache: no target samples were found in expression RDS.", call. = FALSE)
coverage <- do.call(rbind, coverage); missing <- do.call(rbind, missing)
cache <- list(
  schema_version = 2L, status = "complete", created_at = format(Sys.time(), tz = "UTC"),
  builder_md5 = digest::digest(file = .wf_path("01c-ablation03-prepare-biology.R"), algo = "md5"),
  source = list(path = normalizePath(full_path, winslash = "/"), md5 = source_hash),
  signature = list(path = normalizePath(sig_path, winslash = "/"), md5 = signature_hash,
                   config_path = normalizePath(config_path, winslash = "/"),
                   version = anchor_config$version, config_md5 = config_hash),
  sample_key_hash = sample_key_hash, sample_ids = target_ids,
  reference_sample_ids = reference_ids, query_sample_ids = query_ids,
  reference_cohorts = unique(contract$reference$cohort_key),
  sample_contract_md5 = digest::digest(file = file.path(result_dir, "sample-contract.rds"), algo = "md5"),
  duplicate_sample_ids = unique(unlist(lapply(cohorts, function(x) x$sample_id[duplicated(x$sample_id)]))),
  required_genes = required_genes, anchors = anchors, cohorts = cohorts,
  coverage = coverage, missing_genes = missing,
  preprocessing = list(method = "reference_cohort_global_gene_zscore_signature_mean",
                       id_conversion = "exact_match",
                       compression = "gzip", compression_level = 6L)
)
saveRDS(cache, cache_path, compress = "gzip")
rm(list = intersect(c("atlas", "signatures", "raw", "subset_expr", "cache"), ls()))
invisible(gc())
cat(sprintf("cache=%s cohorts=%d samples=%d genes=%d source_md5=%s\n",
            cache_path, index, length(target_ids), length(required_genes), source_hash))

# Structural anchors use exactly the same sample selection and extraction as the
# original stage 04, now performed before any analysis is launched.
structural <- .wf_read("01-representations", "structural-inputs.rds")
prepared <- structural$analysis$prepared
full_d1 <- structural$structural$full_d1
ablation_metadata <- structural$structural$ablation_metadata
biology_cache <- readRDS(cache_path)
source(.ablation03_path("04-ablation03-structural-reproducibility_functions.R"))
output_dir <- out_dir
reference_ids <- Reduce(intersect, list(
  rownames(prepared$reference_direct), rownames(full_d1)
))
external_ids <- Reduce(intersect, list(
  rownames(prepared$query_direct), rownames(full_d1)
))
if (length(intersect(reference_ids, external_ids)) > 0L) {
  stop("structural reproducibility: reference and external samples overlap.", call. = FALSE)
}

# Step 3: Extract the four anchors for every common Direct/d1 sample. The old
# cache supplies frozen signatures and provenance, not the restricted sample set.
cohort_lookup_rows <- unique(ablation_metadata[, c("cohort", "cohort_key")])
if (anyDuplicated(cohort_lookup_rows$cohort)) {
  stop("structural reproducibility: cohort-to-tissue lookup is ambiguous.", call. = FALSE)
}
cohort_key_lookup <- stats::setNames(
  cohort_lookup_rows$cohort_key,
  cohort_lookup_rows$cohort
)
anchor_sample_ids <- sort(unique(c(reference_ids, external_ids)))
anchor_cache_key <- digest::digest(list(
  schema_version = 2L,
  sample_ids = anchor_sample_ids,
  anchors = biology_cache$anchors,
  source_md5 = biology_cache$source$md5,
  cohort_key_lookup = cohort_key_lookup
), algo = "md5")
structural_anchor_path <- file.path(output_dir, "structural-anchor-cache.rds")
anchor_cache <- NULL
if (file.exists(structural_anchor_path)) {
  candidate <- readRDS(structural_anchor_path)
  if (identical(candidate$cache_key, anchor_cache_key)) anchor_cache <- candidate
}
if (is.null(anchor_cache)) {
  if (!file.exists(biology_cache$source$path)) {
    stop("structural reproducibility: complete expression atlas is unavailable.", call. = FALSE)
  }
  expression_atlas <- readRDS(biology_cache$source$path)
  anchor_cache <- .asr_extract_anchor_cache(
    expression_atlas,
    biology_cache$anchors,
    anchor_sample_ids,
    cohort_key_lookup = cohort_key_lookup
  )
  anchor_cache$schema_version <- 1L
  anchor_cache$status <- "complete"
  anchor_cache$cache_key <- anchor_cache_key
  anchor_cache$source <- biology_cache$source
  saveRDS(anchor_cache, structural_anchor_path, compress = "gzip")
  rm(expression_atlas)
  invisible(gc())
}

.wf_receipt("01-biology", "01c-ablation03-prepare-biology.R",
  inputs = c(.wf_output("01-data", "stage-receipt.rds"),
    file.path(result_dir, "stage-receipt.rds"), full_path, sig_path, config_path,
    .ablation03_path("03-ablation03-biology_functions.R"),
    .ablation03_path("04-ablation03-structural-reproducibility_functions.R")),
  outputs = c(cache_path, structural_anchor_path))
