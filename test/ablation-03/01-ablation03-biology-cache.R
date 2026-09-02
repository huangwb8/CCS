#!/usr/bin/env Rscript
# Build the small, auditable biological-anchor cache used by 03-ablation-biology.R.

options(stringsAsFactors = FALSE)

.ablation03_find_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  candidates <- c(
    if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1L])) else character(),
    getwd()
  )
  for (candidate in unique(candidates)) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(candidate, "ablation-03.Rproj"))) return(candidate)
    nested <- file.path(candidate, "test", "ablation-03")
    if (file.exists(file.path(nested, "ablation-03.Rproj"))) {
      return(normalizePath(nested, winslash = "/"))
    }
  }
  stop("ablation-03: cannot locate ablation-03.Rproj.", call. = FALSE)
}

ablation_dir <- .ablation03_find_dir()
root <- normalizePath(file.path(ablation_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(ablation_dir, "tmp", "ablation-biology")
cache_path <- file.path(out_dir, "expression-anchor-cache.rds")
result_dir <- file.path(ablation_dir, "tmp", "ablation-experiment")
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
if (!file.exists(result_dir) || !file.exists(file.path(result_dir, "retrieval.rds"))) {
  stop("ablation-03 cache: retrieval.rds is missing; run the experiment first.", call. = FALSE)
}
if (!file.exists(full_path)) stop("ablation-03 cache: expression RDS not found: ", full_path, call. = FALSE)
if (!file.exists(sig_path)) stop("ablation-03 cache: signature RDS not found: ", sig_path, call. = FALSE)

retrieval <- readRDS(file.path(result_dir, "retrieval.rds"))
neighbours <- retrieval$neighbors[retrieval$neighbors$neighbor_rank <= 15, , drop = FALSE]
target_ids <- sort(unique(c(as.character(neighbours$query_sample), as.character(neighbours$reference_sample))))

# The YAML is intentionally the single frozen entry point.  Its structure is
# small and fixed, so parsing these fields avoids adding a runtime yaml import.
yaml <- readLines(config_path, warn = FALSE, encoding = "UTF-8")
required_yaml_tokens <- c(
  "proliferation:", "immune_tme:", "stromal_tme:", "ifn_il6:",
  "Tumor proliferation rate", "TME_A_Immune2", "TME_B_Stromal2",
  "IFN signaling", "IL6-JAK-STAT3 signaling"
)
if (!all(vapply(required_yaml_tokens, function(token) any(grepl(token, yaml, fixed = TRUE)), logical(1)))) {
  stop("ablation-03 cache: biological-anchors.yml does not match the frozen anchor contract.", call. = FALSE)
}
yaml_value <- function(pattern) {
  line <- grep(pattern, yaml, value = TRUE)
  if (!length(line)) return(NA_character_)
  sub("^[^:]+:\\s*['\"]?", "", sub("['\"]?\\s*$", "", line[1L]))
}
anchor_names <- c("proliferation", "immune_tme", "stromal_tme", "ifn_il6")
signatures <- readRDS(sig_path)
anchors <- list(
  proliferation = signatures[["Conserved-PanCan-TME-subtypes"]][["Tumor proliferation rate"]],
  immune_tme = signatures[["Zeng2021"]][["TME_A_Immune2"]],
  stromal_tme = signatures[["Zeng2021"]][["TME_B_Stromal2"]],
  ifn_il6 = c(signatures[["IFN-IL6"]][[1L]], signatures[["IFN-IL6"]][[2L]])
)
anchors <- lapply(anchors, function(x) sort(unique(as.character(unlist(x)))))
if (any(!nzchar(unlist(anchors)))) stop("ablation-03 cache: frozen signature contains blank genes.", call. = FALSE)
required_genes <- sort(unique(unlist(anchors, use.names = FALSE)))
sample_key_hash <- digest::digest(paste(target_ids, collapse = "\n"), algo = "md5", serialize = FALSE)
source_hash <- digest::digest(file = full_path, algo = "md5")
signature_hash <- digest::digest(list(config = yaml, anchors = anchors), algo = "md5")

atlas <- readRDS(full_path)
cohorts <- list(); coverage <- list(); missing <- list(); index <- 0L
for (tissue in names(atlas)) {
  for (cohort in names(atlas[[tissue]])) {
    raw <- atlas[[tissue]][[cohort]]
    if (is.list(raw) && !is.matrix(raw) && !is.data.frame(raw)) raw <- raw$expr
    ids <- colnames(raw)
    genes <- rownames(raw)
    if (is.null(ids) || is.null(genes)) next
    keep <- ids %in% target_ids
    gene_idx <- match(required_genes, genes)
    if (!any(keep)) next
    # Subset before conversion to avoid materialising a second full atlas matrix.
    subset_expr <- as.matrix(raw[gene_idx[!is.na(gene_idx)], keep, drop = FALSE])
    rownames(subset_expr) <- required_genes[!is.na(gene_idx)]
    key <- paste(tissue, cohort, sep = "/")
    index <- index + 1L
    cohorts[[key]] <- list(
      tissue = tissue, cohort = cohort, cohort_key = key,
      expression = subset_expr, sample_id = ids[keep], gene_id = rownames(subset_expr)
    )
    for (anchor in anchor_names) {
      found <- intersect(anchors[[anchor]], rownames(subset_expr))
      coverage[[length(coverage) + 1L]] <- data.frame(
        tissue = tissue, cohort = cohort, cohort_key = key, anchor = anchor,
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
  schema_version = 1L, status = "complete", created_at = format(Sys.time(), tz = "UTC"),
  source = list(path = normalizePath(full_path, winslash = "/"), md5 = source_hash),
  signature = list(path = normalizePath(sig_path, winslash = "/"), md5 = signature_hash,
                   config_path = normalizePath(config_path, winslash = "/"),
                   version = yaml_value("^version:")),
  sample_key_hash = sample_key_hash, sample_ids = target_ids,
  duplicate_sample_ids = unique(unlist(lapply(cohorts, function(x) x$sample_id[duplicated(x$sample_id)]))),
  required_genes = required_genes, anchors = anchors, cohorts = cohorts,
  coverage = coverage, missing_genes = missing,
  preprocessing = list(method = "cohort_zscore_signature_mean", id_conversion = "exact_match",
                       compression = "gzip", compression_level = 6L)
)
saveRDS(cache, cache_path, compress = "gzip")
rm(list = intersect(c("atlas", "signatures", "raw", "subset_expr", "cache"), ls()))
invisible(gc())
cat(sprintf("cache=%s cohorts=%d samples=%d genes=%d source_md5=%s\n",
            cache_path, index, length(target_ids), length(required_genes), source_hash))
