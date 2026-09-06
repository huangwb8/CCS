#!/usr/bin/env Rscript
# Biological-anchor readout for ablation-03. External RDS inputs are read-only.

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
    if (file.exists(file.path(nested, "ablation-03.Rproj"))) return(normalizePath(nested, winslash = "/"))
  }
  stop("ablation-03: cannot locate ablation-03.Rproj.", call. = FALSE)
}
.ablation03_dir <- .ablation03_find_dir()
root <- normalizePath(file.path(.ablation03_dir, "..", ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(.ablation03_dir, "tmp", "ablation-biology")
fig_dir <- file.path(.ablation03_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

result_dir <- file.path(.ablation03_dir, "tmp", "ablation-experiment")
manifest <- readRDS(file.path(result_dir, "manifest.rds"))
retrieval <- readRDS(file.path(result_dir, "retrieval.rds"))
neighbours <- retrieval$neighbors[retrieval$neighbors$neighbor_rank <= 15, , drop = FALSE]
target_ids <- sort(unique(c(as.character(neighbours$query_sample), as.character(neighbours$reference_sample))))
full_path <- Sys.getenv("CCS_FULL_EXPRESSION_RDS", unset = "E:/Sync/@Analysis/PanCan_Data/Level 1/PanCan_CancerSample_DataListForCCS_GEO+cBioPortal+UCXCXenav20240809.rds")
sig_path <- Sys.getenv("CCS_GENE_SIGNATURE_RDS", unset = "E:/RCloud/database/Signature/report/GeneSignature-HWB.rds")
cache_path <- Sys.getenv("CCS_BIOLOGY_CACHE_RDS", unset = file.path(out_dir, "expression-anchor-cache.rds"))
if (!file.exists(cache_path)) {
  stop("ablation-03 biology: expression-anchor-cache.rds is missing; run 01-ablation03-biology-cache.R first.", call. = FALSE)
}
cache <- readRDS(cache_path)
if (!identical(cache$schema_version, 1L) || !identical(cache$status, "complete")) {
  stop("ablation-03 biology: unsupported or incomplete cache schema.", call. = FALSE)
}
expected_sample_hash <- digest::digest(paste(target_ids, collapse = "\n"), algo = "md5", serialize = FALSE)
if (!identical(cache$sample_key_hash, expected_sample_hash)) {
  stop("ablation-03 biology: cache sample-key hash mismatch; rebuild the cache.", call. = FALSE)
}
if (!identical(normalizePath(sig_path, winslash = "/", mustWork = FALSE), cache$signature$path)) {
  stop("ablation-03 biology: signature source mismatch; rebuild the cache.", call. = FALSE)
}
if (file.exists(full_path)) {
  source_hash <- digest::digest(file = full_path, algo = "md5")
  if (!identical(source_hash, cache$source$md5)) {
    stop("ablation-03 biology: expression source hash mismatch; rebuild the cache.", call. = FALSE)
  }
}
source(file.path(.ablation03_dir, "03-ablation-biology_functions.R"))
anchors <- cache$anchors
coverage <- cache$coverage
coverage$external_query_cohort <- coverage$cohort_key %in% manifest$external_cohorts
# Fit one gene-wise reference scale across all non-query cohorts, then apply
# that fixed transform to both query and reference samples.  Per-cohort z-scoring
# would put every cohort in a different coordinate system and make cross-cohort
# absolute deltas uninterpretable.
reference_keys <- setdiff(names(cache$cohorts), manifest$external_cohorts)
if (!length(reference_keys)) {
  stop("ablation-03 biology: no reference cohorts available for global scaling.", call. = FALSE)
}
global_stats <- lapply(sort(unique(unlist(anchors, use.names = FALSE))), function(gene) {
  values <- unlist(lapply(cache$cohorts[reference_keys], function(cohort) {
    mat <- cohort$expression
    if (!gene %in% rownames(mat)) return(numeric())
    as.numeric(mat[gene, , drop = TRUE])
  }), use.names = FALSE)
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(c(mean = NA_real_, sd = NA_real_, n = length(values)))
  c(mean = mean(values), sd = stats::sd(values), n = length(values))
})
global_stats <- do.call(rbind, global_stats)
rownames(global_stats) <- sort(unique(unlist(anchors, use.names = FALSE)))
global_stats <- as.data.frame(global_stats, stringsAsFactors = FALSE)
global_stats$gene_id <- rownames(global_stats)
global_stats <- global_stats[is.finite(global_stats$mean) & is.finite(global_stats$sd) &
  global_stats$sd > 0, , drop = FALSE]
score_rows <- list(); ii <- 0L
for (cohort in cache$cohorts) {
  mat <- cohort$expression
  ids <- cohort$sample_id
  for (anchor in names(anchors)) {
    genes <- intersect(anchors[[anchor]], intersect(rownames(mat), global_stats$gene_id))
    if (length(genes) < 2L) next
    stats <- global_stats[match(genes, global_stats$gene_id), , drop = FALSE]
    values <- mat[genes, , drop = FALSE]
    z <- sweep(values, 1L, stats$mean, FUN = "-")
    z <- sweep(z, 1L, stats$sd, FUN = "/")
    score <- colMeans(z, na.rm = TRUE)
    ii <- ii + 1L
    score_rows[[ii]] <- data.frame(sample_id = ids, anchor = anchor,
      score = as.numeric(score), stringsAsFactors = FALSE)
  }
}
if (!length(score_rows)) stop("ablation-03 biology: no anchor has at least two cached genes.", call. = FALSE)
scores <- do.call(rbind, score_rows)
scores <- scores[is.finite(scores$score), , drop = FALSE]

boot_ci <- function(x, seed = 20260830L, B = 500L) {
  x <- x[is.finite(x)]; if (!length(x)) return(c(mean = NA_real_, low = NA_real_, high = NA_real_))
  set.seed(seed); draws <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), low = unname(quantile(draws, .025)), high = unname(quantile(draws, .975)))
}
utility_rows <- list(); kk <- 0L
per_query_rows <- list(); pp <- 0L
missing_rows <- list(); mm <- 0L
for (anchor in names(anchors)) {
  sc <- scores[scores$anchor == anchor, c("sample_id", "score")]
  q <- merge(neighbours, sc, by.x = "query_sample", by.y = "sample_id")
  names(q)[names(q) == "score"] <- "query_score"
  q <- merge(q, sc, by.x = "reference_sample", by.y = "sample_id")
  names(q)[names(q) == "score"] <- "reference_score"
  q$abs_delta <- abs(q$query_score - q$reference_score)
  q$utility <- exp(-q$abs_delta)
  per_query <- aggregate(cbind(abs_delta, utility) ~ representation + query_sample + query_cohort,
    data = q, FUN = mean)
  per_query$anchor <- anchor
  pp <- pp + 1L
  per_query_rows[[pp]] <- per_query[, c(
    "anchor", "representation", "query_sample", "query_cohort", "utility"
  )]
  for (rep in unique(per_query$representation)) {
    d <- per_query[per_query$representation == rep, , drop = FALSE]
    ci_u <- boot_ci(d$utility, seed = 20260830L + kk)
    ci_d <- boot_ci(d$abs_delta, seed = 20300830L + kk)
    kk <- kk + 1L
    expected_pairs <- sum(neighbours$representation == rep)
    valid_pairs <- sum(q$representation == rep)
    utility_rows[[kk]] <- data.frame(anchor = anchor, representation = rep,
      query_count = nrow(d), neighbour_pairs = valid_pairs,
      expected_neighbour_pairs = expected_pairs,
      missing_score_pairs = expected_pairs - valid_pairs,
      mean_abs_delta = ci_d[["mean"]], abs_delta_ci_low = ci_d[["low"]],
      abs_delta_ci_high = ci_d[["high"]], utility = ci_u[["mean"]],
      utility_ci_low = ci_u[["low"]], utility_ci_high = ci_u[["high"]])
    mm <- mm + 1L
    missing_rows[[mm]] <- data.frame(
      anchor = anchor, representation = rep,
      expected_neighbour_pairs = expected_pairs,
      valid_neighbour_pairs = valid_pairs,
      missing_score_pairs = expected_pairs - valid_pairs,
      stringsAsFactors = FALSE
    )
  }
}
utility <- do.call(rbind, utility_rows)
per_query_utility <- do.call(rbind, per_query_rows)
missing_pairs <- do.call(rbind, missing_rows)
anchor_inference <- .biology_paired_contrast(
  per_query_utility,
  n_boot = 2000L,
  seed = 20260830L
)
direct <- utility[utility$representation == "Direct-GSClassifier", , drop = FALSE]
d1 <- utility[utility$representation == "Cohort-d1", , drop = FALSE]
contrasts <- merge(d1, direct, by = "anchor", suffixes = c("_d1", "_direct"))
contrasts <- contrasts[, c("anchor", "mean_abs_delta_d1", "mean_abs_delta_direct", "utility_d1", "utility_direct")]
contrasts$d1_minus_direct_abs_delta <- contrasts$mean_abs_delta_d1 - contrasts$mean_abs_delta_direct
contrasts$d1_minus_direct_utility <- contrasts$utility_d1 - contrasts$utility_direct
contrasts$interpretation <- ifelse(contrasts$d1_minus_direct_utility > 0, "d1 higher utility", "Direct higher utility")

write.csv(coverage, file.path(out_dir, "anchor_coverage.csv"), row.names = FALSE)
write.csv(utility, file.path(out_dir, "anchor_utility.csv"), row.names = FALSE)
write.csv(missing_pairs, file.path(out_dir, "anchor_missing_pairs.csv"), row.names = FALSE)
write.csv(contrasts, file.path(out_dir, "anchor_contrasts.csv"), row.names = FALSE)
write.csv(anchor_inference, file.path(out_dir, "anchor_inference.csv"), row.names = FALSE)
saveRDS(list(anchors = anchors, coverage = coverage, utility = utility, contrasts = contrasts,
             inference = anchor_inference,
             retrieval_rows_top15 = nrow(neighbours), source_signature = sig_path,
             cache_path = cache_path, cache_schema_version = cache$schema_version,
             cache_source_md5 = cache$source$md5, cache_sample_key_hash = cache$sample_key_hash,
             scaling = list(method = "reference_cohort_global_gene_zscore",
                            reference_cohort_count = length(reference_keys),
                            reference_cohorts = reference_keys,
                            genes_with_finite_scale = nrow(global_stats))),
        file.path(out_dir, "ablation03-biology.rds"))

# Figures are rendered centrally by 02-ablation03-experiment.Rmd so that the
# HTML, PDF, and preview outputs always share one current data source.  This
# stage writes tabular/RDS products only and does not recreate deprecated
# biological-anchor figures under names that could be mistaken for report plots.
cat(sprintf("anchors=%d coverage_rows=%d utility_rows=%d output=%s\n", length(anchors), nrow(coverage), nrow(utility), out_dir))
