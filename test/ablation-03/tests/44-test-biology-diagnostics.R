# Scientific invariants for continuous biology diagnostics; no real identities.
source(file.path("R", "ablation.R"))
source(file.path("R", "ablation_biology.R"))
config <- yaml::read_yaml("test/ablation-03/config/biology-diagnostics.yml")
raw_config <- config
raw_config$score_definitions <- "common_raw"
raw_config$primary_score_definition <- "common_raw"
raw_config$low_coverage_exploration <- FALSE
# Retain raw-only compatibility alongside the restored rank implementation.
config$score_definitions <- c("common_raw", "common_rank")
config$primary_score_definition <- "common_rank"
config$min_background_genes <- 5000L
set.seed(config$seed)
genes <- paste0("g", seq_len(5002L))
metadata <- function(prefix, cohorts, n) {
  x <- data.frame(sample_id = paste0(prefix, seq_len(cohorts * n)),
    cohort = rep(paste0(prefix, "c", seq_len(cohorts)), each = n),
    cancer_type = "A", tissue = "A", metadata_status = "confirmed",
    evidence_basis = "synthetic_fixture", assay_type = "RNAseq",
    source_system = rep(rep(c("source1", "source2", "source3"), length.out = cohorts), each = n),
    platform_id = "unknown", d1_provenance = if (prefix == "r") "reference" else "external_frozen")
  x$cohort_key <- paste("A", x$cohort, sep = "/")
  x
}
ref <- metadata("r", 6L, 24L)
query <- metadata("q", 4L, 24L)
ids <- c(ref$sample_id, query$sample_id)
expression <- matrix(round(stats::rnorm(length(genes) * length(ids)), 1L),
  nrow = length(genes), dimnames = list(genes, ids))
anchors <- stats::setNames(lapply(seq_along(config$primary_anchors), function(i)
  genes[(i - 1L) * 12L + seq_len(12L)]), config$primary_anchors)
signatures <- list(`IFN-IL6` = stats::setNames(list(anchors$ifn, anchors$il6),
  c("IFNγ signaling", "IL6-JAK-STAT3 signaling")))
mat <- matrix(stats::rnorm(length(ids) * 12L), length(ids),
  dimnames = list(ids, paste0("f", 1:12)))
blocks <- list(`A|r1` = 1:4, `A|r2` = 5:8, `B|r3` = 9:12)
prepared <- list(reference_direct = mat[ref$sample_id, ], query_direct = mat[query$sample_id, ],
  reference_d1 = stats::plogis(mat[ref$sample_id, ]),
  query_d1 = stats::plogis(mat[query$sample_id, ]),
  reference_metadata = ref, query_metadata = query, selected_blocks = blocks,
  module_manifest = list(modules = data.frame(module_id = names(blocks), tissue = c("A", "A", "B"))),
  feature_manifest = list(feature_manifest = data.frame(
    feature = c("g1", "g2:g3", "s1s2"), feature_type = c("single_bin", "gene_pair", "set_pair"))),
  endpoint_eligibility = data.frame(cohort_key = unique(query$cohort_key),
    endpoint = "cancer_retrieval", qualification_status = "estimable"))
cohort_metadata <- rbind(ref, query)
atlas <- list(A = stats::setNames(lapply(split(ids, cohort_metadata$cohort),
  function(x) list(expr = expression[, x, drop = FALSE])), unique(cohort_metadata$cohort)))
# Rebuild named ordering explicitly to avoid split's alphabetical ordering.
atlas$A <- lapply(split(ids, cohort_metadata$cohort), function(x) list(expr = expression[, x, drop = FALSE]))
cache <- list(anchors = anchors)
audit <- .ablation_bio_audit(atlas, cache, prepared, list(genes[1:20]), signatures, config)
stopifnot(length(audit$background) == 5002L, audit$rank_status == "estimable",
  !audit$coverage$estimable[audit$coverage$anchor == "proliferation_disjoint"],
  audit$input_support$complete)
# Sample names in a cohort outside the contract cannot identify its source.
foreign <- atlas
foreign$A$outside_contract <- list(expr = expression[, 1L, drop = FALSE])
foreign_audit <- .ablation_bio_audit(foreign, cache, prepared, list(genes[1:20]),
  signatures, config)
stopifnot(identical(audit, foreign_audit))
changed <- atlas
for (i in seq_along(changed$A)) changed$A[[i]]$expr <- changed$A[[i]]$expr * 7 + 300
other <- .ablation_bio_audit(changed, cache, prepared, list(genes[1:20]), signatures, config)
stopifnot(isTRUE(all.equal(audit$rank_scores, other$rank_scores, tolerance = 1e-12)))
expected <- mean(((rank(expression[, 1L], ties.method = "average") - .5) / nrow(expression))[1:12])
stopifnot(abs(audit$rank_scores[1L, "proliferation"] - expected) < 1e-12)
missing <- atlas
missing$A[[1L]]$expr <- missing$A[[1L]]$expr[-seq_len(7L), ]
limited <- .ablation_bio_audit(missing, cache, prepared, list(genes[1:20]), signatures, config)
stopifnot(!limited$coverage$estimable[limited$coverage$anchor == "proliferation"])
small <- atlas
for (i in seq_along(small$A)) small$A[[i]]$expr <- small$A[[i]]$expr[seq_len(100L), ]
stopifnot(.ablation_bio_audit(small, cache, prepared, list(genes[1:20]), signatures,
  config)$rank_status == "background_below_5000_genes")
# Raw targets need only signature genes; a small background cannot gate them.
raw_audit <- .ablation_bio_audit(small, cache, prepared, list(genes[1:20]),
  signatures, raw_config)
raw_scores <- .ablation_bio_scores(raw_audit)
stopifnot(raw_audit$rank_status == "disabled_by_design", is.null(raw_audit$rank_scores),
  identical(names(raw_scores$values), "common_raw"),
  identical(names(raw_scores$unscaled), "common_raw"),
  all(is.finite(raw_scores$values$common_raw[, raw_config$primary_anchors])),
  isTRUE(all.equal(raw_scores$values$common_raw,
    .ablation_bio_scores(.ablation_bio_audit(atlas, cache, prepared,
      list(genes[1:20]), signatures, raw_config))$values$common_raw)))
overlap <- prepared
overlap$query_metadata$sample_id[1L] <- ref$sample_id[1L]
stopifnot(inherits(try(.ablation_bio_audit(atlas, cache, overlap, list(genes[1:20]),
  signatures, config), silent = TRUE), "try-error"))
eligibility <- .ablation_bio_eligibility(prepared, config, "platform_id")
stopifnot(!any(eligibility$queries$eligible))
restricted <- prepared
restricted$query_metadata$cancer_type[1:24] <- "unsupported"
stopifnot(!any(.ablation_bio_eligibility(restricted, config)$queries$eligible[1:24]))
scores <- .ablation_bio_scores(audit)
constant <- audit
constant$raw_signature[,] <- 1
stopifnot(all(is.na(.ablation_bio_scores(constant)$values$common_raw)))
retrieval <- list()
for (pool in c("all", "same_cancer")) for (arm in c("Direct", "d1")) {
  retrieval[[paste(pool, arm)]] <- .ablation_bio_retrieve(prepared, config, pool, arm = arm)
}
n <- .ablation_bio_bind(lapply(retrieval, `[[`, "neighbors"))
stopifnot(all(table(interaction(n$pool, n$arm, n$query_sample)) == config$k),
  all(n$query_cohort != n$reference_cohort))
readout <- .ablation_bio_readout(prepared, audit, scores, config)
stopifnot(!length(intersect(readout$training_ids, readout$query_ids)),
  length(unique(readout$folds$fold)) == config$readout$folds)
result <- .ablation_bio_inference(prepared, audit, scores, retrieval, readout, config)
stopifnot(nrow(result$inference) > 0L, nrow(result$cohort) > 0L,
  all(result$inference$valid_bootstrap[result$inference$status == "estimable"] >= 1900L))
low <- small
for (i in seq_along(low$A)) low$A[[i]]$expr <-
  low$A[[i]]$expr[setdiff(rownames(low$A[[i]]$expr), c(genes[6:12], genes[17:24], genes[26:36])), ]
# Background and signature gaps affect strict targets, while residual proxies
# retain fixed gene lists and fit reference scales using finite values only.
low$A[[ref$cohort[1L]]]$expr[genes[1L], 1L] <- NA_real_
low_audit <- .ablation_bio_audit(low, cache, prepared, list(genes[1:20]), signatures, config)
low_scores <- .ablation_bio_scores(low_audit)
stopifnot(all(is.na(low_scores$values$common_rank)),
  all(is.na(low_scores$values$common_raw[, c("proliferation", "immune_tme", "stromal_tme")])),
  is.na(low_scores$exploratory$values$common_rank[ref$sample_id[1L], "proliferation"]),
  sum(is.finite(low_scores$exploratory$values$common_rank[, "proliferation"])) == length(ids) - 1L,
  is.na(low_scores$exploratory$values$common_raw[ref$sample_id[1L], "proliferation"]),
  sum(is.finite(low_scores$exploratory$values$common_raw[, "proliferation"])) == length(ids) - 1L,
  low_scores$exploratory$gene_scales$reference_valid_count[1L] == nrow(ref) - 1L)
low_result <- .ablation_bio_exploratory(prepared, low_audit,
  low_scores$exploratory, retrieval, config)
single <- low_result$inference$anchor == "stromal_tme"
stopifnot(all(is.na(low_result$inference$p_value)), all(is.na(low_result$inference$q_value)),
  all(is.na(low_result$inference$ci_low[single])),
  all(low_result$inference$measurement[single] == "single_gene_description"),
  any(low_result$inference$status == "exploratory"),
  all(low_result$availability$sample_count == low_result$availability$valid_count +
    low_result$availability$missing_count))
shift_query <- low_audit
shift_query$raw_signature[, query$sample_id] <- shift_query$raw_signature[, query$sample_id] * 99
shift_query$exploratory_rank[query$sample_id, ] <- .1
stopifnot(identical(low_scores$exploratory$scales,
  .ablation_bio_exploratory_scores(shift_query)$scales),
  identical(low_scores$exploratory$gene_scales,
    .ablation_bio_exploratory_scores(shift_query)$gene_scales))
empty_proxy <- low_scores$exploratory
for (definition in names(empty_proxy$values)) empty_proxy$values[[definition]][, ] <- NA_real_
empty_result <- .ablation_bio_exploratory(prepared, low_audit, empty_proxy, retrieval, config)
stopifnot(all(is.na(empty_result$inference$estimate)),
  all(empty_result$inference$status == "not_estimable"))
stopifnot(all(unique(query$cohort_key) %in% result$experiment_qualification$cohort))
raw_retrieval <- retrieval
raw_retrieval$distance <- .ablation_bio_retrieve(prepared, raw_config,
  "same_cancer", "module_unscaled", "d1")
raw_retrieval$technical_direct <- .ablation_bio_retrieve(prepared, raw_config,
  "same_cancer", "baseline", "Direct", "assay_type")
raw_retrieval$technical_d1 <- .ablation_bio_retrieve(prepared, raw_config,
  "same_cancer", "baseline", "d1", "assay_type")
raw_readout <- .ablation_bio_readout(prepared, raw_audit, raw_scores, raw_config)
raw_result <- .ablation_bio_inference(prepared, raw_audit, raw_scores,
  raw_retrieval, raw_readout, raw_config)
raw_grid <- raw_result$inference[grepl("^grid:", raw_result$inference$comparison), ]
stopifnot(nrow(raw_grid) == 2L * length(config$primary_anchors), all(raw_grid$status == "estimable"),
  all(raw_grid$query_count == nrow(query)),
  !any(grepl("common_rank|rank_pool|all_score|same_score", raw_result$inference$comparison)),
  !any(grepl("common_rank", raw_result$branch_status$experiment)),
  any(grepl("^distance_change:", raw_result$inference$comparison)),
  any(grepl("^technical_change:", raw_result$inference$comparison)),
  all(raw_result$inference$status[grepl("^(distance_change|technical_change):",
    raw_result$inference$comparison)] == "estimable"))
sparse_validation <- retrieval
sparse_validation[[1L]]$validation_lists <- c(list(NULL), sparse_validation[[1L]]$validation_lists)
sparse_result <- .ablation_bio_inference(prepared, audit, scores, sparse_validation, readout, config)
stopifnot(identical(sparse_result$inference, result$inference))
small_query_scores <- scores
small_ids <- query$sample_id[query$cohort_key == query$cohort_key[1L]][-seq_len(5L)]
for (name in names(small_query_scores$values)) small_query_scores$values[[name]][small_ids, ] <- NA_real_
descriptive <- .ablation_bio_inference(prepared, audit, small_query_scores, retrieval, readout, config)
small_grid <- descriptive$cohort[grepl("^grid:", descriptive$cohort$comparison) &
  descriptive$cohort$cohort == query$cohort_key[1L], ]
stopifnot(nrow(small_grid) > 0L, all(small_grid$n == 5L), !any(small_grid$eligible_inference))
for (limited_scores in list(.ablation_bio_scores(constant), .ablation_bio_scores(limited))) {
  unavailable <- .ablation_bio_inference(prepared, audit, limited_scores, retrieval, readout, config)
  stopifnot(nrow(unavailable$inference) > 0L,
    any(unavailable$inference$status == "not_estimable"))
}
all_missing <- scores
for (name in names(all_missing$values)) all_missing$values[[name]][, ] <- NA_real_
unavailable <- .ablation_bio_inference(prepared, audit, all_missing, retrieval, readout, config)
stopifnot(all(unavailable$inference$status == "not_estimable" | grepl("^readout:",
  unavailable$inference$comparison)))
values <- cbind(a = seq_len(4), b = seq_len(4) * 2)
inference <- .ablation_bio_infer_matrix(values, config, "paired-test", matrix(24L, 4, 2))
stopifnot(identical(dim(inference$draws), c(2000L, 4L)),
  abs(inference$summary$ci_low[2L] - 2 * inference$summary$ci_low[1L]) < 1e-12)
same_indices <- .ablation_bio_infer_matrix(values, config, "different-comparison",
  matrix(24L, 4, 2))
stopifnot(identical(inference$draws, same_indices$draws))
values[1L, 2L] <- NA_real_
missing_correlation <- .ablation_bio_infer_matrix(values, config, "missing-correlation",
  matrix(24L, 4, 2), test = FALSE)
stopifnot(missing_correlation$summary$cohort_count[2L] == 3L,
  missing_correlation$summary$valid_bootstrap[2L] == 2000L)
unknown_cancer <- prepared
unknown_cancer$query_metadata$cancer_type[1:24] <- NA_character_
retained <- .ablation_bio_utility(n[n$pool == "all", ], scores, unknown_cancer, config)
stopifnot(all(query$sample_id[1:24] %in% retained$query_sample))
without_support <- .ablation_bio_input_genes(prepared$feature_manifest, list())
stopifnot(!without_support$complete)
cat("continuous biology diagnostic scientific invariants passed\n")
