# Regression cases for the manuscript audit. All inputs are synthetic.
source("R/ccs.R")
source("R/ablation.R")
source("test/ablation-03/03-ablation03-biology_functions.R")
source("test/ablation-03/04-ablation03-structural-reproducibility_functions.R")
failures <- character()
check <- function(name, code) {
  tryCatch({force(code); cat("PASS", name, "\n")}, error = function(e) {
    failures <<- c(failures, name)
    cat("FAIL", name, conditionMessage(e), "\n")
  })
}

check("main readout receives raw query in both arms", {
  ref <- matrix(10 + seq_len(24), 12, dimnames = list(paste0("r", 1:12), c("a", "b")))
  qry <- matrix(20 + seq_len(8), 4, dimnames = list(paste0("q", 1:4), c("a", "b")))
  rm <- data.frame(sample_id = rownames(ref), cohort = rep(c("r1", "r2"), each = 6), cancer_type = "A")
  qm <- data.frame(sample_id = rownames(qry), cohort = rep(c("q1", "q2"), each = 2), cancer_type = "A")
  qm$d1_provenance <- "external_frozen"
  prepared <- list(reference_direct = ref, query_direct = qry, reference_d1 = ref / 100,
    query_d1 = qry / 100, reference_metadata = rm, query_metadata = qm,
    selected_blocks = list(m1 = 1:2),
    query_views = list(cancer_retrieval = list(metadata = qm[1:2, ]), cancer_readout = list(metadata = qm[1:2, ]),
      anchor = list(metadata = qm), technical_excess = list(metadata = qm), learning_curve = list(metadata = qm)))
  env <- new.env(parent = .GlobalEnv)
  env$.ablation_prepare_representation_analysis <- function(...) list(prepared = prepared, anchor = "cancer_type")
  config <- .ablation_representation_default_params()
  config$geometry$search <- "exact"
  config$geometry$k <- 1L
  config$controls$null_rp <- config$controls$null_perm <- FALSE
  env$.ablation_resolve_representation_config <- function(...) config
  env$.ablation_native_geometry_cache_key <- function(...) "fixture"
  env$.ablation_read_native_geometry_cache <- function(...) list()
  retrieval_query_counts <- integer()
  env$.ablation_query_reference_retrieval <- function(reference, query, ...) {
    retrieval_query_counts <<- c(retrieval_query_counts, nrow(query))
    .ablation_query_reference_retrieval(reference, query, ...)
  }
  called <- 0L
  env$.ablation_linear_readout <- function(train, test, ...) {
    called <<- called + 1L
    expected <- if (called == 1L) qry[1:2, , drop = FALSE] else qry[1:2, , drop = FALSE] / 100
    stopifnot(isTRUE(all.equal(test, expected)))
    if (called == 2L) stop(structure(list(message = "captured both arms"), class = c("captured", "error", "condition")))
    list()
  }
  runner <- .ablation_run_representation
  environment(runner) <- env
  tryCatch(runner(methods::new("CCS"), NULL, NULL, tempdir(), list(), 42L, FALSE), captured = function(e) NULL)
  stopifnot(called == 2L)
  stopifnot(identical(retrieval_query_counts, c(2L, 2L, 4L, 4L)))
})

check("biology pairs queries before cohort means", {
  d <- do.call(rbind, lapply(1:3, function(i) data.frame(
    anchor = "a", representation = c("Direct-GSClassifier", "Cohort-d1", "Direct-GSClassifier"),
    query_sample = paste0(i, c("A", "A", "B")), query_cohort = paste0("c", i), utility = c(.8, .8, .2))))
  result <- .biology_paired_contrast(d, n_boot = 100L)
  stopifnot(abs(result$estimate) < 1e-12, result$query_count == 3L)
})

check("constant anchors cannot manufacture high/low states", {
  d <- data.frame(sample_id = paste0("s", 1:30), cohort_key = "T/c", anchor = "a", score = 1, gene_count = 3L)
  result <- tryCatch(.asr_assign_anchor_states(d), error = function(e) NULL)
  stopifnot(is.null(result) || nrow(result) == 0L)
})

check("bank design uses audited tissue identities", {
  manifest <- list(modules = data.frame(module_id = c("Undefined|c1", "Undefined|c2"),
    tissue = "Undefined", cohort = c("c1", "c2")))
  metadata <- data.frame(cohort = c("c1", "c2"), tissue = c("A", "B"))
  resolved <- .ablation_resolve_bank_tissues(manifest, metadata)
  stopifnot(identical(resolved$modules$tissue, c("A", "B")),
    identical(resolved$modules$module_id, manifest$modules$module_id))
})
if (length(failures)) stop(paste("Failed:", paste(failures, collapse = "; ")))
