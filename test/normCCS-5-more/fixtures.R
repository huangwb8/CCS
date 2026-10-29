# Only synthetic data. Restore the caller's RNG, including the no-seed case.
make_norm_fixture <- function(m = 10L, cohorts = 3L, seed = 101L,
                              sparse = FALSE, conflict = TRUE, per_state = 12L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  out <- lapply(seq_len(cohorts), function(c) {
    states <- seq_len(m)
    if (sparse && c > 1L) states <- states[states %% (c + 2L) != 0L]
    n <- sample(seq.int(max(2L, per_state - 3L), per_state + 4L), length(states), replace = TRUE)
    trend <- seq(-1.8, 1.8, length.out = m)
    if (conflict) trend <- trend + stats::rnorm(m, sd = 1.3)
    if (conflict && c %% 3L == 0L) trend <- -trend
    probability <- stats::plogis(trend[states] + (c - 1) / cohorts - 0.5)
    positive <- stats::rbinom(length(states), n, probability)
    # Guarantee both outcomes per cohort, but allow pure and sparse states.
    positive[1L] <- 0L
    positive[length(positive)] <- n[length(n)]
    do.call(rbind, lapply(seq_along(states), function(i) {
      data.frame(Cohort = sprintf("C%02d", c),
                 SampleIDs = sprintf("synthetic_%02d_%03d_%03d", c, states[i], seq_len(n[i])),
                 tumor_type = c("Lung", "Skin", "Gastric")[(c - 1L) %% 3L + 1L],
                 Subtype = sprintf("S%03d", states[i]),
                 response = c(rep(1L, positive[i]), rep(0L, n[i] - positive[i])),
                 stringsAsFactors = FALSE)
    }))
  })
  out <- do.call(rbind, out)
  # Row order matters for OptimizeRank ties; shuffle deterministically.
  out <- out[sample.int(nrow(out)), ]
  rownames(out) <- NULL
  stopifnot(!anyDuplicated(out$SampleIDs), all(out$response %in% 0:1), !anyNA(out))
  out
}

enumerate_cuts <- function(m, k) t(utils::combn(seq_len(m - 1L), k - 1L))

norm_output_dir <- function(path) {
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)) path <- file.path(root, path)
  path <- file.path(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), basename(path))
  candidate <- normalizePath(path, winslash = "/", mustWork = FALSE)
  candidate <- gsub("\\\\", "/", candidate)
  if (!startsWith(tolower(candidate), paste0(tolower(root), "/"))) {
    stop("Output directory must be inside the CCS project.")
  }
  relative <- substring(candidate, nchar(root) + 2L)
  if (!grepl("^(test/normCCS-5-more/products|tmp/tests/[^/]+/(products|reports))$", relative)) {
    stop("Use test/normCCS-5-more/products or an isolated tmp/tests/<run>/{products,reports} directory.")
  }
  dir.create(candidate, recursive = TRUE, showWarnings = FALSE)
  resolved <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  if (!startsWith(tolower(resolved), paste0(tolower(root), "/"))) stop("Output path escapes the project.")
  resolved
}

reference_rank <- function(df, marker = "RR") {
  counts <- plyr::ddply(df, c("Cohort", "Subtype"), function(x) {
    data.frame(size = nrow(x), nResponse = sum(x$response))
  })
  counts$response_rate <- counts$nResponse / counts$size
  cohort_rate <- tapply(counts$nResponse, counts$Cohort, sum) /
    tapply(counts$size, counts$Cohort, sum)
  counts$normalized_response_rate <- (counts$response_rate - cohort_rate[counts$Cohort]) / cohort_rate[counts$Cohort]
  ranked <- plyr::ddply(counts, "Subtype", function(x) {
    data.frame(medianNRR = stats::median(x$normalized_response_rate),
               medianRR = stats::median(x$response_rate))
  })
  switch(marker,
    RR = dplyr::arrange(ranked, dplyr::desc(medianRR)),
    NRR = dplyr::arrange(ranked, dplyr::desc(medianNRR)),
    `RR+NRR` = dplyr::arrange(ranked, dplyr::desc(medianRR), dplyr::desc(medianNRR)),
    `NRR+RR` = dplyr::arrange(ranked, dplyr::desc(medianNRR), dplyr::desc(medianRR)))
}

