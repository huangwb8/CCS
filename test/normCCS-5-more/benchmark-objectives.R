# Existing simple workflow; no cache, scheduler or real patient input.
# Rscript .../benchmark-objectives.R [output-directory] [--config=...] [--section=...]
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
option <- function(name, fallback) {
  values <- args[startsWith(args, paste0("--", name, "="))]
  if (length(values)) sub(paste0("^--", name, "="), "", values[1L]) else fallback
}
config <- normalizePath(option("config", "test/normCCS-5-more/objective-benchmark-config.json"), winslash = "/", mustWork = TRUE)
if (!startsWith(tolower(config), paste0(tolower(normalizePath(".", winslash = "/")), "/"))) stop("Configuration must be inside the project.")
section <- option("section", "all")
if (!section %in% c("all", "quality", "timing", "scaling")) stop("Unknown benchmark section.")
paths <- args[!startsWith(args, "--")]
out <- norm_output_dir(if (length(paths)) paths[1L] else "test/normCCS-5-more/products")
cfg <- jsonlite::fromJSON(config)
source_files <- c("DESCRIPTION", "R/phenotype.R", "R/subtype_norm_dp.R", "R/subtype_norm_objective.R",
                  "test/normCCS-5-more/benchmark-objectives.R",
                  "test/normCCS-5-more/00.Environment.R", "test/normCCS-5-more/fixtures.R",
                  "test/normCCS-5-more/reference/subtype-norm-0.8.5.R",
                  "test/normCCS-5-more/reference/subtype-norm-dp-0.8.5.R", config)
source_identity <- data.frame(path = source_files, md5 = unname(tools::md5sum(source_files)))
identity_file <- file.path(out, "objective-source-identity.csv")
config_file <- file.path(out, "objective-benchmark-config.json")
if (section != "all" && any(file.exists(file.path(out, paste0("objective-", c("quality", "timing", "scaling"), ".csv"))))) {
  if (!file.exists(identity_file) || !file.exists(config_file) ||
      !identical(read.csv(identity_file), source_identity) ||
      !isTRUE(all.equal(jsonlite::fromJSON(config_file), cfg))) {
    stop("Existing benchmark sections have different or missing provenance; rerun all sections together in a fresh directory.")
  }
}
devtools::load_all(quiet = TRUE)
legacy <- new.env(parent = asNamespace("CCS"))
sys.source("test/normCCS-5-more/reference/subtype-norm-0.8.5.R", legacy)
sys.source("test/normCCS-5-more/reference/subtype-norm-dp-0.8.5.R", legacy)
write_product <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE)
product <- paste0("objective-", if (section == "all") c("quality", "timing", "scaling") else section, ".csv")
if (!overwrite && any(file.exists(file.path(out, product)))) stop("Products exist; pass --overwrite explicitly.")

# Independent empirical-quantile integration for exhaustive reference values.
reference_scores <- function(auc, cohorts, objective) {
  if (objective == "median") return(apply(auc, 1L, median))
  real <- auc[, cohorts != "All", drop = FALSE]
  nc <- ncol(real)
  weight <- pmax(0, pmin(seq_len(nc) / nc, 0.25) - (seq_len(nc) - 1) / nc) / 0.25
  lower <- apply(real, 1L, function(x) sum(sort(x) * weight))
  switch(objective, cohort_median = apply(real, 1L, median), mean = rowMeans(real),
         lower_tail = lower, mean_lower_tail = (rowMeans(real) + lower) / 2)
}
best_score <- function(result, objective) {
  if (objective == "median") result$Summary$medianROCAUC[1L] else result$Summary$objectiveScore[1L]
}
describe <- function(result, df, k, method, objective, elapsed) {
  d <- result$SearchDiagnostics
  first <- result$Raw[result$Raw$Info == result$Summary$Info[1L], ]
  real <- first$ROCAUC[first$Cohort != "All"]
  grade <- k - 1L - findInterval(match(df$Subtype, d$ranked_subtypes) - 1L, d$selected_cutpoints)
  group_n <- table(factor(df$Cohort, levels = sort(unique(df$Cohort))), factor(grade, levels = 0:(k - 1L)))
  data.frame(M = length(unique(df$Subtype)), C = length(unique(df$Cohort)), K = k,
             method = method, objective = objective, objective_score = best_score(result, objective),
             cohort_mean = mean(real), cohort_median = median(real),
             cohort_lower_tail = reference_scores(matrix(first$ROCAUC, 1L), first$Cohort, "lower_tail"),
             cohort_auc = paste(sprintf("%s=%.12f", first$Cohort[first$Cohort != "All"], real), collapse = ";"),
             grade_counts = paste(as.integer(table(factor(grade, levels = 0:(k - 1L)))), collapse = ","),
             empty_cohort_grade_cells = sum(group_n == 0), minimum_cohort_grade_n = min(group_n),
             nSample = nrow(df), partitions = choose(d$M - 1L, k - 1L),
             candidates = d$n_unique_candidates, evaluations = d$n_scored_proposals,
             cutpoints = paste(d$selected_cutpoints, collapse = ","),
             mapping = result$Summary$Info[1L],
             elapsed_seconds = elapsed, interval_mib = d$interval_storage_bytes / 1024^2,
             initial_score = if (objective == "median") d$best_initial_median_auc else d$best_initial_objective_score,
             exact_global_optimum = d$exact_global_optimum,
             stringsAsFactors = FALSE)
}

if (section %in% c("all", "quality")) {
  rows <- list()
  for (nc in cfg$quality_c) for (seed in cfg$quality_seeds) for (k in cfg$quality_k) {
    df <- make_norm_fixture(cfg$quality_m, nc, seed, sparse = seed %% 2L == 0L)
    model <- CCS:::norm_prepare_counts(df, reference_rank(df, "NRR"))
    cuts <- enumerate_cuts(model$m, k)
    for (method in cfg$methods) {
      auc <- CCS:::norm_auc_batch(cuts, model, method)
      for (objective in cfg$objectives) {
        cat("Quality", nc, seed, k, method, objective, "\n")
        exact <- max(reference_scores(auc, model$cohorts, objective))
        elapsed <- system.time(result <- CCS:::subtypeNorm(df, k, method, "NRR", 1L, FALSE, "dp", objective))[["elapsed"]]
        row <- describe(result, df, k, method, objective, elapsed)
        row$seed <- seed
        row$exhaustive_best_objective <- exact
        row$gap <- exact - row$objective_score
        stopifnot(row$gap >= -1e-10, row$objective_score >= row$initial_score - 1e-12)
        rows[[length(rows) + 1L]] <- row
      }
    }
  }
  write_product(do.call(rbind, rows), "objective-quality.csv")
}

if (section %in% c("all", "timing")) {
  rows <- list()
  df <- make_norm_fixture(cfg$timing_m, cfg$timing_c, 307L, sparse = TRUE)
  for (k in cfg$timing_k) for (method in cfg$methods) for (rep in seq_len(cfg$timing_reps)) {
    # Reverse order on alternate repetitions to reduce a consistent warmup bias.
    variants <- c("old_dp", "new_dp", "old_bruteforce", "new_bruteforce")
    if (rep %% 2L == 0L) variants <- rev(variants)
    for (variant in variants) {
      cat("Timing", k, method, rep, variant, "\n")
      fn <- if (startsWith(variant, "old")) legacy$subtypeNorm else CCS:::subtypeNorm
      search <- if (endsWith(variant, "bruteforce")) "bruteforce" else "dp"
      gc()
      elapsed <- system.time(result <- suppressWarnings(fn(df, k, method, "NRR", 1L, FALSE, search)))[["elapsed"]]
      d <- result$SearchDiagnostics
      rows[[length(rows) + 1L]] <- data.frame(M = cfg$timing_m, C = cfg$timing_c, K = k,
        method = method, replicate = rep, variant = variant, elapsed_seconds = elapsed,
        candidates = nrow(result$Summary), evaluations = if (is.null(d)) NA_integer_ else d$n_scored_proposals,
        median_auc = result$Summary$medianROCAUC[1L], mapping = result$Summary$Info[1L])
    }
  }
  timings <- do.call(rbind, rows)
  for (part in split(timings, list(timings$K, timings$method, timings$replicate), drop = TRUE)) {
    for (search in c("dp", "bruteforce")) {
      old <- part[part$variant == paste0("old_", search), ]
      new <- part[part$variant == paste0("new_", search), ]
      stopifnot(abs(old$median_auc - new$median_auc) <= 1e-12,
                identical(old$mapping, new$mapping), identical(old$candidates, new$candidates),
                identical(old$evaluations, new$evaluations))
    }
  }
  write_product(timings, "objective-timing.csv")
}

if (section %in% c("all", "scaling")) {
  rows <- list()
  for (m in cfg$scaling_m) for (k in cfg$scaling_k) for (method in cfg$methods) {
    df <- make_norm_fixture(m, cfg$scaling_c, 509L, sparse = TRUE)
    for (objective in cfg$objectives) {
      cat("Scaling", m, k, method, objective, "\n")
      gc()
      elapsed <- system.time(result <- CCS:::subtypeNorm(df, k, method, "NRR", 1L, FALSE, "dp", objective))[["elapsed"]]
      row <- describe(result, df, k, method, objective, elapsed)
      row$seed <- 509L
      rows[[length(rows) + 1L]] <- row
    }
  }
  write_product(do.call(rbind, rows), "objective-scaling.csv")
}
if (!identical(unname(tools::md5sum(source_files)), source_identity$md5)) {
  stop("Benchmark sources changed during execution; rerun affected results before delivery.")
}
file.copy(config, config_file, overwrite = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "objective-session-info.txt"))
write_product(source_identity, "objective-source-identity.csv")
cat("Objective benchmark complete:", section, "\n")
