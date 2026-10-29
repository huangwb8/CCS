# Sequential, explicit benchmark entry. No checkpoint or custom scheduler.
# Rscript test/normCCS-5-more/benchmark.R [output-directory]
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- args[args != "--overwrite"]
config_arg <- args[grepl("^--config=", args)]
args <- args[!grepl("^--config=", args)]
config_file <- if (length(config_arg)) sub("^--config=", "", config_arg[1L]) else
  "test/normCCS-5-more/benchmark-config.json"
config_file <- normalizePath(config_file, winslash = "/", mustWork = TRUE)
if (!startsWith(tolower(config_file), paste0(tolower(normalizePath(".", winslash = "/")), "/"))) {
  stop("Configuration must be inside the project.")
}
cfg <- jsonlite::fromJSON(config_file)
output_dir <- if (length(args)) args[1L] else "test/normCCS-5-more/products"
output_dir <- norm_output_dir(output_dir)
product_names <- c("quality.csv", "timings.csv", "scaling.csv", "group-support.csv", "session-info.txt", "benchmark-config.json")
if (!overwrite && any(file.exists(file.path(output_dir, product_names)))) {
  stop("Benchmark products exist; pass --overwrite explicitly to replace them.")
}
devtools::load_all(quiet = TRUE)
jsonlite::write_json(cfg, file.path(output_dir, "benchmark-config.json"), auto_unbox = TRUE, pretty = TRUE)

# Quality against complete enumeration using the independently validated AUC
# kernel. Timing comparisons below still run the original 0.8.4 code itself.
quality <- list()
for (seed in cfg$quality_seeds) for (k in cfg$quality_k) for (method in c("Raw", "OptimizeRank")) {
  df <- make_norm_fixture(cfg$quality_m, cfg$quality_c, seed, sparse = seed %% 2L == 0L)
  ranked <- reference_rank(df, "NRR")
  model <- CCS:::norm_prepare_counts(df, ranked)
  all_cuts <- enumerate_cuts(cfg$quality_m, k)
  all_auc <- CCS:::norm_auc_batch(all_cuts, model, method)
  exact <- max(apply(all_auc, 1, stats::median))
  elapsed <- system.time(result <- CCS:::subtypeNorm(df, k, method, "NRR", 1L, FALSE, "dp"))[["elapsed"]]
  d <- result$SearchDiagnostics
  quality[[length(quality) + 1L]] <- data.frame(
    seed = seed, M = cfg$quality_m, C = cfg$quality_c, K = k, method = method, rank = "NRR",
    nSample = nrow(df), partitions = nrow(all_cuts),
    exhaustive_best = exact, dp_best = result$Summary$medianROCAUC[1L],
    gap = exact - result$Summary$medianROCAUC[1L],
    refinement_gain = result$Summary$medianROCAUC[1L] - d$best_initial_median_auc,
    elapsed_seconds = elapsed, candidates = nrow(result$Summary),
    cutpoints = paste(d$selected_cutpoints, collapse = ","),
    mapping_sha256 = digest::digest(result$Summary$Info[1L], algo = "sha256", serialize = FALSE))
}
utils::write.csv(do.call(rbind, quality), file.path(output_dir, "quality.csv"), row.names = FALSE)

# Each measured search runs in a fresh R process. Setup/load time is excluded;
# the OS peak includes package loading, and for brute force excludes its worker.
run_one <- function(m, c, k, method, search, seed) {
  callr::r(function(m, c, k, method, search, seed) {
    devtools::load_all(quiet = TRUE)
    source("test/normCCS-5-more/fixtures.R")
    df <- make_norm_fixture(m, c, seed, sparse = TRUE)
    if (search == "legacy") {
      env <- new.env(parent = asNamespace("CCS"))
      sys.source("test/normCCS-5-more/reference/subtype-norm-0.8.4.R", env)
      fn <- env$subtypeNorm
    } else {
      fn <- function(...) CCS:::subtypeNorm(..., norm_search_method = "dp")
    }
    gc()
    elapsed <- system.time(result <- suppressWarnings(fn(df, k, method, "NRR", 1L, FALSE)))[["elapsed"]]
    peak <- ps::ps_memory_info()[["peak_wset"]] / 1024^2
    d <- result$SearchDiagnostics
    info <- result$Summary$Info[1L]
    data.frame(M = m, C = c, K = k, method = method, search = search,
               rank = "NRR", seed = seed, nSample = nrow(df),
               partitions = choose(m - 1L, k - 1L), candidates = nrow(result$Summary),
               elapsed_seconds = elapsed, master_peak_mib = peak,
               interval_mib = if (is.null(d)) NA_real_ else d$interval_storage_bytes / 1024^2,
               median_auc = result$Summary$medianROCAUC[1L],
               cohort_median_auc = stats::median(result$Raw$ROCAUC[result$Raw$Info == info & result$Raw$Cohort != "All"]),
               cutpoints = if (is.null(d)) NA_character_ else paste(d$selected_cutpoints, collapse = ","),
               mapping_sha256 = digest::digest(info, algo = "sha256", serialize = FALSE))
  }, args = list(m, c, k, method, search, seed), libpath = .libPaths())
}
timings <- list()
for (k in cfg$timing_k) for (method in c("Raw", "OptimizeRank")) {
  for (rep in seq_len(cfg$timing_reps)) for (search in c("legacy", "dp")) {
    cat("Timing", k, method, rep, search, "\n")
    row <- run_one(cfg$timing_m, cfg$timing_c, k, method, search, 307L)
    row$replicate <- rep
    timings[[length(timings) + 1L]] <- row
  }
}
utils::write.csv(do.call(rbind, timings), file.path(output_dir, "timings.csv"), row.names = FALSE)

large <- list()
support <- list()
for (m in cfg$scaling_m) for (k in cfg$scaling_k) {
  for (method in c("Raw", "OptimizeRank")) {
    cat("Scale", m, k, method, "\n")
    large[[length(large) + 1L]] <- run_one(m, cfg$scaling_c, k, method, "dp", 509L)
  }
}
utils::write.csv(do.call(rbind, large), file.path(output_dir, "scaling.csv"), row.names = FALSE)

# Group occupancy of the frozen selected mapping: no re-ranking per cohort.
df <- make_norm_fixture(cfg$support_m, cfg$support_c, 509L, sparse = TRUE)
for (k in cfg$support_k) {
  fit <- CCS:::subtypeNorm(df, k, "Raw", "NRR", 1L, FALSE, "dp")
  diag <- fit$SearchDiagnostics
  df$Grade <- k - 1L - findInterval(match(df$Subtype, diag$ranked_subtypes) - 1L, diag$selected_cutpoints)
  size <- stats::aggregate(list(n = df$response, responders = df$response),
                            df[c("Cohort", "Grade")], function(x) length(x))
  positive <- stats::aggregate(df$response, df[c("Cohort", "Grade")], sum)
  size$responders <- positive$x
  grid <- expand.grid(Cohort = sort(unique(df$Cohort)), Grade = 0:(k - 1L), stringsAsFactors = FALSE)
  size <- merge(grid, size, all.x = TRUE)
  size$n[is.na(size$n)] <- 0L
  size$responders[is.na(size$responders)] <- 0L
  size$K <- k
  support[[length(support) + 1L]] <- size
}
utils::write.csv(do.call(rbind, support), file.path(output_dir, "group-support.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
cat("Benchmark complete:", output_dir, "\n")
