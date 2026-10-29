# Investigate every retained quality case above the engineering comparison scale.
source("test/normCCS-5-more/00.Environment.R", encoding = "UTF-8")
devtools::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
paths <- args[!startsWith(args, "--")]
out <- norm_output_dir(if (length(paths)) paths[1L] else "test/normCCS-5-more/products")
input <- file.path(out, "objective-quality.csv")
target <- file.path(out, "objective-gap-diagnostics.csv")
if (file.exists(target) && !"--overwrite" %in% args) stop("Diagnostics exist; pass --overwrite explicitly.")
quality <- read.csv(input)
failed <- quality[quality$gap > 0.005 + 1e-10, ]
rows <- list()
for (i in seq_len(nrow(failed))) {
  x <- failed[i, ]
  df <- make_norm_fixture(x$M, x$C, x$seed, sparse = x$seed %% 2L == 0L)
  ranked <- reference_rank(df, "NRR")
  model <- CCS:::norm_prepare_counts(df, ranked)
  result <- CCS:::subtypeNorm(df, x$K, x$method, "NRR", 1L, FALSE, "dp", x$objective)
  d <- result$SearchDiagnostics
  cuts <- d$selected_cutpoints
  score <- function(cuts) CCS:::norm_objective_score(CCS:::norm_auc_batch(cuts, model, x$method), model$cohorts, x$objective)
  all_cuts <- enumerate_cuts(x$M, x$K)
  values <- score(all_cuts)
  # Count changed boundaries only after matching their positions in the ordered tuple.
  changes <- rowSums(sweep(all_cuts, 2L, cuts, "!="))
  global <- which(values >= max(values) - 1e-12)
  global_cuts <- all_cuts[global[1L], ]
  grade <- x$K - 1L - findInterval(match(df$Subtype, ranked$Subtype) - 1L, global_cuts)
  reference_auc <- vapply(model$rows, function(ids) {
    as.numeric(CCS:::get_perform_markers(df$response[ids], grade[ids], x$K, x$method)$ROCAUC)
  }, numeric(1))
  real <- reference_auc[model$cohorts != "All"]
  nc_real <- length(real)
  quantile_weight <- pmax(0, pmin(seq_len(nc_real) / nc_real, 0.25) -
                           (seq_len(nc_real) - 1) / nc_real) / 0.25
  lower <- sum(sort(real) * quantile_weight)
  reference_objective <- switch(x$objective, median = median(reference_auc),
    cohort_median = median(real), mean = mean(real), lower_tail = lower,
    mean_lower_tail = (mean(real) + lower) / 2)
  stopifnot(max(abs(reference_auc - CCS:::norm_auc_batch(matrix(global_cuts, 1L), model, x$method))) <= 1e-10,
            abs(reference_objective - max(values)) <= 1e-10)
  single <- changes <= 1L
  adjacent <- single
  for (j in seq_len(ncol(all_cuts) - 1L)) {
    outside <- setdiff(seq_len(ncol(all_cuts)), c(j, j + 1L))
    adjacent <- adjacent | rowSums(sweep(all_cuts[, outside, drop = FALSE], 2L, cuts[outside], "!=")) == 0L
  }
  nc <- length(model$cohorts)
  weights <- rbind(rep(1 / nc, nc), c(0, rep(1 / (nc - 1L), nc - 1L)), diag(nc))
  weights <- rbind(weights, -weights)
  loss <- matrix(model$loss, nrow = nc)
  starts <- lapply(seq_len(nrow(weights)), function(j) {
    cost <- matrix(as.vector(crossprod(weights[j, ], loss)), model$m, model$m)
    CCS:::norm_solve_dp(cost, x$K)$cuts
  })
  starts <- unique(c(starts, list(as.integer(floor(seq_len(x$K - 1L) * model$m / x$K)))))
  initial_scores <- score(do.call(rbind, starts))
  best <- if (x$objective == "median") result$Summary$medianROCAUC[1L] else result$Summary$objectiveScore[1L]
  stopifnot(abs(best - x$objective_score) <= 1e-10)
  rows[[i]] <- data.frame(seed = x$seed, C = x$C, K = x$K, method = x$method,
    objective = x$objective, gap = x$gap, cutpoints = paste(cuts, collapse = ","),
    global_cutpoints = paste(all_cuts[global[1L], ], collapse = ","),
    global_reference_objective = reference_objective,
    initial_contains_global = any(initial_scores >= max(values) - 1e-12),
    all_refinements_converged = all(vapply(d$refinements, function(z) z$converged, logical(1))),
    max_passes = max(vapply(d$refinements, function(z) z$passes, integer(1))),
    single_boundary_gain = max(values[single]) - best,
    adjacent_pair_gain = max(values[adjacent]) - best,
    any_pair_gain = max(values[changes <= 2L]) - best,
    minimum_changes_to_global = min(changes[global]))
}
diagnostics <- if (length(rows)) do.call(rbind, rows) else data.frame(seed = integer(), C = integer(), K = integer())
utils::write.csv(diagnostics, target, row.names = FALSE)
paths <- c(input, "test/normCCS-5-more/diagnose-objective-gaps.R", "R/subtype_norm_dp.R", "R/subtype_norm_objective.R")
utils::write.csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))),
                 file.path(out, "objective-gap-diagnostics-sources.csv"), row.names = FALSE)
cat("Diagnosed", nrow(failed), "retained cases.\n")
