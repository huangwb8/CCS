#!/usr/bin/env Rscript

# Frozen-bank design contract for the seven-tissue depth sequence.
source(file.path("R", "ablation.R"))

depth_tissues <- c("ACC", "BRCA", "CRC", "KIRC", "PAAD", "PRAD", "STAD")
capacities <- c(8L, 20L, 41L, 11L, 13L, 16L, 20L)
modules <- do.call(rbind, lapply(seq_along(depth_tissues), function(i) {
  data.frame(
    module_id = sprintf("%s-%02d", depth_tissues[i], seq_len(capacities[i])),
    tissue = depth_tissues[i],
    cohort = sprintf("C%02d", seq_len(capacities[i])),
    stringsAsFactors = FALSE
  )
}))
modules <- rbind(modules, data.frame(
  module_id = paste0("OTHER-", seq_len(13L)),
  tissue = paste0("OTHER-", seq_len(13L)),
  cohort = "C01",
  stringsAsFactors = FALSE
))
module_counts <- c(10L, 25L, 50L, 75L, 100L, 125L, 150L)
baseline <- .ablation_cohort_bank_design(modules, module_counts, 5L, 20260805L)
design <- .ablation_cohort_bank_design(
  modules, module_counts, 5L, 20260805L,
  depth_tissues = depth_tissues, depth_max = 8L
)
repeat_design <- .ablation_cohort_bank_design(
  modules, module_counts, 5L, 20260805L,
  depth_tissues = depth_tissues, depth_max = 8L
)
stopifnot(identical(design$design_hash, repeat_design$design_hash))
score_indices <- .ablation_bank_score_seed_indices(
  modules, module_counts, 5L, 20260805L, design$design,
  depth_tissues = depth_tissues
)
legacy_indices <- match(design$design$design_id, baseline$design$design_id)
stopifnot(identical(score_indices[!is.na(legacy_indices)],
  legacy_indices[!is.na(legacy_indices)]))
stopifnot(all(score_indices[is.na(legacy_indices)] > 500000L))
stopifnot(identical(score_indices,
  .ablation_bank_score_seed_indices(modules, module_counts, 5L,
    20260805L, design$design, depth_tissues = depth_tissues)))

for (family in c("breadth", "matched")) {
  before <- baseline$design[baseline$design$design_family == family, ]
  after <- design$design[design$design$design_family == family, ]
  rownames(before) <- NULL
  rownames(after) <- NULL
  stopifnot(identical(before, after))
}
depth <- design$design[design$design$design_family == "depth", ]
stopifnot(nrow(depth) == 40L)
for (repeat_id in unique(depth$repeat_id)) {
  rows <- depth[depth$repeat_id == repeat_id, ]
  stopifnot(identical(rows$level, seq_len(8L)))
  stopifnot(identical(rows$module_count, 7L * seq_len(8L)))
  stopifnot(all(rows$tissue_count == 7L))
  stopifnot(all(vapply(rows$tissues, function(x) {
    identical(x, sort(depth_tissues))
  }, logical(1L))))
  for (i in 2:8) {
    stopifnot(identical(rows$parent_design_id[i], rows$design_id[i - 1L]))
    stopifnot(all(rows$module_ids[[i - 1L]] %in% rows$module_ids[[i]]))
  }
}

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = conditionMessage)
  stopifnot(is.character(message), !is.na(message), grepl(pattern, message))
}
expect_error(.ablation_cohort_bank_design(modules, module_counts, 2L, 1L,
  depth_tissues = c(depth_tissues, "UNKNOWN"), depth_max = 8L), "absent")
expect_error(.ablation_cohort_bank_design(modules, module_counts, 2L, 1L,
  depth_tissues = depth_tissues, depth_max = 21L), "fewer than")
expect_error(.ablation_cohort_bank_design(modules, module_counts, 2L, 1L,
  depth_tissues = c(depth_tissues, "ACC"), depth_max = 8L), "unique")
duplicate <- modules
duplicate$module_id[2L] <- duplicate$module_id[1L]
expect_error(.ablation_cohort_bank_design(duplicate, module_counts, 2L, 1L,
  depth_tissues = depth_tissues, depth_max = 8L), "module IDs")
expect_error(.ablation_resolve_representation_config(1L,
  list(scaling = list(depth_tissues = depth_tissues))), "set together")

message("seven-tissue depth design passed")
