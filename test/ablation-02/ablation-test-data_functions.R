# Purpose: Keep tissue restoration and data profiling out of the loader script.
# Input: Nested CCS expression data, the audited cohort-to-tissue map, and metadata.
# Parameters: The historical source label defaults to "Undefined".
# Output: Resolved data, provenance-preserving cohort keys, and full audit profiles.

.atd_leaf_expression <- function(leaf) {
  if (is.matrix(leaf) || is.data.frame(leaf)) {
    return(as.matrix(leaf))
  }
  as.matrix(leaf$expr)
}

.atd_cohort_index <- function(data) {
  records <- unlist(lapply(names(data), function(tissue) {
    lapply(names(data[[tissue]]), function(cohort) {
      expr <- .atd_leaf_expression(data[[tissue]][[cohort]])
      data.frame(
        tissue = tissue,
        cohort = cohort,
        cohort_key = paste(tissue, cohort, sep = "/"),
        sample_count = ncol(expr),
        feature_count = nrow(expr),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
  do.call(rbind, records)
}

.atd_validate_tissue_mapping <- function(data, mapping, source_tissue) {
  required_columns <- c("cohort_id", "tissue")
  if (!all(required_columns %in% colnames(mapping))) {
    stop(
      "ablation-test-data: tissue mapping must contain cohort_id and tissue.",
      call. = FALSE
    )
  }
  if (anyDuplicated(mapping$cohort_id) || any(!nzchar(mapping$tissue))) {
    stop(
      "ablation-test-data: tissue mapping contains duplicate cohorts or blank labels.",
      call. = FALSE
    )
  }

  source_cohorts <- names(data[[source_tissue]])
  if (!setequal(source_cohorts, mapping$cohort_id)) {
    stop(
      "ablation-test-data: Undefined cohorts and the audited tissue map differ.",
      call. = FALSE
    )
  }

  target_conflicts <- vapply(seq_len(nrow(mapping)), function(i) {
    mapping$cohort_id[i] %in% names(data[[mapping$tissue[i]]])
  }, logical(1))
  if (any(target_conflicts)) {
    stop(
      "ablation-test-data: resolved tissue already contains cohort(s): ",
      paste(mapping$cohort_id[target_conflicts], collapse = ", "),
      call. = FALSE
    )
  }
  mapping
}

.atd_resolve_tissue_groups <- function(
    data,
    mapping,
    source_tissue = "Undefined"
) {
  mapping <- .atd_validate_tissue_mapping(data, mapping, source_tissue)
  source_data <- data[[source_tissue]]
  audit <- lapply(seq_len(nrow(mapping)), function(i) {
    cohort <- mapping$cohort_id[i]
    resolved_tissue <- mapping$tissue[i]
    expr <- .atd_leaf_expression(source_data[[cohort]])
    data.frame(
      cohort = cohort,
      source_tissue = source_tissue,
      resolved_tissue = resolved_tissue,
      sample_count = ncol(expr),
      feature_count = nrow(expr),
      stringsAsFactors = FALSE
    )
  })
  audit <- do.call(rbind, audit)

  resolved <- data
  for (i in seq_len(nrow(mapping))) {
    cohort <- mapping$cohort_id[i]
    tissue <- mapping$tissue[i]
    resolved[[tissue]][[cohort]] <- source_data[[cohort]]
  }
  resolved[[source_tissue]] <- NULL

  list(data = resolved, audit = audit)
}

.atd_resolve_cohort_keys <- function(
    cohort_keys,
    mapping,
    source_tissue = "Undefined"
) {
  bank_tissue <- sub("/.*$", "", cohort_keys)
  cohort <- sub("^[^/]*/", "", cohort_keys)
  mapped_tissue <- mapping$tissue[match(cohort, mapping$cohort_id)]
  resolved_tissue <- ifelse(
    bank_tissue == source_tissue & !is.na(mapped_tissue),
    mapped_tissue,
    bank_tissue
  )
  paste(resolved_tissue, cohort, sep = "/")
}

.atd_build_data_profile <- function(
    raw_cohort_index,
    resolved_cohort_index,
    metadata,
    tissue_resolution_audit,
    filtered_model_cohorts,
    filtered_cohorts
) {
  cohort_profile <- metadata |>
    dplyr::group_by(.data$tissue, .data$cohort, .data$analysis_set) |>
    dplyr::summarise(
      sample_count = dplyr::n(),
      unique_sample_count = dplyr::n_distinct(.data$sample_id),
      assay_type = paste(sort(unique(.data$assay_type)), collapse = "; "),
      platform_count = dplyr::n_distinct(.data$platform_id),
      source_system = paste(sort(unique(.data$source_system)), collapse = "; "),
      .groups = "drop"
    )

  tissue_profile <- metadata |>
    dplyr::group_by(.data$tissue, .data$analysis_set) |>
    dplyr::summarise(
      sample_count = dplyr::n(),
      unique_sample_count = dplyr::n_distinct(.data$sample_id),
      cohort_count = dplyr::n_distinct(.data$cohort),
      assay_count = dplyr::n_distinct(.data$assay_type),
      source_count = dplyr::n_distinct(.data$source_system),
      .groups = "drop"
    )

  tissue_totals <- metadata |>
    dplyr::group_by(.data$tissue) |>
    dplyr::summarise(
      sample_count = dplyr::n(),
      unique_sample_count = dplyr::n_distinct(.data$sample_id),
      cohort_count = dplyr::n_distinct(.data$cohort),
      median_cohort_size = stats::median(table(.data$cohort)),
      query_sample_count = sum(.data$analysis_set == "external_query"),
      query_fraction = mean(.data$analysis_set == "external_query"),
      .groups = "drop"
    )

  source_profile <- metadata |>
    dplyr::count(.data$tissue, .data$source_system, name = "sample_count") |>
    dplyr::group_by(.data$tissue) |>
    dplyr::mutate(tissue_fraction = .data$sample_count / sum(.data$sample_count)) |>
    dplyr::ungroup()

  resolution_audit <- tissue_resolution_audit |>
    dplyr::mutate(
      cohort_key = paste(.data$resolved_tissue, .data$cohort, sep = "/"),
      analysis_set = ifelse(
        .data$cohort_key %in% filtered_cohorts,
        "external_query",
        "reference_atlas"
      )
    )

  filtered_key_audit <- data.frame(
    bank_cohort_key = filtered_model_cohorts,
    cohort_key = filtered_cohorts,
    remapped = filtered_model_cohorts != filtered_cohorts,
    stringsAsFactors = FALSE
  )

  overview <- data.frame(
    item = c(
      "Raw samples",
      "Raw cohorts",
      "Raw tissue groups",
      "Resolved tissue groups",
      "Remapped cohorts",
      "Remapped samples",
      "Unresolved tissue samples",
      "External filtered cohorts"
    ),
    value = c(
      sum(raw_cohort_index$sample_count),
      nrow(raw_cohort_index),
      dplyr::n_distinct(raw_cohort_index$tissue),
      dplyr::n_distinct(resolved_cohort_index$tissue),
      nrow(resolution_audit),
      sum(resolution_audit$sample_count),
      sum(metadata$tissue == "Undefined"),
      length(filtered_cohorts)
    ),
    stringsAsFactors = FALSE
  )

  list(
    overview = overview,
    raw_cohort_index = raw_cohort_index,
    resolved_cohort_index = resolved_cohort_index,
    resolution_audit = resolution_audit,
    filtered_key_audit = filtered_key_audit,
    cohort_profile = cohort_profile,
    tissue_profile = tissue_profile,
    tissue_totals = tissue_totals,
    source_profile = source_profile
  )
}
