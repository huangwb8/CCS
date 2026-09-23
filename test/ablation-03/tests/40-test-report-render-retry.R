test_retry_render <- function() {
  project_dir <- file.path(getwd(), "test", "ablation-03")
  old_dir <- setwd(project_dir)
  on.exit(setwd(old_dir), add = TRUE)
  source(file.path("renv", "activate.R"))
  source(file.path("targets", "functions.R"))

  report_file <- paste0(basename(tempfile("retry-")), " & report.Rmd")
  sanitized_file <- rmarkdown:::file_name_without_shell_chars(report_file)
  html_file <- sub("\\.Rmd$", ".html", report_file)
  preview_root <- tempfile("report-test-")
  on.exit(unlink(preview_root, recursive = TRUE), add = TRUE)
  on.exit(unlink(file.path(project_dir, c(report_file, sanitized_file, html_file))), add = TRUE)

  writeLines(c(
    "---", "title: Retry report", "output:",
    "  html_document:", "    css: templates/liquid_glass_theme.css",
    "---", "", "Retry succeeded."
  ), report_file, useBytes = TRUE)
  writeLines("leave this existing file alone", sanitized_file, useBytes = TRUE)
  original_collision <- readLines(sanitized_file, warn = FALSE)

  rendered <- .ablation03_render_report(report_file, cache_root = preview_root)
  stopifnot(
    identical(basename(rendered), html_file),
    file.exists(rendered), file.info(rendered)$size > 0L,
    grepl("Retry succeeded", paste(readLines(rendered, warn = FALSE), collapse = "\n"), fixed = TRUE),
    identical(readLines(sanitized_file, warn = FALSE), original_collision),
    length(list.files(".", pattern = "^ablation03-report-.*\\.Rmd$")) == 0L
  )
}
test_retry_render()
cat("report retry rendering passed\n")
