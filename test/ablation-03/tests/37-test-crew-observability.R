project_dir <- file.path(getwd(), "test", "ablation-03")
setwd(project_dir)
source(file.path("renv", "activate.R"))
source(file.path("targets", "functions.R"))

root <- tempfile("ablation03-crew-observability-")
dir.create(root)
on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

config <- .ablation03_observability_config(root)
stopifnot(
  config$backend == "crew_controller_local",
  config$workers >= 1L,
  dir.exists(config$worker_log_dir)
)

controller <- .ablation03_crew_controller(config)
on.exit(controller$terminate(), add = TRUE)
controller$start()
controller$push({
  Sys.sleep(1.2)
  4L
})
controller$wait(seconds_timeout = 20)
result <- controller$pop()
stopifnot(is.data.frame(result), nrow(result) == 1L)
stopifnot(identical(result$result[[1L]], 4L), identical(result$status[[1L]], "success"))
controller$terminate()

worker_logs <- list.files(config$worker_log_dir, full.names = TRUE)
worker_logs <- worker_logs[!dir.exists(worker_logs)]
stopifnot(length(worker_logs) >= 1L)
metrics <- .ablation03_read_resource_metrics(list(observability = config))
stopifnot(is.data.frame(metrics), nrow(metrics) >= 1L)
stopifnot(all(c("pid", "resident") %in% names(metrics)))

pdf <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf)
autometric::log_plot(metrics, metric = "resident")
grDevices::dev.off()
stopifnot(file.exists(pdf), file.info(pdf)$size > 0)
unlink(pdf, force = TRUE)

cat("crew observability tests passed\n")
