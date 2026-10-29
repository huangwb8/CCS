# Render already computed benchmark products; no search is performed here.
args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- args[args != "--overwrite"]
root <- normalizePath("test/normCCS-5-more", winslash = "/", mustWork = TRUE)
output_file <- if (length(args)) normalizePath(args[1L], winslash = "/", mustWork = FALSE) else file.path(root, "report.html")
source(file.path(root, "fixtures.R"), encoding = "UTF-8")
if (dirname(output_file) != root) norm_output_dir(dirname(output_file))
if (file.exists(output_file) && !overwrite) stop("Report exists; use --overwrite explicitly.")
products <- Sys.getenv("BENSZ_NORM_PRODUCTS", file.path(root, "products"))
products <- norm_output_dir(products)
figure_dir <- if (dirname(output_file) == root) file.path(root, "reports/figures") else
  file.path(dirname(output_file), "figures")
if (!rmarkdown::pandoc_available()) {
  candidates <- c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
                  "C:/Program Files/RStudio/bin/pandoc",
                  file.path(Sys.getenv("USERPROFILE"), ".local/bin"))
  found <- candidates[file.exists(file.path(candidates, "pandoc.exe"))]
  if (!length(found)) stop("Pandoc is required; configure RSTUDIO_PANDOC.")
  Sys.setenv(RSTUDIO_PANDOC = found[1L])
}
rmarkdown::render(file.path(root, "report.Rmd"), output_file = basename(output_file),
                  output_dir = dirname(output_file), params = list(products = products, figure_dir = figure_dir),
                  quiet = TRUE, envir = new.env())
