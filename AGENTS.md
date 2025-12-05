# Repository Guidelines

## Project Structure & Module Organization
Core package code lives in `R/`, with `ccs*.R` owning S4 classes, `hyperTuningGS.R` and `importance.R` for modeling, and utilities such as `normalize.R`, `phenotype.R`, and `simulatePA.R` for preprocessing. Generated Rd files remain in `man/`, while the exploratory scripts, saved models, and fixtures under `test/` mirror real workflows (see `test/ccs/project_01/`). Keep prompts, licensing, and README material at the repo root, and add any new auxiliary data next to the script that consumes it.

## Build, Test, and Development Commands
Refresh documentation before committing with `Rscript -e "devtools::document()"`. Package locally with `R CMD build .` followed by `R CMD check ccs_0.7.0.tar.gz` to catch lint, dependency, and vignette issues. During interactive work `Rscript -e "devtools::load_all()"` mirrors `library(CCS)` without reinstalling. Scenario tests run via plain scripts; e.g., `Rscript test/03.test.Classes.R` validates the class definitions and `Rscript test/test.classifier_performance.R` exercises the classifier benchmarks.

## Coding Style & Naming Conventions
Match the existing tidy R style: two-space indent, `<-` for assignment, and `snake_case` for functions and variables while keeping S4 classes and exported methods in `CamelCase`. Keep roxygen headers complete (`@description`, `@param`, `@return`, `@examples`) so `man/` stays in sync after `devtools::document()`. Reuse verbs from the packages declared in `DESCRIPTION`, and update that file before relying on new imports.

## Testing Guidelines
Each script in `test/` targets a reproducible scenario; follow the numbering scheme (`01.test.*`, `02.test.*`, etc.) to signal execution order. Scripts must be runnable through `Rscript` without interactive prompts and should cleanly write their `.rds` artifacts into dedicated subdirectories (e.g., `test/ccs/project_<id>/`). Favor `stopifnot` or lightweight `testthat` expectations to make failures obvious and add sample data snapshots when a change impacts probability outputs or subtype plots.

## Commit & Pull Request Guidelines
Git history favors concise imperative summaries with optional release tags (`v0.7.0`, `Routine update`). Use the same style, referencing issue IDs when relevant. Pull requests should include: purpose, datasets or parameters touched, confirmation that `R CMD check` and the relevant `test/*.R` scripts succeeded, and screenshots or metrics for any plotting/statistics change. Mention large artifacts deliberately excluded from git and coordinate with reviewers before force-pushing rebases.

## Security & Data Handling
Handle genomic `.rds` files as sensitive. Keep PHI and proprietary data out of the repo, scrub metadata in sample fixtures, and document any secure storage requirements in the PR description. Gate remote downloads behind explicit arguments and never embed credentials or tokens in tracked files.
