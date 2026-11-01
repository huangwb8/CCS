<div align="center">
  <h1>CCS · Cohort Congress System</h1>
  <p>
    <a href="https://github.com/huangwb8/CCS"><img src="https://img.shields.io/badge/version-0.8.7-brightgreen.svg" alt="Version 0.8.7"></a>
    <a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-package-276DC3.svg" alt="R package"></a>
    <a href="license.txt"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg" alt="Apache 2.0 license"></a>
  </p>
  <p><a href="README_CN.md">简体中文</a> · <a href="#quick-start">Quick start</a> · <a href="#api-and-documentation">Documentation</a> · <a href="https://github.com/huangwb8/CCS/issues">Issues</a></p>
  <p><strong>A computational framework for personalized pan-cancer genomic classification.</strong></p>
</div>

CCS is an R package that represents a sample through the predictions of cohort-specific models built with **GSClassifier**. It supports subtype discovery from these cohort-probability profiles, prediction for new samples, and analysis of subtype-associated response across cohorts. The downstream **normCCS** workflow groups existing subtype states into ordered grades.

Start with the synthetic example below to explore normCCS without patient data or a pretrained model bank. Training the full genomic classification workflow requires expression data, gene sets, and gene annotations supplied by the user.

## Quick start

### Install

Development and example validation use **R 4.3.1**. [DESCRIPTION](DESCRIPTION) does not specify a minimum R version, and older releases have not been checked here. Enable Bioconductor repositories before installing from GitHub so that dependencies such as `ComplexHeatmap` and those required by `WGCNA` can be resolved:

```r
install.packages(c("remotes", "BiocManager"))
options(repos = BiocManager::repositories())

remotes::install_github("huangwb8/luckyBase", upgrade = "never")
remotes::install_github("huangwb8/GSClassifier", upgrade = "never")
remotes::install_github("huangwb8/CCS", upgrade = "never")

library(CCS)
packageVersion("CCS")
```

The dependency stack includes `xgboost`, `torch`, and `luz`. Deep-learning workflows also require the torch runtime; follow the [R torch installation guide](https://torch.mlverse.org/docs/articles/installation.html) when using that backend. If installation fails, check which dependency failed and its repository or system requirements before retrying. Windows source compilation may require the Rtools version matching your R installation.

### Run a synthetic normCCS analysis

This example creates 480 synthetic observations across three cohorts and eight subtype states, with response rates that increase with state number. It searches for a five-grade mapping on one CPU core.

```r
library(CCS)

df <- expand.grid(
  Cohort = c("A", "B", "C"),
  Subtype = 1:8,
  replicate = 1:20,
  KEEP.OUT.ATTRS = FALSE
)
df$response <- as.integer(df$replicate <= 2 * df$Subtype + 1)
df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
df$tumor_type <- "Lung"

fit <- subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 5, norm_search_method = "dp",
  norm_roc_method = "Raw", norm_objective = "cohort_median",
  numCores = 1, verbose = FALSE
)

normalization <- fit$Data$Normalization
head(normalization$Summary)
grades <- newSubtype(df$Subtype, normalization$Summary$Info[1])
table(grades)
```

`Summary` ranks candidate mappings by their AUC-based scores. Its `Info` column stores each mapping for application to subtype labels with `newSubtype()`. In this example, `grades` assigns all 480 observations to five grades, numbered 0–4. Search details are in `normalization$SearchDiagnostics`, and plots are in `fit$Plot`. The synthetic results demonstrate usage; predictive performance requires patient-data validation.

## From expression data to subtypes

CCS connects cohort-model training, probability-based representation, subtype discovery, and downstream evaluation:

| Task | Main entry points | Inputs or outputs |
| --- | --- | --- |
| Prepare training data | `createCCSData()` | Tissue/cohort-nested expression matrices, gene sets, and gene annotations |
| Train cohort models | `ccs()` | Prepared data → a `CCS` S4 object and saved model files |
| Explore and define subtypes | `dr()`, `cluster()`, `optimizeCluster()` | Cohort-probability representation → embeddings and subtype assignments |
| Fit a subtype classifier and predict | `scaller()`, `predict()` | A clustered CCS object and new expression data |
| Evaluate response and construct grades | `subtypePerformance()`, `newSubtype()` | Subtype labels, cohort labels, and binary response |
| Compare frozen representations | `ablation()` | A fitted CCS object, expression data, and sample metadata |

`ccs()` stores the cohort-probability representation in `object@Data$Probability$d1`. Downstream methods operate on the CCS object; full training is a separate workflow from the synthetic grade-search example. See the [API documentation](#api-and-documentation) for each function's input contract.

For response analysis, supply a complete table containing sample identifiers, subtypes, cohorts, tissues, and response coded as **0/1**. Each cohort must contain both response classes. The example uses the default columns `SampleIDs`, `Subtype`, `Cohort`, and `response`, and explicitly sets `col_tissue = "tumor_type"`.

## Choose a normCCS search

### Search a fixed number of grades

The default `norm_mode = "discrete"` searches only the requested number of grades. The default `norm_search_method = "bruteforce"` retains complete candidate tables; the number of partitions can grow rapidly. Choose `"dp"`, as in the quick start, to use deterministic, single-core dynamic programming with local refinement. DP returns evaluated final candidates and search diagnostics rather than every partition.

DP exactly optimizes a surrogate based on fixed-direction weighted mean Raw AUC. The full search refines candidates against the selected objective and checks them with the reference ROC evaluator. **It does not guarantee a global optimum for the selected public objective.**

### Explore an inherited sequence of grades

Reuse `df` from the quick start to explore K = 1 through K = 8:

```r
continuous_fit <- subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_mode = "continuous", norm_top_n = 10,
  norm_rank_marker = "ShrunkLogOR", norm_objective = "cohort_median",
  norm_roc_method = "Raw", numCores = 1, verbose = FALSE
)

continuous_fit$Data$Normalization$Summary
continuous_fit$Data$Normalization$SearchDiagnostics$layers[[5]]
```

Continuous mode fixes the state order once. K = 1 has a constant score with AUC = 0.5. At each later K, a child preserves every parent boundary and adds one split. Legal single-split extensions are scored with the selected objective before retaining at most `norm_top_n` candidates per layer. This mode always uses inherited DP on one core, regardless of `norm_search_method`.

`Summary` combines retained candidates across K, and the layer diagnostics record parent mappings, cutpoints, and scores. Pruning to a finite number of candidates makes the search approximate; neither global optimality nor improving scores with increasing K is guaranteed.

### Select the scoring objective

| `norm_objective` | Selection criterion |
| --- | --- |
| `"median"` (default) | Legacy median AUC, including pooled `All` |
| `"cohort_median"` | Median AUC across real cohorts |
| `"mean"` | Equal-cohort mean AUC |
| `"lower_tail"` | Mean of the lowest 25% of cohort AUC probability mass, including a fractional boundary cohort |
| `"mean_lower_tail"` | Equal-weight average of the cohort mean and lower-tail mean |

Non-default objectives exclude pooled `All` and require finite AUCs in every real cohort. They rank candidates by `objectiveScore`; `medianROCAUC` retains its legacy definition. Discrete and continuous searches use the chosen objective for candidate selection.

The default `norm_rank_marker = "NRR"` orders states by normalized response rate. Other options are `RR`, `NRR+RR`, `RR+NRR`, `MeanNRR`, `LowerTailNRR`, `ShrunkLogOR`, and `RankFusion`. Definitions and tie-breaking rules are in the [subtypePerformance reference](man/subtypePerformance.Rd).

### Interpret AUC and validation correctly

The legacy default `norm_roc_method = "OptimizeRank"` reorders grades using evaluation labels. The examples explicitly choose `"Raw"` to preserve the grade order, but the ROC evaluator still uses automatic direction selection. **Raw evaluation alone is not a frozen external prediction protocol.**

For external validation, estimate and freeze the state order, cutpoints, mapping, and score direction using training data. Keep external outcomes out of model and mapping selection. An improved training objective does not establish an external predictive gain or clinical benefit.

## API and documentation

| Resource | What it covers |
| --- | --- |
| [Training](man/ccs.Rd) and [data preparation](man/CCS-method.createCCSData.Rd) | Expression inputs, gene sets, annotations, and cohort-model parameters |
| [Subtype classifier](man/CCS-method.scaller.Rd) and [prediction](man/CCS-method.predict.Rd) | Classifier backends and prediction inputs |
| [Response and grade search](man/subtypePerformance.Rd) and [mapping application](man/newSubtype.Rd) | Input columns, objectives, search modes, and result structure |
| [Representation ablation](man/ablation.Rd) | Frozen-representation comparisons and their parameters |
| [External-cache guide](docs/external-cache.md) | Local storage layout and cache handling |
| [Changelog](CHANGELOG.md) | Package changes and development history |

After installation, use R help, for example `?CCS::ccs` or `?CCS::subtypePerformance`. The [online tutorial](https://huangwb8.github.io/ccs.principle/) is intended for step-by-step teaching; the complete, fully usable tutorial is planned for release after formal manuscript publication.

Historical analysis directories referenced by older documentation are not all included in this checkout. Large datasets, trained model banks, and generated reports are not bundled with the README example; reproducing a study requires its inputs and analysis environment.

## Development and reproducibility

Work from the repository root. To load the current source without installing it, use:

```r
devtools::load_all(".")
```

Repository conventions are documented in [AGENTS.md](AGENTS.md), with a [Claude Code adapter](CLAUDE.md). Use R 4.3.1 for Windows development checks. Scenario tests live under `test/` and are run separately from package checks; confirm that the relevant scenario is present before running a command from historical documentation.

For analyses defined by a `targets` graph, use the analysis's `_targets.R` and its actual store for execution, recovery, and monitoring. Configure target-level parallelism through `crew`. Store large data, models, package libraries, and generated reports outside the source tree; the local convention is `D:/cache/ccs/<analysis>/`, with existing analysis cache roots preserved. See the external-cache guide for details.

Human, AI, and tool contributions are recorded with [bensz-auto-contribution](https://github.com/huangwb8/bensz-auto-contribution) in [docs/contribution.bac](docs/contribution.bac). With Python 3.10+ and `bac` installed, verify the ledger with:

```bash
python -m bac --root . --bac-file docs/contribution.bac verify --json
```

Keep patient identifiers, protected health information, and proprietary genomic data out of commits and issue reports. Use de-identified or synthetic examples when reporting a problem.

## Citation, support, and license

The CCS manuscript is currently under review. Citation details will be added when available.

For bugs or usage questions, open a [GitHub issue](https://github.com/huangwb8/CCS/issues) with the CCS version, R session information, and a minimal reproducible example.

**Author and maintainer:** Weibin Huang (黄伟斌) · <hwb2012@qq.com> · [ORCID](https://orcid.org/0000-0003-0159-435X) · [Blog](https://blognas.hwb0307.com/).

CCS is distributed under the **Apache-2.0** license. See [license.txt](license.txt). Dependencies and external datasets retain their respective licenses and access conditions.
