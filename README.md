# CCS: Cohort Congress System

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%203.6.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.7.3-brightgreen.svg)](https://github.com/huangwb8/CCS)

### A computational framework for personalized pan-cancer genomic classification

`CCS` is an R package for personalized pan-cancer genomic classification. It supports cohort-level probability modeling, subtype discovery, and downstream response-oriented analyses, and is designed to work within the `GSClassifier` ecosystem.

## Tutorial

For teaching materials and step-by-step workflow explanations, please visit [Online Tutorial](https://huangwb8.github.io/ccs.principle/).

## Installation

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("huangwb8/luckyBase")
devtools::install_github("huangwb8/GSClassifier")
devtools::install_github("huangwb8/CCS")
```

For local development:

```bash
Rscript -e "devtools::document()"
Rscript -e "devtools::load_all()"
R CMD build .
```

## Citation

The `CCS` manuscript is currently under review.

## Authors

Weibin Huang (黄伟斌); <hwb2012@qq.com>; [Blog](https://blognas.hwb0307.com/); [ORCID: 0000-0003-0159-435X](https://orcid.org/0000-0003-0159-435X)

## License

`CCS` is released under the Apache-2.0 license. See [`license.txt`](license.txt) for details.
