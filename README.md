# CCS: Cohort Congress System

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%203.6.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.8.5-brightgreen.svg)](https://github.com/huangwb8/CCS)

### A computational framework for personalized pan-cancer genomic classification

`CCS` is an R package for personalized pan-cancer genomic classification. It supports cohort-level probability modeling, subtype discovery, and downstream response-oriented analyses, and is designed to work within the `GSClassifier` ecosystem.

## Tutorial

For teaching materials and step-by-step workflow explanations, please visit [Online Tutorial](https://huangwb8.github.io/ccs.principle/). Please note that the tutorial will be made publicly available in its full, fully usable form only after the manuscript is formally published.

## Searching five or more normCCS grades

Use `subtypePerformance(..., n_norm_subtype = 8, norm_search_method = "dp")`
for deterministic, single-core search without enumerating every partition. The
default `"bruteforce"` retains the original complete candidate tables. DP optimizes
fixed-direction mean AUC subproblems and refines the original median AUC objective;
it does not guarantee the global median optimum. Its final candidates are checked
with the original ROC evaluator. `norm_roc_method = "Raw"` preserves the grade
ordering; the legacy `"OptimizeRank"` reorders grades using evaluation labels and
should not be used as a frozen external prediction score.

See the [synthetic benchmark and reproducible report](test/normCCS-5-more/README.md)
for exactness checks, measured speedups, resource usage, and the public API example.

## Citation

The `CCS` manuscript is currently under review.

## AI-Assisted Development

Project-wide AI instructions live in [`AGENTS.md`](AGENTS.md); [`CLAUDE.md`](CLAUDE.md) is the Claude Code adapter. Human, AI, and tool contributions are recorded with [bensz-auto-contribution](https://github.com/huangwb8/bensz-auto-contribution) in [`docs/contribution.bac`](docs/contribution.bac). The ledger requires Python 3.10+ and can be verified with:

```bash
python -m bac --root . --bac-file docs/contribution.bac verify --json
```

## Authors

Weibin Huang (黄伟斌); <hwb2012@qq.com>; [Blog](https://blognas.hwb0307.com/); [ORCID: 0000-0003-0159-435X](https://orcid.org/0000-0003-0159-435X)

## License

`CCS` is released under the Apache-2.0 license. See [`license.txt`](license.txt) for details.
