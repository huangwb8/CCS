# CCS: Cohort Congress System

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%203.6.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.8.6-brightgreen.svg)](https://github.com/huangwb8/CCS)

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

### 选择 normCCS 的目标函数

`norm_objective = "median"` 默认复现包含 pooled `All` 的旧中位数目标。
显式选择 `"cohort_median"`、`"mean"`、`"lower_tail"` 或 `"mean_lower_tail"`
时，仅对真实队列 AUC 评分，队列等权。下尾取最低 25% 概率质量，边界队列按比例计入；
混合目标由均值和下尾均值各占 0.5。穷举和 DP 使用同一目标，DP 整体仍是近似搜索。

```r
fit <- subtypePerformance(df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_search_method = "dp", norm_roc_method = "Raw",
  norm_objective = "mean_lower_tail", numCores = 1)
fit$Data$Normalization$Summary[1, ]
```

新增策略以 `objectiveScore` 择优；`medianROCAUC` 保留旧定义。缺少
`Repeat$norm_objective` 表示使用遗留 `median`。参见
[目标函数报告](test/normCCS-5-more/objective-report.html)和
[复现说明](test/normCCS-5-more/README.md)。新目标的训练内得分改善不等同于外部预测收益。

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
