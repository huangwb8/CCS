# normCCS 5–12 级搜索测试

本目录展示 CCS 0.8.5 的搜索改造：保留默认穷举，新增显式 DP 搜索；所有数据由代码生成，不读取真实患者文件。

建议先阅读 [性能与算法报告](report.html)，源码见 [report.Rmd](report.Rmd)。[验证记录](verification.md) 汇总回归、公开入口和包检查结果。

## 主要改动

`norm_search_method = "dp"` 复用 RR/NRR 排序，把样本压缩为队列与状态的应答计数，用动态规划生成候选，再优化原有的 AUC 中位数。最终候选逐队列通过原 ROC 函数复核。该方法只用单核；没有新增 C++ 或并行依赖。它对加权平均 Raw AUC 子问题精确，对原中位数目标近似，不保证全局最优。

```r
source("test/normCCS-5-more/fixtures.R")
df <- make_norm_fixture(m = 20, cohorts = 4, seed = 42)
# 公开 Raw ROC 沿用原接口的数值状态编码。
df$Subtype <- as.integer(sub("S", "", df$Subtype))
fit <- CCS::subtypePerformance(
  df, col_subtype = "Subtype", col_tissue = "tumor_type",
  n_norm_subtype = 8, norm_search_method = "dp",
  norm_roc_method = "Raw", norm_rank_marker = "NRR",
  numCores = 1
)
fit$Data$Normalization$Summary[1, ]
fit$Data$Normalization$SearchDiagnostics
```

默认 `"bruteforce"` 仍返回所有切分的结果。DP 的 `Raw` 和 `Summary` 只包含最终实际复核的少量候选；`SearchDiagnostics` 记录切点、排序、初值、评分模式、精修轮数、计算时间及参考复核次数。过大的穷举会提前报错，避免分配巨型组合矩阵和结果表。

## 运行方法

从 CCS 仓库根目录，用 `C:/R/R-4.3.1/bin/Rscript.exe`：

```powershell
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/test-search.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/scripts/tests/public-api.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/benchmark.R
& 'C:/R/R-4.3.1/bin/Rscript.exe' --vanilla test/normCCS-5-more/render-report.R
```

重复生成已存在的基准结果或 HTML，需在相应命令末尾显式加 `--overwrite`。全量组件测试也可分成 `test-search.R legacy 2`、`legacy 3`、`legacy 4`、`kernels`、`state` 执行，便于查看每一部分的退出状态。基准按顺序运行，不使用并行调度。已有包工程沿用 simple 流程，没有为此迁移其它分析的 targets 或补建 renv。

`compare-refinement.R` 可复现保存在 `reference/` 中的单切点初版结果；`benchmark-config.json` 列出正式数据规模和种子，隔离测试通过同一入口换用较小配置。

## 文件与证据

| 文件 | 用途 |
|---|---|
| `fixtures.R` | 固定种子合成数据、独立排序与切点枚举辅助函数 |
| `reference/subtype-norm-0.8.4.R` | 修改前的原始搜索函数，用作回归和性能基线 |
| `test-search.R` | 旧结果等价、DP 精确子问题、pROC 方向和并列、调用者状态测试 |
| `scripts/tests/public-api.R` | 实际运行公开函数及 `n_norm_subtype=NULL` 分支 |
| `scripts/tests/smoke-test.R` | 隔离产品目录的计算与报告轻量测试 |
| `benchmark.R` | 正确性差距、三次直接计时、高维搜索与分级人数检查 |
| `products/quality.csv` | 每个合成问题的完整穷举差距与精修收益 |
| `products/timings.csv` | 原版和 DP 的每次直接计时 |
| `products/scaling.csv` | 107／216 状态、20 队列、5–12 级的实际运行结果 |
| `products/group-support.csv` | 各队列各等级人数与应答人数 |
| `products/session-info.txt` | 实际 R 与依赖版本 |

固定合成数据上的加速和最优值匹配只说明工程行为。真实患者全量验证未运行；更高 K 的临床收益、跨队列稳定性和校准需要另行研究。
