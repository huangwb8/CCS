# ablation-03：生物学锚点与表示性质

## 运行约定

推荐在 RStudio 中打开本目录下的 `ablation-03.Rproj`。所有脚本会自动识别当前目录或仓库根目录，不需要手动 `setwd()`：

- `01-ablation03-test-data.R`：准备数据审计缓存；
- `02-ablation03-experiment.R`：运行完整消融实验；
- `02-ablation03-cohort-scaling.R`：单独重算 cohort scaling；
- `01-ablation03-biology-cache.R`：一次性构建表达矩阵/anchor 子集缓存；
- `03-ablation-biology.R`：运行生物学锚点评估（也可使用同名 Python 脚本）。

中间结果统一写入本目录的 `tmp/`，图表写入本目录的 `figures/`。外部数据路径可通过 `CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、`CCS_FULL_EXPRESSION_RDS` 和 `CCS_GENE_SIGNATURE_RDS` 覆盖，便于在不同机器上复现或排查。

本目录是 ablation-02 的独立分析副本：完整纳入 02 的几何、retrieval、readout、decoder、learning-curve、technical-source 与 cohort-scaling 分析，并在其上新增 biological-anchor readout；不覆盖或写入 `test/ablation-02/`。

## 执行顺序

1. `02-ablation03-experiment.Rmd` 读取并完整展示 `test/ablation-02/tmp/ablation-experiment/` 的 02 基线结果。
2. `01-ablation03-biology-cache.R` 在一次性准备阶段读取完整表达矩阵，按冻结 manifest 提取目标样本/基因并写入 `expression-anchor-cache.rds`；脚本完成后释放大对象。
3. `03-ablation-biology.R` 只消费带 source/signature/sample hash 校验的缓存，不会回退读取完整表达矩阵。
4. 以 cohort 内 z-score signature mean 计算增殖、免疫 TME、基质 TME、IFN/IL6 四类独立 anchor。
5. 在相同 top-15 邻居边界上配对比较 Direct-GSClassifier 与 Cohort-d1，输出覆盖审计、utility、bootstrap 区间和图形；结果作为 HTML 的新增章节。

## 主要产物

- `tmp/ablation-biology/anchor_coverage.csv`
- `tmp/ablation-biology/anchor_utility.csv`
- `tmp/ablation-biology/anchor_contrasts.csv`
- `tmp/ablation-biology/expression-anchor-cache.rds`（一次性缓存，含 schema、来源哈希、样本键哈希与覆盖审计）
- `figures/figure-01-anchor-coverage.pdf` 与 `.jpg`
- `figures/figure-02-biological-utility.pdf` 与 `.jpg`
- `02-ablation03-experiment.html`

外部表达矩阵和 signature RDS 均为只读输入；当前报告不把 cancer_type 当作 biological utility，也不执行 PAM50/CMS 的临时反推。
