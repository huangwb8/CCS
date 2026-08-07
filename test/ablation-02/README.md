# Ablation 02 执行顺序

该目录按文件名前缀表达推荐执行顺序。主分析分为数据审计与消融实验两个阶段，契约测试从 `10` 开始，避免与主分析阶段混淆。

## 主分析

| 顺序 | 文件 | 用途 | 主要产物 |
|---|---|---|---|
| `00` | `00.Environment.R` | 统一加载依赖与绘图环境 | 当前 R 会话环境 |
| `01` | `01-ablation-test-data.R` | 加载完整数据、恢复 tissue 标签并构建未筛选审计画像 | `tmp/ablation-experiment/data-profile.rds` |
| `01` | `01-ablation-test-data.Rmd` | 展示数据构成、来源异质性与审计边界 | `01-ablation-test-data.html`、数据审计 PDF 图 |
| `02` | `02-ablation-experiment.R` | 运行 d1 与 Direct-GSClassifier 的完整消融计算 | `tmp/ablation-experiment/` 下的 RDS/CSV |
| `02` | `02-ablation-experiment.Rmd` | 展示消融结果、统计证据与结论边界 | `02-ablation-experiment.html`、消融结果 PDF 图 |

同一阶段的 `_functions.R` 文件只保存该阶段专用函数，不作为独立执行入口。
两个阶段复用原有 `tmp/ablation-experiment/` 计算缓存，避免仅因入口编号变化复制大体积 RDS；执行顺序由代码与报告文件名表达。

## 契约测试

按以下顺序运行可先覆盖轻量内部合同，再进入依赖真实数据的集成合同：

```powershell
Rscript test/ablation-02/10-test-ablation-data-contract.R
Rscript test/ablation-02/11-test-ablation-representation.R
Rscript test/ablation-02/12-test-ablation-readout.R
Rscript test/ablation-02/13-test-ablation-decoder.R
Rscript test/ablation-02/14-test-ablation-real-encoder.R
Rscript test/ablation-02/15-test-ablation-public-representation.R
Rscript test/ablation-02/16-test-ablation-public-readout.R
```

`10`、`14`、`15`、`16` 需要项目约定的真实数据与冻结模型路径；可通过 `CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、`CCS_ABLATION_MODEL_ROOT` 等环境变量覆盖机器默认路径。
