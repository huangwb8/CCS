# Ablation 02 执行顺序

该目录按文件名前缀表达推荐执行顺序。主分析分为数据审计与消融实验两个阶段；契约测试集中在 `tests/` 中，从 `10` 开始，避免与主分析阶段混淆。

## 主分析

| 顺序 | 文件 | 用途 | 主要产物 |
|---|---|---|---|
| `00` | `00.Environment.R` | 统一加载依赖与绘图环境 | 当前 R 会话环境 |
| `01` | `01-ablation-test-data.R` | 加载完整数据、恢复 tissue 标签并构建未筛选审计画像 | `tmp/ablation-experiment/data-profile.rds` |
| `01` | `01-ablation-test-data.Rmd` | 展示数据构成、来源异质性与审计边界 | `01-ablation-test-data.html`、数据审计 PDF 图 |
| `02` | `02-ablation-experiment.R` | 运行 d1 与 Direct-GSClassifier 的完整消融及二维 cohort-evidence scaling | `tmp/ablation-experiment/` 下的 RDS/CSV |
| `02` | `02-ablation-cohort-scaling.R` | 复用已有输入，单独重算 breadth/depth/matched-size scaling | `tmp/ablation-experiment/cohort_scaling.rds` schema v2 |
| `02` | `02-ablation-experiment.Rmd` | 展示消融、二维 scaling、matched-size 对照与证据边界 | `02-ablation-experiment.html`、消融结果 PDF 图 |

同一阶段的 `_functions.R` 文件只保存该阶段专用函数，不作为独立执行入口。
两个阶段复用原有 `tmp/ablation-experiment/` 计算缓存，避免仅因入口编号变化复制大体积 RDS；执行顺序由代码与报告文件名表达。
已有完整消融产物时，可单独运行 `Rscript test/ablation-02/02-ablation-cohort-scaling.R`，无需重算 retrieval、learning curve 与 decoder。新 schema 的 `design` 保存 bank 组成与排除原因，`metrics` 保存逐 bank 长表，`contrasts` 保存 repeat-level slope、交互和 matched-size 差，`diagnostics` 分开记录完整 Direct 主合同与 `Direct-GSClassifier-TSP` 敏感性合同；缺失独立生物锚点时，`reasons` 明确返回 `not_evaluated`。

## 契约测试

`tests/` 中的脚本按以下顺序运行，可先覆盖轻量内部合同，再进入依赖真实数据的集成合同：

```powershell
Rscript test/ablation-02/tests/10-test-ablation-data-contract.R
Rscript test/ablation-02/tests/11-test-ablation-representation.R
Rscript test/ablation-02/tests/12-test-ablation-readout.R
Rscript test/ablation-02/tests/13-test-ablation-decoder.R
Rscript test/ablation-02/tests/14-test-ablation-real-encoder.R
Rscript test/ablation-02/tests/15-test-ablation-public-representation.R
Rscript test/ablation-02/tests/16-test-ablation-public-readout.R
Rscript test/ablation-02/tests/17-test-ablation-cohort-scaling.R
Rscript test/ablation-02/tests/18-test-ablation-report-figure-02.R
```

`10`、`14`、`15`、`16` 需要项目约定的真实数据与冻结模型路径；`17` 使用不平衡的合成 tissue bank，锁定 breadth/depth 嵌套性、matched-size 配对、block 等权、双 Direct 合同、指标角色与 repeat bootstrap。真实数据路径可通过 `CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、`CCS_ABLATION_MODEL_ROOT` 等环境变量覆盖机器默认值。
