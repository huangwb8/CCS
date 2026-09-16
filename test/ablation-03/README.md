# ablation-03：生物学锚点与表示性质

## 运行约定

推荐在 RStudio 中打开本目录下的 `ablation-03.Rproj`。所有脚本会自动识别当前目录或仓库根目录，不需要手动 `setwd()`：

- `01-ablation03-test-data.R`：准备数据审计缓存；
- `02-ablation03-experiment.R`：运行完整消融实验；
- `tests/26-test-endpoint-eligibility.R`：检查 external candidate 保留与癌种相关 endpoint 的分层资格合同；
- `02-ablation03-cohort-scaling.R`：单独重算 cohort scaling；
- `01-ablation03-biology-cache.R`：一次性构建表达矩阵/anchor 子集缓存；
- `03-ablation-biology.R`：运行生物学锚点评估。
- `04-ablation03-structural-reproducibility.R`：比较 Direct 与 d1 的独立队列间生物状态几何可复现性。

中间结果统一写入本目录的 `tmp/`，图表写入本目录的 `figures/`。外部数据路径可通过 `CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、`CCS_FULL_EXPRESSION_RDS` 和 `CCS_GENE_SIGNATURE_RDS` 覆盖，便于在不同机器上复现或排查。

本目录是 ablation-02 的独立分析副本：完整纳入 02 的几何、retrieval、readout、decoder、learning-curve、technical-source 与 cohort-scaling 分析，并在其上新增 biological-anchor readout；不覆盖或写入 `test/ablation-02/`。

## 执行顺序

1. 运行 `02-ablation03-experiment.R`：自动准备数据，复用 CCS 中已有 d1，重建或读取严格校验的 Direct 缓存，生成本目录的完整实验结果。
2. 运行 `01-ablation03-biology-cache.R`：校验主实验来源，按 `sample-contract.rds` 固定的完整 reference/query 范围提取表达锚点。signature 按名称选择，来源与配置均记录哈希。
3. 运行 `03-ablation-biology.R`：完整 reference 拟合共同基因尺度，按 `anchor_retrieval.rds` 中所有候选 query 的 top-15 邻居评价；两臂完整后先按 query 配对，再以 cohort 为独立单位推断。
4. 运行 `04-ablation03-structural-reproducibility.R`：计算双向 module-bank 投影及共同 tissue 的匹配敏感性；并列分位边界不强制切成高低状态。
5. 渲染 `01-ablation03-test-data.Rmd` 和 `02-ablation03-experiment.Rmd`。主报告读取本目录 `tmp/` 的本轮产物，渲染前验证三个阶段的来源指纹。

癌种检索、technical excess、readout 和 learning curve 使用其预先声明的可估计子集；连续 anchor 与结构分析保留全部候选外部 cohort，实际有效覆盖另行审计。readout 两边均从原始表示进入函数，由每个训练折拟合尺度，避免 query 重复标准化。

任何受追踪输入/代码改变后，按依赖顺序重跑；旧产物不能与新产物拼接后渲染。仅重跑 scaling 时可使用 `02-ablation03-cohort-scaling.R`，它会验证完整输入和配置，并更新主实验指纹；后续生物学与结构阶段需重新生成。

## 主要产物

- `tmp/ablation-experiment/direct-feature-cache.rds`（持久化 Direct-GSClassifier 特征缓存；缓存键校验表达矩阵、样本顺序、冻结模型元数据与 feature manifest，命中后跳过 Direct 特征重建）

- `tmp/ablation-biology/anchor_coverage.csv`
- `tmp/ablation-biology/anchor_utility.csv`
- `tmp/ablation-biology/anchor_contrasts.csv`
- `tmp/ablation-biology/anchor_inference.csv`（cohort-level 效应量、95% CI、P 值和 BH 校正 P 值）
- `tmp/ablation-biology/expression-anchor-cache.rds`（一次性缓存，含 schema、来源哈希、样本键哈希与覆盖审计）
- `tmp/ablation-biology/anchor_missing_pairs.csv`（每个 anchor/表示的有效与缺失邻居对数）
- `tmp/ablation-structural-reproducibility/structural_entity_audit.csv`
- `tmp/ablation-structural-reproducibility/structural_pair_comparisons.csv`
- `tmp/ablation-structural-reproducibility/structural_summary.csv`
- `tmp/ablation-structural-reproducibility/structural_similarity_matrices.rds`
- `tmp/ablation-structural-reproducibility/structural_directional_summary.csv`
- `tmp/ablation-structural-reproducibility/structural_directional_pair_comparisons.csv`
- `tmp/ablation-structural-reproducibility/structural_directional_sample_audit.csv`
- `tmp/ablation-structural-reproducibility/structural_module_bank_audit.csv`
- `tmp/ablation-structural-reproducibility/structural_matched_bank_design.csv`
- `tmp/ablation-structural-reproducibility/structural_matched_bank_repeats.csv`
- `tmp/ablation-structural-reproducibility/structural_matched_bank_summary.csv`
- `tmp/ablation-structural-reproducibility/ablation03-structural-reproducibility.rds`
- `figures/figure-01-native-geometry.pdf` 至 `figures/figure-14-structural-matched-bank-sensitivity.pdf`（出版矢量图；编号按报告模块保留）
- `tmp/plot-previews/02-ablation-experiment-*/figure-*.jpg`（按运行时间戳保存的 200 dpi 预览）
- `02-ablation03-experiment.html`

外部表达矩阵和 signature RDS 均为只读输入；当前报告不把 cancer_type 当作 biological utility，也不执行 PAM50/CMS 的临时反推。

共享准备层会保留通过样本对齐、无 fit/target 重叠和矩阵可用性检查的全部 external candidate cohorts。`endpoint_eligibility.csv` 同时记录 `candidate`、`estimable` 与 `not_estimable`；癌种 retrieval、technical excess、linear readout 和 learning curve 使用预先声明的 reference cancer-cohort 支持门槛，geometry、连续 anchor 与结构复现不再被该门槛静默截断。报告中的 `query_cohort_count` 是癌种 readout 的可估计集合；完整候选数使用 `external_candidate_cohort_count`。

结构可复现性模块与现有 anchor readout 回答不同问题：前者检查生物状态原型之间的几何能否跨独立队列复现，后者检查单个 query 的局部邻居是否保持连续锚点接近。这里的 externality 由“目标样本是否参与当前 module-bank 的训练”定义，因此 reference/external 两组可以交换 bank 与 target 角色；反向投影用于检验互惠性和 bank 依赖，不替代生产模型的正向外部验证。队列对共享 cohort 节点，故报告使用 cohort-node bootstrap 区间；22 对 22 重复匹配报告的是设计敏感性范围，不是置信区间。
