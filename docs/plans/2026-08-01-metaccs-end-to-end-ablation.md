# CCS 端到端 metaCCS 表征消融扩展实施计划

## 问题是什么

- **当前现象：** `CCS::ablation()` 的 cohort experiment 把 Direct-GSClassifier 与 Cohort d1 比较到共同秩表示和 linear probe；tissue-first experiment 则固定完整 d1，只比较 Two-stage 与 One-stage。两者都没有直接回答“cohort bank 的变化能否穿过实际的 two-stage UWOT、d3 和 DBSCAN，最终改善 metaCCS”。
- **影响：** 目前只能确认 Cohort d1 小幅改变了局部纯度和 cohort mixing，不能判断这些变化在完整聚类流程中是有效信号、被下游放大后的工程差异，还是没有任务价值的几何扰动。
- **已知实现事实：** 现用参数表中目标 `paramMD5` 对应 `dimension = c(5, 2)`、`n_neighbors = 30`、`min_dist = 0.01`、`spread = 0.75`、`set_op_mix_ratio = 1`、`eps = 0.02`、`minPts = 20`，与当前 ablation 默认值一致；实际流程还显式使用 Euclidean metric 和运行线程数。
- **关键公平性问题：** Cohort d1 的列天然具有 `tissue|cohort|feature` 分块，而 Direct-GSClassifier 没有同构的 tissue block。直接把全局 GSClassifier 输入送入现有 two-stage helper 会把“是否有 tissue-first 结构”与“是否有 cohort transformation”混在一起。
- **`record` 的实际边界：** 参数表保存的是 DR/DBSCAN 数值超参数；`record` 在现有代码中把一次 DBSCAN 解的原始 cluster ID 合并为 0–3 四类 `normSubtype`。原始 cluster ID 会随 arm、seed 和参数变化，不能把同一个 `record` 字符串直接套到 Direct 结果。

## 要达到什么目标

- **完成后的变化：** `ablation()` 新增独立的 `experiment = "metaccs"`，在同一批样本、同一随机计划、同一 two-stage UWOT 预算和同一 DBSCAN 预算下，直接比较 `Cohort-d1` 与 `Direct-GSClassifier` 两条端到端流程。
- **可回答的问题：** 输出应能区分生物结构、跨 cohort 混合、降维保真、聚类稳定性和聚类工程特征，形成 `Cohort-d1 - Direct-GSClassifier` 的配对差值，而不是用“图更成簇”代替任务价值。
- **兼容性：** 保持现有 `cohort`、`scaling`、`tissue_first` 行为与结果文件不变；新实验不受 Gate 1 阻断，因为它本身就是验证 cohort representation 是否有下游价值的直接证据。
- **可追溯性：** 每个结果都记录样本、Direct feature block、DR 参数集、DBSCAN 参数集、resample、UMAP seed 和输入 hash，并保存逐样本 cluster assignment。

## 不在本次处理范围

- 不重新训练 cohort submodels，也不改变 `CCS::dr()` 或 `CCS::cluster()` 的公共行为。
- 不把 raw DBSCAN cluster 数量、低 noise fraction 或更好看的二维图单独定义为成功。
- 不默认把现有 `record` 作为 Direct arm 的四类 metaCCS 真值或后处理规则。
- 不在这一改动中解决 d1 是否严格 out-of-fold 生成的问题；manifest 和文档继续披露该限制。
- 不把 transductive UMAP/DBSCAN 结果写成未见 cohort 的归纳预测能力；若需要该主张，应另行增加可 transform 的 held-out-cohort 设计。

## 改进方向

### 定义不混杂 tissue-first 的 Direct 对照

新增 tissue-level GSClassifier feature manifest。对每个 tissue，汇总其冻结 cohort models 实际使用过的单基因分箱、普通基因对和 gene-set 对比，并由 GSClassifier 原生预处理函数统一重建，形成 `Direct-GSClassifier` tissue blocks；该 arm 不读取任何 d1 数值，也不使用 biology 标签。

两条主路径定义为：

```text
Direct-GSClassifier:
各 tissue 冻结特征集上的原始 GSClassifier 输入
→ tissue-specific UWOT d2
→ global UWOT d3
→ column scaling
→ DBSCAN

Cohort-d1:
冻结 cohort models 的 d1，按 tissue 分块
→ tissue-specific UWOT d2
→ global UWOT d3
→ column scaling
→ DBSCAN
```

这种 Direct baseline 保留与 Cohort arm 相同的冻结特征支持和 tissue-first 架构，只消除 cohort model transformation。特征在多个 tissues 中出现时允许进入多个 tissue blocks，但必须在 feature manifest 中显式记录；不能静默复制完整全局输入到每个 tissue。

第一阶段的 tissue 顺序、随机种子、目标维数和 `n_neighbors` 应联合确定。若任一 arm 的某 tissue 因列数、去重后样本数或有效秩不足而需要降档，两组都使用共同可实现的较小值，避免一组获得更多 d2 维数或更大邻域预算。全局第二阶段采用同样的共同约束。

### 增加固定参数与等预算两种比较模式

固定参数模式回答“当前实际 CCS 配置是否优于 Direct”。两组共享目标参数表中的 UWOT 与 DBSCAN 参数；包函数只接收已经解析的参数，不读取项目外部的参数 RDS。建议把 `metric = "euclidean"` 纳入 `params$dr`，把 `n_threads` 作为已审计的运行参数，而不是默认为固定机器上的 24。

仅使用当前参数仍可能偏向原 Cohort 流程，因为这些参数最初在 Cohort 几何上选定。因此增加等预算 grid 模式作为必要敏感性分析：

- 两组评估完全相同的 DR parameter sets 和 DBSCAN parameter sets；
- 每个组合拥有稳定的 `dr_param_set_id` 和 `cluster_param_set_id`；
- 相同 d3 可复用于多个 DBSCAN 参数，避免重复计算；
- 结果按相同参数集配对，不在包内事后挑选最有利于某一 arm 的单个组合；
- 如以后需要自动选参，必须另行预注册共同的无监督目标、tie-break 和 tuning/evaluation 切分。

公共配置建议在现有 `params$dr`、`params$cluster` 基础上增加一个小型 `params$metaccs` 节点：

- `resample_seeds`、`umap_seeds` 和 `subsample_fraction`；
- `parameter_mode = "fixed"` 或 `"grid"`；
- 可选的 `dr_grid`、`cluster_grid`；
- `direct_feature_mode = "tissue_model_union"`，第一版只支持这一种有明确定义的主对照；
- `retain_assignments = TRUE`，默认保留轻量 cluster assignment，不默认保存所有 d3 矩阵。

配置验证应在重型计算前检查：两个阶段的维数、UWOT 参数、DBSCAN 参数、seed 非空、subsample 范围、grid 列名和 Direct tissue feature coverage。

### 新增 Experiment 4 的执行与结果合同

新增 `.ablation_experiment_metaccs()` 作为短入口，复用现有样本准备、分层抽样、two-stage reduction、DBSCAN、bootstrap 和稳定性 helpers。每个 `resample_id × umap_seed × parameter_set` 先生成一次共同样本索引，再依次运行两个 arms，保证严格配对。

新实验至少返回：

- `metrics`：逐 arm、resample、seed、参数集的长表；
- `summary` 与 `contrasts`：先在同一 resample 内聚合 UMAP seeds，再以 resample 为统计单位计算 `Cohort-d1 - Direct-GSClassifier`；只有全样本 seed 重复时应明确标为算法波动，而不是生物抽样 CI；
- `stratified`：逐 tissue 的邻域保真、purity 和 mixing，避免总体结果掩盖小 tissue 损失；
- `stability` 与 `stability_contrasts`：各 arm 跨 seed/resample 的 neighbor Jaccard、ARI 和双向 cluster Jaccard及其配对差值；
- `cross_arm_agreement`：相同 repeat 下两组聚类的 ARI/cluster Jaccard，只描述两条流程改变了多少，不把高一致性自动解释为更好；
- `assignments`：`sample_id`、arm、resample、seed、参数集、raw cluster 和 noise 标记；
- `parameter_manifest` 与 `audit`：完整参数、共同维数/邻域降档、Direct feature hash、输入 hash 和 seed 衍生关系。

结果保存为 `experiment-04-metaccs.rds`，并挂到 `CCSAblation$experiments$metaccs`。现有 `audit.csv` 使用列并集扩展，不改变旧实验字段含义。

### 用任务相关指标判定是否真正增益

复用现有指标：trustworthiness、continuity、tissue kNN retention、biology purity、cohort mixing、cluster count、cluster-size entropy 和 noise fraction。增加 raw cluster 对独立 biology 标签的 label-invariant 指标，至少包括 non-noise 样本上的 ARI、NMI或等价信息量指标，以及加权 cluster purity；同时报告 non-noise coverage，防止通过把困难样本变成 noise 虚增一致性。

端到端增益不能由单一指标决定。支持 Cohort-based 的结果应同时满足：

- biology concordance 或独立临床/pathway anchor 有预注册的实质提高；
- 跨 seed/resample cluster stability 不下降；
- noise、几何保真和小 tissue 指标没有不可接受的损失；
- mixing 的提高不是由 biology purity 或同 subtype 跨 cohort 邻域损失换来的；
- 固定参数结论在共同参数 grid 中不是孤立的单点优势。

新实验第一版只输出各方向的配对证据，不生成未经预注册权重的综合分数，也不自动宣布 winner。

### 把 `record` 放在正确的后处理层

现有 `record` 可作为当前 Cohort reference solution 的版本化注释写入分析 manifest，但不得按 raw cluster ID直接应用到 Direct 或其他 seed。当前 Experiment 4 的主要聚类比较使用 ARI、NMI和 cluster Jaccard等对标签编号不敏感的指标。

若论文主终点必须是最终四类 normalized metaCCS，需要后续单独实现 arm-neutral normalization：在独立 tuning data 上冻结“raw cluster → 4类”的生物学规则或 reference transfer，再对两个 arms 的 evaluation data 使用同一规则。该规则未被编码之前，报告应把 Experiment 4 称为 raw metaCCS clustering comparison，而不是四类 normalized metaCCS 的最终比较。

## 实施范围与顺序

1. 先稳定当前正在进行的 `test/ablation/ablation.R → R/ablation.R` 迁移，确认测试不再 source 已删除文件；保留现有未提交修改，不覆盖用户工作。
2. 扩展冻结 feature manifest，使其同时提供完整 GSClassifier feature union、TSP 子集、module-to-feature 和 tissue-to-feature 映射，并在 manifest 中加入对应 hash 与 `breakVec`。
3. 抽出成对 two-stage reduction helper，实现共同 tissue 顺序、共同可实现维数、共同 `n_neighbors` 和相同派生 seed；在其上建立 `Direct-GSClassifier` 与 `Cohort-d1` arms。
4. 接入 `experiment = "metaccs"`、独立结果文件、assignments、parameter manifest 和扩展 audit；该分支不依赖 Gate 1。
5. 先完成 frozen fixed-parameter 路径，再增加相同 parameter grid 的批量评估与缓存，避免在核心流程尚未对齐时引入选参复杂度。
6. 增加 cluster biology、稳定性和配对差值指标，明确 resample 与 UMAP seed 的统计层级。
7. 更新 roxygen/Rd、ablation 方法说明、测试和 CHANGELOG；由于新增向后兼容的公共能力，实施时按项目 SemVer 规则评估从 0.8.x 提升到下一个 minor 版本。

## 受影响的组件

- `R/ablation.R`：公共 experiment 注册、配置、feature manifest、paired pipeline、指标、结果和 audit。
- `test/ablation/test-ablation.R`：合成数据下的结构、公平性、确定性、边界和回归测试。
- `test/ablation/test-ablation-real-layered.R`：真实冻结模型上的最小端到端 smoke test。
- `test/ablation/ablation.md`：Experiment 4 的问题、路径、统计单位、`record` 边界和解释规则。
- `man/ablation.Rd`、`NAMESPACE`：由 roxygen 重新生成；若只扩展现有导出函数，NAMESPACE 不应产生新的公开符号。
- `DESCRIPTION`、`CHANGELOG.md`：实施完成后的版本与变更记录。

## 如何确认完成

- 使用目标 `paramMD5` 的参数行运行 Cohort arm时，在相同样本、seed 和无动态降档条件下，d3 与 `CCS::dr(..., cover = TRUE)`、cluster 与 `CCS::cluster(eps = 0.02, minPts = 20)` 达到预期的确定性等价。
- 合成测试证明 Direct arm 只读取重建的 GSClassifier 输入和 tissue feature manifest，不读取 d1 数值或 biology 标签；人为改变 d1 只影响 Cohort arm。
- 两组每个 repeat 的 sample IDs、tissue 顺序、共同第一阶段维数、共同邻域上限、UMAP seeds、最终维数和 DBSCAN 参数完全相同，并由 audit 断言。
- 固定参数模式产生两个 arms、全部核心指标、配对 contrasts、stability、cross-arm agreement 和逐样本 assignments。
- Grid 模式对每个 parameter set 两组行数一致，d3 对 DBSCAN grid 被复用，不把参数组合或 UMAP seeds误当作独立 cohort 样本。
- 当 cluster IDs 被任意重编号时，ARI/NMI/cluster Jaccard不变；noise-only、单 cluster、缺少 biology 类别和小 tissue 均有清晰结果或显式错误。
- 旧三个 experiments 的合成测试和真实 smoke test继续通过，已有结果对象字段保持兼容。
- 执行 `devtools::document()`、核心 ablation scripts、`R CMD build` 和 `R CMD check`；大型真实数据只写入测试临时目录，不进入 Git。

## 风险与待确认事项

- **Direct baseline 的估计对象：** `tissue_model_union` 保留冻结模型使用过的全部 GSClassifier 输入特征，因此比较的是“在相同冻结特征支持下，cohort model transformation 是否有增益”，不是“完全原始、未选择的所有表达特征是否更好”。
- **参数选择偏倚：** 当前数值参数来自 Cohort 流程，fixed 模式可能偏向 Cohort；共同 grid 是等机会结论的必要敏感性分析。
- **四类 normalized metaCCS：** `record` 依赖原始 cluster IDs且目前没有 arm-neutral 生成规则。它不阻塞 raw DBSCAN solution 的端到端比较，但阻塞“最终四类 metaCCS 已被公平比较”的强表述。
- **统计单位：** UMAP seeds只是算法重复；只有独立 resamples、cohort resamples或外部数据才支持更强的不确定性推断。
- **计算量：** two-stage UWOT 按 tissue、arm、resample、seed 和 DR grid 成倍增长，需要缓存 d3并限制默认 grid；不得通过减少某一 arm 的预算来提速。
- **当前工作树：** `R/ablation.R`、文档、NAMESPACE、DESCRIPTION 和测试迁移已有未提交修改。实施者必须基于这些修改增量工作，并在开始前确认删除的旧测试源码已被正确替代。
