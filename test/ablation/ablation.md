# `ablation()` 函数设计说明

## 设计目标

`ablation()` 是挂接在现有 CCS 体系旁边的实验入口。它不重新训练 cohort submodels，也不修改 `CCS` 对象，而是读取冻结的 TSP 特征、完整 d1 和样本注释，自动执行 CCS 分层与端到端消融计划中的四个核心实验：

- cohort representation 是否有必要；
- cohort module 数量是否存在 scaling；
- tissue-first 两步降维是否有必要；
- cohort transformation 能否穿过 two-stage UWOT 和 DBSCAN，形成更好的 raw metaCCS solution。

之所以采用“冻结模型后的评估器”，而不是新建一套训练管线，是因为本轮需要隔离的是 CCS 各层的作用。若在每个消融组中重新训练模型，组间差异会同时混入训练随机性、模型容量和调参机会，很难再归因于被删除的层。冻结已有 150 个 cohort modules 后，实验只改变表示路径、模块子集或降维层，其余条件可以保持一致。

这里的目标是评价 cohort bank 是否把 TSP 重表达成更稳定、更容易解码生物结构的空间，而不是声称 d1 创造了 TSP 中原本不存在的信息。d1 是 TSP 的确定性函数；有意义的证据必须同时包括生物效用、cohort 混合、生物结构保护和零模型对照，不能只看 d1 与 TSP 是否不同。

## 总体架构

主函数放在文件最前面，只负责组织流程；所有实验和通用能力都放在后续私有函数中。当前调用关系可以概括为：

```text
ablation()
├─ validate and merge configuration
├─ prepare frozen inputs
│  ├─ recover module blocks
│  ├─ recover model-used TSP features
│  ├─ align expression, d1, and metadata
│  └─ create manifest and hashes
├─ Experiment 1: cohort representation
│  ├─ Direct / Cohort / Null-RP / Null-Perm
│  ├─ paired fold summaries and contrasts
│  └─ Gate 1
├─ Experiment 2: cohort-axis scaling       [records Gate 1; enforcement is optional]
│  ├─ nested module sequences
│  ├─ d1 scaling curves and saturation fit
│  └─ representative d3/cluster stability
├─ Experiment 3: tissue-first              [independent of Gate 1]
│  ├─ Two-stage embedding
│  └─ One-stage embedding
├─ Experiment 4: end-to-end metaCCS        [independent of Gate 1]
│  ├─ tissue-model-union Direct-TSP / Cohort-d1
│  ├─ paired two-stage reduction and DBSCAN
│  ├─ fixed parameters / shared parameter grid
│  └─ assignments, stability and cross-arm agreement
└─ combine audit rows and save CCSAblation
```

这个结构有三个目的：

- 公共 API 保持稳定，实验细节不会不断膨胀到 `ablation()` 中；
- 每个实验入口都可以单独测试，通用统计函数也能被多个实验复用；
- 将来新增实验时，主要增加一个 `.ablation_experiment_*()` 模块，而不需要改动 CCS 包现有函数。

## 公共接口

`ablation()` 只保留八个公共参数：

| 参数 | 作用 | 这样设计的原因 |
|---|---|---|
| `object` | 提供冻结 d1、module 信息和模型位置的 `CCS` 对象 | 复用 CCS 当前数据结构，不再定义平行对象体系 |
| `data` | 原始 RNA 表达量或 tissue/cohort 嵌套列表 | 实验一必须从真实表达量重建模型使用的 TSP，才能建立 Direct 对照 |
| `metadata` | 可选样本注释 | 嵌套数据可自行派生注释；单矩阵输入则必须显式提供 |
| `experiment` | 选择 `cohort`、`scaling`、`tissue_first`、`metaccs` | 允许分阶段运行，避免每次都执行全部高成本实验 |
| `output.dir` | 独立输出目录 | 避免覆盖现有 CCS 产物，并让一次运行的所有文件集中归档 |
| `params` | 递归覆盖 `.ablation_default_params()` 指定叶节点的命名嵌套列表 | 可只修改某个实验或 Gate cutoff，而不必复制整套配置 |
| `seed` | 主随机种子 | 从一个入口派生所有随机过程，保证整次实验可复现 |
| `verbose` | 控制进度信息 | 兼容交互运行和脚本批处理 |

`.ablation_merge_lists()` 对 `params` 做递归合并。因此用户可以只覆盖一个叶节点，例如 `list(gate1 = list(enforce = TRUE, min_gain = 0.02))`，不必重复提供整套 Gate 配置。这种接口比不断给主函数增加参数更容易维护，也能让以后新增实验参数时保持向后兼容。Gate cutoff 的具体含义与使用建议见下文“Gate 1”。

## 主函数的执行顺序

`ablation()` 按以下顺序运行：

1. 检查 `object` 是否为 `CCS`，解析实验名称，合并并验证配置。
2. 检查输出目录。目录非空时默认停止，只有 `params$cover = TRUE` 才允许覆盖同名结果。
3. 调用 `.ablation_prepare_input()` 对齐 d1、RNA/TSP 和 metadata。
4. 在任何耗时实验前保存 `manifest.rds` 和 `config.rds`。
5. 按依赖关系运行实验。`scaling` 会自动先运行实验一并计算 Gate 1；默认继续探索，只有显式设置 `params$gate1$enforce = TRUE` 且 Gate 未通过时才停止。`tissue_first` 与 `metaccs` 独立运行。
6. 合并各实验审计表，构造 `CCSAblation` 对象，并保存 CSV/RDS 结果。

先保存 manifest 和配置，是为了即使长时间计算中途失败，也能知道该次运行使用了什么样本、特征、模块和参数。实验结果分别保存，而不是只在最后保存一个大对象，也降低了单个阶段失败后排查问题的难度。

## 输入准备

### 冻结 module 边界

`.ablation_module_manifest()` 从 d1 列名的 `tissue|cohort|feature` 约定恢复 module block。每个 `tissue|cohort` 对应一个完整 block，记录其列范围和宽度。

block 必须保持完整，因为一个 cohort model 的多列输出共同表示该模型。Null-Perm 和 scaling 都以完整 block 为单位；若拆成单列置换或抽样，就会同时破坏 module 内部结构，消融因素不再单一。

### 模型实际使用的 TSP

`.ablation_extract_tsp_features()` 优先从 `object@Model` 读取所有重复和类别模型使用的 gene pairs；模型未内嵌时，再通过 `.ablation_model_path_map()` 读取冻结的 `modelFit.rds`。Direct、Null-RP 和机制分析只使用这些特征的并集。

这样做是为了保持输入公平。若 Direct 路径额外获得 cohort bank 从未使用过的 gene pairs，它与 Cohort 的差异就不再只来自 cohort transformation。

### 样本对齐与重复样本

`.ablation_flatten_expression()` 将单矩阵或嵌套列表统一成 expression matrix 和 metadata。嵌套列表合并前先取共同基因集合。随后 `.ablation_prepare_input()` 取 d1 行名、表达矩阵列名和 metadata sample ID 的交集。

重复 sample ID 会整组排除，而不是任意保留第一列。原因是重复 ID 无法保证 d1、表达量和标签之间的一一对应；任意选择可能形成静默错配。排除数量会写入 manifest。

若表达矩阵行名不能直接覆盖模型基因，函数复用 `GSClassifier::geneMatch()` 和 CCS 对象已有的 gene annotation 规则，而不是另写一套 ID 转换。最终 TSP 用 `geneA >= geneB` 编码为 1。

### metadata 规范化

`.ablation_prepare_metadata()` 识别常见别名并统一为：

- `sample_id`；
- `cohort`；
- `tissue`；
- `biology`。

缺失 `tissue` 时尝试从 `object@Data$CancerType` 补充；缺失 `biology` 时默认使用 `tissue`。默认值保证基础流程可以运行，但正式生物学结论最好提供独立、样本级且跨多个 cohort 分布的 `biology` 标签，否则 tissue 与 cohort 的混杂可能限制解释强度。

## 配置与随机性

`.ablation_default_params()` 集中管理四个实验的默认参数。不同随机过程使用显式命名、带偏移的 seed 集合：

- `rp_seeds`：Null-RP；
- `permutation_seeds`：Null-Perm；
- `scaling_embedding_seeds`：实验二的下游嵌入；
- `tissue_seeds`：实验三的抽样、降维和聚类重复；
- `metaccs$resample_seeds`：实验四的 tissue × cohort 分层样本构成重复；
- `metaccs$umap_seeds`：固定样本集上的算法重复。

这种显式拆分不是统计要求，而是审计设计。看到一个异常结果时，可以直接定位是哪类随机过程和哪个 seed；新增随机操作时也可以为其分配新的派生规则，避免无意改变已有实验。

大样本上需要限量时，`.ablation_stratified_sample()` 在 tissue × cohort 分层之间轮转抽样，而不是简单全局抽样。它用于 `max_samples`、几何评估、机制分析和嵌入保真度子集，目的是避免大 cohort 垄断计算子集。

## 实验一：cohort representation

### 科学问题

实验一回答：在相同输入和相同下游表示预算下，冻结的 cohort bank 是否产生了比直接压缩 TSP 更有用的表示。

它包含四组：

| 组 | 当前实现 | 排除的替代解释 |
|---|---|---|
| `Direct` | 训练 fold 标准化 TSP，再拟合固定秩 PCA | 不经过 cohort bank 时的能力 |
| `Cohort` | 对冻结 d1 做训练 fold 标准化和同秩 PCA | 完整 cohort-informed representation |
| `Null-RP` | 标准化 TSP 乘以稀疏 Achlioptas 随机投影矩阵 | 收益是否只是投影或维度预算造成 |
| `Null-Perm` | 在每个 cohort 内独立置换完整 d1 blocks，再做同秩 PCA | 收益是否只来自高维边际分布，而不依赖正确样本对应 |

Null-RP 不交换样本，只把一个样本自己的 TSP 随机组合成等维特征。Null-Perm 则保留每个 block 的数值和内部多列结构，但破坏不同 blocks 属于同一样本的对应关系。两种 Null 回答不同问题，因此不能合并。

### grouped validation 与共同秩

`.ablation_grouped_folds()` 以完整 cohort 为验证单位，并尽量平衡各 fold 的样本量。同一 cohort 不会同时进入训练和测试集合。`.ablation_scale_train_apply()` 和 `.ablation_fit_pca()` 只用训练 fold 估计均值、标准差和 PCA rotation，再原样应用到测试 cohort。

Direct、Cohort 和 Null 使用共同目标秩 `rank_q`。实际秩会根据特征数、训练样本数和各 fold 的可实现秩动态截断。这样比较的是相同表示预算下的质量，而不是谁拥有更多列。默认主秩为 50，并在 25、50、100 上做可行范围内的敏感性分析。

### 指标体系

实验一不使用单一总分，而是保留互补指标：

- `linear_cka`：整体几何是否相似；
- `distance_spearman`：样本远近顺序是否一致；
- `knn_jaccard`：局部近邻集合是否一致；
- `effective_rank`：表征实际使用了多少维度；
- `cohort_mixing`：同一 biology 内近邻来自不同 cohort 的比例；
- `biology_purity`：近邻保持相同 biology 的比例；
- `macro_auroc` 和 `balanced_accuracy`：低容量 linear probe 的跨 cohort 解码能力；
- `discordant_rank_change`：biology 不同的 TSP 近邻在候选空间中是否被推远；
- `concordant_cross_cohort_rank_change`：同 biology、跨 cohort 的 TSP 近邻如何变化。

`cohort_mixing` 必须在相同 biology 内计算，否则组织差异会被误写成 cohort 分离。linear probe 只保留至少出现在两个 cohort 的标签，避免模型直接利用 cohort 身份猜标签。

`.ablation_dimension_free_geometry()` 另外直接比较原始 TSP 和完整 d1。CKA、距离相关和 kNN Jaccard 不要求列数相同，因此这部分用于证明“表示确实发生了什么变化”；固定秩实验再回答“这种变化是否更有用”。

### 统计汇总

`.ablation_metric_summary()` 先对同一 fold 内的随机 seed 重复取平均，再以 fold 为统计单位执行 1,000 次 bootstrap。这样不会把 20 个随机投影或置换 seed 误当成 20 个独立 cohort。

`.ablation_paired_contrasts()` 在相同 fold 上计算：

- `Cohort - Direct`；
- `Cohort - Null-RP`；
- `Cohort - Null-Perm`。

差值为正表示 Cohort 的该指标更高。是否“更高更好”仍取决于指标含义；例如 CKA 只描述相似程度，effective rank 也不是性能指标，不能单独用于放行。

### Gate 1

Gate 1 是实验一与实验二之间的证据判定规则，不是新的消融组或性能指标。它回答的是：在研究“增加多少 cohort modules 才足够”之前，完整 cohort representation 是否已经显示出相对 Direct 和零模型的可信增益，同时没有明显破坏生物结构或 cohort 混合。

设计 Gate 1 的原因是实验二只衡量不同 module 子集对完整 d1 的逼近、饱和与下游稳定性，不能证明完整 d1 本身有用。即使 d1 没有生物学价值，增加 modules 也可能产生平滑的饱和曲线。Gate 因此用于防止把“更好地重建一个无效表示”误解为 cohort bank 具有可扩展的生物学收益；启用强制模式时，它也可以避免在前置证据不足时继续投入高成本计算。

`.ablation_gate_one()` 默认以 `balanced_accuracy` 为 `primary_metric`，使用以下三个配对差值 95% CI 的下界，而不是只看点估计：

```text
Cohort - Direct
Cohort - Null-RP
Cohort - Null-Perm
```

三个下界都必须不低于 `min_gain`。同时，`biology_purity` 和 `cohort_mixing` 的 `Cohort - Direct` 下界必须分别不低于 `-purity_tolerance` 和 `-mixing_tolerance`。若小样本或类别不足导致 `balanced_accuracy` 无法估计，当前实现回退到 `biology_purity`。

默认 `gate1$enforce = FALSE`：Gate 仍会计算并保存在 scaling 结果中，但无论通过与否都会继续实验二。这一默认值适合方法开发、敏感性分析和失败诊断，可以保留更丰富的探索性结果。探索性 scaling 结果必须与 Gate 状态一起报告；Gate 未通过时，不能把饱和曲线表述为实验一已经支持 cohort representation。

建议仅在以下场景设置 `gate1$enforce = TRUE`：

- 确认性分析中，需要把实验一作为预注册的 go/no-go 条件；
- scaling 或下游 embedding 计算成本很高，需要在前置证据不足时提前停止；
- 多批数据使用同一冻结决策规则，且 cutoff 已在查看本批结果前确定。

cutoff 应按指标原始尺度解释：

- `min_gain` 是三个主要对比共同要求的最小 CI 下界。`0` 表示要求 Cohort 的增益在当前置信水平下不低于零；正值（如 `0.02`）表示要求至少两个百分点的最小可信增益。
- `purity_tolerance` 和 `mixing_tolerance` 是允许的最大退化量。`0` 表示不接受可信退化；`0.01` 表示相应指标的 `Cohort - Direct` CI 下界最低可到 `-0.01`，即容忍最多一个百分点的下降。
- cutoff 应依据领域中的最小重要差异、历史重复实验的自然波动或独立 pilot 数据预先确定，不能在查看当前结果后为了通过 Gate 而调整。

例如，下面的配置要求主要指标至少有 0.02 的可信增益，同时允许 biology purity 和 cohort mixing 各下降最多 0.01：

```r
params = list(
  gate1 = list(
    enforce = TRUE,
    primary_metric = "balanced_accuracy",
    min_gain = 0.02,
    purity_tolerance = 0.01,
    mixing_tolerance = 0.01
  )
)
```

## 实验二：cohort-axis scaling

### 为什么使用冻结 blocks 的嵌套序列

实验二不训练 10、25 或 50 cohort 的新模型，而是从完整 d1 中抽取冻结 module blocks。这样模块数量是唯一主要变化因素，也避免不同规模拥有不同训练预算。

`.ablation_nested_module_sequences()` 先按 tissue 分层，在 tissue 内随机排列 modules，再用归一化次序交错各 tissue。每条序列都满足：较大规模严格包含较小规模。默认评估 10、25、50、75、100、125、150 modules，共 100 条序列。

嵌套关系非常重要。若 25-module 和 50-module 组来自两个无关随机子集，二者差值会混入 module 组成差异；当前设计在同一 sequence × fold 内比较相邻规模，可以把差值解释为增加一批 modules 的边际变化。

### d1 scaling 曲线

每个 module 子集都与完整 d1 在同一 fold、同一 `rank_q` 下拟合 PCA。当前指标包括：

- 对完整 d1 的 CKA 和 kNN 恢复；
- effective rank；
- cohort mixing 和 biology purity；
- linear probe 的 AUROC 和 balanced accuracy；
- module 数和累计 d1 列数。

同时报告 module 数与累计维数，是因为不同 cohort block 的宽度可能不同。只看列数会把“某些 module 输出更宽”误写成“cohort 数更多”。

`.ablation_scaling_summary()` 在相同 sequence × fold 内计算相邻规模增量并 bootstrap CI。`.ablation_saturation_fit()` 拟合：

$$
S(m)=S_{\infty}-a\exp(-m/\tau).
$$

代码对 `tau` 做网格搜索，并在线性条件下求 `S_inf` 和 `a`，比直接非线性优化更稳定。较小的 `tau` 表示更早接近平台。曲线只描述当前冻结 bank 内的 retrospective scaling，不能代替真正加入第 151 个独立 cohort 后的外推验证。

### 下游 d3 与聚类

完整 d1 的几何恢复不等于下游 metaCCS 更稳定。因此 `.ablation_scaling_embedding()` 默认在 25、50、100、150 modules 上选取 10 条序列，每条使用 10 个 seed，重新运行当前 Two-stage d3 和 DBSCAN。

同一 seed 在不同 module 数下使用相同的 tissue × cohort 分层样本子集。输出包括 d1→d3 的保真度、组织内近邻保留、混合/纯度、聚类数量、簇大小熵和噪声比例。`.ablation_scaling_embedding_stability()` 在相同 sequence 和 module 数内两两比较 seed，报告近邻 Jaccard、ARI 和 cluster Jaccard。

聚类数量和 entropy 只作描述。更多 clusters 也可能来自噪声或碎片化；只有邻域、生物结构和重采样稳定性同时改善，才能把变化解释为分辨率提高。

## 实验三：tissue-first

### 配对路径

实验三固定完整 d1，只比较：

| 组 | 路径 |
|---|---|
| `Two-stage` | `d1 → tissue-specific d2 → global d3` |
| `One-stage` | `d1 → direct global d3` |

每个 seed 下，两组共享完全相同的分层样本子集、最终维数、UMAP 参数和 DBSCAN 参数。因此组间唯一主要差异是是否保留 tissue-first 层。

`.ablation_two_stage_embedding()` 从 d1 列名提取 tissue reference，先调用 `.ablation_reduce_by_reference()` 对每个 tissue block 独立降到第一阶段维数，再拼接并做全局降维。默认维数为 `c(5, 2)`。小 tissue 或低秩 block 会动态降低第一阶段维数，并限制 `n_neighbors`，避免要求不可实现的秩。

`.ablation_one_stage_embedding()` 跳过 tissue 分块，直接将完整 d1 降到最终二维。两条路径都通过 `getFromNamespace()` 复用 CCS 的 `drCCSProbability`、`data_for_dr`、`CORE_DR` 和 `repairCCS`，目的是最大限度保持与现有 CCS 降维行为一致，而不是在消融脚本中复制一套近似实现。

### 指标与稳定性

DBSCAN 在逐维标准化的 d3 上运行，标签 0 作为噪声。`.ablation_embedding_metrics()` 计算：

- trustworthiness：惩罚进入低维近邻的高维远邻；
- continuity：惩罚在低维中丢失的高维近邻；
- tissue kNN retention；
- cohort mixing 和 biology purity；
- cluster count、cluster-size entropy 和 noise fraction。

`.ablation_tissue_stratified_metrics()` 逐 tissue 计算保真与混合指标，并将样本数处于下四分位的 tissue 标记为 `small_tissue`。这样可以识别“总体均值改善，但小 tissue 被过度压缩或分割”的副作用。

`.ablation_embedding_stability()` 在 Two-stage 和 One-stage 内分别两两比较 seeds。它先对齐共同样本，再计算 d3 近邻 Jaccard、adjusted Rand index 和双向 cluster Jaccard。最终 `.ablation_two_group_contrasts()` 在相同 seed 上计算 `Two-stage - One-stage` 的配对 CI。

## 实验四：端到端 raw metaCCS

### 公平的 Direct-TSP 对照

实验四比较两条共享 tissue-first 架构的路径：

```text
Direct-TSP:
各 tissue 冻结模型实际使用过的 TSP 并集
→ tissue-specific d2
→ global d3
→ column scaling
→ DBSCAN

Cohort-d1:
冻结 cohort model d1，按 tissue 分块
→ tissue-specific d2
→ global d3
→ column scaling
→ DBSCAN
```

`.ablation_frozen_feature_manifest()` 同时恢复 module-to-feature 和 tissue-to-feature 映射。Direct-TSP 只读取 TSP 与该映射，不读取 d1 数值或 biology 标签；同一个 TSP 可因被多个 tissue banks 使用而进入多个 tissue blocks，这种重复会记录在 tissue feature manifest 与 hash 中。

两组每个 tissue 的第一阶段维数、`n_neighbors`、tissue 顺序和派生 seed 联合确定。若任一组因列数、有效秩或去重后样本数无法实现目标预算，两组共同使用较小值。全局第二阶段采用相同规则。共同降档写入 `parameter_audit`，因此不会出现一组在低秩数据上获得更多维度或邻域预算的情况。

### 固定参数与共同 grid

`params$metaccs$parameter_mode = "fixed"` 使用 `params$dr` 和 `params$cluster` 的当前 CCS 配置，回答现有流程下的端到端差异。`"grid"` 模式用 `dr_grid` 与 `cluster_grid` 覆盖同一套基础参数，两组评估完全相同的笛卡尔积。

每个 DR 和 DBSCAN 参数集根据完整配置生成稳定 ID。相同 `resample × UMAP seed × DR set` 的 d3 只计算一次，再复用于对应的 DBSCAN grid；`dr_cache_key` 记录这一复用边界。函数不根据结果自动挑选参数，也不生成未经预注册的综合 winner 分数。

### 指标、配对与统计单位

实验四复用 trustworthiness、continuity、tissue kNN retention、cohort mixing、biology purity、cluster count、cluster-size entropy 和 noise fraction，并增加非噪声样本上的 label-invariant 指标：

- cluster 与独立 biology 标签的 ARI；
- normalized mutual information；
- 按 cluster 大小加权的 purity；
- non-noise coverage。

`metrics` 和 `stratified` 保留逐 arm、resample、UMAP seed 与参数集结果。`summary` 与 `contrasts` 先在同一 resample 内平均 UMAP seeds，再以 resample 为单位计算 `Cohort-d1 - Direct-TSP`。当 `subsample_fraction = 1` 时，即使提供多个 resample seeds，也只能标为 `algorithm_variation_only`，不能解释成生物抽样置信区间。

`stability` 在同一 arm 与参数组合内比较不同运行，报告共同样本上的近邻 Jaccard、ARI 和 cluster Jaccard；`stability_contrasts` 再计算两组配对差值。`cross_arm_agreement` 只描述同一运行中两条路径改变了多少，不把高一致性或低一致性自动解释为性能优劣。

### raw cluster 与 `record` 的边界

`assignments` 保存逐样本 raw DBSCAN cluster 与 noise 标记。现有 `record` 依赖某一次 Cohort solution 的原始 cluster ID，不会应用到 Direct-TSP、其它 seed 或 grid 参数。Experiment 4 因此是 raw metaCCS clustering comparison，不等价于四类 normalized metaCCS 的最终比较。

若后续需要比较最终四类标签，必须在独立 tuning data 上冻结 arm-neutral 的 normalization 或 reference-transfer 规则，再原样应用到两组 evaluation data。

## 与外部泛化计划的衔接

新增消融计划中的 d1 单样本确定性、filtered cohorts 外部泛化、小样本 learning curve、原始表达扰动和 model-bank 删除，需要从原始 RNA 重新执行冻结 model bank，并管理预测缓存与模型子集版本。当前 `ablation()` 只接收已冻结 d1，不会用占位结果假装完成这些实验。

本轮已把可由现有输入直接支持的样本构成敏感性纳入 `metaccs$resample_seeds`、逐样本 assignments 与跨运行 stability。外部样本重编码等实验应在明确冻结预测合同、cache key 和 filtered cohort 训练边界后另行接入。

## 输出与审计

每次运行都会生成 manifest、config、audit 和总结果；四个实验文件则只在请求相应实验时生成：

| 文件 | 内容 |
|---|---|
| `manifest.rds` | 样本数、TSP 特征、module blocks、版本和输入哈希 |
| `config.rds` | 完整合并后的实际配置 |
| `experiment-01-cohort.rds` | 实验一指标、汇总、配对差值和 Gate 所需信息 |
| `experiment-02-scaling.rds` | 实验二曲线、饱和拟合和下游 embedding 结果 |
| `experiment-03-tissue-first.rds` | 实验三整体、分 tissue 和稳定性结果 |
| `experiment-04-metaccs.rds` | 实验四指标、配对差值、参数审计、稳定性和逐样本 raw cluster |
| `audit.csv` | 可人工查看的逐指标审计长表 |
| `ablation-result.rds` | 完整 `CCSAblation` 返回对象 |

`CCSAblation` 包含调用、manifest、配置、各实验结果、统一 audit 和规范化输出路径。audit 保留 experiment/group/fold、随机种子、实际秩、k、module sequence、module count、输入哈希以及投影、置换、UMAP seed。其目的不是只保证“同一个 seed 能重跑”，还要能回答某一行结果来自哪批样本、哪组特征和哪条 module 序列。

## 子函数的职责分层

| 层次 | 关键函数 | 职责 |
|---|---|---|
| 配置 | `.ablation_default_params()`、`.ablation_merge_lists()` | 默认值与局部覆盖 |
| 输入 | `.ablation_prepare_input()`、`.ablation_module_manifest()`、`.ablation_extract_tsp_features()` | 冻结特征恢复与样本对齐 |
| 变换 | `.ablation_fit_pca()`、`.ablation_random_projection()`、`.ablation_permute_blocks()` | 公平表示与零模型 |
| 通用指标 | `.ablation_linear_cka()`、`.ablation_mixing_purity()`、`.ablation_probe()` | 几何、生物结构与解码能力 |
| 实验一 | `.ablation_experiment_cohort()`、`.ablation_cohort_rank_metrics()`、`.ablation_gate_one()` | cohort representation 与启动门槛 |
| 实验二 | `.ablation_experiment_scaling()`、`.ablation_scaling_summary()`、`.ablation_scaling_embedding()` | module scaling 与下游稳定性 |
| 实验三 | `.ablation_experiment_tissue_first()`、`.ablation_tissue_embeddings()` | Two-stage/One-stage 对照 |
| 实验四 | `.ablation_experiment_metaccs()`、`.ablation_paired_two_stage_embeddings()` | Direct-TSP/Cohort-d1 端到端公平对照 |
| 实验四参数 | `.ablation_metaccs_parameter_manifest()`、`.ablation_parameter_sets()` | fixed/grid 等预算参数与稳定 ID |
| 稳定性 | `.ablation_trust_continuity()`、`.ablation_embedding_stability()`、`.ablation_ari()` | 降维和聚类复现性 |
| 审计 | `.ablation_build_manifest()`、`.ablation_add_audit_hashes()`、`.ablation_rbind()` | 跨实验追溯与表结构统一 |

这种分层遵循“实验入口只组织本实验，通用 helper 不知道具体实验”的原则。例如 ARI、kNN Jaccard 和 paired contrast 不绑定某个组名，因此实验二和实验三可以直接复用。

## 如何增加新的消融实验

当前结构是显式注册，而不是自动扫描函数。新增实验时建议：

1. 在 `experiment` 的可选值中增加新名称；
2. 在 `.ablation_default_params()` 中增加该实验真正需要的配置；
3. 新建 `.ablation_experiment_<name>()`，让它只返回该实验的 metrics、summary、audit 和状态；
4. 优先复用现有输入、fold、指标、bootstrap 和审计 helper；
5. 在 `ablation()` 中增加一个短分支和独立结果文件；
6. 在合成测试中验证结构、随机性和边界，再用真实数据做最小运行验证。

显式注册需要改动少量主函数代码，但调用顺序、依赖和输出文件一目了然。对于只有四个到少量实验的脚本，这比引入复杂插件注册器更符合 KISS；当实验数量明显增多时，再考虑用实验注册表替代分支。

## 当前实现的边界

这份文档描述的是当前 `ablation.R` 的实际行为，以下边界需要在正式分析前知晓：

- 函数可以保证 grouped validation，但不能仅凭对象自动判断每个 d1 值是否来自严格 out-of-fold 预测；若样本参与过对应 cohort module 训练，相关结果需要标为描述性分析。
- 默认 `biology = tissue` 只是可运行的回退值，不等价于独立生物学验证。强结论应提供跨 cohort 的 outcome-independent sample-level anchors。
- 当前近邻和距离实现固定为 Euclidean。对二值 TSP，Euclidean 距离与 Hamming 距离的排序单调等价；当前没有实现 Jaccard 距离分支，`distance` 主要作为审计字段。
- module sequence 当前只按 tissue 分层，没有进一步按平台或训练样本量平衡；manifest 记录 block 宽度和累计维数，但不包含这些扩展属性。
- grouped folds 和样本集合在各组间一致，但部分逐组距离抽样和 probe 使用派生 seed，并不保证每一组使用完全相同的随机样本对或优化轨迹。正式主分析若要求最严格的随机配对，应进一步锁定并共享这些索引。
- 当前统计输出以 bootstrap CI 为主，没有实现计划中次要指标的 P 值或 FDR 校正。
- Gate 默认会计算但不强制阻断（`enforce = FALSE`），且 `min_gain = 0` 尚未编码项目特异的最小重要差异。确认性运行若启用 `enforce = TRUE`，应根据完整模型的自然波动或独立 pilot 数据预注册 cutoff，而不是看完结果后调整。
- tissue-first 通过 CCS 的非导出内部函数保持实现一致，因此比复制代码更忠实，但也会对 CCS 内部函数签名变化更敏感。
- metaccs 同样复用 CCS 的 `data_for_dr`、`CORE_DR` 和 `repairCCS`；共同降档保证组间公平，但内部函数签名变化仍需定向回归测试。
- Experiment 4 的 UMAP 与 DBSCAN 是 transductive 流程，输出是当前参与样本集合上的 raw clustering；不能写成未见样本的归纳分类性能。
- fixed 参数来自当前 Cohort 流程，可能偏向 Cohort-d1；共同 grid 是必要敏感性分析，但函数不会事后自动选取最有利参数。
- 当前没有 arm-neutral 的四类 cluster normalization，`record` 只可作为现有 Cohort reference solution 的版本化注释。

这些边界不妨碍脚本用于当前测试和方法开发，但它们决定了最终论文中可以使用多强的因果或泛化表述。

## 测试覆盖

当前测试分为三层：

| 测试 | 主要覆盖 |
|---|---|
| `test-ablation.R` | 合成 CCS 对象、module manifest、TSP 提取、grouped folds、Null 可复现性、四个实验的返回结构、Direct 隔离、共同预算、cluster 标签置换不变性和 grid 合同 |
| `test-ablation-real-data.R` | 真实 150 modules、496 TSP 特征、固定 400 样本以及实验一重复运行一致性 |
| `test-ablation-real-layered.R` | 真实数据上的 scaling embedding、Two-stage/One-stage 与 Direct-TSP/Cohort-d1 端到端路径 |

这些测试证明代码路径、返回结构和随机复现机制能够工作，但不替代对真实科学结论、标签独立性和最小重要差异的人工审查。
