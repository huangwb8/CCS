# `R/ablation.R` 设计与使用说明

`R/ablation.R` 是 CCS 的消融实验入口。它不重新训练 CCS 的 cohort 模型，而是把已经冻结的 `CCS` 对象拆成可比较的表示（Direct-GSClassifier 与 Cohort-d1），在相同样本、相同随机种子和相同评估预算下回答以下问题：

- `Cohort-d1` 是否保留了 Direct 表示的几何结构与邻域关系？
- 在跨 cohort 的外部查询上，哪一种表示更能检索到相同生物学标签的参考样本？
- d1 的模块数量、组织优先降维、完整 metaCCS 流程会如何改变结果？
- 观察到的收益是否超过随机投影（Null-RP）或块内置换（Null-Perm）控制？

脚本约 7,500 行，采用“公共入口 + 共享准备上下文 + 实验分支 + 统一落盘”的结构。本文按执行顺序解释它，而不是逐行罗列实现。

ablation-03 的新 targets 入口使用包导出的
`ablation_prepare_representation_inputs()`、
`ablation_run_representation_nodes()`、
`ablation_make_learning_curve_jobs()` 和
`ablation_make_scaling_jobs()`。这些接口默认以 `cache_dir = NULL` 在内存中
准备输入；targets 负责持久化、失效和恢复，包函数不创建 `SUCCESS`、stage
receipt 或 targets store。正式分析必须加载已安装的 `CCS` 包，不能在
`_targets.R`、Rmd 或 helper 中 `source()` 仓库内的 `R/ablation.R`。

## 先看整体结构

公共函数位于 `R/ablation.R:265`，内部函数均以 `.ablation_` 开头。可以把一次运行理解为下面的管线：

```text
CCS object + raw expression + metadata
                │
                ├─ experiment = "representation"
                │      └─ reference/external split
                │         → native geometry
                │         → cross-cohort retrieval
                │         → null controls
                │         → supervised readout + learning curve
                │         → optional scaling/decoder
                │
                └─ layered experiment vector
                       ∈ {cohort, scaling, tissue_first, metaccs}
                       └─ shared frozen context
                          → Experiment 1: cohort representation
                          → Gate 1
                          → Experiment 2: cohort-axis scaling
                          → Experiment 3: tissue-first reduction
                          → Experiment 4: end-to-end metaCCS
```

默认 `experiment = "representation"`。单独传入旧值 `"cohort"` 时，它只作为兼容别名映射到 representation，并发出弃用警告；`"cohort"` 出现在包含其它 layered 分支的实验向量中时，仍表示 layered Experiment 1。`representation` 不能与 layered 实验混用（`R/ablation.R:414-465`）。

## 公共入口如何分派

公共入口只做四件事：规范化实验名称、拒绝未知实验、处理兼容别名、把执行交给专用 orchestrator。核心代码如下：

```r
experiment <- unique(as.character(experiment))
layered_experiments <- c("scaling", "tissue_first", "metaccs")

if (identical(experiment, "cohort")) {
  warning("experiment = 'cohort' is deprecated; using 'representation'.")
  experiment <- "representation"
}

if (identical(experiment, "representation")) {
  return(.ablation_run_representation(
    object, data, metadata, output.dir, params, seed, verbose
  ))
}

.ablation_run_layered(
  object, data, metadata, experiment, output.dir, params, seed, verbose
)
```

这层不进行模型计算，因此扩展新实验时，优先修改分派和对应的实验函数，不要把新逻辑塞进公共入口。

## 输入契约与样本对齐

### `object`

必须是 `CCS` S4 对象，至少包含：

- `object@Data$Probability$d1`：已有的 d1 概率矩阵，行名是样本 ID；
- `object$Model`：按组织/模块保存的冻结模型；
- `object@Data$filtered.cohort`：默认用于划分外部 query cohort；
- `object$Repeat`：GSClassifier 的 feature、gene annotation 与模型元数据。

两个 orchestrator 都会在计算前检查 `methods::is(object, "CCS")`；公共入口负责分派，并由选中的 orchestrator 执行该检查（`R/ablation.R:474-528`、`5077-5089`）。

### `data`

支持两种形式：

1. 一个表达矩阵或 data.frame：行是基因，列是样本；必须另传 `metadata`。
2. 推荐的嵌套列表：`tissue -> cohort -> leaf`。`leaf` 可以直接是表达矩阵，也可以是 `list(expr = ..., subtype = ...)`。

```r
data <- list(
  Lung = list(
    TCGA = list(expr = lung_tcga, subtype = lung_label),
    GEO  = list(expr = lung_geo,  subtype = lung_label_geo)
  ),
  Breast = list(
    TCGA = list(expr = breast_tcga, subtype = breast_label)
  )
)
```

`.ablation_flatten_expression()`（`R/ablation.R:1222`）会完成以下工作：

- 将所有 leaf 的基因做并集；缺失基因填 `NA`；
- 生成 `sample_id/cohort/tissue/biology` 元数据；
- 发现重复样本 ID 时整行排除，并记录到 `excluded_duplicate_samples`；
- 保留输入列表顺序，随后由 metadata 对齐矩阵行。

### `metadata`

脚本会在 `.ablation_prepare_metadata()`（`R/ablation.R:1289`）中识别常见列名别名并统一为：

| 统一列 | 用途 |
| --- | --- |
| `sample_id` | 与表达矩阵列名、d1 行名对齐 |
| `cohort` | 分组折叠、跨 cohort 检索、训练/测试划分 |
| `tissue` | 模块覆盖、组织分层抽样与 tissue-first |
| `biology` | 外部生物学邻域一致性（可选） |
| `cancer_type` | representation 的默认主 anchor；必须作为独立列提供 |

矩阵输入必须显式提供 metadata；只有嵌套列表输入能够从列表层级派生基础 metadata。规范化过程保证 `sample_id/cohort/tissue/biology` 四个核心列：缺少 `tissue` 时才尝试从 CCS 对象派生，缺少 `biology` 时使用 `tissue`。`cancer_type` 不属于自动生成的核心列；representation 默认以它作为主 anchor，并且 endpoint 资格表也依赖该列，因此调用方应同时提供独立的 `tissue` 与 `cancer_type`，避免只有 `cancer_type` 时它被识别为 `tissue` 别名。所有关键矩阵都按 `sample_id` 重新排序，避免“矩阵位置相同但样本含义不同”的隐性错误（`R/ablation.R:1289-1324`、`4927-4938`）。

## 两套参数 schema

脚本保留两套相互隔离的配置：

### Layered schema

`.ablation_default_params()`（`R/ablation.R:708`）按科学问题分组：

| 分组 | 关键字段 | 默认含义 |
| --- | --- | --- |
| `general` | `rank`, `k`, `n_folds`, `bootstrap` | PCA rank、邻域大小、整 cohort 折数、bootstrap 次数 |
| `general$dr` | `method`, `dimension`, `n_neighbors` | UWOT/降维参数；维度默认 `c(5, 2)` |
| `general$cluster` | `eps`, `minPts` | DBSCAN 参数 |
| `cohort` | `rank_sensitivity`, `rp_seeds`, `permutation_seeds` | Experiment 1 的 rank 敏感性与 null 重复 |
| `scaling` | `counts`, `sequences`, `embedding_counts` | d1 模块扩展曲线与下游 embedding 子集 |
| `scaling$gate` | `enforce`, `primary_metric`, `min_gain` | Gate 1 是否阻断 scaling |
| `tissue_first` | `seeds`, `subsample_fraction` | paired 组织分层抽样 |
| `metaccs` | `resample_seeds`, `umap_seeds`, `parameter_mode` | 样本组成与算法随机性的双层重复 |

参数先由 `.ablation_normalize_params()` 兼容旧的 flat 写法，再由 `.ablation_merge_lists()` 递归合并；未知字段、重复指定同一路径和不合法值会在昂贵计算前失败（`R/ablation.R:777-939`）。例如：

```r
params <- list(
  general = list(
    rank = 50L,
    n_folds = Inf,
    dr = list(method = "UWOT", dimension = c(5L, 2L)),
    cluster = list(eps = 0.02, minPts = 20L)
  ),
  cohort = list(
    rp_seeds = 20260727 + 1:20,
    permutation_seeds = 20261727 + 1:20
  ),
  scaling = list(gate = list(enforce = TRUE, min_gain = 0.02))
)
```

### Representation schema

representation 使用独立的 `.ablation_representation_default_params()`（`R/ablation.R:4446`），避免 layered 参数意外影响默认主分析：

| 分组 | 关键字段 | 作用 |
| --- | --- | --- |
| `comparison` | `module_ids`, `direct_group`, `cohort_group` | 选择模块与输出组名 |
| `provenance` | `external_cohorts`, `max_reference_samples`, `require_external` | 定义 reference/query 来源与样本上限 |
| `anchors` | `primary`, `primary_role`, `technical` | 主标签、独立性声明、技术 anchor |
| `geometry` | `k`, `search`, `n_trees`, `search_k` | exact/Annoy 检索与近似搜索验证 |
| `validation` | `learning_fractions`, `repeats`, `lambda`, `inner_folds`, `numCores`, `workers` | grouped readout 与确定性 learning-curve job 调度；总线程预算约为 `numCores × workers` |
| `scaling` | `enabled`, `module_counts`, `biology_anchors` | 可选 cohort bank 扩展 |
| `controls` | `null_rp`, `null_perm`, `null_rp_seeds` | paired null 对照 |
| `tradeoffs` | `decoder`, `decoder_rank`, `decoder_lambda` | 从 d1 解码 Direct feature 的辅助诊断 |
| `output` | `cover` | 是否允许覆盖非空输出目录 |

## Layered 工作流的生命周期

`.ablation_run_layered()` 将运行拆成三个阶段（`R/ablation.R:474-506`）：

```r
context <- .ablation_prepare_layered_context(...)
execution <- .ablation_execute_layered_experiments(context, seed, verbose)
.ablation_finalize_layered_run(context, execution, output.dir, match.call(), verbose)
```

### 共享准备上下文

`.ablation_prepare_layered_context()` 只做一次昂贵的输入准备：检查对象、解析参数、拒绝非空输出目录、调用 `.ablation_prepare_input()`、生成 manifest，并保存 `manifest.rds` 与 `config.rds`（`R/ablation.R:512-566`）。这样四个实验不会各自重新编码表达矩阵，也保证所有结果引用同一份样本/feature hash。

### Experiment 1：cohort representation

入口 `.ablation_experiment_cohort()`（`R/ablation.R:1892`）的设计要点：

1. 以 cohort 为单位做 grouped folds，避免同一 cohort 同时出现在训练和测试；
2. 对每个 rank 复用相同的 permutation 结果，保证 rank 比较是 paired 的；
3. `Direct`、`Cohort`、`Null-RP`、`Null-Perm` 使用同一 fold 和 test sample；
4. 主 rank 用于 Gate 1，其余 rank 只进入 sensitivity；
5. 同时计算 dimension-free geometry（CKA、距离排序相关、kNN Jaccard、effective rank）。

```r
folds <- .ablation_grouped_folds(
  prepared$metadata$cohort,
  n_folds = config$general$n_folds,
  seed = seed
)

permuted <- lapply(config$cohort$permutation_seeds, function(seed_i) {
  .ablation_permute_blocks(
    prepared$d1,
    prepared$module_manifest$blocks,
    prepared$metadata$cohort,
    seed_i
  )
})
```

结果包含 `metrics`（长表）、`summary`、`contrasts`、`rank_sensitivity`、`dimension_free_geometry` 和 `audit`。

### Experiment 2：cohort-axis scaling

`.ablation_experiment_scaling()`（`R/ablation.R:2379`）生成 tissue-balanced 的 nested module sequences。每条 sequence 的小模块集合是大集合的前缀，因此相邻 module count 可以做增量比较：

```r
for (sequence_id in seq_along(sequences)) {
  order <- sequences[[sequence_id]]
  for (module_count in counts) {
    module_ids <- order[seq_len(module_count)]
    candidate <- prepared$d1[, columns, drop = FALSE]
    # full 与 candidate 使用相同 rank_q、fold 和 test 样本
  }
}
```

每个 fold × sequence × module count 记录：

- 与完整 d1 的线性 CKA、kNN Jaccard、effective rank；
- cohort mixing 与 biology purity；
- 可选线性 XGBoost probe 的 macro AUROC 与 balanced accuracy；
- 在代表性 module count 上继续做两阶段降维、DBSCAN 与稳定性比较。

### Gate 1

当请求 `scaling` 时，执行器先运行 Experiment 1，再调用 `.ablation_gate_one()`。Gate 1 默认只记录不阻断（`enforce = FALSE`）；设置为 `TRUE` 时，若主指标的 paired CI、purity 或 mixing 不满足阈值，Experiment 2 会返回 `status = "stopped_by_gate_one"`，而不是伪造一个空的 scaling 结果（`R/ablation.R:599-631`）。

### Experiment 3：tissue-first reduction

`.ablation_experiment_tissue_first()`（`R/ablation.R:2772`）只改变降维顺序：

- **Two-stage**：每个 tissue 先从 d1 降到 d2，再拼接并降到 d3；
- **One-stage**：完整 d1 直接降到 d3。

两条 arm 在每个 seed 使用同一批 tissue × cohort 抽样，随后统一做 DBSCAN、embedding fidelity、cluster size entropy、tissue stratified retention 和稳定性指标。脚本还标记最小样本量四分位的 tissue，检查 tissue-first 是否只对大 tissue 有利。

### Experiment 4：end-to-end metaCCS

`.ablation_experiment_metaccs()`（`R/ablation.R:3508`）把样本组成重复（`resample_seeds`）、降维随机性（`umap_seeds`）和可选参数网格组合起来。每个 run 同时计算 Direct-GSClassifier 与 Cohort-d1：

```r
paired <- .ablation_paired_two_stage_embeddings(
  direct_blocks = direct_i,
  cohort_d1 = cohort_d1,
  dr_config = parameter_manifest$dr[[dr_row]],
  seed = umap_seed
)
```

结果中的 `inference_scope` 会区分：

- `resample_variation`：`subsample_fraction < 1` 且有多个样本组成重复；
- `algorithm_variation_only`：只改变算法 seed，不能解释为样本组成不确定性。

## Representation 主工作流（默认路径）

### 参考集、外部 query 与 d1 provenance

`.ablation_prepare_representation_input()`（`R/ablation.R:4632`）先读取冻结模块与 feature manifest，再将 cohort 分为：

- **reference**：不属于 `filtered.cohort`，且必须已有 d1 行；
- **external query**：属于 filtered cohort，可以没有 d1 行。

reference 与 query 都只消费 `object@Data$Probability$d1` 中已有的行。query 缺少预计算 d1 时，函数发出警告并将该样本排除；如果排除后没有可用 query，则直接报错。保留的 reference/query 分别记录为 `d1_provenance = "in_sample"` 与 `"external_frozen"`。

```r
precomputed_query_ids <- intersect(query_ids, rownames(d1))
excluded_query_d1_ids <- setdiff(query_ids, precomputed_query_ids)
query_ids <- precomputed_query_ids
query_d1 <- d1[query_ids, expected_columns, drop = FALSE]
```

这个边界把 d1 准备责任留给调用方，避免下游分析静默预测新的 d1。若 `anchors$primary_role = "independent"` 且 query provenance 全部为 `external_frozen` 或 `out_of_fold`，最终 evidence level 才会是 `confirmatory`；否则是 `descriptive`（`R/ablation.R:4418-4443`）。

### Phase 1：native geometry

Direct 和 d1 不直接拼接比较，而是在各自 native scale 上转换：

- Direct：按 reference 的均值/标准差做 standardized Euclidean；
- d1：先按列标准化，再按 module block 加权，使每个模块拥有相同的总平方距离权重。

```r
direct_scaled <- .ablation_scale_train_apply(
  prepared$reference_direct,
  prepared$query_direct
)
d1_scaled <- .ablation_module_balanced_transform(
  prepared$reference_d1,
  prepared$query_d1,
  prepared$selected_blocks
)
```

`.ablation_native_geometry()`（`R/ablation.R:4791`）在 reference 上报告 linear CKA、distance Spearman、kNN Jaccard、Direct/d1 effective rank，并检查每个 d1 module block 是否近似 simplex（行和约等于 1）。

### Phase 2：cross-cohort retrieval

`.ablation_query_reference_retrieval()`（`R/ablation.R:4088`）严格执行“query 对 reference”的检索：对每个 query 排除同 cohort reference，再取 `k = c(5, 15, 30)` 邻居。它同时支持 exact 和 Annoy：

- exact：逐 query 计算欧氏距离；
- Annoy：先取候选邻居，再过滤同 cohort；
- Annoy 模式还会在同一 query/reference 任务上与 exact 对照，计算 recall，低于 `min_annoy_recall` 时终止。

主要输出指标：`top1_label_match`、`top_k_label_rate`、MRR，以及可选技术列的 observed match、expected match 和 match excess。Direct 与 Cohort-d1 的 per-sample 结果通过 sample ID 配对，避免不同 query 构成造成偏差。

### Phase 3：Null controls

主检索结果与 null 结果分开保存：

- **Null-RP**：在 Direct 空间上使用稀疏 Achlioptas 随机投影，改变维度但保留随机几何；
- **Null-Perm**：在 cohort 内按整个 module block 置换。若每个 query cohort 的 anchor 都是常量，则没有 label-level power，脚本会返回 `not_eligible` 而不是运行无意义的检验。

所有 null 使用独立、显式的 seed，并沿用主检索的 k、技术列和跨 cohort 约束。

### Phase 4：supervised readout 与 learning curve

若 `validation$enabled = TRUE` 且有可估计的 query cohort，`.ablation_linear_readout()` 会：

1. 以 cohort 为组做 inner CV；
2. 从 `lambda` 候选中选择 balanced accuracy 最好的值；
3. 在固定 lambda 下训练线性 XGBoost readout；
4. 返回 overall、by-cohort、逐样本 predictions 与概率矩阵。

Direct 与 d1 使用同一训练 cohort、同一 query view、同一随机种子；learning curve 再在固定 test set 上改变训练 cohort fraction，并对每个 fraction 做 paired repeats。

### Phase 5：scaling 与 decoder 辅助分析

representation 内的 `scaling$enabled` 与 layered 的 scaling 不同：它固定 Direct feature contract，对冻结 cohort bank 进行 breadth/depth/matched-size 设计，重点关注非冗余、技术稳健性、外部 biology 和邻域稳定性。

`tradeoffs$decoder = TRUE` 时，`.ablation_decode_direct_features()` 在 reference 上拟合 d1 → Direct 的 ridge decoder，再在 query 上按 feature type 使用原生损失：

| feature type | 主要指标 |
| --- | --- |
| `gene_pair` | balanced accuracy、Brier |
| `single_bin` | Spearman、MAE |
| `set_pair` | Spearman、RMSE |

它是解释性/机制诊断，不会替换主 retrieval endpoint。

## 关键内部组件

### Frozen feature manifest

`.ablation_module_manifest()`、`.ablation_frozen_feature_manifest()`（`R/ablation.R:960-1082`）从 CCS 模型中提取：模块 ID、模块所属 tissue、d1 block 列、Direct feature、TSP feature、break vector 和 feature type。所有后续矩阵都依赖这个 manifest，避免从列名猜测模块边界。

### Direct-GSClassifier 特征重建

`.ablation_gsclassifier_matrix()`（`R/ablation.R:1351`）依据 frozen feature manifest 从输入表达矩阵重建 Direct-GSClassifier 特征，并严格保持 manifest 中的列顺序。该矩阵是 representation 路径的主要一次性限速步骤：默认写入 `output.dir/direct-feature-cache.rds`，缓存键同时绑定表达矩阵、样本顺序、冻结模型元数据和 feature manifest；后续运行只有在全部校验一致时才复用，失配或缓存损坏则安全重建并原子替换。representation 路径不会用这些特征补算 d1；d1 始终来自 `object@Data$Probability$d1`。

原生几何中的精确高维 kNN 与有效秩另行写入 `output.dir/native-geometry-cache.rds`。缓存键绑定 Direct/d1 输入、样本顺序、完整 geometry 参数与 seed；旧版 `native_geometry.rds` 不再自动提升为新缓存，因为旧 manifest 不能证明完整输入内容一致。缓存命中只跳过完全相同的精确计算，不会切换为近似近邻算法。

representation 的 retrieval（含 anchor retrieval 与 controls）、readout、learning curve 和 decoder 使用彼此独立的多版本节点 checkpoint。内容键绑定实际输入指纹、节点相关参数、seed、schema、算法修订、递归代码依赖和运行库版本；只修改 learning curve 设计不会牵连 retrieval/readout。缓存值与 state 分离并原子提交，只有 state、key、文件 MD5 与 value hash 全部一致的 `complete` 条目才可命中；`running`、截断或契约不匹配均安全重算。正式 manifest 的 `node_cache` 保留命中/未命中、写入状态和失效原因，恢复缓存仍与面向报告的正式 RDS 分离。

### ablation-03 审核修复后的输入与结果契约

readout 和 learning curve 的 reference/query 均传入原始表示，由训练折内部拟合变换，并应用于该折测试集；最终模型只在完整 reference 上拟合尺度。这里的原始 d1 指 CCS 中已有的概率矩阵，不重新预测 cohort 模型。

`anchor_retrieval.rds` 为所有候选外部样本提供邻居，独立于癌种资格限定的 `retrieval.rds`；`sample-contract.rds` 保存固定 reference/query 身份。连续 anchor 用完整 reference 拟合基因尺度，要求两臂的 top-15 邻居分数完整，再先按 query 配对、后按 cohort 推断。IFN/IL6 signature 由名称解析，并校验 signature、配置和表达来源内容。结构状态若在分位边界并列则不强行切分。

bank scaling 使用经审计的 tissue 标签，同时保留原始模块 ID；缓存绑定 reference/query 内容、metadata、feature manifest 和验证参数。manifest 的 `input_key` 独立记录完整输入，数值未改变的原生几何仍可使用自身严格缓存键。主实验、生物学、结构三个阶段各写入 `stage-receipt.rds`，下游和 Rmd 渲染核验输入/产物内容指纹，拒绝混合批次。ablation-03 的 biology 阶段仅使用 R 实现。

### ablation-03 的准备、分析与报告分层

ablation-03 的运行入口已按业务分层，见 [ablation-03 运行说明](../test/ablation-03/README.md)。`01a/01b/01c` 顺序准备审计数据、Direct/d1 表示与表达 anchor 缓存；完整样本契约在 `01b` 生成，`01c` 不再依赖主实验结果。结构表达提取前移，尺度拟合和统计比较仍留在分析阶段。

`02-ablation03-representation.R`、`03-ablation03-biology.R` 与 `04-ablation03-structural-reproducibility.R` 都是同一个实验的模块入口，读取准备好的输入，不隐式运行准备脚本。每个模块有同名 Rmd，渲染为同名 HTML；数据准备另有 `01-ablation03-data-overview.Rmd`。报告保留原有图表与摘要计算，表示报告不再等待后两个模块。

包内 `.ablation_run_representation()` 保留公共 `ablation()` 的原始输入准备路径，并委托 `.ablation_run_prepared_representation()` 执行抽出的原计算体；targets 入口通过上面的导出适配器消费同一计算体。提取不调整科学计算、默认参数或随机种子。准备缓存的来源凭据及其上游依赖会递归核验，旧凭据不能直接作为新流程的完成证明。已有 Direct、精确几何和 scaling 缓存仍按原内容键判断复用，本次迁移只完成代码与串行依赖骨架，未运行完整数据准备、实验或报告渲染。

### 抽样与可重复性

`.ablation_stratified_sample()` 对 tissue × cohort 分层；每个实验拥有独立 seed 区间：

```text
base seed                 主实验
base + 10,000             layered scaling
base + 20,000             tissue-first
base + 30,000             metaCCS
base + 40,000 / 50,000    representation decoder / readout 辅助阶段
```

脚本还保存 module sequence hash、sample hash、feature manifest hash 和 config hash，使“同一实验是否真的使用同一输入”可以复核。

## 返回对象与落盘文件

所有路径最终返回 class 为 `CCSAblation` 的 list。layered 结果由 `.ablation_finalize_layered_run()` 组装（`R/ablation.R:677-703`）：

```r
result <- structure(
  list(
    call = call,
    manifest = context$manifest,
    config = context$config,
    experiments = execution$experiments,
    audit = audit,
    output.dir = normalizePath(output.dir, winslash = "/", mustWork = TRUE)
  ),
  class = "CCSAblation"
)
```

### Representation 输出

| 文件 | 内容 |
| --- | --- |
| `manifest.rds` | 样本、cohort、feature、anchor、evidence level 与 config hash |
| `native_geometry.rds` | Direct/d1 几何和 d1 block 诊断 |
| `retrieval.rds` | neighbors、per-sample、summary、paired 与 search validation |
| `readout.rds` | supervised readout、predictions、by-cohort 指标 |
| `learning_curve.rds` | 不同训练 cohort fraction 的 paired 曲线 |
| `cohort_scaling.rds` | 可选 frozen bank scaling |
| `tradeoffs.rds` | feature type 统计、simplex 检查、decoder |
| `endpoint_eligibility.rds/csv` | candidate 与 estimable endpoint 资格表 |
| `excluded-query-d1.csv` | 因缺少预计算 d1 而被排除的 query 样本及原因 |
| `audit.csv` | 面向审阅者的紧凑摘要 |

### Layered 输出

```text
manifest.rds
config.rds
experiment-01-cohort.rds
experiment-02-scaling.rds
experiment-03-tissue-first.rds
experiment-04-metaccs.rds
audit.csv
ablation-result.rds
```

只请求某个实验时，不会创建未请求分支的伪结果；scaling 因 Gate 1 停止时会明确记录 `status = "stopped_by_gate_one"`。

## 最小使用示例

### 默认 representation

```r
result <- ablation(
  object = ccs_fit,
  data = query_data,
  metadata = query_metadata,
  experiment = "representation",
  output.dir = "results/ablation-representation",
  params = list(
    geometry = list(search = "exact", k = c(5L, 15L, 30L)),
    validation = list(enabled = TRUE),
    controls = list(null_rp = TRUE, null_perm = TRUE)
  ),
  seed = 20260727
)

result$evidence_level
result$retrieval$summary
result$readout$overall
```

representation 参数位于 `params$geometry`、`params$validation` 等嵌套分组中；不要把 layered 的 `general$...` 参数混进这条路径。

### Layered 实验

```r
result <- ablation(
  object = ccs_fit,
  data = nested_data,
  metadata = metadata,
  experiment = c("scaling", "tissue_first", "metaccs"),
  output.dir = "results/ablation-layered",
  params = list(
    general = list(rank = 50L, bootstrap = 1000L),
    scaling = list(gate = list(enforce = TRUE, min_gain = 0.02)),
    metaccs = list(parameter_mode = "fixed")
  )
)

result$experiments$cohort$contrasts
result$experiments$scaling$summary
result$experiments$tissue_first$contrasts
result$experiments$metaccs$contrasts
```

请求 `scaling` 时会自动先运行 layered Experiment 1，因此结果仍包含 `result$experiments$cohort`，无需在实验向量中额外写入已弃用的单值入口 `"cohort"`。

### 内嵌 smoke fixture

脚本顶部保留一个默认关闭的 deterministic fixture。仅在显式设置环境变量时运行：

```powershell
$env:CCS_ABLATION_RUN_SMOKE = "true"
Rscript -e "library(CCS); ablation(object = NULL, data = NULL)"
```

上面的调用适用于已安装包。开发 CCS 包本身时可以使用 `pkgload::load_all()`
做源码测试，但这不属于 ablation-03 的正式运行方式。

它会构造一个小型 CCS 对象，直接调用 layered orchestrator，并断言 `CCSAblation`、四个 group 和 audit 表均存在（`R/ablation.R:275-412`）。这不是生产数据测试，但适合快速检查入口、参数、输出结构是否破坏。

## 解读结果时的边界

1. **先看 `endpoint_eligibility` 和 `evidence_level`。** candidate query 不等于可用于 cancer-labelled endpoint 的 estimable query。
2. **把 primary、diagnostic、mechanism 分开。** cancer readout 和 lineage readout 被标为 diagnostic 时，不能替代外部 biology 的主证据。
3. **不要把算法重复当作样本重复。** metaCCS 的 `umap_seeds` 只反映算法变化；只有多个 resample 才能支持 `resample_variation` 推断。
4. **比较 Direct 与 d1 时使用 paired 字段。** `delta_*` 是同一 query、同一 fold 或同一 resample 内的差值。
5. **Annoy 结果必须检查 search validation。** 若 recall 低于阈值，脚本会要求增加 `n_trees/search_k` 或改用 exact。
6. **输出目录默认不可覆盖。** 只有显式设置 `params$output$cover = TRUE`（representation）或 `params$general$cover = TRUE`（layered）才允许写入非空目录。

## 维护与扩展建议

- 新增指标：优先在对应 metric helper 中返回命名向量，再由长表汇总器统一 bootstrap；不要在 orchestrator 中拼接散乱列。
- 新增实验：实现 `.ablation_experiment_<name>()`，在 layered dispatch 中声明依赖关系，并为 audit 补齐 sample/feature/config hash。
- 改变输入契约：同步修改 flatten、metadata alias、manifest 和 embedded smoke fixture。
- 改变输出 schema：更新本文档的 result/file 表，并保留 `schema_version` 或 status 字段，避免下游静默误读。
- 任何随机步骤都应从显式 seed 派生；不要在循环内依赖不可见的全局 RNG 状态。

换句话说，`R/ablation.R` 的核心不是“运行很多指标”，而是把冻结表示、独立 query、paired 对照、资格判定和审计证据绑定在同一条可复现流水线上。

## d1 输入边界（当前契约）

表示层 `ablation()` 只消费 `object@Data$Probability$d1` 中已有的样本行。
query 缺少 d1 时，函数发出警告并取已有样本交集，不在下游分析中读取冻结模型
补算 d1；排除的样本 ID 与原因写入 `excluded-query-d1.csv`。
