# `R/ablation.R` 业务逻辑与函数地图

`R/ablation.R` 是 CCS 的冻结表示消融分析入口。它不重新训练 CCS 的 cohort 模型，而是从已经拟合好的 `CCS` S4 对象中构造两类可比较表示：

- `Direct-GSClassifier`：依据冻结 GSClassifier 模型的完整 feature contract，从 RNA 表达重建原生输入；
- `Cohort-d1`：直接复用 `object@Data$Probability$d1` 中已有的模块概率。

它用配对、跨 cohort、可审计的方式回答：

- d1 是否保留了 Direct 的几何结构、生物学邻域和技术稳健性；
- cohort 模块数量、tissue-first 降维和完整 metaCCS 流程如何影响结果；
- 观察到的收益是否超过随机投影、模块内置换或单阶段降维对照。

因此，这不是一个“删除变量后重跑模型”的小函数，而是一条从输入对齐、冻结表示重建、资格判定、配对评估、缓存到审计输出的完整流水线。

## 一张图看懂入口

~~~mermaid
flowchart TD
    A[CCS object] --> B[ablation()]
    X[expression matrix 或 tissue/cohort 嵌套 data] --> B
    M[metadata] --> B
    B --> D{experiment}
    D -->|representation 默认| R[Representation orchestrator]
    D -->|cohort 单值旧别名| R
    D -->|cohort/scaling/tissue_first/metaccs| L[Layered orchestrator]
    R --> R1[冻结 module/feature manifest]
    R1 --> R2[reference/query 划分]
    R2 --> R3[Direct 重建 + d1 复用]
    R3 --> R4[native geometry]
    R3 --> R5[cross-cohort retrieval]
    R3 --> R6[Null-RP / Null-Perm]
    R3 --> R7[readout + learning curve]
    R3 --> R8[可选 bank scaling + decoder]
    L --> C[共享 frozen context]
    C --> E1[Experiment 1 cohort]
    E1 --> G{Gate 1}
    G -->|通过或不强制| E2[Experiment 2 scaling]
    G -->|强制且失败| STOP[scaling stopped_by_gate_one]
    C --> E3[Experiment 3 tissue-first]
    C --> E4[Experiment 4 metaCCS]
    R --> O[CCSAblation + RDS/CSV]
    L --> O
~~~

公共入口只负责分派；真正的业务逻辑在后面的 `.ablation_*` 函数组完成。

## 公共入口与生命周期

### `ablation()` 的实验分派

~~~r
ablation(
  object = ccs_fit,
  data = expression_data,
  metadata = sample_metadata,
  experiment = "representation",
  output.dir = "results/ablation",
  params = list(),
  seed = 20260727,
  step = "all"
)
~~~

`experiment` 支持：

- `"representation"`：当前默认的独立表示比较流程；
- `"cohort"`：单独传入时是 deprecated alias，发出 warning 后映射到 `"representation"`；放在多元素实验向量中时仍表示 layered Experiment 1；
- `"scaling"`、`"tissue_first"`、`"metaccs"`：layered 实验，可传一个或多个；
- `representation` 不能和 layered 实验混用，未知名称在计算前直接报错。

`step` 支持以下阶段：

| step | 输入 | 动作 | 返回类 |
| --- | --- | --- | --- |
| `all` | object/data/metadata | 完整运行并写正式产品 | `CCSAblation` |
| `context` | object/data/metadata，或含 object/data 的 input | 准备输入、manifest、配置和缓存布局 | `CCSAblationContext` |
| `plan` | context | 生成节点与 job 描述，不执行计算 | `CCSAblationPlan` |
| `run` | plan；也接受 context 并自动建 plan | 执行 representation runner | `CCSAblationRun` |
| `result` | run | 复制 runner 产品并写 reviewer-facing 结果 | `CCSAblation` |

staged execution 当前只支持 `experiment = "representation"`。阶段对象通过 `input` 传递；非 `all` 阶段不会把中间对象伪装成正式结果。

`cache.root` 只指定中间缓存根目录，默认是项目下被忽略的 `.ccs-cache/ablation/`；`output.dir` 保存正式产品，两者不能相同。非空输出目录默认拒绝写入：representation 使用 `params$output$cover = TRUE`，layered 使用 `params$general$cover = TRUE`。

## 输入怎样变成可比较表示

| 函数 | 功能 | 被谁消费 |
| --- | --- | --- |
| `.ablation_module_manifest()` | 从 d1 列名和模型恢复 tissue、module、d1 block 边界 | 所有模块级实验 |
| `.ablation_frozen_feature_manifest()` | 提取 Direct feature、TSP、break vector、feature type | Direct、scaling、decoder |
| `.ablation_flatten_expression()` | 统一矩阵或 `tissue -> cohort -> leaf` 嵌套列表 | 两条主流程 |
| `.ablation_prepare_metadata()` | 识别别名并统一 sample/cohort/tissue/biology | 所有评估 |
| `.ablation_gsclassifier_matrix()` | 按冻结 feature contract 重建 Direct 输入，不拟合新模型 | Direct 表示 |
| `.ablation_prepare_input()` | layered 的一次性全体输入准备 | 四个 layered 实验 |
| `.ablation_prepare_representation_input()` | representation 的 reference/query 划分、d1 过滤和 cache key | 默认流程 |
| `.ablation_build_manifest()` | 保存样本、feature、模块、版本和 config hash | audit 和复核 |

`data` 可以是表达矩阵/data.frame，也可以是：

~~~r
data <- list(
  Lung = list(
    TCGA = list(expr = lung_tcga, subtype = lung_label),
    GEO  = list(expr = lung_geo,  subtype = lung_label_geo)
  ),
  Breast = list(
    TCGA = list(expr = breast_tcga, subtype = breast_label)
  )
)
~~~

嵌套列表的叶节点可以是矩阵，或 `list(expr = ..., subtype = ...)`。脚本按基因并集对齐，缺失基因填 `NA`；重复 sample ID 整组排除并记录。矩阵输入没有列表层级，因此必须另传 `metadata`。

metadata 至少需要 `sample_id` 和 `cohort`。常见别名会被识别；`tissue` 缺失时才尝试从 CCS 派生，`biology` 缺失时退回 `tissue`。representation 的默认主 anchor 是 `cancer_type`，它不是自动生成的核心列，调用方应同时提供独立的 `tissue` 和 `cancer_type`。

### 两条准备边界

layered 通过 `.ablation_prepare_input()` 只保留同时出现在 d1、表达矩阵和 metadata 中的样本。

representation 通过 `.ablation_prepare_representation_input()` 将样本分成：

- reference：不属于 `filtered.cohort`，且必须已有 d1 行；
- external query：属于 `filtered.cohort`，先作为候选 query；
- query 若缺少预计算 d1，会警告、写入 `excluded_query_d1_ids` 并排除；不会在下游静默补算 d1。

保留样本标记为 `d1_provenance = "in_sample"` 或 `"external_frozen"`。只有主 anchor 为 independent 且两侧 provenance 都满足 external/out-of-fold 时，`evidence_level` 才是 `confirmatory`；否则为 `descriptive`。

## 两套参数 schema

representation 和 layered 使用两套隔离配置。不要把 layered 的 `general$...` 混入 representation，也不要把 representation 的 `geometry$...` 混入 layered。

### representation schema

默认值来自 `.ablation_representation_default_params()`：

| 组 | 关键字段 | 业务含义 |
| --- | --- | --- |
| `comparison` | module_ids, direct_group, cohort_group | 选择模块和组名 |
| `provenance` | external_cohorts, max_reference_samples, max_query_samples, require_external | reference/query 来源和样本上限 |
| `anchors` | primary, primary_role, bank_aligned, technical, min_reference_cohorts, endpoint_min_reference_cohorts | 主标签、独立性、技术列和 endpoint 资格 |
| `geometry` | k, search, n_trees, search_k, exact_validation_queries, min_annoy_recall, geometry_samples, distance_pairs | native geometry 与 exact/Annoy retrieval |
| `validation` | enabled, learning_fractions, repeats, inner_folds, lambda, nrounds, min_class_n, numCores, workers | readout 和 learning curve |
| `scaling` | enabled, module_counts, sequences, direct_feature_type, sensitivity_feature_type, biology_anchors, score_reference_samples, score_query_samples, lambda, bootstrap | representation 内的 frozen bank scaling |
| `controls` | null_rp, null_rp_rank, null_rp_seeds, null_perm | paired null |
| `tradeoffs` | decoder, decoder_rank, decoder_lambda, decoder_max_reference_samples, decoder_max_query_samples | d1→Direct decoder |
| `output` | cover, cache_direct | 覆盖输出目录和 Direct cache |

### layered schema

默认值来自 `.ablation_default_params()`：

| 组 | 关键字段 | 业务含义 |
| --- | --- | --- |
| `general` | rank, k, distance, n_folds, bootstrap, max_samples, probe, probe_label, probe_nrounds, numCores, fidelity_samples, cover | 共享验证和输出设置 |
| `general$dr` | method, dimension, n_neighbors, min_dist, spread, set_op_mix_ratio, metric, n_threads | 降维后端 |
| `general$cluster` | eps, minPts | DBSCAN |
| `cohort` | rank_sensitivity, geometry_samples, distance_pairs, mechanism_samples, rp_density, rp_seeds, permutation_seeds | Experiment 1 |
| `scaling` | counts, sequences, embedding_counts, embedding_sequences, embedding_seeds, subsample_fraction, gate | Experiment 2 和 Gate 1 |
| `tissue_first` | seeds, subsample_fraction | Experiment 3 |
| `metaccs` | resample_seeds, umap_seeds, subsample_fraction, parameter_mode, dr_grid, cluster_grid, direct_feature_mode, retain_assignments | Experiment 4 |

`.ablation_normalize_params()` 兼容旧 flat 名称；`.ablation_merge_lists()` 递归合并；`.ablation_validate_config()` 在昂贵计算前拒绝未知字段、非法 seed、不可行维度和错误的 metaCCS grid。

## Representation：默认主流程

执行链是 `.ablation_run_representation()` → `.ablation_prepare_representation_analysis()` → `.ablation_run_prepared_representation()`。

~~~mermaid
flowchart LR
    P[prepared reference/query] --> S[Direct 标准化]
    P --> B[d1 module-balanced transform]
    S --> G[native geometry]
    B --> G
    S --> Q[query-reference retrieval]
    B --> Q
    Q --> N[Null-RP / Null-Perm]
    S --> L[linear readout]
    B --> L
    L --> LC[paired learning curve]
    B --> SC[可选 cohort bank scaling]
    B --> D[d1 -> Direct decoder]
    G --> F[manifest + audit]
    Q --> F
    LC --> F
    SC --> F
    D --> F
~~~

### 表示变换、资格和 retrieval

`.ablation_scale_train_apply()` 只在 reference 上估计 Direct 的中心和尺度，再应用到 query；`.ablation_module_balanced_transform()` 在 reference 边界标准化 d1 block 并重新加权。

`.ablation_native_geometry()` 报告 linear CKA、距离排序 Spearman、kNN Jaccard、Direct/d1 effective rank，以及 d1 block 的概率范围、行和和 simplex 诊断。

`.ablation_endpoint_eligibility()` 按 query cohort 统计 reference 对各 cancer type 的支持，为 `cancer_retrieval`、`technical_excess`、`cancer_readout`、`learning_curve` 标记 `estimable`/`not_estimable`；`.ablation_endpoint_view()` 返回各 endpoint 的 query 子集。

`.ablation_query_reference_retrieval()` 对每个 query 只在不同 cohort 的 reference 中找邻居，支持 exact 和 Annoy。Annoy 由 `.ablation_validate_neighbor_search()` 与 exact 对照计算 recall，低于 `min_annoy_recall` 终止。指标包括 `top1_label_match`、`top_k_label_rate`、MRR，以及技术列的 observed/expected match 和 match excess。`.ablation_bind_retrieval()` 按 sample ID 配对 Direct/d1 并生成 `delta_*`。

### Null、readout 和 learning curve

`.ablation_projection_matrix()`/`.ablation_random_projection()` 实现稀疏 Achlioptas Null-RP；`.ablation_permute_blocks()` 在 cohort 内按完整 module block 置换；`.ablation_null_perm_eligibility()` 防止对恒定 anchor 运行无解释力的 Null-Perm。

`.ablation_linear_readout()` 按 cohort 做 grouped inner CV，选择 lambda 后用 xgboost `gblinear` 拟合 readout。`.ablation_sample_training_cohorts()`、`.ablation_learning_curve_job()`、`.ablation_learning_curve_execute_job()`、`.ablation_learning_curve_worker()` 和 `.ablation_learning_curve()` 共同完成 fraction × repeat × representation 的串行/PSOCK 调度和 job checkpoint。

### representation scaling 和 decoder

`params$scaling$enabled = TRUE` 时，`.ablation_representation_scaling()` 构造 tissue-balanced nested cohort bank，比较 breadth、within-tissue depth 和 matched-size，并报告非冗余、技术匹配、外部 biology 邻域一致性和稳定性；癌种 readout 仅为 diagnostic lineage。

`params$tradeoffs$decoder = TRUE` 时，`.ablation_decode_direct_features()` 拟合 d1→Direct 的 ridge decoder：gene_pair 用 balanced accuracy/Brier，single_bin 用 Spearman/MAE，set_pair 用 Spearman/RMSE。decoder 是机制诊断，不替换 retrieval endpoint。

## Layered：共享 context 和四个实验

~~~mermaid
sequenceDiagram
    participant U as ablation()
    participant C as prepare_layered_context
    participant E as execute_layered_experiments
    participant F as finalize_layered_run
    U->>C: CCS + data + metadata + params
    C->>C: resolve config / prepare input / build manifest
    C-->>E: immutable context
    E->>E: Experiment 1 cohort
    E->>E: Gate 1
    E->>E: Experiment 2 scaling
    E->>E: Experiment 3 tissue-first
    E->>E: Experiment 4 metaCCS
    E-->>F: experiments + audit_parts
    F-->>U: CCSAblation + formal files
~~~

`.ablation_prepare_layered_context()` 只准备一次 prepared、manifest、config 和 cache layout；`.ablation_execute_layered_experiments()` 决定分支和依赖；`.ablation_finalize_layered_run()` 合并 audit 并写正式产品。

### Experiment 1：cohort representation

入口：`.ablation_experiment_cohort()`。

1. `.ablation_grouped_folds()` 以完整 cohort 为验证单位，避免 leakage；
2. `.ablation_cohort_rank_metrics()` 在每个 fold 比较 Direct、Cohort、Null-RP、Null-Perm；
3. 四组共享 fold、test sample 和 rank；
4. `.ablation_metric_rows()` 收集 CKA、距离 Spearman、kNN Jaccard、effective rank、cohort mixing、biology purity、probe 和 selective reconstruction mechanism；
5. `.ablation_metric_summary()` 对 fold bootstrap，`.ablation_paired_contrasts()` 计算 paired CI；
6. `.ablation_dimension_free_geometry()` 在 native 空间补充几何比较。

### Gate 1

`.ablation_gate_one()` 检查 Cohort 相对 Direct、Null-RP、Null-Perm 的主指标 CI 是否达到 `min_gain`，并检查 purity/mixing 的下降是否在容差内；`balanced_accuracy` 不可估计时回退到 `biology_purity`。请求 `scaling` 时，`gate$enforce = FALSE` 只记录不阻断；`TRUE` 且失败时 Experiment 2 返回 `status = "stopped_by_gate_one"`。

### Experiment 2：cohort-axis scaling

入口：`.ablation_experiment_scaling()`。`.ablation_nested_module_sequences()` 生成 tissue-balanced、可嵌套模块顺序；同一 sequence 的前 `m` 个 module 构成 `m`-module 子集，支持 paired marginal gain。

每个 sequence × module count × grouped fold 计算与完整 d1 的 CKA/kNN Jaccard/effective rank、cohort mixing、biology purity、可选 tissue probe，以及代表性 module count 的两阶段 embedding、DBSCAN 和 seed stability。`.ablation_scaling_summary()` 汇总曲线、相邻增量和 saturation fit；`.ablation_scaling_embedding()`/`.ablation_scaling_embedding_stability()` 检查 d3/cluster 稳定性。

### Experiment 3：tissue-first reduction

入口：`.ablation_experiment_tissue_first()`。两条 arm 使用同一 tissue × cohort 子样本和 seed：

- Two-stage：`.ablation_reduce_by_reference()` 对每个 tissue 得到 d2，再由 `.ablation_two_stage_embedding()` 合并为 d3；
- One-stage：`.ablation_one_stage_embedding()` 将完整 d1 直接降到 d3。

`.ablation_embedding_metrics()` 计算 trustworthiness、continuity、tissue kNN retention、cohort mixing、biology purity、cluster count、cluster size entropy 和 noise fraction。`.ablation_tissue_stratified_metrics()` 防止大 tissue 掩盖小 tissue；`.ablation_embedding_stability()` 比较不同 seed 的 neighborhood Jaccard、ARI 和 cluster Jaccard。

### Experiment 4：end-to-end metaCCS

入口：`.ablation_experiment_metaccs()`。`.ablation_metaccs_parameter_manifest()` 展开 fixed 或 DR × DBSCAN grid；`.ablation_tissue_subsample()` 生成 resample；`.ablation_paired_two_stage_embeddings()` 用共同 block、可行维度、邻居预算和 seed 生成 Direct/d1 两条 arm；`.ablation_dbscan()` 聚类（0 为 noise）；`.ablation_cluster_biology()` 计算 cluster-biology ARI、NMI、weighted purity 和 non-noise coverage；三个 `.ablation_metaccs_*` 汇总函数处理稳定性和 paired 差异。

`umap_seeds` 只反映算法随机性；只有多个 `resample_seeds` 且 `subsample_fraction < 1` 时，`inference_scope` 才是 `resample_variation`。

## `.ablation_*` 子函数职责地图

下表按业务职责组织，而不是按源码出现顺序罗列。

| 层 | 关键函数 | 业务职责 | 连接 |
| --- | --- | --- | --- |
| 入口 | ablation, .ablation_dispatch_step | 实验/stage 分派 | 两个 orchestrator |
| 配置 | .ablation_default_params, .ablation_representation_default_params, .ablation_normalize_params, .ablation_merge_lists, .ablation_validate_config, .ablation_resolve_*_config | 默认值、flat 兼容、合并、校验 | 所有分支 |
| manifest | .ablation_module_manifest, .ablation_frozen_feature_manifest, .ablation_model_features, .ablation_model_break_vectors, .ablation_feature_type | 冻结 module/feature contract | 输入、scaling、decoder |
| 输入 | .ablation_flatten_expression, .ablation_prepare_metadata, .ablation_prepare_input | 表达合并、metadata 和样本对齐 | layered context |
| representation 输入 | .ablation_prepare_representation_input, .ablation_endpoint_eligibility, .ablation_endpoint_view | reference/query、provenance、资格 | representation runner |
| Direct 重建 | .ablation_gsclassifier_matrix, .ablation_direct_expression_fingerprint, .ablation_direct_feature_cache_key | 重建和缓存 Direct | 两条主流程 |
| 变换/几何 | .ablation_scale_train_apply, .ablation_module_balanced_transform, .ablation_fit_pca, .ablation_linear_cka, .ablation_distance_spearman, .ablation_knn, .ablation_knn_jaccard, .ablation_effective_rank | 防泄漏尺度、PCA 和结构指标 | geometry、cohort、readout |
| retrieval | .ablation_query_reference_retrieval, .ablation_validate_neighbor_search, .ablation_bind_retrieval | query→reference 检索和配对差异 | representation |
| null | .ablation_projection_matrix, .ablation_random_projection, .ablation_permute_blocks, .ablation_null_perm_eligibility | 随机/置换控制 | representation、Experiment 1 |
| 监督评估 | .ablation_xgb_linear_predict, .ablation_linear_readout, .ablation_probe, .ablation_classification_metrics, .ablation_binary_auc | readout、probe、AUROC/BA | readout、curve、scaling |
| learning curve | .ablation_sample_training_cohorts, .ablation_learning_curve_job, .ablation_learning_curve_execute_job, .ablation_learning_curve_worker, .ablation_learning_curve | fraction/repeat 调度和 job checkpoint | representation |
| layered orchestration | .ablation_prepare_layered_context, .ablation_execute_layered_experiments, .ablation_finalize_layered_run | 共享准备、Gate 依赖、统一输出 | 四个 layered 实验 |
| scaling | .ablation_nested_module_sequences, .ablation_experiment_scaling, .ablation_scaling_summary, .ablation_scaling_embedding | module bank 扩展和下游稳定性 | Experiment 2 |
| tissue-first | .ablation_tissue_embeddings, .ablation_two_stage_embedding, .ablation_one_stage_embedding, .ablation_embedding_metrics | 两种降维顺序比较 | Experiment 3 |
| metaCCS | .ablation_metaccs_parameter_manifest, .ablation_parameter_sets, .ablation_paired_two_stage_embeddings, .ablation_experiment_metaccs | Direct/d1 端到端降维聚类 | Experiment 4 |
| 统计 | .ablation_metric_summary, .ablation_paired_contrasts, .ablation_two_group_contrasts, .ablation_bootstrap_mean, .ablation_rbind | 长表、paired CI、audit 合并 | 所有实验 |
| cache/audit | .ablation_resolve_cache_layout, .ablation_node_cache_key, .ablation_cached_node, .ablation_atomic_save_rds, .ablation_atomic_write_csv, .ablation_build_manifest | 内容寻址缓存、原子写入、审计 | 所有正式结果 |

## 缓存、随机性和审计

`.ablation_resolve_cache_layout()` 将缓存分为 `context/`、`plan/`、`preparation/`、`nodes/`、`jobs/`、`runner/`、`state/`。`.ablation_node_cache_key()` 绑定实际输入、样本/feature hash、节点参数、seed、schema、algorithm revision、递归代码依赖和运行库版本。

`.ablation_cached_node()` 只有在 state、key、文件 MD5 和 value hash 均一致且为 `complete` 时才命中；中断、`running`、损坏或契约不匹配都会安全重算。Direct、native geometry、retrieval、readout、learning curve job、scaling fit 和 decoder 各有独立缓存边界。

每个随机步骤从显式 seed 派生。layered 默认在 base、`+10000`、`+20000`、`+30000` 附近分别分配主实验、scaling、tissue-first、metaCCS；representation 的 readout、learning curve、scaling、decoder 另有偏移。sample hash、feature hash、config hash、module sequence hash 和 input key 用于复核配对是否真实成立。

## 返回对象和正式文件

两条主流程最终返回 class 为 `CCSAblation` 的 list。

### representation 输出

| 文件 | 内容 |
| --- | --- |
| manifest.rds | 样本、cohort、feature、endpoint 资格、provenance、evidence level、cache/node 状态 |
| native_geometry.rds | Direct/d1 几何和 d1 block 诊断 |
| retrieval.rds | neighbors、per-sample、summary、paired、search validation |
| anchor_retrieval.rds | 不受 cancer endpoint 资格裁剪的完整候选 query retrieval |
| sample-contract.rds | reference/query 身份 |
| readout.rds | supervised readout、predictions、by-cohort |
| learning_curve.rds | 训练 cohort fraction 的 paired 曲线 |
| cohort_scaling.rds | 可选 frozen bank scaling |
| tradeoffs.rds | feature type、simplex、decoder |
| endpoint_eligibility.rds/csv | candidate/estimable/not_estimable 资格表 |
| excluded-query-d1.csv | 缺少预计算 d1 而排除的 query 及原因 |
| audit.csv | retrieval summary 加 evidence 和输入计数 |

### layered 输出

~~~text
manifest.rds
config.rds
experiment-01-cohort.rds
experiment-02-scaling.rds
experiment-03-tissue-first.rds
experiment-04-metaccs.rds
audit.csv
ablation-result.rds
~~~

只请求的分支才执行。请求 `scaling` 会自动先执行 layered Experiment 1；Gate 1 强制失败时，Experiment 2 明确保存 `status = "stopped_by_gate_one"`。

## 最小使用示例

### 默认 representation

~~~r
result <- ablation(
  object = ccs_fit,
  data = expression_data,
  metadata = metadata,
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
~~~

### Layered 实验

~~~r
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
~~~

### staged 调用

~~~r
context <- ablation(
  object = ccs_fit, data = nested_data, metadata = metadata,
  experiment = "representation", step = "context",
  cache.root = ".ccs-cache/ablation"
)
plan <- ablation(step = "plan", input = context)
run <- ablation(step = "run", input = plan)
result <- ablation(
  step = "result", input = run,
  output.dir = "results/ablation-representation"
)
~~~

### 内嵌 smoke fixture

脚本顶部保留默认关闭的 deterministic fixture。显式设置环境变量后，它会构造小型 CCS 对象，调用 layered cohort orchestrator，并断言 `CCSAblation`、四组比较和 audit 表存在：

~~~powershell
$env:CCS_ABLATION_RUN_SMOKE = "true"
Rscript -e "library(CCS); ablation(object = NULL, data = NULL)"
~~~

## 结果解读边界

1. 先看 `endpoint_eligibility` 和 `evidence_level`；candidate query 不等于可用于 cancer-labelled endpoint 的 estimable query。
2. 把 primary、diagnostic 和 mechanism 分开；lineage readout 或 decoder 不能替代外部 biology 的主证据。
3. `umap_seeds` 是算法重复，不等于样本组成重复；只有多个 resample 才支持 `resample_variation`。
4. Direct/d1 的 `delta_*` 必须在同一 query、fold、sequence 或 resample 内解释。
5. Annoy 必须检查 search validation；recall 低于阈值时增加 `n_trees/search_k` 或改用 exact。
6. d1 缺失 query 会被记录和排除；先检查 `excluded-query-d1.csv`。
7. 输出目录默认不可覆盖；确认目录和 `cover` 设置后再运行。

## 维护原则

- 新增指标：在对应 metric helper 返回命名向量或长表，再由统一汇总器 bootstrap；
- 新增实验：实现 `.ablation_experiment_<name>()`，在 `.ablation_execute_layered_experiments()` 中声明依赖和 audit；
- 改变输入契约：同步更新 flatten、metadata alias、manifest、资格表和 smoke fixture；
- 改变输出 schema：同步更新本文档的结果对象与文件表，并保留 `status`/`schema_version`；
- 改变缓存边界：同步更新 cache key、state 校验和 manifest 的 node 状态；
- 不把 `.ablation_*` 函数当作独立脚本；新逻辑应挂在“入口 → 准备 → 实验 → 汇总 → 落盘”层次上。
