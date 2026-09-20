# `R/ablation.R` 业务逻辑与函数地图

`R/ablation.R` 是 CCS 的冻结表示消融分析入口。它不重新训练 CCS 的 cohort 模型，而是从已经拟合好的 `CCS` S4 对象中构造两类可比较表示：

- `Direct-GSClassifier`：依据冻结 GSClassifier 模型的完整 feature contract，从 RNA 表达重建原生输入；
- `Cohort-d1`：直接复用 `object@Data$Probability$d1` 中已有的模块概率。

它用配对、跨 cohort、可审计的方式回答：

- d1 是否保留了 Direct 的几何结构、生物学邻域和技术稳健性；
- 观察到的收益是否超过随机投影、模块内置换等配对 null 对照；
- 扩大 cohort 模块库（frozen bank scaling）能否继续改善检索与 readout。

因此，这不是一个“删除变量后重跑模型”的小函数，而是一条从输入对齐、冻结表示重建、资格判定、配对评估、缓存到审计输出的完整流水线。

## 一张图看懂入口

~~~mermaid
flowchart TD
    A[CCS object] --> B[ablation]
    X[expression matrix 或 tissue/cohort 嵌套 data] --> B
    M[metadata] --> B
    B --> S{step}
    S -->|all 默认| R[Representation orchestrator]
    S -->|context/plan/run/result| L[staged lifecycle 复用同一流程]
    R --> R1[冻结 module/feature manifest]
    R1 --> R2[reference/query 划分]
    R2 --> R3[Direct 重建 + d1 复用]
    R3 --> R4[native geometry]
    R3 --> R5[cross-cohort retrieval]
    R3 --> R6[Null-RP / Null-Perm]
    R3 --> R7[readout + learning curve]
    R3 --> R8[可选 bank scaling + decoder]
    R4 --> O[CCSAblation + RDS/CSV]
    R5 --> O
    R7 --> O
    R8 --> O
    L --> O
~~~

公共入口只负责分派；真正的业务逻辑在后面的 `.ablation_*` 函数组完成。

## 公共入口与生命周期

### `ablation()` 的调用形式

~~~r
ablation(
  object = ccs_fit,
  data = expression_data,
  metadata = sample_metadata,
  output.dir = "results/ablation",
  params = list(),
  seed = 20260727,
  step = "all",
  input = NULL,
  cache.root = NULL
)
~~~

`ablation()` 只保留 representation 工作流；历史的 `experiment` 参数（含 `"cohort"` 兼容别名与 layered 实验）已在 0.8.3 后移除，不再接受。

`step` 支持以下阶段：

| step | 输入 | 动作 | 返回类 |
| --- | --- | --- | --- |
| `all` | object/data/metadata | 完整运行并写正式产品 | `CCSAblation` |
| `context` | object/data/metadata，或含 object/data 的 input | 准备输入、manifest、配置和缓存布局 | `CCSAblationContext` |
| `plan` | context | 生成节点与 job 描述，不执行计算 | `CCSAblationPlan` |
| `run` | plan；也接受 context 并自动建 plan | 执行 representation runner | `CCSAblationRun` |
| `result` | run | 复制 runner 产品并写 reviewer-facing 结果 | `CCSAblation` |

阶段对象通过 `input` 传递；非 `all` 阶段不会把中间对象伪装成正式结果。

`cache.root` 只指定中间缓存根目录；`output.dir` 保存正式产品，两者不能相同。未提供
`cache.root` 时使用本次 R 会话的临时目录，进程结束后不保证保留，也不会在当前工作目录
创建 `.ccs-cache`。需要恢复、审计或复用缓存的分析必须显式提供 `cache.root`。非空输出
目录默认拒绝写入，需要 `params$output$cover = TRUE` 显式放开。

运行时可用 `CCS_ABLATION_CORES`、`CCS_ABLATION_WORKERS` 和
`CCS_ABLATION_MEMORY_GB` 声明总线程、PSOCK worker 与内存预算。配置在准备完成后按
四个只读表示矩阵的三倍体积估计单 worker 内存；预算不足以容纳一个 worker 时直接
报错，足够时只下调 worker 并重新分配每个 worker 的 XGBoost 线程，不改变 seed、job
顺序或科学参数。

## 输入怎样变成可比较表示

| 函数 | 功能 |
| --- | --- |
| `.ablation_module_manifest()` | 从 d1 列名和模型恢复 tissue、module、d1 block 边界 |
| `.ablation_frozen_feature_manifest()` | 提取 Direct feature、TSP、break vector、feature type |
| `.ablation_flatten_expression()` | 统一矩阵或 `tissue -> cohort -> leaf` 嵌套列表 |
| `.ablation_prepare_metadata()` | 识别别名并统一 sample/cohort/tissue/biology |
| `.ablation_gsclassifier_matrix()` | 按冻结 feature contract 重建 Direct 输入，不拟合新模型 |
| `.ablation_prepare_representation_input()` | reference/query 划分、d1 过滤和 cache key |
| `.ablation_build_manifest()` | 保存样本、feature、模块、版本和 config hash |

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

metadata 至少需要 `sample_id` 和 `cohort`。常见别名会被识别；`tissue` 缺失时才尝试从 CCS 派生，`biology` 缺失时退回 `tissue`。默认主 anchor 是 `cancer_type`，它不是自动生成的核心列，调用方应同时提供独立的 `tissue` 和 `cancer_type`。

### 准备边界：reference 与 external query

`.ablation_prepare_representation_input()` 将样本分成：

- reference：不属于 external cohort，且必须已有 d1 行；
- external query：属于 external cohort，先作为候选 query；
- query 若缺少预计算 d1，会警告、写入 manifest 的 `excluded_query_d1_*` 字段并排除；不会在下游静默补算 d1。

external cohort 默认取 `object@Data$filtered.cohort`，可用 `params$provenance$external_cohorts` 覆盖；`require_external = TRUE` 时没有任何可用 external cohort 会直接报错。

保留样本标记为 `d1_provenance = "in_sample"`（reference）或 `"external_frozen"`（query）。只有主 anchor 为 independent 且两侧 provenance 都满足 external/out-of-fold 时，`evidence_level` 才是 `confirmatory`；否则为 `descriptive`，并在 manifest 记录具体原因。

## 参数 schema

默认值来自 `.ablation_representation_default_params()`，配置由 `.ablation_resolve_representation_config()` 逐字段校验，非法取值在昂贵计算前直接报错：

| 组 | 关键字段 | 业务含义 |
| --- | --- | --- |
| `comparison` | module_ids, direct_group, cohort_group | 选择模块和组名 |
| `provenance` | external_cohorts, max_reference_samples, max_query_samples, require_external | reference/query 来源和样本上限 |
| `anchors` | primary, primary_role, bank_aligned, technical, min_reference_cohorts, endpoint_min_reference_cohorts | 主标签、独立性、技术列和 endpoint 资格 |
| `geometry` | k, search, n_trees, search_k, exact_validation_queries, min_annoy_recall, geometry_samples, distance_pairs | native geometry 与 exact/Annoy retrieval |
| `validation` | enabled, learning_fractions, repeats, inner_folds, lambda, nrounds, min_class_n, numCores, workers | readout 和 learning curve |
| `scaling` | enabled, module_counts, sequences, direct_feature_type, sensitivity_feature_type, biology_anchors, score_reference_samples, score_query_samples, lambda, bootstrap | frozen cohort bank scaling |
| `controls` | null_rp, null_rp_rank, null_rp_seeds, null_perm | paired null |
| `tradeoffs` | decoder, decoder_rank, decoder_lambda, decoder_max_reference_samples, decoder_max_query_samples | d1→Direct decoder |
| `output` | cover, cache_direct | 覆盖输出目录和 Direct cache |

## 主流程：representation 比较

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

`.ablation_linear_readout()` 用 `.ablation_grouped_folds()` 按 cohort 做分组 inner CV，选择 lambda 后用 xgboost `gblinear` 拟合 readout。`.ablation_sample_training_cohorts()`、`.ablation_learning_curve_job()`、`.ablation_learning_curve_execute_job()`、`.ablation_learning_curve_worker()` 和 `.ablation_learning_curve()` 共同完成 fraction × repeat × representation 的串行/PSOCK 调度和 job checkpoint。

### scaling 和 decoder

`params$scaling$enabled = TRUE` 时，`.ablation_representation_scaling()` 通过 `.ablation_cohort_bank_design()` 构造 tissue-balanced nested cohort bank，比较 breadth、within-tissue depth 和 matched-size，并报告非冗余、技术匹配、外部 biology 邻域一致性和稳定性；癌种 readout 仅为 diagnostic lineage。

`params$tradeoffs$decoder = TRUE` 时，`.ablation_decode_direct_features()` 拟合 d1→Direct 的 ridge decoder：`gene_pair` 报告 balanced accuracy/Brier，`single_bin` 和 `set_pair` 报告 Spearman 相关。decoder 是机制诊断，不替换 retrieval endpoint。

## `.ablation_*` 子函数职责地图

下表按业务职责组织，而不是按源码出现顺序罗列。

| 层 | 关键函数 | 业务职责 |
| --- | --- | --- |
| 入口 | ablation, .ablation_dispatch_step | step 分派与阶段调度 |
| 配置 | .ablation_representation_default_params, .ablation_resolve_representation_config | 默认值、合并和逐字段校验 |
| manifest | .ablation_module_manifest, .ablation_frozen_feature_manifest, .ablation_model_features, .ablation_model_break_vectors, .ablation_feature_type, .ablation_extract_tsp_features, .ablation_extract_direct_features, .ablation_model_path_map | 冻结 module/feature contract |
| 输入 | .ablation_flatten_expression, .ablation_prepare_metadata, .ablation_limit_metadata, .ablation_stratified_sample | 表达合并、metadata 别名和样本对齐 |
| representation 输入 | .ablation_prepare_representation_input, .ablation_endpoint_eligibility, .ablation_endpoint_view | reference/query、provenance、资格 |
| Direct 重建 | .ablation_gsclassifier_matrix, .ablation_direct_expression_fingerprint, .ablation_direct_feature_cache_key, .ablation_read_direct_feature_cache | 重建和缓存 Direct 输入 |
| 变换/几何 | .ablation_scale_train_apply, .ablation_module_balanced_transform, .ablation_fit_pca, .ablation_linear_cka, .ablation_distance_spearman, .ablation_knn, .ablation_knn_jaccard, .ablation_effective_rank, .ablation_native_geometry, .ablation_mixing_purity | 防泄漏尺度、PCA 和结构指标 |
| retrieval | .ablation_query_reference_retrieval, .ablation_validate_neighbor_search, .ablation_bind_retrieval, .ablation_neighbor_index_jaccard | query→reference 检索和配对差异 |
| null | .ablation_projection_matrix, .ablation_random_projection, .ablation_permute_blocks, .ablation_null_perm_eligibility | 随机/置换控制 |
| 监督评估 | .ablation_readout_transform, .ablation_xgb_linear_predict, .ablation_linear_readout, .ablation_classification_metrics, .ablation_binary_auc, .ablation_probe, .ablation_eligible_probe_labels | readout、probe、AUROC/BA |
| learning curve | .ablation_sample_training_cohorts, .ablation_learning_curve_job, .ablation_learning_curve_execute_job, .ablation_learning_curve_worker, .ablation_learning_curve, .ablation_make_learning_curve_jobs | fraction/repeat 调度和 job checkpoint |
| scaling | .ablation_representation_scaling, .ablation_representation_scaling_v1, .ablation_cohort_bank_design, .ablation_cohort_bank_matrices, .ablation_module_score_matrix, .ablation_scaling_metric_row, .ablation_score_bank_geometry, .ablation_scaling_coverage, .ablation_score_bank_retrieval, .ablation_scaling_direct_contracts, .ablation_representation_scaling_summary, .ablation_resolve_bank_tissues, .ablation_make_scaling_jobs | frozen bank 构造、评分和汇总 |
| 机制诊断 | .ablation_decode_direct_features, .ablation_selective_reconstruction, .ablation_predict_module_from_direct | d1→Direct decoder 与模块重建 |
| staged 生命周期 | .ablation_make_representation_context, .ablation_make_representation_plan, .ablation_run_representation_plan, .ablation_finalize_representation_stage, .ablation_stage_value | context/plan/run/result 阶段对象 |
| 证据分级 | .ablation_evidence_level | confirmatory/descriptive 判定 |
| cache/audit | .ablation_resolve_cache_layout, .ablation_node_code_identity, .ablation_node_cache_key, .ablation_cached_node, .ablation_read_node_cache, .ablation_read_fit_cache, .ablation_inspect_node_cache, .ablation_promote_legacy_native_geometry, .ablation_atomic_save_rds, .ablation_atomic_write_csv, .ablation_build_manifest | 内容寻址缓存、原子写入、审计 |

## 缓存、随机性和审计

`.ablation_resolve_cache_layout()` 将缓存分为 `context/`（含 `preparation/`）、`plan/`、`nodes/`、`jobs/`、`runner/`、`state/`。`.ablation_node_cache_key()` 绑定实际输入、样本/feature hash、节点参数、seed、schema、algorithm revision、递归代码依赖和运行库版本。

`.ablation_cached_node()` 只有在 state、key、文件 MD5 和 value hash 均一致且为
`complete` 时才命中。state 记录 `run_id`、PID、主机、开始/更新时间、job/参数摘要、
耗时、结果大小和可用的峰值工作集；错误或中断原子写为 `failed`。同主机 owner PID
仍存活时拒绝接管，PID 已退出时先记录 `stale` 与恢复来源再重算；跨主机状态不会仅凭
超时被擅自接管。Direct、native geometry、retrieval、readout、learning curve job、
scaling fit 和 decoder 各有独立缓存边界。

每个随机步骤从显式 seed 派生：readout/learning curve、bank scaling 和 decoder 分别使用独立偏移（如 `+1000`、`+20000`、`+30000` 附近），null 控制直接由主 seed 派生。sample hash、feature hash、config hash、module sequence hash 和 input key 用于复核配对是否真实成立。

## 返回对象和正式文件

主流程最终返回 class 为 `CCSAblation` 的 list，并写出：

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
| audit.csv | retrieval summary 加 evidence 和输入计数 |
| ablation-result.rds | 汇总结果对象 |

Direct 特征与 native geometry 的中间缓存（`direct-feature-cache.rds`、`native-geometry-cache.rds`、`cohort-scaling-fit-cache.rds`）只写入 `cache.root`，不进入正式输出目录。

## 最小使用示例

### 默认整跑

~~~r
result <- ablation(
  object = ccs_fit,
  data = expression_data,
  metadata = metadata,
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

### staged 调用

~~~r
context <- ablation(
  object = ccs_fit, data = nested_data, metadata = metadata,
  step = "context",
  cache.root = "D:/cache/ccs/my-ablation"
)
plan <- ablation(step = "plan", input = context)
run <- ablation(step = "run", input = plan)
result <- ablation(
  step = "result", input = run,
  output.dir = "results/ablation-representation"
)
~~~

## 结果解读边界

1. 先看 `endpoint_eligibility` 和 `evidence_level`；candidate query 不等于可用于 cancer-labelled endpoint 的 estimable query。
2. 把 primary、diagnostic 和 mechanism 分开；lineage readout 或 decoder 不能替代外部 biology 的主证据。
3. Direct/d1 的 `delta_*` 必须在同一 query、fold、sequence 或 resample 内解释。
4. Annoy 必须检查 search validation；recall 低于阈值时增加 `n_trees/search_k` 或改用 exact。
5. d1 缺失 query 会被记录和排除；先检查 manifest 的 `excluded_query_d1_*` 字段。
6. 输出目录默认不可覆盖；确认目录和 `cover` 设置后再运行。

## 维护原则

- 新增指标：在对应 metric helper 返回命名向量或长表，再由统一汇总器 bootstrap；
- 新增 endpoint：挂在“资格判定 → endpoint view → runner → finalize”链上，并同步资格表与 audit；
- 改变输入契约：同步更新 flatten、metadata alias、manifest、资格表和示例；
- 改变输出 schema：同步更新本文档的结果对象与文件表，并保留 `status`/`schema_version`；
- 改变缓存边界：同步更新 cache key、state 校验和 manifest 的 node 状态；
- 改变公共 API：同步更新 `R/ablation.R` 的 roxygen 与 `man/`，再 `devtools::document()`；
- 不把 `.ablation_*` 函数当作独立脚本；新逻辑应挂在“入口 → 准备 → 评估 → 汇总 → 落盘”层次上。
