# `ablation()` 业务逻辑

`R/ablation.R` 只有一个公共入口：`ablation()`。它针对已经冻结的 CCS 对象，构造 Direct-GSClassifier 与 Cohort-d1 两种表示，并在同一 reference/query 边界上完成几何、检索、readout、learning curve、可选 cohort-bank scaling 与 decoder 诊断。

## 公共调用

```r
result <- ablation(
  object = ccs_fit,
  data = expression_data,
  metadata = sample_metadata,
  output.dir = "results/ablation",
  params = list(
    geometry = list(search = "exact"),
    validation = list(enabled = TRUE),
    scaling = list(enabled = FALSE)
  ),
  seed = 20260727
)
```

`ablation()` 不再接受 `experiment` 参数。`experiment = "cohort"` 不是兼容入口，也不再映射到 representation；调用者直接使用默认工作流。旧 layered/cohort 编排、Gate 1 和对应参数已经删除，避免保留与正式分析脱节的历史分支。

## 生命周期阶段

需要断点恢复时，使用同一个入口的 `step`：

```r
context <- ablation(
  object = ccs_fit, data = expression_data, metadata = sample_metadata,
  step = "context", cache.root = ".ccs-cache/ablation"
)
plan <- ablation(step = "plan", input = context)
run <- ablation(step = "run", input = plan)
result <- ablation(
  step = "result", input = run,
  output.dir = "results/ablation"
)
```

阶段依次返回 `CCSAblationContext`、`CCSAblationPlan`、`CCSAblationRun` 和 `CCSAblation`。`cache.root` 只保存可复用中间产物，`output.dir` 保存正式结果。

## 参数与输出

参数默认值来自 `.ablation_representation_default_params(seed)`，按 `comparison`、`provenance`、`anchors`、`geometry`、`validation`、`controls`、`tradeoffs`、`scaling` 和 `output` 分组。未知字段在计算前报错；不再支持旧 flat/layered 参数 schema。

主要输出包括：

- `manifest.rds`：输入样本、feature contract、cohort 边界、证据等级和配置哈希；
- `native_geometry.rds`、`retrieval.rds`、`readout.rds`、`learning_curve.rds`：表示比较结果；
- `cohort_scaling.rds`：仅在 `params$scaling$enabled = TRUE` 时生成；
- `tradeoffs.rds`：Direct feature decoder 诊断；
- `audit.csv` 与 `ablation-result.rds`：统一审计表和最终对象。

## 输入边界

`data` 可以是表达矩阵，也可以是 `tissue -> cohort -> list(expr = ...)` 的嵌套结构。函数只消费 CCS 对象中已有的 d1；缺失 query d1 的样本会记录到 `excluded-query-d1.csv` 并排除，不在下游自动重算。
