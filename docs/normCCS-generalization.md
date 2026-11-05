# normCCS 冻结映射与泛化验证

新接口在训练数据中完成状态排序、等级划分和参数选择，再用冻结映射预测新患者。预测不读取结局；评价始终把 response=1 作为阳性，并规定等级越高越可能应答。原 `subtypePerformance()`、`newSubtype()` 和旧 ROC 语义保持兼容。

这些接口及精确搜索适配函数统一放在 [R/phenotype.R](../R/phenotype.R) 底部，由 `R/RcppExports.R` 连接 C++ 内核。

这套流程解决的是给定状态后的映射学习与评价。若上游 metaCCS 状态、特征筛选或 caller 曾使用留出患者的结局，仍须在上游另行隔离；本接口不会重新训练模型库或发现状态。

## 拟合、应用与评价

`fitNormCCS()` 接受 `Subtype`、`Cohort`、`response` 和 `SampleIDs`；列名可配置。`col_source` 指定独立研究来源，默认使用 cohort。一个来源可包含多个 cohort，同一 cohort 不能属于多个来源。患者标识必须在整份输入中唯一；来源别名、重复患者和关联患者须事先处理。相同 ID 会报错，程序不能识别不同 ID 下的同一患者。

```r
fit <- fitNormCCS(
  training,
  k = 3,
  col_source = "study_id",
  rank_marker = "ShrunkLogOR",
  prior_strength = 10,
  min_grade_n = 15,
  min_grade_sources = 2,
  min_source_n = 3,
  min_grade_positive = 2,
  min_grade_negative = 2,
  regularization = 0.5
)
grades <- predict(fit, external)  # 不读取 external$response
evaluation <- evaluateNormCCS(fit, external)
evaluation$by_source
```

这些参数数值用于展示接口，不是经过临床验证的推荐值。拟合对象保存状态顺序、切点、映射、超参数、等级支持、训练来源和患者 ID 的 SHA-256 摘要，不保存逐患者结局。对象可用 `saveRDS()` 保存，用 `readRDS()` 恢复后直接预测。

训练中未见的状态或缺失状态默认预测为 `NA`，也可用 `predict(fit, states, unseen = "error")` 拒绝。评价需要完整状态标签，但允许它们未出现在训练中；结果同时报告总人数、映射人数、覆盖率和仅可映射患者的 AUC。未映射患者不会被临时补到某一级。

`evaluateNormCCS()` 始终拒绝与训练相同的患者 ID，默认也拒绝训练中出现过的来源。仅当研究目标是同来源内的新患者预测时，显式使用 `allow_seen_sources = TRUE`。它不会修改模型、调整映射、翻转方向或重新排序。单一结局类别、没有可映射患者的来源保留 NA AUC 和原因，不从来源均值中悄悄删除。

## 排序收缩与等级约束

新流程默认使用既有 `ShrunkLogOR` 排序，不再使用评价标签重排等级。`prior_strength = 0` 保留已有 Jeffreys 平滑规则；正值在每个训练折把小样本状态的概率向该 cohort 的基线收缩：

\[
p_{0c}=\frac{Y_c+0.5}{N_c+1},\qquad
\widetilde p_{sc}=\frac{Y_{sc}+0.5+\tau p_{0c}}{N_{sc}+1+\tau}.
\]

按各 cohort 的 `logit(tilde_p) - logit(p0)` 等权平均排序。状态未在某 cohort 出现时不增加虚构观测。正的收缩强度只支持 `ShrunkLogOR`，其他既有排序规则在强度为零时可用。收缩可能稳定稀疏状态，也可能减弱真实稀有信号，需要在训练侧验证。

每个候选等级同时满足以下可配置条件：

| 参数 | 含义 | 单次拟合默认值 |
| --- | --- | --- |
| `min_grade_n` | 合并患者数底线 | 1 |
| `min_grade_positive` | 合并应答患者数底线 | 0 |
| `min_grade_negative` | 合并未应答患者数底线 | 0 |
| `min_grade_sources` | 支持该等级的独立来源数 | 1 |
| `min_source_n` | 某来源至少有多少患者进入该级才算支持 | 1 |
| `regularization` | `mean(1/sqrt(n_g))` 的惩罚乘数 | 0 |

事件数底线作用于合并患者，不代表每个来源都有足够事件。合法区间先按训练数据固定，再交给 C++ 完整搜索，人数惩罚同时参与择优。返回结果再用逐记录 midrank AUC 和支持人数独立复核，容差为 `1e-10`。

K 支持 2–8，状态最多 512 个。`max_candidates` 默认一亿，超过预算会报错，不能把不完整搜索当成精确最优；调用方可显式提高预算。内存保护和原生中断继续生效。精确性只针对固定训练排序、K、约束和惩罚后的目标。

`objective` 沿用已有五种聚合定义，但新接口的真实评分坐标是独立来源。`median` 仍包含 pooled `All`，其他目标只使用真实来源。排序仍先计算 cohort 内基线；多个 cohort 属于同一来源时，搜索评价按来源合并。全部 AUC 使用固定方向。

## 嵌套分组验证

`validateNormCCS()` 把每个完整来源分配到一个折。外层估计对未见来源的表现；内层只利用外层训练来源选择配置。每个内层训练折都会重新计算状态顺序、分界和支持约束。候选不能借用全数据排序，也不能因为某个留出来源表现差而把它删去。

```r
grid <- expand.grid(
  k = 2:4,
  prior_strength = c(0, 10),
  regularization = c(0, 0.5),
  min_grade_n = 10,
  min_grade_sources = 2
)
result <- validateNormCCS(
  development,
  grid = grid,
  col_source = "study_id",
  outer_folds = 5,
  inner_folds = 3,
  selection = "one_se",
  min_coverage = 0.95,
  permutations = 199,
  seed = 489
)
result$evaluation$by_source
result$stability$pairwise
final_fit <- result$fit
```

至少需要三个独立来源，并且每个外层训练部分还须能分出两个内层来源。请求折数最多使用现有来源数；不会改成患者随机划分。来源和事件不足、约束过强或训练状态覆盖不足时可能无法选型。

候选须在所有内层折拟合成功，每个留出来源达到 `min_coverage`，并且 AUC 可估计，才能参与比较。选择指标是所有内层留出来源等权的平均固定方向 AUC。`selection = "best"` 选择最高分；`"one_se"` 接受与最高分相差不超过其来源级标准误的方案，依次优先较少等级、较强来源及人数/事件约束、较强惩罚和排序收缩。同分规则预先固定，不参考外层成绩。这个标准误是选型启发式，不能解释为泛化误差的置信区间。

候选全部失败会抛出 `normCCSSelectionError`；错误对象保留 `summary`、`folds` 和每个候选各折的 `evidence`，可用 `tryCatch(..., normCCSSelectionError = identity)` 检查。失败不会被替换为一个临时可运行的配置。

外层每个来源只被评价一次，返回逐来源、逐 cohort 的 AUC 和覆盖率。不同折可能选出不同 K 和映射，因此不报告把不同模型分数混排的 pooled AUC。外层覆盖未达门槛时显式标记 `coverage_failed`。全部开发来源还会独立执行一次分组内层选型，训练最终部署对象；部署对象不是外层验证模型中成绩最好的那一个。

## 稳定性与置换检验

`stability$per_state` 报告各状态在多少个外层训练折中可见，以及归一化等级 `grade/(K-1)` 的均值、标准差和极差。`pairwise` 在两个模型的共同状态上比较排序 Spearman 相关和等级绝对变化；只有 K 相同时才给出精确等级一致率。每折原始状态、等级和排序位置也保留。不同折的切点位置不直接比较，因为它们的状态顺序可能不同。

`permutations > 0` 时，在每个来源内置换结局，固定同一套外层和内层来源划分，每次重新执行内层选型和外层拟合。检验统计量为来源等权外层 AUC，单侧 p 值为 `(1 + 超过或等于观察值的次数)/(B + 1)`。任一置换失败或覆盖不达门槛时保留失败明细和有效次数，p 值为 NA；零次置换返回 `not_run`。

来源内置换要求来源内结局可交换。如果患者相关、治疗分层或其他设计因素破坏交换性，应先调整分析设计；当前检验不自动建模这些因素。它评估给定状态下完整映射学习流程的空假设表现，不能修复上游曾使用测试结局的问题。函数局部管理随机种子，成功或失败都恢复调用者的 RNG 状态。

## 执行与数据管理

新接口不自动写盘，不改变工作目录或全局并行后端。`num_cores` 默认 1，仅控制单次原生搜索线程，避免嵌套验证中的线程超配。实际分析若已有 targets DAG，应把拟合、验证和交付接入原 `_targets.R`；不另设正式运行器。

`result$predictions` 和评价输出包含患者 ID 与结局，应按敏感数据管理。模型与批量结果写入分析自身的外置缓存；包开发验证使用 `D:/cache/ccs/normCCS-generalization/`。测试入口为 `tmp/normCCS-generalization/run.ps1`，包括独立小空间穷举、冻结预测、来源隔离、置换和旧接口回归。编译与包检查另由同目录脚本执行，实际结果见[验证记录](../tmp/normCCS-generalization/README.md)。
