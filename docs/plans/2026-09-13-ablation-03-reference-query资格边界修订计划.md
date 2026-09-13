# ablation-03 reference/query 资格边界修订实施计划

## 通俗解释：究竟发生了什么

- **一句话说明：** 当前 ablation-03 先用“癌种在 reference 中至少有两个 cohort”筛掉 query，再把剩下的 query 交给几何、技术来源、表达锚点和 readout 等所有分析；这个门槛只对部分癌种标签分析有必要，却被误用成了所有分析的统一入场规则。
- **具体场景：** 这像一座医院先要求每个病例所属疾病必须在旧医院有两个独立病例，才允许病例进入任何检查。这个要求对“比较疾病诊断”有帮助，但对测量影像几何或检查设备偏差并不必要；否则本来可以使用的病例会被提前丢掉。
- **对应到本问题：** external cohort 是待观察的新病例，reference cohorts 是用于建立 bank 和背景的旧病例；`cancer_type` 的跨 cohort 支持是癌种 readout 的资格条件，不是所有表示分析的外部性定义。
- **改变前后：** 现在 43 个 external 候选 cohort 在共享准备阶段只剩 17 个；改进后共享层保留全部可用候选，癌种 readout 仍可只对可估计的子集做主要推断，而 geometry、独立表达 anchor 和结构复现按各自质量条件使用更完整的目标集合。

## 专业判断：问题在哪里

- **当前现象：** `R/ablation.R` 在共享 representation preparation 阶段依据 reference 中每个 `cancer_type` 的 cohort 数过滤 query。当前配置为 `min_reference_cohorts = 2`，因此 43 个 external 候选中仅 17 个进入后续对象；同一过滤结果随后被 geometry、retrieval、technical inference、biological anchor 和 readout 共同消费。
- **已知原因：** 该门槛原本用于避免 probe 通过 cohort 身份猜测 label；当某癌种在 reference 中只有一个 cohort 时，癌种与 cohort/platform 可能完全混杂，癌种 readout 的跨队列解释确实不稳健。
- **核心缺陷：** 这是 endpoint-specific 的“可估计性/反混杂”规则，却被实现成 global query exclusion。它不改变 fit/target 是否重叠，也不是 externality 的定义；对不使用癌种标签的 geometry、连续表达 anchor 和结构复现，提前删除 query 会造成不必要的样本和癌种谱选择。
- **数据边界：** 当前完整 external 候选为 43 个 cohort；主分析中 17 个是通过癌种支持门槛的 estimable query，而不是全部 external query。150 个 reference modules 与 147 个可用 reference cohorts 的差异来自重复 sample ID 清理，属于独立的数据可用性问题，本计划不把两者混为同一资格规则。

## 要达到什么目标

### 完成后的变化

- 建立“候选 query → endpoint-specific qualification → estimable/not_estimable”三层数据合同，不在共享数据层按癌种支持数整体删除 external cohorts。
- 癌种 readout 和依赖癌种背景的 technical endpoint 继续保留明确的跨 cohort 支持门槛，但只在对应 endpoint 内应用。
- geometry、独立表达 anchor、结构复现及其他不依赖癌种标签的分析，使用全部满足自身样本、anchor、距离和无泄漏条件的 external candidates。
- 所有报告和机器可读产物同时记录候选数、endpoint 可估计数、排除数及原因，使“没有结果”与“没有纳入样本”可以区分。
- 保持 Direct 与 d1 在每个 endpoint 使用同一目标样本集合；不重新训练冻结 module，不引入新数据，不修改项目版本号。

### 不在本次处理范围

- 不改变 `min_reference_cohorts = 2` 作为癌种 readout 主要推断的保守门槛；是否提升为 3 或更高，留作后续敏感性分析，不在本次同时决定。
- 不把单 reference cohort 的癌种结果包装成确认性跨队列证据；这些结果可作为描述性或显式 `not_estimable` 记录。
- 不修改 full d1 的 module 数、模块训练过程、外部数据划分来源或样本重复清理规则。
- 不把本次 endpoint-specific 资格修订扩展到 d2、d3、metaCCS、临床结局或新生物标签开发。

## 改进方向

### 分离共享候选层与 endpoint 资格层

共享数据层保留所有通过样本对齐、无 fit/target 重叠和矩阵可用性检查的 external cohort。每个 cohort 生成稳定的资格表，至少区分：`candidate`、`estimable`、`descriptive_only` 和 `not_estimable`，并记录 endpoint、原因、分母、支持 cohort 数、样本数和 hash。

共享层不再直接改写 `prepared$query_metadata` 为癌种可估计子集。各分析函数接收完整 candidate target 与 endpoint mask，或接收经过明确命名的 endpoint-specific view，避免后续阶段误把一个 endpoint 的筛选当成全局事实。

对普通读者而言，这相当于先保留所有能安全检查的病例，再由每项检查说明自己需要哪些资格，而不是用一张诊断规则决定所有检查能否进行。

### 将癌种 readout 的门槛限定在癌种相关 endpoint

癌种 retrieval、癌种 linear readout 及其 learning-curve 继续要求目标癌种在 reference 中至少有两个独立 cohort 支持，作为主要外部推断的最低复制条件。支持不足的 cohort 不从共享数据删除，而是在该 endpoint 标记为 `not_estimable` 或 `descriptive_only`；汇总推断只使用预先定义的 estimable 集合，并报告候选集合与估计集合的差异。

如果某癌种只在一个 reference cohort 出现，仍可保留其样本用于邻居或其他诊断性计算，但不得把该癌种的 readout 结果解释为独立跨 cohort 泛化。若目标类别在 reference 完全不存在，则癌种 readout 直接标记为类别不可预测，并保留原因而不是静默丢弃。

对普通读者而言，这会把“这项病例能否用于癌种诊断验证”与“这项病例能否用于其他检查”分开。

### 让不依赖癌种标签的分析使用完整 external candidates

target-only geometry、Direct/d1 的距离和邻域比较、四个连续表达 anchor、结构可复现性以及不以癌种作为标签的技术诊断，按自身的样本数、anchor 覆盖、距离方差、kNN 可形成性和 fit/target 外部性判断资格，不继承癌种支持门槛。

technical excess 若继续采用“同癌种 reference 背景”这一条件化 estimand，则只对能够建立该背景的 endpoint/cohort 计算；无法建立背景时标记 `not_estimable`。这不应阻止同一 cohort 参与几何或连续 anchor 分析。若未来增加不条件化的技术来源指标，应作为单独 endpoint 命名和校正，不能与现有 cancer-conditioned excess 混合。

对普通读者而言，即使一个病例不能用于公平比较癌种诊断，也仍可能用于判断地图是否变形、连续生物信号是否保留或设备痕迹是否改变。

### 统一报告、缓存与回归审计

更新 manifest 和报告口径，明确展示 `external_candidate_cohort_count`、各 endpoint 的 `estimable_cohort_count`、`not_estimable_count`、`descriptive_only_count` 及原因。正向与反向验证分别记录 fit bank、target candidates 和 endpoint mask；不能只保存最终 17 个 query 的子集而丢失 43 个候选的审计信息。

缓存 hash 必须包含候选集合、endpoint mask、标签支持表、样本键、module manifest 和配置。新增回归测试覆盖：某癌种仅一个 reference cohort、癌种在 reference 完全缺失、同一 external cohort 在不同 endpoint 资格不同、Direct/d1 目标样本错配、fit/target 泄漏和行顺序改变。

## 实施范围与顺序

1. 先冻结 43 个 external candidate 与 147 个可用 reference cohort 的共享数据合同，新增按 endpoint 记录资格和原因的审计表；不再在共享 preparation 阶段删除癌种支持不足的 query。
2. 将癌种支持规则移入 cancer-type retrieval/readout 及癌种条件化 technical endpoint，生成明确的 estimable mask，并验证现有 17 个主要 query 的结果在合同一致时可复现。
3. 将 geometry、连续 anchor 和结构复现切换到完整 candidate target，保留各自的最小样本/实体/距离资格；对不可估计项返回显式状态，不降低门槛强行补数。
4. 更新 ablation-03 的函数、测试、README、Rmd 及 reciprocal-validation 产物 schema；完成定向测试、既有 Stage 04 回归和报告渲染后，再决定是否运行 22 对 22 matched-bank 敏感性分析。

## 如何确认完成

- 共享审计显示 external candidates 未因 cancer-type 支持数被全局删除；当前数据应能同时追溯 43 个候选和各 endpoint 的 estimable 子集。
- cancer-type readout 的主要推断仍只使用预先定义的跨 cohort 支持集合；不支持的类别或 cohort 保留 `not_estimable`/`descriptive_only` 原因、分母和样本数。
- geometry、四个连续 anchor 和结构复现的 Direct/d1 目标样本集合不再被癌种 readout mask 意外截断；两种表示在每个 endpoint 完全配对。
- technical-conditioned endpoint 在缺少同癌种背景时显式不可估计，不把该 cohort 从其他 endpoint 删除。
- 两个方向的 fit/target sample 与 cohort overlap 均为零；所有中心化、尺度、anchor 统计量、CV 和阈值仍只由各自 bank 侧拟合。
- 合成测试覆盖单 reference cohort、无 reference 类别、endpoint 间不同资格、缺失 anchor、单 cohort、零方差、目标样本错配和行顺序变化。
- 使用 `C:\R\R-4.3.1` 运行相关 `test/ablation-03/tests/`、Stage 05 reciprocal-validation、现有 Stage 04 回归及 `02-ablation03-experiment.Rmd` 渲染；检查 HTML 中候选数、估计数和方向标签动态来自机器可读产物。
- `git diff --check` 和 `python -m bac --root . --bac-file docs/contribution.bac verify --json` 通过；实现阶段如行为发生变化，再在 `CHANGELOG.md` 的 `[Unreleased]` 记录，不修改项目版本号。

## 技术补充（按需阅读）

### 建议的数据合同字段

共享方向审计至少保留：`direction`、`bank_role`、`target_role`、`cohort_key`、`candidate_status`、`endpoint`、`qualification_status`、`qualification_reason`、`reference_support_cohort_count`、`target_sample_count`、`fit_target_overlap`、`sample_hash`、`module_manifest_hash` 和 `config_hash`。

manifest 不应再用单一 `query_cohort_count` 同时表达候选数和癌种 readout 可估计数。建议保留旧字段的兼容含义，并新增候选与 endpoint-specific 汇总；分析产物 schema 变更时同步更新读取端和回归测试。

### 主要受影响位置

- `R/ablation.R`：拆分共享 preparation 与 endpoint-specific eligibility，保留完整 candidate matrices，并在 manifest/audit 中写入分层计数。
- `test/ablation-03/02-ablation03-experiment_functions.R`：将当前单一 `min_reference_cohorts` 配置改为癌种 endpoint 的显式配置入口。
- `test/ablation-03/04-ablation03-structural-reproducibility.R` 及其 helper：使用完整 candidate target，结构实体资格单独审计。
- `test/ablation-03/02-ablation03-experiment.Rmd`、`test/ablation-03/README.md`：并列展示 candidate、estimable 和 not-estimable 口径，避免把 17 误写成全部 external cohorts。
- `test/ablation-03/tests/`：新增 endpoint-specific 资格合同测试，并保留现有数据、表示、readout 和结构复现回归。

## 风险与待确认事项

- 采用全部 43 个候选后，geometry 和 anchor 结果的 cohort 异质性可能增大；这是信息恢复而非失败，应通过 cohort-level 区间和覆盖审计表达。
- 癌种 readout 的 2-cohort 门槛仍不能消除平台与癌种混杂；它只提供最低复制条件，不能宣称因果去混杂。
- 某些 endpoint 的候选数与 estimable 数不同，报告和图表必须同时显示分母，避免读者把 `not_estimable` 误解为零效应。
- 反向 external-bank 分析的类别覆盖可能不足；应沿用同一 endpoint-specific 规则，但不能因为癌种 readout 不可估计而删除反向结构或 anchor target。
- 本计划属于已观察结果后的方法学强化；实施后应在报告中标记为 exploratory/sensitivity extension，不把资格边界修订包装成事前确认性检验。
