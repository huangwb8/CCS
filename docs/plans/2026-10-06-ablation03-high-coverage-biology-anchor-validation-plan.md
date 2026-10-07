# ablation-03 高覆盖生物锚点验证优化计划

## 实施约定

2026-10-07 修订：负责人在查看增殖选定合同的 4,848 个背景基因后，明确授权将共同秩背景门槛由 5,000 降为 3,000。当前执行以[门槛修订计划](2026-10-07-ablation03-background-threshold-3000.md)及配置为准；下文 5,000 门槛与禁止自行降低的规则记录原设计，不阻止此次明确授权的修订。

按负责人本轮要求，本计划作为增量分析实施：广队列基线及低覆盖机制诊断作为不同证据层保留，全部按当前五锚点重新计算，高覆盖结果独立保存并在报告中并列解释。测量合同、覆盖率档和推断规则先于新效应冻结；不按结果调整签名或选队列。

新增设计使用 `test/ablation-03/config/biology-high-coverage.yml`，避免改写只读原始输入。分析专用函数位于 `test/ablation-03/R/biology_high_coverage.R`，由既有 `_targets.R` 发现，I/O 在 `targets/biology_high_coverage.R`；计算复用已安装 CCS，不改公共包源码或版本。正式结果只写入现有 cache root 的 `biology-high-coverage/`。

队列搜索采用冻结的可用性种子规则：按原签名基因可用频率生成全体、单队列及两队列交集种子；每个种子纳入全部兼容且有癌种 reference 支持的队列，再对完整有限值求实际共同签名并稳定收敛。候选按 query 队列数、已知来源数、query 数、共同基因数、reference 队列数依次比较，相同者保留固定顺序首项。这是可审计的有界搜索，不证明全局最大广度。所有档位保留，主档按本计划的 8 队列／3 来源目标选择；不足时保留有限验证或敏感性层级。

同一队列可以包含不同平台的样本。平台、assay 和癌种标签保留在样本层，技术限制按样本判断；资格汇总表列出队列内实际出现的全部标签。队列角色与来源组保持可追溯，不能因混合平台静默排除整个有效队列。

| 估计对象 | 人群与单位 | 区间与检验 | 结果去向 |
|---|---|---|---|
| 高覆盖邻居效用 d1−Direct | 同一冻结合同和评分资格下的配对 query；先队列内平均，再队列等权 | 2,000 次 cohort bootstrap 95% CI；预设双侧 cohort sign-flip 零差检验；主五锚点固定 BH 分母 5 | primary、grid、cohort |
| 连续读取器标准化 MAE d1−Direct | reference-only 训练，外部 query 配对；队列等权 | 同一 cohort 条件区间；各评分／读取器五锚点 BH；正值更差 | readout_metrics、readout_inference |
| 距离／技术限制差中差 | 先在同一患者匹配原／修改分支，再队列等权 | 条件 cohort 区间；敏感性不做新零假设检验 | sensitivity |
| C_core、留一队列／来源 | 保持冻结签名；C_core 取五个秩同癌种主比较实际通过人数门槛的共同队列，合同交集另存为描述 | 条件 cohort 区间，未做来源层随机抽样外推；不发布 p/q | core、leave_one |
| 覆盖、来源与基因数 | 冻结输入描述 | 无抽样目的，不给 CI/p | frontier、eligibility、gene_lists |

原尺度与秩资格可能因非有限背景而不同，各格明确保存实际人数；跨评分／候选池的机制差须在共同患者上配对，不能直接相减不同人群均值。没有科学依据的非劣效界值时只做估计，不宣称等价或保留。

## 当前方案与科学问题

IFNγ 与 IL6-JAK-STAT3 分别使用原签名，当前共五个锚点。原合并通路方案及其结果已弃用；本文件不保留旧结果快照。查看原合并结果后的修订时点、正式重算与验收记录见[通路拆分重分析计划](2026-10-06-ablation03-IFNγ与IL6锚点拆分重分析计划.md)，最新结果见[生物锚点报告](../../test/ablation-03/02.02.00.%20生物锚点分析.html)。

高覆盖验证要解决的是测量对象是否仍能代表原签名，以及在这个测量合同下，d1 与 Direct 的近邻效用和连续信息可读取性是否存在差异。广队列分析可能混合跨队列测量定义和癌种组成的影响；低覆盖残余基因仅作机制探索。每项锚点独立选择可测队列，避免为统一版式再取最低公分母。

五锚点 × 四覆盖档 × 两评分 × 两候选池形成 80 个预设格子。广队列基线、严格诊断、低覆盖探索和高覆盖验证全部按五项定义更新；任何不可估计项保留状态，不降低门槛，也不从五项 BH 分母移除。测量合同依据可用性选择，不读取两臂效应。

## 要达到什么目标

### 完成后的变化

1. 为每个预设 biological anchor 建立一个**完全由 measurement eligibility 决定**的 High-Coverage Biological Anchor Validation Set。
2. 在不查看 d1−Direct biological effect 的前提下冻结：
   - query cohorts；
   - reference cohorts；
   - 实际使用的 signature genes；
   - score definition；
   - candidate-pool rule；
   - inference rule；
   - 随机种子和 bootstrap 规则。
3. 以“**高覆盖完整 signature + rank-based score + same-cancer candidate pool**”作为首选的 measurement-controlled 主比较；仅当预设 rank 背景门禁通过时发布该主结果。
4. 同时保留 raw/rank × all-cancer/same-cancer 的完整预定对照，用于区分尺度效应与 lineage composition。
5. 增加 reference-only continuous readout，使结果能够区分：
   - 表示中 biological information 本身难以恢复；
   - biological information 仍在，但默认 d1 距离没有形成相应的邻居结构。
6. 保留原 43-cohort baseline，作为 breadth / heterogeneous deployment 条件下的现实表现，不覆盖、不重写、不删除。
7. 最终论文能够明确回答“在什么 measurement contract 下 biology 被保留到什么程度”，而不是只给一个模糊的 “preserved / lost” 二元结论。

### 不在本次处理范围

- 不重新训练 CCS/d1 bank。
- 不更换或根据结果挑选新的 biological signatures。
- 不删除原始 43-cohort baseline。
- 不根据观察到的 d1−Direct 方向决定 cohort、基因、阈值或 platform。
- 不把 expression anchor 当作独立病理、蛋白、治疗反应或临床结局验证。
- 不把 CI 跨 0 自动解释为等价或无损。
- 不为了得到更好结果而搜索“最优”距离、最优候选池或最优 cohort 子集。
- 不将 low-coverage exploratory proxy 升级成完整 anchor 的确认性证据。

## 改进方向

### 建立 outcome-blind 的高覆盖 cohort 资格合同

当前设计把全部契约 cohort 一次性求共同基因交集，优先保证 cohort breadth，代价是 signature identity 被严重稀释。新的资格合同应反过来：先保护 biological construct，再在满足 construct validity 的前提下最大化 cohort breadth。

对每个 anchor 独立建立二元基因可用性矩阵：

`cohort × signature gene`

cohort 的选择只能读取以下信息：

- cohort role（reference/query）；
- signature gene 是否存在；
- 样本中是否为有限值；
- query 数；
- cancer label 与证据状态；
- reference cancer support；
- assay/source/platform 等预先存在的 measurement metadata；
- d1 external-frozen provenance。

**禁止读取：**

- Direct utility；
- d1 utility；
- d1−Direct delta；
- readout MAE；
- 任何根据 d1“表现好不好”产生的派生结果。

这样可以保证后续 subset 不是 outcome-driven cherry-picking。

对普通读者而言，这一步意味着：先决定“哪些医院真的做了足够完整的同一种检查”，再比较两个算法，而不是先看算法在哪些医院表现好，再挑医院。

### 用 coverage–breadth frontier 决定 cohort 数量，而不是固定追求 43 个 cohort

预先定义若干 measurement coverage 层级，例如：

- 60%
- 70%
- 80%
- 90%

每个层级都只计算资格，不计算 d1−Direct biological effect，并输出：

- 共同 signature gene 数；
- 可纳入 query cohort 数；
- 可纳入 query 数；
- 独立 source-system 数；
- 支持这些 query 的 reference cohort 数；
- rank background 的共同基因数。

推荐使用以下**结果无关的冻结规则**选择 primary coverage tier：

> 对每个 anchor，选择能够保留至少 8 个独立 query cohorts、至少 3 个独立 source groups、每个 query cohort 至少 20 个有效 query，并满足癌种 reference 支持的最高 coverage tier。

这里的 “8 个 cohort / 3 个 source groups” 是 publication-grade 的设计目标，不是数学定理。如果实际数据无法达到：

- 5–7 个 query cohorts：保留为有限外部验证，明确 breadth 边界；
- 少于 5 个：不升级为 primary preservation 证据，只作为 high-coverage sensitivity；
- 当前统计代码的 `min_inference_cohorts = 3` 可继续作为“能否计算区间”的技术下限，但不能等同于“足够有说服力”的论文标准。

所有预设 coverage tiers 都应保留在结果中形成一条 coverage–breadth frontier，避免只展示最终挑中的一档。

### 每个 anchor 使用自己的高覆盖 cohort 集合

不再要求 proliferation、immune、stromal、IFNγ、IL6-JAK-STAT3 必须使用完全相同的 query cohorts。

分别定义：

- `C_proliferation`
- `C_immune`
- `C_stromal`
- `C_ifn`
- `C_il6`

原因是五个 signatures 的基因数和平台可测性差异很大。为了统一图形而取五者交集，会再次把分析拖回“最小公分母”。

另定义：

`C_core = C_proliferation ∩ C_immune ∩ C_stromal ∩ C_ifn ∩ C_il6`

只作为五锚点横向可比的 sensitivity analysis。

主结论按 anchor-specific validation set 报告；`C_core` 用于回答“如果强制同一 cohort 人群，方向是否一致”。

### 对每个 anchor 冻结唯一的共同基因集合

即使 cohort 被判定为“高覆盖”，也不能让不同 cohort 各用不同的 signature 子集计算同一个 anchor。

对每个 anchor，在最终冻结的 query + reference measurement set 上定义：

`G_anchor* = 所有纳入 cohort 共同拥有的原 signature genes`

之后：

- 所有 reference 样本；
- 所有 query 样本；
- Direct/d1 两臂；
- raw/rank 两种评分；

必须使用完全相同的 `G_anchor*`。

必须保存：

- 原 signature gene count；
- `G_anchor*` gene count；
- coverage fraction；
- 缺失的原 signature genes；
- gene-list hash。

任何样本级缺失不得通过“这个样本少用一个基因”静默解决。主分析要求对冻结 `G_anchor*` 完整可测；不满足者作为明确缺失处理。

这一步的意义是：确保所有队列里的 “proliferation score” 真的是同一个数学对象。

### 将 full-signature rank + same-cancer 设为首选主比较

新的 high-coverage 分析继续保留 2×2 设计：

| Score | Candidate pool |
|---|---|
| high-coverage raw | all cancers |
| high-coverage raw | same cancer |
| high-coverage rank | all cancers |
| high-coverage rank | same cancer |

首选 primary endpoint：

> **high-coverage rank score × same-cancer candidate pool**

理由：

- high coverage 最大程度恢复原 biological signature；
- rank score 降低跨平台绝对表达尺度和单调变换差异；
- same-cancer pool 减少 lineage/cancer composition 对邻居关系的贡献；
- 仍使用冻结的外部 query，不参与模型训练或调参。

rank score 只有在冻结 cohort 集合的共同背景基因达到当前预设 `min_background_genes = 5000` 时才可作为 primary。如果门禁不通过，不降低阈值以追求结果；该 anchor 的 rank-primary 状态应为 `not_estimable`，同时报告 high-coverage raw 和 measurement limitation。

其余三格属于预设机制对照，不根据结果决定是否展示。

### 对 query 与 reference 同时实施 measurement eligibility

只筛 query cohort 不够。

如果 query 的 anchor score 很完整，而 top-15 reference neighbor 来自 signature 严重缺失的 cohort，utility 仍然不是统一测量。

因此每个 anchor 需要建立独立的**score-eligible reference atlas**：

- reference cohort 满足同一个 `G_anchor*`；
- reference cohort 与 query cohort 不重叠；
- same-cancer 分支至少有 2 个独立 reference cohorts；
- 候选 reference sample 数不少于 `k = 15`；
- 所有 candidate samples 对 `G_anchor*` 有完整有限值。

Direct 和 d1 必须使用相同的 reference candidate pool。

不能先用全 reference atlas 检索，再把没有 anchor score 的 neighbor 丢掉；这样会让两臂因缺失模式不同而产生 selection bias。

### 把“邻居效用”与“信息可读取性”分开

单独看 `ΔU` 无法判断信息丢失发生在哪一层。

建议在同一 high-coverage anchor target 上同步运行 reference-only continuous readout：

- Direct → biological score
- d1 → biological score

优先保留现有 ridge readout；只有满足冻结的 reference-CV 触发规则时才运行 XGBoost，两臂同时运行。

外部 query 只用于最终评估，不参与：

- lambda 选择；
- feature scaling 拟合；
- nonlinear-reader 触发；
- threshold 调整。

重点联合解释两类指标：

| 邻居 utility | Continuous readout | 更支持的解释 |
|---|---|---|
| d1 ≈ Direct | d1 ≈ Direct | biology 在当前任务下基本保留 |
| d1 < Direct | d1 ≈ Direct | 信息仍可读取，但默认 d1 几何/距离没有形成相应局部邻居 |
| d1 < Direct | d1 < Direct | 更支持表示本身丢失或压缩了该 continuous biological information |
| d1 ≈ Direct | d1 < Direct | 需要检查局部邻居任务与全局 readout 是否测到不同结构，不能简单归因 |

这一步比单纯追求 `ΔU ≈ 0` 更重要，因为它直接区分“表示内容”和“表示几何”。

### 不用“CI 跨 0”证明 biology preserved

如果论文要使用 “preserved”“no meaningful loss” 或 “non-inferior” 一类措辞，需要预先给出一个 biologically meaningful loss margin：

`δ_anchor > 0`

保留门禁：

`lower 95% CI(ΔU) > -δ_anchor`

才能支持“在预设容忍损失范围内未见有意义下降”。

`δ_anchor` 不能在看到 high-coverage d1−Direct 结果后反推。优先选择以下依据之一：

1. 外部文献中同一 signature 的重复测量/平台重现误差；
2. reference atlas 中预先设计的 split-half / resampling reliability，对 anchor measurement 自身噪声进行校准；
3. 独立的临床/生物学效应尺度，说明多大的 utility 变化才足以改变实际 biological-neighbor interpretation。

如果无法给出有科学依据的 margin，本轮仍然可以非常有价值地报告 estimate + CI，但只写：

> “the disadvantage was substantially attenuated under a high-coverage, lineage-controlled measurement contract”

而不写：

> “biology was proven to be preserved.”

### 保留 breadth baseline，建立三层证据结构

最终报告按以下层级组织：

#### Primary：High-Coverage Biological Anchor Validation

回答：

> 在 measurement-controlled 条件下，d1 是否保留局部 continuous biology？

核心条件：

- anchor-specific high-coverage cohort set；
- frozen common signature；
- rank score（若背景门禁通过）；
- same-cancer candidate pool；
- cohort-equal paired inference。

#### Secondary：Broad-Cohort Robustness

保留原 43-cohort / 14,060-query baseline，回答：

> 在真实 heterogeneous multi-cohort deployment 条件下，d1 与 Direct 的 biological-neighbor utility 有什么差异？

不覆盖当前正式结果，不改变其定义。

#### Exploratory：Mechanism Diagnostics

保留：

- low-coverage proxy；
- raw ↔ rank；
- all-cancer ↔ same-cancer；
- distance changes；
- technical restriction；
- module diagnostics；
- leave-one sensitivity。

它们解释“为什么会出现差异”，但不承担完整 biology preservation 的确认性结论。

### 增加 coverage threshold 的稳定性曲线

不要只输出最终 high-coverage threshold 的一个点。

对预设的 60/70/80/90% tiers，全部保存 measurement eligibility；在满足最小推断人数的 tiers 上，按同一冻结分析规则计算 effect，形成：

- x 轴：signature coverage threshold；
- y 轴：d1−Direct utility；
- 同时显示有效 cohort 数和 query 数。

这张图回答一个非常重要的问题：

> 当 biological construct 越来越完整时，d1−Direct 的结论是否稳定？

可能出现的模式：

- coverage 越高，`ΔU` 越接近 0：支持 measurement artifact 参与原负差；
- coverage 越高，`ΔU` 仍稳定约 -0.05：更支持真实 local biological geometry trade-off；
- effect 随 threshold 大幅震荡：说明结论高度依赖 measurement population，不能用单一全局结论概括。

所有 tier 必须预先固定并完整报告，禁止只显示方向最有利的一档。

### 增加 source-level 与 leave-one 稳健性

主分析仍以 cohort 为主要独立单位，但如果多个 cohort 来自同一 source system，不能自动把它们看成完全独立的泛化重复。

在 high-coverage validation set 中同步报告：

- leave-one-cohort；
- leave-one-source；
- 每个 source 的 cohort 数；
- source-level direction consistency。

如果有效 source groups 足够，可增加 source-block bootstrap 作为补充；如果不足，则明确说明推断主要是 cohort-level conditional inference。

不因为某个 cohort 对结果不利而单独删除。任何排除必须由预先冻结的 measurement/QC 规则触发。

## 实施范围与顺序

1. **先冻结 measurement-only 设计，不运行新的 d1−Direct biological effect。**  
   从当前正式 biology cache 构建每个 anchor 的 cohort × gene availability matrix、coverage–breadth frontier、query/reference 资格和 source 分布。此阶段的产物不得包含 utility delta 或 readout performance。

2. **按预设规则生成 anchor-specific validation contracts。**  
   对每个 anchor 固定 coverage tier、query cohorts、reference cohorts、`G_anchor*`、rank background、same-cancer reference support 和所有 hash；生成独立 contract 文件，正式 effect target 只消费 contract，不自行再次筛选 cohort。

3. **在冻结 contract 上重算高覆盖 biological scores 与近邻。**  
   Direct/d1 使用相同 query、相同 score-eligible reference atlas、相同 `G_anchor*` 和相同 k；分别产生 raw/rank × all/same-cancer 四格，不从原 low-coverage proxy 继承基因集合。

4. **运行配对 query → cohort-equal 推断。**  
   先在同一 query 上形成 d1−Direct difference，再在 cohort 内求均值，最后 cohort 等权；沿用 2,000 次 cohort bootstrap、有效次数门禁和预设 sign-flip/BH 框架。不同 anchor 可有不同 cohort 集合，分母必须显式报告。

5. **同步运行 continuous readout。**  
   使用相同高覆盖 biological targets，reference-only 训练，external query 只评估；将 readout 与 neighbor utility 联合解释。

6. **完成预设 sensitivity。**  
   包括 coverage tiers、`C_core`、raw/rank、all/same-cancer、leave-one-cohort、leave-one-source，以及已有 distance/technical diagnostics 中与新 contract 可合法复用的部分。

7. **更新报告证据层级。**  
   生物锚点报告首先展示 measurement contract 和 coverage–breadth frontier，然后展示 high-coverage primary，之后才展示 broad baseline 与 exploratory diagnostics，避免读者先看到显著性再寻找测量定义。

8. **保持原科学结果不可变。**  
   原 `anchor-inference.csv`、`anchor-cohort-deltas.csv` 和现有 low-coverage diagnostic 作为既有证据保留；新结果写入独立 `high_coverage` 产物，不能覆盖旧文件。

## 建议的数据与产物边界

沿用现有 `_targets.R`，由分析专用 `test/ablation-03/R/biology_high_coverage.R` 调用已安装 CCS 的既有计算函数；不建立第二套 runner，不修改包源码。

建议新增独立配置，例如：

`test/ablation-03/config/biology-high-coverage.yml`

配置只保存冻结设计参数，例如：

- coverage tiers；
- primary tier selection rule；
- minimum query cohort target；
- minimum source-group target；
- minimum query per cohort；
- minimum reference cohorts；
- rank-background gate；
- primary score/pool；
- bootstrap 与 seed；
- 是否启用 non-inferiority margin。

建议增加语义清晰的 targets：

- `biology_high_coverage_measurement_frontier`
- `biology_high_coverage_contracts`
- `biology_high_coverage_scores`
- `biology_high_coverage_retrieval`
- `biology_high_coverage_inference`
- `biology_high_coverage_readout`
- `biology_high_coverage_sensitivity`

具体名称可以调整，但必须保持：

> measurement selection target 不依赖任何 effect target。

正式报告建议导出至少以下聚合产物：

- high-coverage measurement frontier；
- per-anchor frozen cohort contract；
- per-anchor frozen gene list/hash；
- query/reference eligibility；
- primary inference；
- raw/rank × all/same-cancer grid；
- cohort-level deltas；
- coverage-tier sensitivity；
- continuous-readout metrics；
- leave-one-cohort/source；
- branch status / not-estimable reason。

患者级表达、患者级邻居清单和模型对象继续只保留在受控 cache，不提交仓库。

## 建议图表

### Figure A：Coverage–Breadth Frontier

每个 anchor 单独展示：

- signature coverage threshold；
- eligible query cohort count；
- query count；
- source-group count；
- common gene count。

作用：证明 cohort reduction 是 measurement-driven，而不是 result-driven。

### Figure B：High-Coverage Primary Forest Plot

展示：

- d1−Direct utility；
- 95% CI；
- cohort count；
- query count；
- common/total genes；
- coverage fraction。

如果存在预设 non-inferiority margin，同时画 `-δ_anchor`，不要只画 0。

### Figure C：Coverage-Tier Stability

展示 60/70/80/90% 各 tier 的 effect 与 CI，并在副轴或标签中标记 cohort/query 数。

作用：直接观察“提高 construct fidelity 后结论怎样变化”。

### Figure D：Utility × Readout Mechanism Map

横轴：

`Δ neighbor utility`

纵轴：

`Δ readout error` 或统一为“d1 相对 Direct 的 biological performance”。

五个 anchor 各一个点/区间，帮助区分：

- information loss；
- geometry mismatch；
- preserved information。

### Figure E：Cohort Heterogeneity / Leave-One

保留 cohort-level delta 与 leave-one-source/leave-one-cohort 结果，避免总体均值掩盖极端 cohort。

## 如何确认完成

### 科学资格

- [ ] 每个 primary anchor 都有独立的 frozen measurement contract。
- [ ] cohort 选择代码不读取 Direct/d1 utility、delta、readout performance。
- [ ] 每个 anchor 的 query 与 reference 使用同一个 `G_anchor*`。
- [ ] `G_anchor*` 属于原预设 signature，不引入根据结果挑选的新基因。
- [ ] 所有 primary query cohort 达到冻结的 coverage、query-count 和 provenance 门禁。
- [ ] same-cancer primary 中，每个 query 的 candidate pool 满足至少两个独立 reference cohorts 和完整 top-15。
- [ ] Direct 与 d1 的 query、candidate pool、anchor score 和缺失规则完全配对。
- [ ] rank-primary 只有在共同背景达到冻结门槛时才发布；不足时保留 not-estimable，不事后放宽。

### 统计推断

- [ ] 先在 query 级配对，再 cohort 内汇总，再 cohort 等权。
- [ ] 不把数千个患者当成数千个独立跨 cohort 重复。
- [ ] 95% CI 明确条件于冻结 bank、measurement contract、reference atlas 与 query cohorts。
- [ ] 五个 anchor 的预设检验族完整保留，并按既定规则做 multiplicity control。
- [ ] CI 跨 0 不写成 equivalence。
- [ ] 只有在 `δ_anchor` 于 effect 计算前有独立依据并冻结时，才进行 non-inferiority / preservation claim。
- [ ] leave-one-cohort 和 leave-one-source 不按结果方向筛选。

### 防止结果驱动选择

建议增加一个自动化不变量测试：

> 将现有 effect/delta 列删除、打乱或替换为随机数后，measurement contract 输出必须逐字节一致。

同时检查：

- [ ] cohort list hash 不受 effect 文件变化影响；
- [ ] gene list hash 不受 effect 文件变化影响；
- [ ] primary coverage tier 不受 effect 文件变化影响；
- [ ] threshold selection target 的依赖图中不存在 biological inference/readout target。

这是本轮最重要的防 cherry-picking 证据之一。

### 回归与可复现性

- [ ] 当前五锚点广队列 baseline 与诊断回算的点估计、人数在 1e-10 容差内一致，不复用旧合并效应。
- [ ] 高覆盖分支不覆盖当前广队列科学文件；本次拆分不改变原表示分析文件哈希。
- [ ] 低覆盖 exploratory 结果保持独立，不被新高覆盖结果静默替换。
- [ ] 合成 fixture 覆盖：不同 cohort 缺不同 signature genes、样本级非有限值、reference 支持不足、rank background 不足、source 重复、癌种标签未知。
- [ ] 同一 frozen contract 重跑得到相同 cohort/gene list 与相同 seed 派生结果。
- [ ] Rmd 只读取 targets 产物，不在报告阶段重新筛 cohort 或重新计算 effect。
- [ ] 正式 store、cache root、source-code identity 和 installed CCS package identity 继续通过现有 provenance 门禁。

## 结果出来后如何判读

### 情形 A：高覆盖 rank + same-cancer 后仍约 `ΔU ≈ -0.05`

更支持：

> d1 相比 Direct 存在真实且跨 cohort 稳定的局部 continuous biological geometry trade-off。

如果 continuous readout 同样变差，进一步支持 information compression/loss；如果 readout 接近 Direct，则更偏向 geometry/distance 问题。

### 情形 B：高覆盖后 `ΔU` 明显缩小并接近 0

更支持：

> 原 broad-cohort negative anchor effect 至少部分来自 cohort-dependent biological measurement 和/或 lineage composition，而不能全部解释为 d1 biological information loss。

若有预设 `δ_anchor` 且 CI 完整落在 non-inferiority boundary 内，才可进一步讨论“无有意义损失”。

### 情形 C：neighbor utility 仍差，但 d1 readout 与 Direct 接近

这是非常有价值的机制结果：

> biology 仍编码在 d1 中，但当前欧氏距离/模块权重不能把这种信息组织成局部邻居。

下一步应优化 metric/observer，而不是直接修改 representation。

### 情形 D：high-coverage rank、raw、same-cancer、all-cancer 全部一致变差，readout 也变差

这是最强的 representation-loss 证据。

此时可以更有底气把 biology preservation 写成 CCS 当前版本的真实 trade-off，并讨论未来如何加入 biology-preserving constraint。

### 情形 E：结果随 coverage tier 剧烈变化

说明：

> biological preservation 不是一个与 measurement contract 无关的单一常数。

这本身也是重要方法学发现，应把 measurement invariance 作为 CCS 跨 cohort 表示评估的一部分，而不是只追求一个“最好看”的 effect。

## 风险与待确认事项

- 高覆盖 subset 可能显著减少 cohort breadth。减少 cohort 本身不是失败，只要 measurement contract 更可信；但需要诚实区分“construct validity 提高”和“external breadth 减少”。
- IFNγ 与 IL6-JAK-STAT3 的签名大小、可测性分别检查，达到 80–90% 全局共同覆盖的支持可能不同。不要为了图形一致降低其它 anchor 的标准；每个 anchor 独立选择最高可行 tier。
- 当前许多 cohort 的 `expression_unit` / preprocessing evidence 仍为 unknown。rank scoring 可以降低但不能完全消除这一风险；不能因此宣称 platform invariance 已被证明。
- same-cancer pool 会改变科学问题：它回答“癌种内部 continuous biological state 是否保留”，不是“泛癌所有邻居关系是否保留”。因此 all-cancer 分支仍需保留。
- subset 选择即使只基于 gene availability，也是在已看过原 broad result 后提出的新设计，应在论文中标记为后续冻结的 validation/sensitivity，而不要伪装成最初未见数据的完全 confirmatory study。
- 如果某 anchor 最终只能保留很少 source groups，即使患者数很大，也不能把 query 数量误当作跨来源泛化证据。
- 非劣效界值是最容易被 reviewer 质疑的部分；没有可辩护的 biological/measurement 依据时，宁可只做 estimation，也不要事后选一个能让结论通过的 margin。

## 推荐的最终证据叙事

本轮优化完成后，ablation-03 的生物学部分建议形成以下逻辑链：

1. **Broad observation**  
   报告当前五锚点在广队列中的各自效应、区间与有效人数，并说明测量和癌种组成限制。

2. **Measurement audit**  
   强制全 cohort 共同测量使 signature coverage 严重下降，说明 broad result 混合了 representation 与 biological measurement heterogeneity。

3. **High-coverage validation**  
   在 outcome-blind、high-coverage、anchor-specific validation sets 中，用完整度明显更高的统一 signature 重新定义 biological state。

4. **Lineage-controlled primary**  
   以 high-coverage rank + same-cancer 检验癌种内部 continuous biological geometry。

5. **Mechanism separation**  
   用 continuous readout 判断 biology 是真的难以恢复，还是仍存在于 d1 但默认 metric 没有利用。

6. **Breadth–fidelity synthesis**  
   将 high-coverage validation 与原 43-cohort baseline 并列，而不是互相取代，最终回答：

   > CCS/d1 在追求跨 cohort 稳定表示时，哪些 biological structures 能在测量一致条件下保留，哪些局部关系会发生真实 trade-off，以及观察到的 trade-off 中有多少依赖于 cohort-specific measurement system。

这比单独追求一个 `ΔU ≈ 0` 更具有方法学价值，也更符合 CCS “cohort effect 不只是 nuisance，而是观察 biology 的一部分”这一核心研究方向。
