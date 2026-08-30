# CCS `ablation-02` 方法学审计与重设计摘要

> 项目：CCS（Cohort Congress System）  
> 仓库：<https://github.com/huangwb8/CCS>  
> 审查重点：<https://github.com/huangwb8/CCS/tree/master/test/ablation-02>  
> 审查日期：2026-08-29

## Executive assessment

当前 `ablation-02` 已经是一套相当成熟的 **representation forensics / evidence-boundary audit**：它认真处理了 frozen bank、外部 query 编码、query-to-reference retrieval、同 cohort 邻居排除、cohort-level 统计推断、近似检索 recall 门禁、bank composition 审计以及 contract tests。

但从方法学论文的审稿标准看，它目前主要证明了：

- Cohort-d1 相比 Direct-GSClassifier **确实改变了表示几何**；
- cancer-type 和 technical-source 的局部可恢复性有所下降；
- d1 保留了部分 Direct 信息，但当前 decoder 无法完整恢复 Direct。

它尚未证明：

- d1 去除的是 nuisance variation，而不是无选择地抹平所有信息；
- d1 保留或增强了真正重要的 biological-state information；
- 表现差异来自 CCS 的 cohort-congress mechanism，而不是普通压缩、模块等权、外部监督、特征重复计权或一般性的非线性重编码；
- breadth 与 depth 的效应已被因果或准因果地识别；
- 增加 cohort 数量带来了可迁移的边际生物学信息。

因此，当前最准确的定位是：

| 证据类型 | 当前状态 |
|---|---|
| Evidence of transformation | 强 |
| Evidence of compression / bottleneck | 部分，且 decoder 设计存在混杂 |
| Evidence of invariance | 提示性，只有 nuisance 一侧证据 |
| Evidence of useful biological representation | 尚未建立 |
| Evidence that cohort-congress mechanism is necessary | 尚未建立 |
| Evidence that breadth/depth effects are causal | 尚未建立 |

最核心的判断是：

> 当前 CCS 最大的问题不是 ablation 数量太少，而是缺少一个同时包含 **positive biological targets** 和 **competing mechanistic controls** 的决定性实验。

---

## CCS 的算法—表示—假设—观测链条

### 算法结构

CCS 首先在每个 cohort 内，根据 gene modules 的表达结构构造 cohort-relative pseudo-subtypes。随后，每个 cohort 独立训练一个 GSClassifier：

- 输入为 Direct-GSClassifier 特征，包括 single-bin、gene-pair/TSP 和 set-pair 特征；
- 预测目标为该 cohort 内生成的 pseudo-subtype；
- 每个 cohort model 输出对若干 pseudo-subtypes 的概率。

对一个新样本 $x$，若 Direct 表示为 $\phi(x)$，第 $m$ 个 cohort model 为 $F_m$，则 Cohort-d1 可以写成：

$$
d1(x) = [F_1(\phi(x)), F_2(\phi(x)), \ldots, F_M(\phi(x))].
$$

因此，d1 更接近一组冻结的 **cohort-conditioned probes / experts** 的响应拼接，而不是简单的 cohort voting。

### 合理的理论价值

在 cohort bank 固定后，d1 是 Direct 的确定性函数。它不可能创造 Direct 中原本不存在的信息。它可能产生价值的途径包括：

- 提供更有利的有限样本归纳偏置；
- 降低表示的有效复杂度；
- 抑制 cohort-specific shortcuts；
- 把跨 cohort 一致的状态转化为更容易被低容量模型读出的边界；
- 将外部 cohort bank 学到的先验迁移到新任务。

### 正确的逻辑图

| 层级 | 正确表述 |
|---|---|
| Algorithm | 在每个 cohort 内构造模块对齐的伪亚型，并训练 cohort-specific GSClassifier |
| Representation | 将多个 cohort model 的概率输出拼接为 d1 |
| Hypothesized benefit | 外部 cohort prior 将 Direct 特征重组为更容易跨队列读取的状态证据 |
| Observable consequence | 技术域信息下降，但预定义生物状态被保留或增强；外部泛化优于公平压缩对照 |
| Decisive experiment | cross-fitted d1 + biological anchors + matched compression/metric controls + cohort-structure randomization |

---

## CCS 最合理的 claims

| Claim 类型 | 最合理的表述 | 当前地位 | 当前证据状态 |
|---|---|---|---|
| Algorithmic claim | d1 是多个冻结 cohort-specific GSClassifier 输出的结构化拼接 | 核心 | 支持 |
| Representation claim | cohort bank 实质性改变了 Direct 的几何与可读信息分布 | 核心/支撑 | 支持 |
| Generalization claim | frozen bank 能编码未见 query cohorts | 核心 | query 侧支持，reference 侧存在不对称 |
| Robustness claim | d1 降低技术来源相关的邻域结构 | 支撑 | 部分支持 |
| Biological claim | d1 保留或增强 pathway/TME/subtype 等生物状态 | 潜在核心 | 当前不应 claim |
| Representation improvement | d1 是更好的 biological representation | 潜在核心 | 当前不应 claim |
| Compression claim | d1 是更低维、更高效的表示 | 支撑 | 当前不应按列数 claim |
| Sample-efficiency claim | d1 在较少训练 cohorts 下更易读出重要任务 | 支撑 | cancer-type 任务不支持；其他任务未知 |
| Cohort-evidence claim | 增加独立 cohort evidence 带来可测量的边际收益 | 机制核心 | 探索性 |
| Breadth/depth claim | 新 tissue breadth 与同 tissue depth 有不同信息价值 | 机制核心 | 尚未被识别 |
| Mechanism-specificity claim | 真正 cohort structure 优于随机分组和普通压缩 | 核心机制 | 未检验 |

一个尤其重要的语言边界是：

> “representation changed” 不能写成 “representation improved”。

同样：

> “cancer-type information decreased” 不能写成 “biological-state information increased”。

---

## Claim–evidence map

| Claim | 当前证据 | 能否支持 | 最大漏洞 | 需要的额外证据 |
|---|---|---|---|---|
| d1 是 cohort-conditioned transformation | 源码、输出 block、几何差异 | 可以 | 无主要科学漏洞 | 保持实现与报告可追溯 |
| d1 改变了 geometry | CKA、distance correlation、kNN overlap | 可以，只能支持 changed | 缺少独立 bank/sample 重复 | 多个 bank/sample composition repeats |
| d1 降低 lineage recoverability | retrieval 与 linear readout 均下降 | 可以 | cancer 与 cohort/tissue/technical source 高度绑定 | within-cancer、matched-domain 分析 |
| d1 降低 technical-source recoverability | same-technology neighbor excess 下降 | 部分可以 | 可能是整体信息抹平 | biological anchor 保留 + domain-held-out probe |
| d1 改善 biological representation | 无 positive biological endpoint | 不可以 | 没有正向生物学目标 | external pathway/TME/subtype/outcome anchors |
| d1 是有益压缩 | rank-50 decoder 可恢复部分 Direct | 很弱 | decoder 自身引入额外 PCA 瓶颈 | full-rank decoder、PCA/RP ceiling controls |
| d1 优于一般降维 | 少量 rank-100 RP seeds | 不可以 | rank、有效容量、seed 数不匹配 | PCA/RP rank grid，20–50 seeds |
| 优势来自 cohort organization | 无随机 cohort structure 对照 | 不可以 | 任意监督非线性 ensemble 都可能得到类似结果 | true cohorts vs matched random pseudo-cohorts |
| breadth/depth 价值不同 | internal scaling metrics | 探索性 | composition、sample size、platform、quality 同时变化 | paired marginal cohort addition |
| 增加 cohort evidence 改善 utility | effective rank、redundancy、stability | 不可以 | internal geometry 不是 utility | biological utility 的 marginal gain |
| d1 可泛化到未见 cohort | query external encoding | 部分可以 | reference own-module 是 in-sample | own-block masking / OOF |
| d1 更高效 | 未测 runtime/storage/label efficiency | 不可以 | d1 并不天然更低维或更便宜 | runtime、memory、performance-vs-label-budget |

---

## 当前 `ablation-02` 的主要优点

- query 通过 frozen bank 外部编码；
- retrieval 使用 query-to-reference 而非 query 内自检；
- 排除同 cohort 邻居；
- 主要统计推断以 cohort 而非 sample 为单位；
- 对近似检索实现 recall 门禁；
- 对 bank composition 和数据边界进行了较系统审计；
- 在没有 biological anchors 时明确返回 `not_evaluated`；
- representation、readout、decoder、real encoder 和 scaling 均有 contract tests。

这些设计明显优于把上万个样本直接当作独立生物重复的常见泛癌分析。

但是，统计单位正确并不等于 estimand 正确。特别是 cancer type 通常在每个 query cohort 内恒定，因此当前所谓 cohort-level balanced accuracy 实际上更接近：

> cohort-macro correct-class recall / cohort-equal accuracy。

建议在论文中明确重命名，而不要与全体样本上的 global balanced accuracy 混用。

---

## 当前实验逐项审计

| 实验 | 真正检验的问题 | 最大风险 | 显著时最强结论 | 建议 |
|---|---|---|---|---|
| Native geometry | d1 是否改变 Direct 几何 | 单个 frozen estimate，无天然好坏方向 | 表示发生 transformation | 保留，描述性 |
| Query-to-reference retrieval | d1 是否改变癌种局部可恢复性 | reference own-block in-sample；cancer 是粗标签 | lineage 对邻域的支配减弱 | 保留，改名 lineage diagnostic |
| Cancer linear readout | 低容量模型读出癌种能力是否改变 | information budget 不公平；cohort 内单一癌种 | 癌种线性可读性下降 | 修改 estimand 名称 |
| Learning curves | 差异是否只是小样本效应 | repeats 太少；源码与报告不一致 | cancer readout 差异不是小样本偶然 | 补跑或降 supplement |
| Technical-source excess | d1 是否减少技术同源邻居 | 可能只是普遍信息损失 | technical-domain association 下降 | 保留但必须搭配 biology endpoint |
| Module variance concentration | 变化是否集中在少数 modules | module 相关，variance 受校准影响 | 仅描述 variation distribution | Supplement |
| Decoder | d1 中多少 Direct 信息可线性恢复 | d1 先被压到 rank 50；固定 ridge | rank-50 linear decodability | 重写解释 |
| Null-RP | 随机投影能否复现 lineage 改写 | 仅 3 seeds，rank 不匹配 | 只排除少数特定 RP | 重做后保留 |
| Null-Perm | 标签关系是否偶然 | cohort 内 cancer label 恒定，置换无效 | 无有效零分布 | 删除 |
| TSP sensitivity | 结果是否由 gene-pair 子集决定 | 不是 dimension-matched control | 仅说明 feature-type sensitivity | Supplement |
| Breadth sequence | 新 tissue 是否增加非冗余证据 | tissue、quality、platform、sample size 同时变化 | 只能描述特定加入路径 | 大幅修改 |
| Depth sequence | 同 tissue 增加 cohort 是否有收益 | 每一级同时给多个 tissue 加 cohort | 只能描述当前 eligible tissues 的路径 | 大幅修改 |
| Matched-size breadth/depth | 相同 module 数下 tissue diversity 是否有价值 | 只匹配 module count，未匹配质量与信息容量 | 不足以回答异质 cohort 是否更有价值 | 降 supplement |
| Neighborhood stability | 不同 bank compositions 是否得到相似邻域 | 稳定不等于正确 | compositional reproducibility | 保留为 secondary metric |

---

## 最薄弱的三个环节

### 缺少 positive biological anchor

目前最强的“正向”证据是技术来源聚集下降，但没有证明：

- pathway/TME/state 得以保留；
- molecular subtype 得以保留；
- downstream utility 得以保留；
- 跨 cohort 生物一致性提高。

因此无法区分：

- selective invariance；
- indiscriminate smoothing。

这是当前最可能导致 Reviewer #2 拒稿的缺口。

### 缺少能识别 cohort mechanism 的公平 comparator

Direct 与 d1 同时差异在：

- 非线性表示函数；
- 外部 pretrained cohort prior；
- feature grouping；
- feature reuse；
- module-level weighting；
- 表示维数与有效秩。

因此，即使 d1 将来优于 Direct，也不能自动归因于 cohort-congress mechanism。

### reference 与 query 的表示 provenance 不对称

query d1 是完全 external frozen encoding，而 reference sample 的自身 cohort block 来自一个训练时看过该样本的模型。即使 retrieval 排除了同 cohort 邻居，这种 fingerprint 仍然保留在 reference vector 中。

最低成本敏感性分析：

- mask reference sample 的 own-cohort block；
- 对 query 和 reference 使用相同剩余 blocks；
- 重跑 geometry、retrieval、readout 和 technical-source analyses。

决定性版本：

- 对 reference sample 的自身 block 生成 K-fold OOF prediction。

---

## Dimensionality 与 compression 的关键问题

当前实验中：

- Direct：529 features；
- Cohort-d1：583 probability columns。

因此 d1 不能被描述为普通意义上的“降维表示”。更准确的语言是：

- structured probability re-encoding；
- information bottleneck；
- effective-rank compression；
- task-aligned nonlinear representation。

当前 decoder 先将 d1 压到 rank 50，再用固定 $\lambda=1$ 的 ridge 恢复 Direct。因此，它最多说明：

> Direct 在 rank-50 linear summary of d1 中的可恢复程度。

它不能证明：

> 完整 583 列 d1 本身不可逆或已经丢失了对应信息。

建议的 comparator matrix：

- Direct-full；
- Direct-blocked/module-balanced；
- Direct-reuse-weighted；
- PCA(Direct)；
- RandomProjection(Direct)，20–50 seeds；
- d1-full；
- PCA(d1)；
- d1-own-block-masked；
- d1-own-block-OOF；
- 可选 Autoencoder(Direct)。

建议比较整条 **utility–effective-rank Pareto frontier**，而不是只比较一个 rank。

### 不同结果的解释

| 结果 | 最合理解释 |
|---|---|
| d1 ≈ PCA/RP(Direct) | 主要是 generic compression / regularization |
| d1 ≈ Direct-blocked | 主要是 module weighting 与 feature reuse |
| d1 > PCA/RP/blocked Direct，但 random pseudo-cohorts ≈ d1 | 主要是 generic supervised nonlinear ensemble |
| true cohorts > random groups，但 shuffled subtype ≈ true cohorts | cohort grouping 有用，pseudo-target 语义未被证明 |
| true cohorts > shuffled target 与 random groups | 对 cohort organization + cohort-relative target 的机制证据最强 |
| d1 只在 query tissue 已存在于 bank 时有效 | 主要是 tissue matching，不是 pan-cancer transfer |
| d1 低标签量更好、全数据相近 | 支持 finite-sample inductive bias |
| d1 所有标签量都更差 | 可能丢失任务相关信息或 comparator 不公平 |

---

## Cohort-evidence scaling 的核心问题

### 当前 breadth/depth 是否被识别？

当前设计在构造层面区分了：

- breadth：加入此前未覆盖 tissue 的首个 cohort；
- depth：固定 tissue set 后增加同 tissue cohort；
- matched-size：比较相同 module 数的 breadth-heavy 与 depth-heavy bank。

但这不等于因果效应被识别。

breadth 同时改变：

- tissue identity；
- cohort identity；
- platform/source；
- cohort sample size；
- heterogeneity；
- pseudo-label entropy；
- model calibration；
- feature coverage。

matched-size 目前只匹配 module count，没有匹配：

- 总训练样本数；
- module 输出列数；
- pseudo-subtype class balance；
- submodel quality；
- platform/source；
- tumor composition；
- cohort heterogeneity；
- effective rank。

因此，它不能回答：

> 一个新的 heterogeneous cohort 是否比同 tissue 下增加一个 cohort 更有价值？

它只能回答：

> 按照当前构造的两种 bank composition，在相同 module 数下，若干内部几何指标是否不同？

### 更有解释力的 marginal cohort contribution

建议把 scaling 的中心问题改为：

> 固定 bank $B$，加入一个候选 cohort $c$ 后，对固定 query set 的边际影响是什么？

定义：

$$
\Delta U_c(B)=U(B\cup c)-U(B),
$$

$$
\Delta T_c(B)=T(B\cup c)-T(B),
$$

$$
\Delta G_c(B)=G(B\cup c)-G(B),
$$

$$
\Delta R_c(B)=R(B\cup c)-R(B).
$$

其中：

- $\Delta U$：biological utility；
- $\Delta T$：technical leakage；
- $\Delta G$：geometry / stability；
- $\Delta R$：redundancy。

候选 cohort 配对比较：

- same-tissue vs new-tissue；
- technically similar vs technically different；
- biologically similar vs biologically different；
- high-quality vs low-quality submodel；
- high-heterogeneity vs low-heterogeneity cohort。

候选之间尽量匹配：

- sample size；
- assay/platform；
- output block width；
- pseudo-label entropy；
- submodel CV quality；
- missing-gene rate；
- tumor/normal composition。

统计单位应为：

> base-bank × candidate-cohort addition

而不是 query sample。

这套设计比当前 breadth/depth curve 更接近可解释的 intervention。

---

## Biological-anchor strategy

当前最优先的 biological anchors 不是机械增加几十个 pathways，而是选择少量、互补、可跨 cohort 使用的目标。

### Proliferation

使用冻结的外部 signatures，例如 Hallmark E2F targets 和 G2M checkpoint。

价值：

- 几乎所有癌种内部都有连续变异；
- 跨 cohort 可比较；
- 不依赖稀疏临床 metadata；
- 更接近 tumor-intrinsic biological state。

### Immune/TME state

优先考虑：

- ESTIMATE immune/stromal scores；
- 或冻结的 IFN-γ / inflammatory signatures。

它可以检验 d1 是否保留与微环境和细胞组成相关的连续状态。

### BRCA PAM50

PAM50 是明确的乳腺癌 intrinsic subtype system，适合在多个独立 BRCA cohorts 上做 leave-one-cohort-out subtype readout。

### CRC CMS

CMS 是结直肠癌中广泛使用的 transcriptomic subtype system，适合检验 d1 是否保留疾病内而不是疾病间的生物结构。

### 避免新的 circularity

因为 CCS pseudo-subtypes 本身来自 gene-set/module expression，biological anchor 不应完全由同一组输入 genes 构造。

建议：

- 使用冻结的外部 signature；
- 尽量排除 CCS 的 32 个输入 genes；
- 最好再排除 Direct gene-pairs 涉及的 genes；
- 同时报告 full-signature 与 disjoint-gene signature；
- anchor 不参与 bank selection、rank selection 或 metric selection。

### 推荐成功判据

d1 不必在所有 biological endpoints 上显著优于 Direct。更合理的 selective-invariance criterion 是：

1. biological-anchor performance 对 Direct 非劣效；
2. technical-source recoverability 显著下降；
3. worst-cohort performance 不恶化；
4. 至少一个低标签量或跨平台场景优于公平压缩对照。

| Biology | Technical leakage | 解释 |
|---|---|---|
| 上升 | 下降 | 最强 selective representation |
| 基本不变 | 下降 | 支持 selective invariance |
| 上升 | 基本不变 | 支持 task-aligned reorganization |
| 下降 | 下降 | indiscriminate smoothing / trade-off |
| 下降 | 基本不变 | 纯信息损失 |
| 基本不变 | 基本不变 | 当前任务上无实用增益 |

---

## Redesigned ablation hierarchy

### Tier 1 — 必做

| 实验 | Scientific question | Comparator | Statistical unit | Primary metric | Reviewer impact |
|---|---|---|---|---|---|
| Reference own-block cross-fitting | 当前结果是否受训练内 fingerprint 影响？ | original、masked、OOF d1 | query cohort/tissue | 所有核心 endpoint 的 paired delta | 极高 |
| External biological-anchor benchmark | d1 是否选择性保留/增强生物状态？ | Direct、blocked Direct、PCA、RP、d1 | cohort/tissue | within-cancer AUROC、Spearman、worst-cohort | 极高 |
| Matched compression/metric controls | 差异是否只是 compression/weighting？ | Direct、blocked Direct、PCA、RP、d1 | cohort/tissue，RP seed 仅作算法重复 | external anchor utility | 极高 |
| True cohort structure vs random structure | cohort organization 是否必要？ | 真 cohort、matched random pseudo-cohort、shuffled subtype labels | bank repeat/cohort/tissue | biological utility 与 technical leakage | 极高 |

### Tier 2 — 高价值

| 实验 | Scientific question | 价值 |
|---|---|---|
| Paired marginal cohort contribution | 新 tissue、同 tissue depth、技术差异各自带来多少边际信息？ | 直接重构 breadth/depth 机制证据 |
| Leave-one-platform/source-out | d1 是否真正提高技术域外泛化？ | 将 technical clustering 下降转化为外部 utility |
| Leave-one-tissue-out bank | d1 是否依赖 query tissue 已存在于 bank？ | 区分 pan-cancer transfer 与 tissue lookup |
| Decoder ceiling analysis | 当前 reconstruction failure 来自 d1 还是 decoder？ | 修正 compression/information-loss claim |

### Tier 3 — Nice-to-have

- Autoencoder-compressed Direct；
- module dropout / module importance；
- random TSP / gene-pair controls；
- d2/d3/metaCCS stability；
- runtime、memory、storage benchmark。

---

## Reviewer #2 simulation

按可能导致拒稿的严重程度排序：

1. **论文把 cancer-type 和 technical-source recoverability 下降等同于 biological representation 改善，但没有任何独立生物状态或下游任务。** — 完全未解决。
2. **reference atlas 包含 own-cohort model 的 in-sample d1 block，而 query 完全 out-of-sample，表示 provenance 不对称。** — 已承认但未解决。
3. **Direct 与 d1 不共享相同 information/inductive-bias budget；d1 使用 150 个 pretrained models 和 module-balanced weighting。** — 完全未解决。
4. **没有随机 pseudo-cohort、shuffled pseudo-subtype 或 non-cohort ensemble，因此 cohort organization 的必要性没有被证明。** — 完全未解决。
5. **breadth/depth 分析没有真正隔离 breadth 或 depth；matched-size 只匹配 module count。** — 部分意识到但未解决。
6. **information-loss claim 基于 rank-50 PCA + 固定 ridge，无法归因于 d1 本身。** — 部分解决但解释过强。
7. **cancer type 被当作 independent confirmatory anchor，但它在 cohort 内恒定并与 tissue/cohort composition 绑定。** — 未解决。
8. **same-platform neighborhood excess 下降不等于 batch correction，也可能是局部结构普遍丢失。** — 部分解决。
9. **learning-curve repeats 在当前源码与报告中不一致，artifact provenance 不清楚。** — 未解决。
10. **尚未证明对 unseen tissue、unseen platform、downstream d2/d3/metaCCS 或效率的优势。** — 部分或完全未解决。

---

## 实验信息增益排序

| Rank | Experiment | Scientific information gain | Reviewer impact | Cost | Priority |
|---:|---|---|---|---|---|
| 1 | External biological-anchor benchmark + Direct/PCA/RP/blocked-Direct controls | 极高 | 极高 | 中 | 立即 |
| 2 | Reference own-block masking + OOF sensitivity | 高 | 极高 | 低到中高 | 立即 |
| 3 | True cohorts vs matched random pseudo-cohorts / shuffled labels | 极高 | 极高 | 高 | 核心机制优先 |
| 4 | Paired marginal cohort contribution | 极高 | 很高 | 中高 | anchors 后执行 |
| 5 | Leave-one-platform/source-out | 高 | 很高 | 中 | 高 |
| 6 | Full-rank decoder + PCA/RP ceilings | 中高 | 中高 | 低中 | 高 |
| 7 | Leave-one-tissue-out bank | 中高 | 高 | 中 | 高 |
| 8 | Autoencoder-compressed Direct | 中 | 中低 | 高 | 延后 |
| 9 | Random TSP/gene-pair controls | 中低 | 中低 | 中 | Supplement |
| 10 | Module dropout | 中 | 中 | 中高 | 核心结论成立后 |

### 如果只能做一个实验

做一个多臂、锁定式的：

> **Cross-fitted selective-invariance benchmark**

最小 arms：

- Direct；
- module-balanced/blocked Direct；
- PCA(Direct)；
- RP(Direct)，多个 seeds；
- d1-own-block-masked；
- 最好再加 d1-own-block-OOF。

主要 endpoints：

- proliferation；
- immune/TME；
- PAM50 或 CMS；
- secondary technical-domain recoverability。

这一个实验能同时回答：

- d1 有没有 biological utility；
- 是否只是 compression；
- 是否只是 weighting；
- 是否具有 selective invariance。

### 如果只能做三个实验

- external biological-anchor + matched-control benchmark；
- own-module masking / OOF；
- true cohort bank vs matched random pseudo-cohort bank。

它们分别回答：

- 有没有用；
- 现有证据是否干净；
- 是否真由 cohort congress 导致。

### 如果只能做五个实验

在上述三个基础上增加：

- paired marginal cohort contribution；
- leave-one-platform/source-out selective-invariance。

### 最可能改变 Reviewer #2 判断的实验

> External biological-anchor benchmark。

### 即使结果对 CCS 不利，也最值得做的实验

> True cohort structure vs matched random pseudo-cohort structure。

如果随机分组与真 cohort bank 相同，这会否定 “cohort congress” 作为独特机制，但能及时防止项目继续在错误机制上堆实验，并帮助重新定义论文贡献。

---

## If I were the PI, I would do next…

### 先冻结证据版本

- 固定 commit；
- 清空或严格校验旧缓存；
- 统一 learning repeats；
- 写入 Git SHA、config hash、data/model manifest hash、RDS 时间戳；
- 将 cancer type 从 `independent confirmatory` 改为 `diagnostic_lineage`；
- 将整体 evidence level 暂时设为 `exploratory_mechanistic`。

### 立即运行 own-module-masked sensitivity

这是成本最低、最能快速判断现有 retrieval/readout 证据是否受 reference fingerprint 影响的实验。

### 接入三个 biological anchors

优先：

- proliferation；
- immune/TME；
- PAM50 或 CMS。

不要先扩展几十个高度相关的 pathway signatures。

### 补齐公平 comparator

- Direct-blocked/module-balanced；
- PCA rank grid；
- RP rank grid；
- d1-PCA rank grid。

Autoencoder 可以暂缓。

### 设置 Gate 1

只有满足以下条件，才继续投入新的 scaling：

- d1 在至少一个外部生物 anchor 上优于公平 controls，或在多个 anchors 上非劣效；
- technical leakage 同时下降；
- own-module-masked/OOF 结论稳定；
- 结果不是只发生在 bank 已覆盖的同 tissue query 中。

### Gate 1 通过后，再做 marginal cohort addition

当前 breadth/depth curves 可以保留为 exploratory supplement，但不应继续承担主要机制 claim。

---

## Minimal publishable ablation set

一个正规 bioinformatics / computational biology journal 可接受的最低集合，建议包括：

- 可追溯的数据与版本 manifest；
- patient/sample overlap 审计；
- query external encoding；
- reference own-block-masked sensitivity，最好关键结果有 OOF；
- Direct-full、Direct-blocked、PCA、RP、d1 的公平比较；
- proliferation、immune/TME、PAM50 或 CMS；
- biology 与 technical leakage 的联合 selective-invariance analysis；
- cohort-level / tissue-level inference；
- leave-one-cohort-out；
- worst-cohort performance；
- full-rank / rank-grid decoder；
- geometry、cancer retrieval、module variance 和旧 breadth/depth 降为机制性或 supplement 结果。

最低可支持的中心结论应是：

> d1 所削弱的 lineage/technical information 对预定义的 biological-state readout 并非必要，并且这种选择性在未见 cohort 中可复现。

---

## Strongest-possible ablation set

最强版本可进一步包括：

- 所有 reference own-block 的完整 cross-fitting；
- patient-level overlap 与镜像数据去重；
- Direct、blocked Direct、PCA、RP、AE、d1 的 rank–utility Pareto curves；
- true cohorts、random pseudo-cohorts、shuffled pseudo-labels；
- 多个 pan-cancer continuous anchors；
- 多疾病 molecular subtype anchors；
- mutation-associated transcriptional states；
- 至少一个外部 clinical/outcome task；
- leave-one-cohort-out；
- leave-one-platform-out；
- leave-one-source-out；
- leave-one-tissue-out bank；
- paired marginal cohort addition；
- candidate biological/technical novelty 的 factorial analysis；
- 基于 marginal utility 的 saturation analysis；
- d2/d3/metaCCS 在不同 bank、seed 和 query composition 下的稳定性；
- 完全锁定的 external validation collection；
- runtime、memory、storage 和 label-efficiency 分析。

---

## Bottom line

- 当前 `ablation-02` 最强的结论是 **d1 进行了结构化重编码**，不是 **d1 更好**。
- cancer-type recoverability 下降不能解释为 biological-state information 增加。
- technical-source clustering 下降不能单独解释为 batch correction。
- 当前最危险的缺口是没有 positive biological anchor。
- reference own-cohort block 的 in-sample provenance 必须通过 masking / OOF 处理。
- Direct 与 d1 的 information/inductive-bias budget 不公平，必须加入 blocked Direct、PCA 和多 seed RP。
- d1 为 583 列、Direct 为 529 维，不能 claim literal dimensionality reduction。
- rank-50 decoder 只能支持 rank-50 linear decodability，不能证明完整 d1 的信息不可逆。
- 当前 breadth/depth matched-size 只匹配 module 数，不能回答 heterogeneous cohort 的边际价值。
- 下一步最高信息增益实验是：**cross-fitted/own-block-masked d1，在 proliferation、immune/TME 和 PAM50/CMS 上，与 Direct、blocked Direct、PCA 和 RP 做严格的外部 cohort-level 对照。**
