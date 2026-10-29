# CCS normCCS 搜索算法优化计划：保留穷举与新增动态规划

> 目标仓库：[`huangwb8/CCS`](https://github.com/huangwb8/CCS)，默认分支 `master`  
> 核心修改位置：[`R/phenotype.R`](https://github.com/huangwb8/CCS/blob/master/R/phenotype.R)、[`R/ccs_base.R`](https://github.com/huangwb8/CCS/blob/master/R/ccs_base.R)  
> 制订日期：2026-10-29  
> 状态：实施计划，**未修改仓库代码、未运行 CCS 实际队列数据**。

## 研究目标与基本判断

当前 CCS 通过大量 `metaCCS` 分子状态构建固定的 `normCCS` 分级。现有论文只在 2、3、4 级候选中选择最优方案；这一结论不能推出四级是更大分级空间的最优解。增加分级数，可能减少分数并列、保留更多 ICI 应答排序信息，但不保证独立队列 AUC 一定提升。

本次工程改造只解决一个清晰的问题：**给现有 `subtypePerformance()` 增加算法选择参数，使旧穷举行为不变；实现一个以动态规划为核心、能高效探索 `n_norm_subtype >= 5` 的新搜索方法，尽可能逼近旧方法定义的优化目标。** 不改变 CCS 主体、32 基因面板、metaCCS 生成算法，也不自动修改论文的既有四级映射。

预期实现原则：正确性优先、旧 API 向后兼容、固定排序与同一批候选映射、输出可追溯、运行时间及内存可控。动态规划算法要清楚标示**哪里是精确最优、哪里只是近似优化**。

## 源码核查与实际瓶颈

已读取仓库当前源码，关键事实如下。

| 位置 | 已确认行为 | 改造含义 |
|---|---|---|
| `R/phenotype.R:40–74` | 对外导出的 `subtypePerformance()` 提供 `n_norm_subtype`、`norm_roc_method`、`norm_rank_marker`、`numCores` 等参数 | 在这里增加一个公共算法选择参数，向下透传 |
| `R/phenotype.R:217–230` | 仅当 `n_norm_subtype != NULL` 时调用内部 `subtypeNorm()` | 保持 `NULL` 语义不变 |
| `R/phenotype.R:846–889` | `subtypeNorm()` 计算各原始状态跨队列 `medianNRR` / `medianRR`，按 `norm_rank_marker` 固定排序 | **不改排序规则**；DP 与穷举必须基于同一排序 |
| `R/phenotype.R:890` | `combn(1:(M-1), K-1)` 预先生成全部连续切点 | 主要的组合爆炸与内存风险 |
| `R/phenotype.R:895–980` | `foreach` + `doParallel` 遍历切点；每个组合重新映射全部患者、重新计算总体和队列 AUC | 主要的重复计算成本 |
| `R/phenotype.R:988–1015` | 对所有候选以 `Info` 分组，计算 `median(ROCAUC, na.rm=TRUE)`，按中位数降序返回；`Raw` 含 `All` 与逐队列行 | **真实遗留目标是 AUC 中位数，而且包含 `All` 汇总行** |
| `R/ccs_base.R:516–585` | `multiXClass_roc()` 按每一类别的 one-vs-rest ROC AUC 对类别重新排序，再计算整体 ROC | `OptimizeRank` 与固定等级的原始 AUC 不等价 |
| `R/ccs_base.R:635–682` | `get_perform_markers()` 中二级采用 `compareRealPred2()`；三级及以上依据 `Raw` / `OptimizeRank` 调用不同性能计算 | DP 必须维持这些现存统计分支的边界 |

**候选空间的性质很重要**：现有穷举并非允许任意 metaCCS 状态随意合并，而是在预定 RR/NRR 排序上选择 `K-1` 个切点，得到 `K` 个**连续区间**。这恰好是动态规划最擅长处理的问题。新 DP 只要保持排序，就不需要缩小现有合法候选空间。

组合数为 `choose(M-1, K-1)`（`M` 为推导数据中实际出现并进入排序的原始状态数，不应直接使用全部 216 作为实际搜索维度）：

| M | K=4 | K=5 | K=6 | K=8 | K=10 |
|---:|---:|---:|---:|---:|---:|
| 107 | 192,920 | 4,967,690 | 101,340,876 | 24,370,067,800 | 3,283,866,636,050 |
| 216 | 1,633,355 | 86,567,815 | 3,653,161,793 | 3,817,554,073,685 | 2,282,897,336,063,630 |

上述计数是组合数学结果，不是实测执行时间。当前实现会先物化 `combn()` 矩阵，再对每个候选生成字符串 `Info` 和多行 AUC 结果，计算与内存成本均随候选数急剧增长。

需要特别注意的两个统计边界：

- 旧 `medianROCAUC` 计算同时纳入 `All` 与各 `Cohort` 行。这是遗留目标的精确定义，不等于对 20 个队列单独求中位数。新算法应先复现这个目标；另在诊断中报告“仅真实队列”的中位数，以便后续独立决定是否修订目标，**不得静默改变**。
- 默认 `norm_roc_method='OptimizeRank'`：`multiXClass_roc()` 使用评价人群自身的应答标签重排类别顺序。若在留出或外部队列上同样重排，就不是固定分型顺序的无泄漏预测评估。此次工程保持遗留模式的可比性，但真正用于论文主要外部验证时应以**锁定等级顺序的 Raw ROC AUC**为主要指标；`OptimizeRank` 只作为明确标注的探索性遗留评分模式。

## 公共 API：只增加一个必要参数

建议新增：

```r
subtypePerformance(
  ..., 
  n_norm_subtype = 2,
  norm_search_method = c("bruteforce", "dp"),
  ...
)
```

并向内部透传：

```r
subtypeNorm(
  df,
  n_norm_subtype = 2,
  norm_roc_method = "Raw",
  norm_rank_marker = "RR",
  norm_search_method = c("bruteforce", "dp"),
  numCores = NULL,
  verbose = TRUE
)
```

内部用 `match.arg(norm_search_method)` 解析。**为了避免破坏依赖位置参数的旧脚本，实施时建议将新参数附加到现有形参列表末尾，而不是机械照抄上面的展示顺序。**

- `"bruteforce"`：默认；走遗留算法。现有成功调用、表结构、`Info` 字符串、排序依据原则上不变。
- `"dp"`：执行“预统计 → 精确均值 AUC 动态规划 → 原始中位数目标局部精修 → 精确复核”。在 K≥5 主要使用，K=2–4 用于一致性基准测试。
- `"bruteforce"` 运行前仅计算组合数量（用 `lchoose()` 等避免数值溢出），对预计产生灾难性内存分配的组合应提前提示并安全停止，明确建议改用 `"dp"`；小规模兼容测试不受影响。
- `n_norm_subtype=NULL`：继续不执行 normalization 搜索，不因新增参数而改变行为。
- 新模式不自动根据 K 替用户选择算法，不静默回退到暴力搜索。
- 输入校验：`K` 必须为整数，`2 <= K <= M`；验证二分类应答标签、状态编号、队列内有效性，以及 `numCores` 正整数、缺失和极端情况。确认不合法时早期报错，不能让 `combn()` 或 `pROC` 在深层报难以理解的异常。

使用示例：

```r
# 旧流程：不传新参数，保留原行为
res_old <- subtypePerformance(data = df, n_norm_subtype = 4)

# 显式启用原暴力算法
res_bf <- subtypePerformance(
  data = df, n_norm_subtype = 4,
  norm_search_method = "bruteforce"
)

# 搜索 8 级：启用新算法，并建议使用固定等级的原始 AUC
res_dp <- subtypePerformance(
  data = df, n_norm_subtype = 8,
  norm_search_method = "dp",
  norm_roc_method = "Raw",
  norm_rank_marker = "NRR",
  verbose = TRUE
)
```

## 算法总体设计：DP 负责全局高效搜索，局部精修逼近遗留中位数最优

### 固定排序与充分统计量

首先严格复用旧 `subtypeNorm()` 计算出的排序 `s_1,...,s_M`，从**低应答倾向到高应答倾向**统一编码（旧代码实际先降序排列，再把首段编码为最高等级；反转存储时须保持最终映射一致）。

只需对每个队列、每个排序状态预聚合两项计数：

- `P[c,i]`：队列 `c` 中处于状态 `i` 的 responder 人数；
- `N[c,i]`：队列 `c` 中处于状态 `i` 的 non-responder 人数。

另外建立 `All` 的 pooled 计数（注意这不是独立队列）。预计算前缀和、有效队列掩码和每个队列的 `P_total`、`N_total`。仅有一种结局的队列原样按遗留评估口径设为不可估计/NA，不错误算为 AUC=0.5。可通过计数计算任何区间的应答/非应答数与固定等级 AUC；无需为每个切点重新复制 1,177 行患者数据。

### 精确 DP 的目标：固定排序下的**加权平均 Raw AUC**

对于某个队列，定义状态 `i` 的人数 `p_i,n_i`，所有应答与非应答总数为 `P,N`。把每个状态都作为单独等级时的 AUC 记为 `AUC_fine`。

如果将连续区间 `[a,b]` 合并为同一等级，相比细粒度排序，AUC 的变化量可用区间内的正负样本配对直接计算：

`loss_c(a,b) = [sum_{a <= i < j <= b} (n_i p_j - p_i n_j)] / (2 P_c N_c)`。

因此，对一个由 `K` 个连续区间组成的分组，

`AUC_c(K blocks) = AUC_fine,c - sum_blocks loss_c(block)`。

这里 `loss` 可以为负数：如果两个状态的实际应答率与先验排序相反，把它们合并反而可能改善 AUC。因此**不得假定增加 K 就必然提高性能**。

使用固定队列权重 `w_c`（第一版采用有效队列等权，`All` 不进入 DP 的权重）定义区间代价：

`C(a,b) = sum_c w_c * loss_c(a,b)`。

把区间代价预计算后，进行标准区间动态规划：

```text
DP[k, j] = min_{i=k-1,...,j-1} {
    DP[k-1, i] + C(i+1, j)
}
```

`DP[1,j] = C(1,j)`；记录 `backpointer`，回溯得到 `K-1` 个切点。每个区间至少包含一个实际排序状态。计算复杂度约为 `O(C*M^2 + K*M^2)`，空间约 `O(C*M^2 + K*M)`（能进一步优化）。**对于“固定排序 + 固定加权平均 Raw AUC”这个代理目标，DP 找到的是全局最优连续分段。**

重要说明：上式是严格的排序 AUC 公式，但实际代码 `pROC::multiclass.roc()` 的自动 `direction` 与二级分组特殊分支可能存在口径差异。实施前必须在随机和真实样本子集上验证“计数公式 = 旧函数 Raw AUC”，并检查 `direction`、factor levels、NA、单类别队列、并列值；不一致时先定位原因，不能悄悄使用近似公式冒充旧结果。

### 中位数不是可加目标：不能直接声称 DP 精确复现暴力最优

遗留穷举的最终依据是：

`medianROCAUC = median(c(AUC_All, AUC_cohort_1, ..., AUC_cohort_C), na.rm=TRUE)`。

中位数是多个队列 AUC 的非线性函数；无法一般性地分解为独立区间代价 `C(a,b)`。因此标准 `DP[k,j]` **不能直接保证找到遗留目标的全局最优**，即使 DP 已经在均值 AUC 上获得严格全局最优。

建议采用一个工程上简单、质量可检验的两阶段方法：

**第一阶段：生成一组高质量 DP 初始候选。**

1. 同一份预聚合数据上，分别运行：仅有效队列等权均值、适当加入 pooled AUC 权重的均值、强调当前较弱队列的确定性 reweight 方案等少数固定权重组合。
2. 每次 DP 均返回一组切点，去重后得到候选池。第一版可采用 5–10 个固定权重方案，不引入额外公共参数。
3. 保留一个“分组尽可能均衡”的参照切点作为诊断初值；不要依靠它替代 DP。

**第二阶段：在原始目标上做局部精修。**

1. 为候选切点快速计算逐队列固定顺序 AUC，并按遗留定义加入 `All` 行，计算精确的 `medianROCAUC`。
2. 对每个候选，逐个移动某一个边界到其他**合法**位置（保持边界有序、所有组非空）；评估其原始中位数目标，接受严格改善，迭代到无单边界改善。
3. 为减少局部最优影响，再测试有限的双边界协同移动，或采用小规模 beam search 保留若干高分解；作为后续增强选项，第一版先实现单边界全位置搜索 + 多初值。
4. 多起点结果按**原始** `medianROCAUC` 排序；若中位数相同，使用预设固定次级规则（例如 cohort-only mean Raw AUC、排序切点字典序）保证可复现。
   - 每个候选的快速评分应基于事先聚合的 cohort × grade 响应计数；仅对进入最终候选集的映射运行完整的既有 pROC 函数。若快速评分与旧实现不一致，应回退使用精确评分函数，并将差异作为测试失败调查。
5. 对排在前面的候选，调用仓库现有 `get_perform_markers()` 复核 `All` 和各队列 AUC，确保最终 `Raw` / `Summary` 不是只由代理公式生成。

`OptimizeRank` 的特殊处理：既然它会对每个候选、每个队列使用真实结局进行重新排序，不能把它当成固定顺序 AUC 的区间可加目标。**允许 DP 在 Raw AUC 上寻找初值，再在 `OptimizeRank` 的遗留评分函数上做候选精修和最终复核**；输出清楚说明“代理搜索 + `OptimizeRank` 目标近似精修”。如精修开销过大，优先限制候选数，而不是假称得到精确全局最优。正式外部效果的评价不能在目标队列上重新执行 `OptimizeRank`。

不建议第一版直接引入遗传算法、强化学习、贝叶斯优化或外部优化包。现在的状态空间有明显的连续分段结构，先利用这个结构更合适。

## 内部代码组织建议

仅把新增逻辑放入 `R/phenotype.R` 中的私有 helper（如内容明显膨胀，可独立成 `R/subtype_norm_dp.R`；R 包会自动加载），不新增不必要的导出函数。

```text
subtypePerformance(..., norm_search_method="bruteforce")
  └─ subtypeNorm(..., norm_search_method=...)
       ├─ prepare_norm_rank()          # 共用：保持现有排序和标签语义
       ├─ if bruteforce:
       │    └─ subtype_norm_bruteforce() # 封装旧搜索，尽量不改变算法
       └─ if dp:
            ├─ prepare_norm_counts()    # Cohort × 排序状态的 P/N 计数
            ├─ build_auc_interval_cost()# Raw AUC 区间代价，含一致性校验
            ├─ solve_norm_dp()          # 精确求解加权均值目标
            ├─ refine_norm_median()     # 多初值 + 遗留中位数目标精修
            └─ finalize_norm_candidates() # 恢复 Info、Raw、Summary
```

其中 `prepare_norm_rank()` 必须保证与原来 `norm_rank_marker` 的全部四种模式一致；排序完全相等时保留稳定的确定性次序；`Info` 映射分组也必须继续输出旧代码的 `0..K-1` 标签和每组状态列表。

为保持 API 兼容，建议返回：

```r
list(
  Raw = ...,      # 与旧流程相同的列和 Cohort/Info 语义；DP 只含最终候选
  Summary = ...,  # 至少含 Info、medianROCAUC；最佳候选位于首行
  SearchDiagnostics = list(
    method = "dp",
    objective = "legacy_median_auc",
    surrogate = "weighted_mean_raw_auc",
    exact_global_optimum = FALSE,
    M = ..., K = ...,
    n_dp_starts = ..., n_unique_candidates = ...,
    selected_cutpoints = ...,
    elapsed_seconds = ...
  )
)
```

旧 `bruteforce` 分支可完全保留 `Raw`、`Summary` 两个元素，或仅在向后兼容测试允许时添加无侵入 diagnostics。**DP 不应伪造所有从未枚举过的候选：**它的 `Raw` 只含已实际检查的最终候选，`Summary` 只对这些候选排序，文档要明确与暴力搜索的“完整候选表”不同。

为了控制复杂度，第一版允许只返回 DP 多起点及精修后的独立候选（例如前若干个），在结果中记录候选数。不应在 `K >= 5` 时继续创建巨型 `combn()` 或保存全体候选的字符串 `Info`。

## 必做的正确性与科学一致性测试

### 源码等价性测试：不影响原算法

- **不传参数**与显式 `norm_search_method="bruteforce"`：在相同输入下 `Raw` / `Summary`、最佳 `Info`、AUC 和组号完全一致（允许纯运行时/日志元数据差异）。
- 原 2/3/4 级运行在固定种子、相同队列、相同 `norm_rank_marker`、`Raw`/`OptimizeRank` 下，结果可与修改前快照比较。
- 检查 `n_norm_subtype=NULL`、无效 K、只有一个结局的队列、缺失标签、分型样本很少、RR/NRR 排序并列等情况。
- 对 `numCores`、并行后端关闭与 Windows 环境测试；既有 `stopImplicitCluster()` 异常路径要确保不残留 worker。

### 数学正确性：小数据全枚举验证

构造 M=6–12 个排序状态、C=3–5 个队列、K=2–5 的合成数据，并显式构造：

- 某些状态有很少的 responder，某些状态全 responder / 全 non-responder；
- 局部状态排序反转，区间损失可以为负；
- 队列间排序冲突；
- 队列不存在部分状态、`All` 与 cohort-only 中位数排序不一致；
- 分组数量等于状态数；
- 各分组出现并列 ROC 分数。

判定门槛：

1. **DP 的加权平均 Raw AUC 代理目标**必须与同一代理目标的穷举最优值在浮点容差内相等（这里可以严格证明/检验全局最优）。
2. 前缀计数计算的各候选 Raw AUC 必须与现有 `get_perform_markers()` 在定义对齐后数值一致（建议误差 `<=1e-10`，如 pROC 自动方向导致差异必须先解决）。
3. **DP + 精修的遗留中位数目标**应与原始穷举 `medianROCAUC` 比较差距；它不一定完全相等，必须报告 `gap = exhaustive_best - dp_best`，不得把近似写成精确。
4. 基于固定病例生成多种随机合成测试，输出匹配率、最优目标差值、精修提升及耗时，用于选择合理的初值和局部搜索策略。

### 真正的跨队列性能评价

算法开发、排序、切点与超参数选择只使用原有 20 个 ICI 推导队列。在推导数据上以**队列为单位**做留出/嵌套评估：每次只在训练队列估计状态排序与分组，然后在未参加排序的队列上计算**冻结映射、固定顺序** Raw AUC；不要使用单队列评价标签进行 `OptimizeRank` 重排。

在未参与参数选择的独立数据上对 4、5、6、8、10、12 级做比较，记录每队列 AUC 及其差异，而不是仅报告最好的合并 AUC。当前 Kim2018、Zachary2024、Binghua2024 和 QL1706 已被用来观察四级结果，因此新的探索应明确为后续方法优化，不能反过来把这些队列当成全新的、未被观察过的确认性数据；如果新模型要成为论文正式主版本，最好增加真正新的外部队列。

必须同时输出每级样本量/应答者数量、出现了多少空等级、稀疏等级和 K 增加时的校准及稳定性；不能只因为开发集 AUC 增高就断言更适合临床使用。

## 性能基准与验收门槛

基准表建议保存为 CSV/Markdown：`M`、`C`、`K`、`norm_rank_marker`、`norm_roc_method`、搜索方法、候选数、实际运行时间、峰值 RAM、最佳遗留中位数 AUC、最佳 cohort-only AUC、与穷举最优值的差、映射哈希。

- **正确性门槛**：小规模的均值 DP 最优值与均值穷举严格一致；遗留中位数代理差距显式可见；原 `bruteforce` 快照不变。
- **质量目标（可调，不是数学保证）**：在真实可完成穷举的 K=2–4 上，DP+精修的最佳 legacy median AUC 与暴力最优差尽量达到 `<=0.005`，记录每一组未达标情形；不允许把该工程阈值解释为统计显著性。
- **性能目标**：在有代表性的固定硬件上，K=4 的 DP 应显著快于穷举，优先以 `>=10x` 作为工程目标；K=5–12 不应预先生成全组合矩阵，峰值内存应随 `C*M^2 + K*M` 规模增长。具体耗时由基准实测，不能预先承诺。
- **可解释性**：同一个 `df`、`K`、排名模式和算法参数必须产生相同的 `Info`、切点与最优值；`SearchDiagnostics` 明确说明是否有全局最优保证。

## 实施步骤与修改顺序

### 阶段 A：冻结旧行为，建立基线

1. 新增小型独立场景测试，例如 `test/norm-search/01-test-legacy-contract.R`，记录 `subtypePerformance()` / `subtypeNorm()` 在 K=2–4 的旧结果快照。
2. 提取排序与候选评估中能够共享的代码前，先验证原始默认调用一致性。
3. 建立足够小、可全枚举的合成数据夹具，不包含真实患者隐私数据。

**完成标准**：旧版 `Raw`、`Summary` 的关键列及最优映射有可比较的基线。

### 阶段 B：新增算法开关，原方法保持原样

1. 给 `subtypePerformance()`、`subtypeNorm()` 增加末尾参数 `norm_search_method`，文档默认 `bruteforce`。
2. 将既有 `combn` + `foreach` 逻辑封装为内部暴力分支，勿在此阶段重写评分目标。
3. 检查 `Repeat` 中保存实际算法选择，输出日志包含所用算法。

**完成标准**：所有旧调用不变，新参数只决定进入哪个分支。

### 阶段 C：实现 DP 的精确均值子问题

1. 复用旧排序，建立 P/N 计数、前缀和、区间损失；分开有效 cohort 和 `All`。
2. 实现 O(KM²) 的 DP 与 backpointer，记录切点和精确均值目标。
3. 用小数据穷举验证代理目标精确性和 AUC 口径。

**完成标准**：仅对“固定排序 + 加权平均 Raw AUC”声称精确；运行不依赖 `combn()`。

### 阶段 D：近似求解遗留中位数目标

1. 基于少数确定性权重生成多 DP 初值。
2. 对候选进行单边界全位置局部精修，使用**相同的 legacy median 评分定义**。
3. 对获胜候选调用原性能函数复核；生成与既有 `Info`、`Raw`、`Summary` 兼容的输出。
4. 将代理目标、遗留目标、候选数、收敛原因、近似性质记录在 `SearchDiagnostics`。

**完成标准**：在 K≤4 的可穷举场景，质量与性能均有量化结果；K≥5 能运行。

### 阶段 E：验收与文档

1. 跑合成测试、旧版回归测试、`devtools::document()`、`R CMD build`、`R CMD check` 和相关场景脚本；按 `AGENTS.md` 在 Windows 使用项目指定的 R 环境。
2. 更新 `man/subtypePerformance.Rd`（由 roxygen2 生成）、`CHANGELOG.md` 的 `[Unreleased]`、必要的 README 或用例，不擅自更改版本号。
3. 若在 CCS 仓库实际实施代码修改，按项目 `AGENTS.md` 使用 BAC 正规工具记录贡献和验证，而不是直接改 `.bac` 文件。
4. 用推导队列开展 K=2–12 的跨队列敏感性实验；**这属于算法实现后的科学分析，不是本计划已经完成的结果**。

**完成标准**：旧算法保留可复现，新算法具有明确质量边界，并可交付实际研究使用。

## 风险与不做的事情

- **风险：DP 与暴力最优不一致。** 原因是原目标为中位数，而非可加均值。措施是明确区分目标、用多初值与精修降低差距、公开实际 gap。
- **风险：`OptimizeRank` 评价偏乐观。** 原因是用评价标签重新确定排序。措施是保留历史兼容模式，但把冻结分级顺序的 Raw AUC 作为独立泛化评估主要指标。
- **风险：K 增大导致小等级不稳定。** 措施是同时审计分型规模、不同队列覆盖与应答者数；先不在核心 DP 中引入未被批准的复杂样本数约束，若要加入须另行明确统计定义。
- **风险：默认参数或返回结构改变旧脚本。** 措施是新增参数附加末尾、默认暴力、仅 DP 模式添加诊断字段，并跑旧版快照测试。
- **风险：单纯用更高外部 AUC 反向选择 K。** 措施是所有 K 的选择在推导集内部完成，独立验证集仅作评估，避免多重试验后的选择性报告。

本次**不做**：修改 32 基因面板、改变 CCS/d1 构建、实现任意无序状态分配、引入复杂黑箱优化依赖、替换既有手稿的四级 normCCS 映射、直接声称新 K 的临床效用或 AUC 已经提升。

## 给 Codex 的简明执行指令

> 仔细阅读 `R/phenotype.R` 内 `subtypePerformance()` 与 `subtypeNorm()`，以及 `R/ccs_base.R` 内 `get_perform_markers()`、`multiXClass_roc()`。严格遵循本计划分阶段实施。公共接口只新增 `norm_search_method = c("bruteforce", "dp")`，默认行为必须与旧版完全一致。DP 必须利用原 RR/NRR 排序的连续分段结构，先精确优化队列加权平均 Raw AUC，再用多起点与局部精修逼近原有包括 `All` 的 medianROCAUC；不允许把代理目标误称为原目标的精确全局最优。所有高 K 路径不得调用 `combn()` 物化所有组合。先用可穷举的小数据严密验证 AUC 口径、切点与目标差值，再进行真实数据 K=5–12 实验；完成文档、回归测试及包检查，不修改版本号或发表结论。

## 相关源码链接

- [`R/phenotype.R`](https://github.com/huangwb8/CCS/blob/master/R/phenotype.R)：`subtypePerformance()` / `subtypeNorm()`。
- [`R/ccs_base.R`](https://github.com/huangwb8/CCS/blob/master/R/ccs_base.R)：`multiXClass_roc()` / `get_perform_markers()`。
- [`man/subtypePerformance.Rd`](https://github.com/huangwb8/CCS/blob/master/man/subtypePerformance.Rd)：现有 API 文档。
- [`CCS/paper/extraTex/back/supplementary-methods.tex`](https://github.com/huangwb8/Manuscripts/blob/main/CCS/paper/extraTex/back/supplementary-methods.tex)：论文中现有 normCCS 推导与独立验证边界。
