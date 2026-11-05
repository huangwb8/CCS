# normCCS K≤8 高性能精确组合搜索实施计划

> 适用项目：[`huangwb8/CCS`](https://github.com/huangwb8/CCS)  
> 计划日期：2026-11-04  
> 计划性质：可直接交给 Codex 在**本地工作区**分阶段实施的工程与验证规范；**本文不表示代码已修改或性能已实测**。  
> 主要目标：在固定 metaCCS 状态库与固定排序下，为 **2≤K≤8** 构建高性能、尽可能无约束、可核验最优性的 normCCS 搜索器。典型规模 M≈107；对 M=108、缺失状态和不同队列数亦应正确工作。

## 一句话结论

**先以 R/C++ 实现与原始评分严格一致的流式穷举内核，争取 K≤5–6 的无约束全局精确最优；再通过可证明有效的分支定界（B&B）和混合整数规划（MILP）挑战 K=7–8；将“保留 K4 切点的精确细分”作为有价值但有约束的对照，而非冒充无约束最优。**

搜索器负责解决“训练目标是否确实找到了最优组合”，独立外部验证负责解决“最优训练组合是否有临床可迁移价值”。两者严格分开。

## 范围、非目标与不可突破的边界

- 只针对已有 metaCCS → normCCS 的连续分段搜索，K 为 2–8 的整数，实际 M 从训练输入统计，不硬编码 107。
- 保持 d1、metaCCS、32 基因模型、上游模型、样本资格和历史论文 K4 映射不变；当前任务不负责重新学习 metaCCS 或重写外部验证。
- 不把既有 DP、连续继承模式或论文分析结果删掉；先做**旁路新引擎 + 完整基线回归**，验证后才考虑接入公开接口。
- 绝不依据外部 ICI 队列的反应标签优化切点、排序、评分方向、K 或候选选择；历史外部队列已用于方法探索时，结论标记为探索性，不能重称独立确认。
- 对每一个结果精确标记：`optimal_proven`、`feasible_not_proven`、`constrained_optimal` 或 `heuristic`，绝不因为使用了“exact”命名就自动宣称全局最优。
- 不承诺 K=8 必然能在固定时限内无约束证最优；必须报告已探索节点数、当前最优解、可信上界、最优性差距及时间/内存限制。
- K=8 是本项目的**预设实用上限**；不是“八级以上绝对没有临床价值”的经验或理论结论。

## 源码现状与必须先审计的事实

以 2026-11-04 核对的 GitHub `master` 版（`DESCRIPTION` 为 0.8.7）为起点；Codex **必须优先检查本地工作区最新代码、尚未提交的修改及实际测试目录**，若本地不同，以本地真实状态为准，记录差异，不覆盖现有工作。

| 位置 | 当前观察与实施关系 |
|---|---|
| `R/phenotype.R` → `subtypePerformance()` | 对外参数含 `n_norm_subtype`、`norm_roc_method`、`norm_rank_marker`、`norm_search_method`、`norm_mode`、`norm_objective`、`norm_top_n`；不要破坏位置参数与默认返回契约 |
| `R/phenotype.R` → `subtypeNorm()` | 离散穷举用 `choose()` 估计数量，超过 100 万候选或结果内存预算则停止；用 `combn()` 整体物化切点，再用 `foreach`/`ddply` 对每个方案计算 AUC，存在大量重复开销 |
| `R/phenotype.R` → `norm_prepare_counts()` | 已预计算各队列各排序状态的 responder/non-responder 计数、区间损失、前缀和及含 `All` 的 pooled 统计；是新引擎应优先复用的数学基础 |
| `R/phenotype.R` → `norm_auc_batch()` | 已支持快速批量评分 `Raw`、`OptimizeRank`，并对 ROC 自动方向、类别重排、部分并列行为做兼容处理；应作为**首选 R 快速参考实现**，但不能不经核对就假设 C++ 逐位一致 |
| `R/phenotype.R` → `subtype_norm_dp()` | 多 DP 初值 + 单切点/相邻双切点精修；DP 仅对其固定方向、加权均值 Raw AUC 代理子问题可证明精确，并不保证公开目标全局最优 |
| `R/phenotype.R` → `subtypeNorm_continuous()` | 已有保留父代切点逐层细分模式；它是**受继承限制的搜索**，有限 beam 不保证无约束最优 |
| `R/ccs_base.R` → `get_perform_markers()` | K=2、K>2 的 `Raw`/`OptimizeRank` 可能分别走不同底层实现；这是所有 fast scorer 必须精确复核的真值入口 |
| `AGENTS.md` | 约定 R 4.3.1/Windows、`targets`/`crew`、外置缓存、测试位置、包检查、BAC 贡献审计与禁止未经授权升版 |
| `docs/plans/2026-11-01-normCCS-paper-goal-refactor.md` | 现行论文分析的目标与训练/外部隔离约定；历史“参数优化”计划已注明退役，不应混用为当前执行规范 |

**第一步应保存当前源码指纹、当前四级参考切点、排序、全部必要的逐队列 AUC、当前 Raw/Summary/Repeat 结构和时间/内存基线。** 若本地大数据仅可通过正式 targets 读取，应在原正式 DAG 内读取；不要另写替代运行器生成正式结果。

## 搜索空间与为什么值得专用优化

固定 M 个已排序 metaCCS 状态，K 个非空连续等级由 K−1 个严格递增切点决定：

$$
1\le t_1<\cdots<t_{K-1}\le M-1,\qquad
N(M,K)=\binom{M-1}{K-1}.
$$

以 M=107 为例，组合数如下（精确整数，非运行时间）：

| K | 无约束切分数 | 首选处理 |
|---:|---:|---|
| 2 | 106 | 传统穷举与新引擎逐项等价 |
| 3 | 5,565 | 两者全量等价 |
| 4 | 192,920 | 全量等价 + 实际速度基准 |
| 5 | 4,967,690 | C++ 流式全量精确穷举 |
| 6 | 101,340,876 | C++ 并行穷举 + B&B 试验 |
| 7 | 1,705,904,746 | 优先 B&B / MILP，留流式穷举可恢复路径 |
| 8 | 24,370,067,800 | 重点攻克剪枝/优化器；不预物化全体结果 |

受限对照：如果 K8 必须保留固定 K4 的三个切点，则只需增加四个切点，候选数为

$$
\binom{106-3}{7-3}=\binom{103}{4}=4,421,275,
$$

比无约束 K8 少约 5,512 倍。**但是该结果只能证明“给定 K4 三切点约束下”最优**。若强制继承，最优解可能在原 K4 切点之外，不可写 `optimal_proven_unconstrained=TRUE`。

## 先定义“精确”：评分语义合同

“搜索算法精确”至少有两层含义，必须分别证明：

1. **评分精确**：对一个固定的切分，快速计算器与当前 `get_perform_markers()`/`norm_objective_score()` 的逐队列和汇总得分一致（约定浮点容差与并列规则）。
2. **优化精确**：在当前定义的合法分段全集中，能证明返回的候选得分最高；只有穷尽所有合法方案、或使用被证明有效的上界安全剪枝且已完成证明，或 MILP 得到且复核通过证书/求解状态，才允许宣称。

### 必须显式识别评分模式

| 搜索评分 | 搜索是否能直接写为区间可加 | 精确性要求 |
|---|---|---|
| 固定方向 `Raw` AUC（推荐作为可证明优化的主研究问题） | **每队列 AUC 可加** | 适合区间代价、均值 DP、median MILP、B&B |
| 当前 legacy `Raw`（自动方向） | 区间可加的基础 AUC + **按候选决定的方向规则** | 需要精确复现方向选择；不能默默用固定方向替代 |
| 当前 `OptimizeRank` | **不可直接视为固定等级方向的区间可加目标** | 采用原语义 fast scorer 的完整穷举，或先证明新上界/数学模型正确，才能叫 exact |
| `median`（含 pooled `All`） | 汇总目标非可加 | 不能将均值 DP 的最优作为 median 的全局最优 |
| `cohort_median`、`mean`、`lower_tail`、`mean_lower_tail` | 需分别建模 | 固定方向下部分可 MILP/DP 精确求解；若方向/重排改变则重新证明 |

遗留 `median` 具体是 `median(c(AUC_All, AUC_1, …, AUC_C), na.rm=TRUE)`；`All` 和真实队列有重叠患者，不能当作额外独立验证队列，也不能误以为 pooled AUC 等于各队列 AUC 的算术平均。真实训练队列若有 20 个，含 `All` 就是 21 项，median 是第 11 大值；如果仅做真实 20 队列中位数，则是第 10 与 11 大值的平均。

**排序规则、各队列响应因子水平、ROC 方向、K=2 特例、稀疏状态、类别不存在、Tie、NA、`All` 的统计口径都是科学算法合同的一部分。**

在不改变旧接口默认语义的前提下，新引擎建议显式记录：

```text
rank_marker       = RR | NRR | NRR+RR | RR+NRR | ...已有合法值
roc_semantics      = legacy_raw_auto | legacy_optimize_rank | fixed_raw
objective          = median | cohort_median | mean | lower_tail | mean_lower_tail
search_backend     = stream_exact | bnb_exact | milp_exact | dp_baseline
constraint_mode    = unrestricted | inherited_k4
```

注意：`fixed_raw` 是需要独立定义、独立验证的新增统计语义；它不等同于简单地把旧 `norm_roc_method = 'Raw'` 改个名字。

## 数学内核：一次聚合，反复评分

固定从**高分/高应答倾向到低分**的排序状态 $1,\dots,M$；实际需要先核对现有 `norm_rank_subtypes()` 的顺序与 `norm_info()` 的标签方向，避免编码反向。

对 cohort $c$、状态 $i$，预统计 responder 数 $p_{ci}$ 和 non-responder 数 $n_{ci}$，总数分别为 $P_c,N_c$，均必须非零。对 `All` 独立合并患者计数。

在高分优先方向下，以所有 M 个状态保持不同等级时的 AUC 为：

$$
A_c^{\mathrm{fine}}=
\frac{1}{P_c N_c}\sum_{i=1}^{M}p_{ci}
\left(\sum_{j>i}n_{cj}+\frac{n_{ci}}2\right).
$$

将区间 $[a,b]$ 合并成同一等级，使该区间内正负样本对由按状态排序变为并列；固定方向下损失为

$$
L_c(a,b)=\frac{1}{2P_c N_c}
\sum_{a\le i<j\le b}
\left(p_{ci}n_{cj}-n_{ci}p_{cj}\right).
$$

注意 $L_c(a,b)$ **可以为负值**。排序局部违反真实应答率时，合并反而可能提高 AUC；实现中绝不能把负损失截断为零，也不能假设增加 K 必然改善 AUC。

任意 K 个连续非空区间构成的切分 $\mathcal P$ 满足：

$$
\boxed{A_c(\mathcal P)=A_c^{\mathrm{fine}}-\sum_{[a,b]\in\mathcal P}L_c(a,b)}
$$

预计算 `loss[c,a,b]`、前缀计数、必要的方向/OptimizeRank 辅助信息后，单个候选只需处理约 `C × K` 规模数据，而不是重建全部患者数据帧。严格复用当前 `norm_prepare_counts()` 及 `norm_auc_batch()` 的定义；先用 R 做 golden reference，再把**已证明一致的热点循环**移到 C++。

## 算法路线 A：C++ 流式全穷举（第一优先级）

### 核心思路

**搜索空间完全不缩减，只减少每个候选的固定开销和内存。** 在固定 M、K、排序、评分函数与目标下，扫描全部 $\binom{M-1}{K-1}$ 个切点组合，精确计分，仅保留最优与 Top N；这是评价其他高 K 搜索算法的黄金基准。

### 必须实现的关键技术

- 使用**字典序组合迭代器**，边生成边打分；禁止 `combn()` 全量预分配、逐候选 `data.frame`/`group_by`/`paste(Info)`、逐候选 `pROC`。
- R 一次性预聚合为连续内存数组，C++ 使用预计算 `loss`、`base_auc`、计数及方向辅助量。
- `K−1≤7`，可把切点放在紧凑固定数组/小向量；评分阶段不进行大对象分配。
- 每批例如 4,096–65,536 个方案，可按实际性能调优；先提供单核可靠基线，再选择最优批量。
- 用 Top-N 小顶堆或定长有序容器，只保存少数切点与评分；完整 `Raw` 全组合表仅在 K≤4 的兼容验证模式中可选输出。
- 组合计数、已访问数量、任务偏移与评估次数用至少 **64 位整数**；不得用 R 32 位 integer 存储 243.7 亿计数。
- 排序 `Info` 字符串和传统性能表只对 Top N 回填；最终对每个返回方案跑原 `get_perform_markers()` 逐队列复核。
- 暴力模式中新实现的输出规模与旧模式不同，应明示 `result_policy='top_n_only'`，不能伪造所有未存储方案的 `Raw` 或 `Summary`。

### 迭代器与确定性

```text
输入 M,K、固定排序、评分语义、预聚合数据、top_n
初始化 cuts = [1,2,...,K-1]
while cuts 有效:
    计算每队列 AUC（包括 All）
    计算目标评分
    按固定比较器维护 Top N
    访问量 +1
    使用组合后继算法找到下一个 cuts
完成后：验证 visited == choose(M-1,K-1)
         对 Top N 重新调用 R 原版 ROC 与 objective 评分复核
         exact_global_optimum = TRUE（仅当无约束且全集覆盖）
```

Top N 同分时先执行现行目标分支的并列排序逻辑（旧默认 vs 新目标可能不同），再固定切点字典序；不得让线程完成顺序影响获胜者。计算容差只能用于正确性判断；若引入“近似并列”决策，必须定义并测试它是否改变原函数的最优语义。

### 原生代码选择

优先 **C++17 + Rcpp**：本 R 包已是 R 生态，适合新增 `src/norm_exact.cpp`、必要的 `src/Makevars.win`/`src/Makevars`、`R/RcppExports.R`，R 层适配函数统一追加到 `R/phenotype.R` 底部（按 2026-11-05 的源码整理要求，不再单独托管）。先核对当前 Windows R 4.3.1 编译工具链和包依赖，按 R 包标准把 `Rcpp` 加入 `Imports`、`LinkingTo` 并运行 `Rcpp::compileAttributes()`；不得把 Rcpp 当作已经存在的仓库依赖。若本地编译环境或构建约束有明确障碍，可评估纯 C++ `.Call`；**不建议第一版同时维护 C++ 与 Rust 两套等价内核**。

## 算法路线 B：可证明安全的 Branch-and-Bound

### 目标

在**不漏掉任何可能优于当前解的方案**的前提下，尽早淘汰不可能产生最优结果的整棵切点搜索子树。与 DP 代理不同，B&B 必须用**原目标**选优，并用数学上成立的乐观上界剪枝。

### 最有希望的严格上界：每队列乐观完成值

把搜索节点定义为：已经完成前 $b$ 个状态的 $g$ 个分段，尚需用 $r=K-g$ 段覆盖后缀 $[b+1,M]$。设固定前缀已经造成 cohort $c$ 的损失 $D_c^{\mathrm{prefix}}$。用区间 DP 计算该 cohort 后缀在恰好 $r$ 段情况下**可能达到的最小总损失** $D_{c,\min}^{\mathrm{suffix}}(b,r)$，则

$$
U_c=A_c^{\mathrm{fine}}-D_c^{\mathrm{prefix}}-
D_{c,\min}^{\mathrm{suffix}}(b,r)
$$

是任何该节点合法完成方案的**固定方向 AUC 上界**。每个 cohort 的最乐观后缀方案可以不同，所得到的向量虽然可能无法同时实现，但仍是有效的逐坐标上界。

中位数、队列均值、下尾均值均对每个 AUC 单调不减，因此将 $U_c$ 代入同一个目标函数得到的值，是该搜索节点的**合法乐观上界**。对含 `All` 的目标，把 pooled `All` 作为独立的评分坐标计算上界，但不要称它独立队列。

更紧的实现可以加入共享切点可行域约束，先不必实现。利用当前 DP 获得高质量初始解（incumbent）能够更早剪枝，**DP 负责找到初值，不负责提供未经证明的剪枝结论**。

对于 `legacy_raw_auto`，如果评分确实只可能取正向 $A_c$ 或反向 $1-A_c$，需要同时计算可行区间 $[L_c,U_c]$，并使用 $\max(U_c,1-L_c)$ 作为不依赖方向判定的乐观上界，再核对实际 pROC 规则。不可直接把固定方向 $U_c$ 当作自动方向分数上界。

对于 `OptimizeRank`，**不能直接复用上述固定方向上界**，因为不同候选和队列会改变类别顺序。未获得严格证明的上界时，只能退回 $1$ 这一安全但通常无用的上界，或对候选做完整穷举；禁止“貌似合理”的启发式剪枝却标记全局精确。

### B&B 实施与证明要求

- 先在固定方向 `Raw` 上实现深度优先或 best-bound 搜索。
- 候选前缀可行性：剩余 $M-b$ 个状态必须至少容纳 $K-g$ 个非空区间。
- 有可行初解、可信的 upper bound；只在 `bound < incumbent`（含明示的浮点安全裕量）时剪枝。需要枚举全部并列 Top N 时剪枝条件相应调整。
- 完成时给出 `best_score`、`global_upper_bound`、`gap`、`nodes_expanded`、`nodes_pruned`、`n_candidates_evaluated`、`termination_status`。
- 用 M≤12 的可穷举随机测试逐节点检查：**任何被剪枝节点的所有合法后代分数都不得超过记录的上界**；与全穷举结果严格一致。
- 数值稳定性：上界使用保守的向上误差余量，不能因为 1e−12 的浮点误差误杀真正最优解。对边界与同分情况单独测试。

## 算法路线 C：MILP / 数学规划（K=7–8 的重点攻关）

### 将分段编码为 DAG 上的 K 条弧

构建节点 $0,1,\ldots,M$，对每个 $0\le i<j\le M$ 建有向弧 $(i,j)$，表示把状态 $i+1,\ldots,j$ 合并成一个等级。

设 $x_{ij}\in\{0,1\}$ 表示是否选弧。约束从 0 到 M 有一条完整连续路径：

- 源点净出流为 1；终点净入流为 1；其他节点流量守恒。
- $\sum_{i<j}x_{ij}=K$，因此恰好 K 个非空等级。
- 图有向无环，因此上述单位流整数解对应完整连续切分，不允许乱序或跳漏状态。

M=107 时弧数为 $107\times108/2=5,778$，与枚举 243.7 亿组合相比，建模规模具有探索价值，但**变量少不等于求解一定快**。

每个 cohort 的固定方向 AUC 变为线性表达式：

$$
A_c(x)=A_c^{\mathrm{fine}}-
\sum_{0\le i<j\le M}L_c(i+1,j)x_{ij}.
$$

### 中位数的严格 MILP 表示

假设 $n$ 个纳入目标的 AUC 都在 $[0,1]$。

**奇数 $n=2h+1$**：中位数为降序第 $h+1$ 个值。引入连续变量 $t\in[0,1]$、二进制 $z_c$：

$$
\max t,\quad
\sum_c z_c\ge h+1,\quad
 t\le A_c(x)+1\cdot(1-z_c)\quad \forall c.
$$

当至少 $h+1$ 个 $z_c=1$ 时，阈值 $t$ 不得高于这些 cohort 的 AUC；最大化 $t$ 即得到精确的中位数。这里 1 是基于 AUC∈[0,1] 的有效 big-M；可以通过预先计算更紧界改善求解。

**偶数 $n=2h$**：降序第 $h$ 与 $h+1$ 个值的均值才是常规定义的 median。定义 $t_h,t_{h+1}$，分别要求至少 $h$ 与 $h+1$ 个 AUC 达到对应阈值，最大化 $(t_h+t_{h+1})/2$，不能错误地只优化其中一个顺序统计量。

对于 `median`（含 All）和 `cohort_median`（不含 All），必须按各自的 n 创建模型。对 `mean` 直接线性优化；对 lower-tail 可研究如下辅助变量表示（$\alpha=0.25$）：

$$
\max_{x,\eta,u}\quad \eta-\frac{1}{\alpha C}\sum_c u_c,
\qquad u_c\ge \eta-A_c(x),\quad u_c\ge0.
$$

`mean_lower_tail` 再加上 $0.5\,\mathrm{mean}(A_c)$ 与 $0.5$ 倍上述下尾项。必须先与当前 `norm_objective_score()` 的 C=1,2,3,4,5,8 边界和分数精确对齐。

**MILP 适用性警告：以上模型只对固定方向、可加 AUC 严格成立。** 现有 `legacy_raw_auto` 的方向翻转及 `OptimizeRank` 的结局驱动重排不能直接纳入同一套线性模型，除非另行建立完整等价的逻辑约束并通过小规模穷举验证。否则结果只能叫 `milp_fixed_raw` 全局最优，不能标成 `legacy_optimize_rank` 最优。

### 求解器与实现建议

- 首选单独实验适配器（在本地许可、安装与版本确认后）试用 **SCIP**，有条件再测试 **Gurobi**；不在 CCS R 包核心强制绑定某个商用求解器。
- 初版优先由 R 生成稀疏 MIP 数据，导出 LP/MPS 到外置缓存目录，用本地 solver 可执行文件或官方接口求解；解析求解状态、primal 候选、dual bound、gap、求解日志。
- 将 DP/stream/继承式搜索获得的候选转为 MIP feasible start（solver 支持时），并与无 warm start 的求解比较。
- 对固定方向 median 的 small-M 案例，必须与 C++ 全穷举最优值、切点和 AUC **完全核对**。
- 维护多个优秀解时，对完整路径上的 K 条选中弧增加 no-good cut：$\sum_{e\in \text{selected path}}x_e\le K-1$，在同一模型上反复求次优；注意只有每轮再次完成最优性证明才叫精确 Top N。
- solver 在时间限制处退出时保存 `best_feasible`、`best_bound`、`gap`、`status=time_limit`；不写“证明全局最优”。

参考资料：
- SCIP 官方：<https://www.scipopt.org/>
- SCIP 接口与文档：<https://scipopt.org/doc/html/>
- Gurobi R 接口：<https://docs.gurobi.com/projects/optimizer/en/current/reference/r/overview.html>
- Rcpp package attributes：<https://search.r-project.org/CRAN/refmans/Rcpp/html/compileAttributes.html>

## 可选补充：受约束 K4→K8 精确细分

在论文 K4 切点固定的前提下，新增 4 个切点做 442 万个方案的完整穷举，作为：

1. K=8 能够实际完成的高质量临床基线。
2. 新评分内核在 K8 规模下的端到端性能基准。
3. B&B、MILP 的可行初解来源。
4. 检查继承分级是否比自由重新分界更稳定的分析对象。

其 `constraint_mode='inherited_k4'`，`optimality_scope='within_constraint'`；即使真的遍历完，也只声明约束最优，不得与自由 K8 最优混为一谈。将来放宽 1/2/3 个固定切点，只能视为不同受限子空间，不能说已经覆盖全集。

## 工程结构：R 负责管线，C++ 负责热点，优化器负责证明

建议**先在现有代码附近新增模块，不大量移动旧 `phenotype.R`**：

```text
R/phenotype.R（底部辅助函数）     # 参数校验、预统计适配、策略分发、返回转换
R/norm_exact_objectives.R       # 评分契约和目标一致性测试辅助
src/norm_exact.cpp              # 组合迭代、原生评分、Top-N、可选 B&B
R/norm_exact_milp.R             # MIP 模型构造、solver 适配（可选依赖）
test/normCCS-5-more/tests/      # 小规模精确性、兼容与边界测试（依项目实际目录）
test/normCCS-exact/             # 新的实验场景/基准代码及说明（若本地尚无则新建）
docs/plans/                    # 本实施计划与补充记录
```

正式分析及所有大型 benchmarking 数据在 `D:/cache/ccs/normCCS-exact/`，遵守 `AGENTS.md` 的外置规则。若接入已有 `test/normCCS-optimization` 正式 DAG，则 **沿用它当前 `_targets.R` 和既有 targets store**，不要另建替代入口，也不要重算或覆盖不相关的旧分析。资源目录、临时目录和监控均按既有配置统一。

### 拟议 API（仅为实现指导，不必一步到位）

```r
# 新内部入口，避免在第一阶段改变公共函数默认行为
norm_exact_search <- function(
  df,
  k,
  rank_marker,
  roc_semantics = "legacy_raw_auto",
  objective = "median",
  backend = c("stream_exact", "bnb_exact", "milp_exact"),
  constraint_mode = c("unrestricted", "inherited_k4"),
  inherited_cuts = NULL,
  top_n = 20L,
  time_limit_sec = Inf,
  memory_limit_gib = 8,
  num_cores = 1L,
  verbose = TRUE
) { ... }
```

当 `backend='milp_exact'` 且请求 `legacy_optimize_rank` 等尚无精确模型的评分模式时，**明确报错**或明确标为不支持；不得静默换成固定方向代理目标并返回一个看似同目标的答案。

返回对象建议（可复用现有 `Raw`、`Summary` 与 `SearchDiagnostics`，但应明示仅 Top N）：

```r
list(
  Raw = reference_verified_top_candidates,
  Summary = top_candidates_summary,
  SearchDiagnostics = list(
    M = M, K = K,
    backend = "stream_exact",
    score_semantics = "legacy_raw_auto",
    objective = "median",
    constraint_mode = "unrestricted",
    optimality_scope = "unrestricted_global",
    status = "optimal_proven",
    optimality_proven = TRUE,
    score_best = best_score,
    score_upper_bound = best_score,
    absolute_gap = 0,
    relative_gap = 0,
    n_possible_candidates = possible_count_as_double_or_string,
    n_candidates_evaluated = visited_as_double_or_string,
    top_n = 20L,
    score_reference_tolerance = 1e-10,
    input_hash = "...", rank_hash = "...", code_hash = "...",
    elapsed_seconds = NA_real_, peak_memory_mb = NA_real_,
    solver_version = NA_character_
  )
)
```

若截断/超时，`optimality_proven=FALSE`，`status='time_limit'` 或 `memory_limit`，上界和 gap 不可求则填 `NA` 而不是 0；允许返回已经找到并复核的最佳候选。

## 并行、缓存与可恢复性

### CPU 并行的正确拆分

- 组合枚举可按**第一个切点的固定值**划分互不重叠的任务；每任务继续枚举剩余 K−2 个切点。第一切点 $t_1=s$ 的规模是 $\binom{M-1-s}{K-2}$；使用组合数估算做均衡调度，不要简单按切点数等分工作量。
- 新 C++ 内核先保证单线程/单分区确定性，随后选择进程级分区或原生线程；若正式分析已启用 `targets`/`crew`，优先由 `targets` 分支并行，单任务内原生线程数量严格受控，**避免嵌套并行和 BLAS/OpenMP 过度订阅**。
- 每个分区保存独立 Top N、访问量、完成标记、输入签名；最终通过确定性全局排序归并 Top N。
- 不要并发向同一个大 RDS/CSV 追加写；分区临时输出使用原子重命名和校验哈希。

### 中断续跑

- 独立“分区任务”是最可靠的续跑单位；未完成分区重算，已完成且哈希一致的分区复用。
- 由原正式 `_targets.R` 追踪有效计算及缓存，不发明第二套正式调度/恢复系统。
- `input_hash` 应覆盖状态排序、队列名单、聚合计数、响应因子水平、目标、ROC 模式、约束、切点锚点及源码版本；任何关键项改变都使旧结果失效。
- checkpoint 或日志不记录患者级数据，敏感原始矩阵和大型产物仍保留本地指定缓存与备份。

## 测试设计：没有这些测试不能称 exact

### 金标准等价测试（第一阶段必过）

针对 M=5–12、C=1–8 的合成数据，枚举合法 K=2–min(8,M)，包含：

- RR、NRR 与其组合排序，排序分数完全并列、状态编号非连续、某队列缺失部分状态。
- 不平衡队列、极低/极高应答率、正负样本稀疏但每个有效队列两类均存在。
- AUC 接近 0.5、方向可能改变、等级并列、负区间损失、$K=M$ 和 K2 特例。
- `All` 与各真实队列的 AUC 排序冲突；20 队列奇偶中位数、真实队列尾部 25% 权重边界。
- `Raw` 与 `OptimizeRank` 对同一切分的不同评价；不同因子水平、自动方向、浮点并列顺序。
- 原版 `get_perform_markers()` 与原生 fast scorer 对**所有候选**逐队列对照，`max(abs(auc_fast - auc_reference)) <= 1e-10`；不允许只查赢家。
- 以相同目标对比原版穷举 vs 新全穷举，最优值差 ≤1e−10，最佳候选/同分规则一致（若同最优不同切点，必须能解释是并列策略而不是数值错误）。

### 结构性证明测试

- 穷举组合访问量等于 $\binom{M-1}{K-1}$，没有重复、没有漏项、无非法切点，K≤8 且 M=107 的数量运算不溢出。
- B&B 的上界用独立小规模全枚举检查，验证每个被剪枝节点所有后代的分数都 ≤ 该节点上界；任何违反立即阻断提交。
- MILP 的每组返回切分满足恰好 K 个连续非空区间；各目标的求解分数与独立手算/小规模穷举一致；MIP `optimal` 不能直接替代 R 评分复核。
- 受限 K4 细分的最优值，与同约束全量枚举一致；但不能误标无约束最优。
- 时间限制/中断情形产生合法可复核候选且 `optimality_proven=FALSE`；完成时证明标志必须与实际求解状态一致。

### 旧行为与软件测试

- 不传新参数、旧 `bruteforce`、旧 `dp`、旧 `continuous` 的默认 `Raw`/`Summary`/`Repeat`、`Info` 字符串和参考 K4 结果保持回归一致。
- public API 位置参数与命名参数均正常；非法 K、缺失状态、非法 cores、无可用双分类队列、求解器未安装明确报错。
- Windows R 4.3.1 源码编译、`Rcpp::compileAttributes()`、`devtools::document()`、包构建/检查、动态库加载卸载都要通过。
- 使用 `test/normCCS-5-more/run-unit-tests.R`（如果在本地存在）运行单元测试，按本地 `AGENTS.md` 要求追加必要测试；公开 GitHub 的文件树不一定包含全部本地 `test/` 内容，不能凭不存在的文件名直接假定测试完成。
- 并行结果在 1/2/4 核及中断恢复后同分解、同切点、同评分；无遗留 worker；内存保护确实生效。

## 性能基准与阶段门禁

在同一硬件、相同输入、相同评分模式、相同目标、相同并行数及固定 BLAS 线程设置下记录：

`M, C, N_patients, K, candidate_count, evaluated_count, pruned_count, per_core_throughput, total_seconds, peak_RSS, preprocessing_seconds, scoring_seconds, merging_seconds, best_score, upper_bound, gap, optimality_status, best_cuts, input_hash, software_versions`。

### 重要的性能推算（仅数学量级，不是实测承诺）

M=107、K=8 需要约 243.7 亿个候选。如果整体有效吞吐仅有 100 万候选/秒，完整穷举仍需约 **6.77 小时**；1000 万候选/秒也需要约 **40.6 分钟**。真实吞吐取决于 cohort 数、目标函数、自动方向、CPU 核数、内存带宽及编译优化，必须先按小规模实测推断，而非承诺“C++ 一定数分钟跑完”。

### 阶段门禁

| 阶段 | 工作 | 必须通过才继续 |
|---|---|---|
| P0：审计/冻结 | 读取本地最新 `AGENTS.md`、目标函数、现有测试/正式 DAG；留存 K4 参考结果 | 评分合同、数据来源与现有行为可复现 |
| P1：预统计/快速评分 | 固化 R golden scorer；抽离/优化批量计数接口 | 所有小规模候选与原 ROC 全量等价 |
| P2：C++ 流式全穷举 | 从 K2–4 回归到 K5，再跑 K6 资源试验；Top-N、确定性、恢复分区 | K≤5 对应无约束全局最优，K6 是否可在预算内完成有实测结论 |
| P3：受限 K8 | 固定论文 K4 切点，442 万候选全部完成 | 有约束精确最优、速度基准和完整外部隔离 |
| P4：B&B | 实现经过验证的逐队列乐观上界，基于 fixed_raw 证明安全剪枝 | 每个小规模案例与全穷举一致，并记录 K6–8 节点压缩率 |
| P5：MILP | 路径整数规划 + 精确 median/order statistic 建模；接入 solver 状态和 bounds | 小规模固定方向最优与全穷举一致，K7/8 有实际求解日志 |
| P6：全局目标评估 | 比较 unrestrained K2–8 的 exact/未证最优与继承方案、DP 候选 | 每个 K 都有状态、分数与界限；无条件性或无证明结论不冒充最优 |
| P7：论文/外部评价 | 冻结训练所得映射，再运行现有正式评价管线 | 所有外部指标只是评估，不参与重新搜索；保留完整队列结果 |

如果 P4 或 P5 不能在合理资源内完成全局证明，仍应输出 **Best-so-far + Valid Upper Bound + Gap**，如未能获得有效上界则明确为未知，而非推迟交付或改称“精确”。

## 临床使用性和科学判读

- K 从 2 到 8 的训练 AUC 最优值不必单调：切分数固定且新增切点可改变得分，某些切分反而损失应答排序信息。
- 除训练 AUC 外，至少记录每级患者数、应答者数、非应答者数、跨队列等级覆盖、空等级及最低支持度。
- 外部队列统一使用训练确定的排序、映射和方向；`OptimizeRank` 若依靠目标队列结局重排类别，只能作为明确标记的描述性/事后评价，不能冒充临床可部署的前瞻性预测 AUC。
- 进行 cohort-level 留出或重新抽样时，**必须在每个训练折内重算排序与所有搜索**，不能先用全部队列定排序/切点后再计算留出 AUC。
- 明确展示与论文 K4 的比较：训练目标、独立外部表现、波动与不确定性，以及新增等级所换取的解释复杂度。最终推荐 K 并不必然是 8。
- 如果多个 K 的训练得分非常接近，优先关注模型稳定、等级支持与对临床决策的增量价值，而非只按训练最优 AUC 选择。

## 本次任务的改动权限与交付规范

1. **先交付 P0–P2，再开启下一阶段**。Codex 可以在用户授权的一次任务内持续推进各阶段，但每过一个门禁必须先把证据写入本地验证文档；失败时不跨过门禁编造完成状态。
2. 先给精确/近似搜索并存的新内部入口；**未经证明默认输出兼容，不替换旧 `norm_search_method="bruteforce"` 默认行为**。新入口成熟后，可提出附加 `stream_exact`/`bnb_exact`/`milp_exact` 选项的兼容改动，遵守位置参数规则。
3. 不擅自修改 `DESCRIPTION` 版本号、发布标签、GitHub Release、当前论文的固定四级映射、正式训练/外部样本、归档结果；涉及需要更新的依赖、NAMESPACE、roxygen 和 `CHANGELOG [Unreleased]` 按 `AGENTS.md` 正常维护。
4. 有正式 `targets` 的研究仅修改/执行其原正式 DAG。小型纯算法单元测试可以独立运行；它们不得冒充正式训练或验证分析。
5. 大文件写入 `D:/cache/ccs/normCCS-exact/` 或现有正式缓存根目录，不回退至坚果云同步目录；内存/线程/缓存/输入哈希/恢复边界必须可检查。
6. 依仓库 `AGENTS.md` 用 `bac` 工具记录、验证贡献链，不直接编辑 `docs/contribution.bac`。
7. 最终报告必须清楚区分：**已实施、已验证、仅设计、资源失败、理论上可行但尚未实现**，不能用“预计加速”替代 benchmark。

## Codex 最终应交付的文件

- 新搜索器的 R/C++ 源码、编译接入文件、必要的新依赖声明。
- 小数据全候选黄金测试、随机/对抗用例、旧代码回归测试、并行/中断测试。
- K2–8 可复现 benchmark 配置与每配置实际结果（包含不能完成、超时和未证最优的记录）。
- 一个 `README.md`：说明每种求解器**在哪些评分语义下**能证明全局最优；一张算法能力矩阵尤其重要。
- 一个 `verification.md`：完整命令、环境、编译器、R 包、样本来源、通过/失败项、源码哈希、最优性证据和测量时间。
- `CHANGELOG.md` 的 `[Unreleased]` 记录、必要的 roxygen/Rd 文档，以及经校验的 BAC 记录。
- 当新引擎在当前固定方向或遗留目标下不能证明无约束 K8 最优时，交付**实测负结果及明确下一步**，不得悄悄改用 DP 假装完成。

## 可直接粘贴到 Codex 的实施指令

```text
请在当前本地 CCS 仓库实施 docs/plans/ 下的《normCCS K≤8 高性能精确组合搜索实施计划》（我会把 Markdown 文件放进该目录）。

优先阅读 AGENTS.md、R/phenotype.R、R/ccs_base.R、当前 normCCS 的最新分析入口、全部相关测试与现有工作区改动。先做 P0 审计并保存基线；不要覆盖用户未提交的修改，不要修改论文固定 K4、metaCCS、训练与外部数据名单、版本号、正式 targets store。

按 P1→P2→P3→P4→P5 的先后顺序执行；第一优先级是复现 legacy ROC/AUC/objective 语义的快速计分内核和 C++ 流式精确穷举。小规模必须对每一个候选的逐队列 AUC 与原始 R 实现做独立核对，通过后才能继续 K5/K6。使用 Rcpp/C++17 优先，必要时按包规范增加依赖。所有大文件外置 D:/cache/ccs/，所有正式分析沿原 _targets.R、crew 和 tar_make 工作，不另造运行器。

B&B 的任何剪枝必须有严格合法上界；MILP 首先仅优化经过验证的固定方向 AUC，不准将求解的代理目标伪称旧 OptimizeRank 或 auto-direction 全局最优。每个配置需记录 score_semantics、objective、constraint_mode、best_score、best_bound、gap、visited、runtime、peak RAM、optimality_status。对于没完成证明的 K7/K8 返回 feasible_not_proven 和真实界限，绝不能写全局最优。

每阶段写验证证据与性能基准，遇到关键回归不通过必须修复后再推进。保持现有默认行为和 API 兼容；维护 CHANGELOG、文档与 BAC。最终提供已实施/未实施清单、验证命令、性能表、剩余风险、可复现的交付路径。不要夸大加速效果或外部泛化结果。
```

---

### 参考实现与设计依据

- CCS 当前源码：<https://github.com/huangwb8/CCS/blob/master/R/phenotype.R>
- ROC 实现：<https://github.com/huangwb8/CCS/blob/master/R/ccs_base.R>
- 仓库工程约束：<https://github.com/huangwb8/CCS/blob/master/AGENTS.md>
- 现行论文搜索方案：<https://github.com/huangwb8/CCS/blob/master/docs/plans/2026-11-01-normCCS-paper-goal-refactor.md>
- 历史 DP 搜索讨论：<https://github.com/huangwb8/CCS/blob/master/docs/plans/CCS_normCCS_DP_search_optimization_plan_2026-10-29.md>

> 本计划的数学公式说明的是**在对应评分定义下的计算等价与可能的最优性证明**，不是对 K=8 实际可解速度或外部预测增益的保证。先验证语义，再优化速度，最后评价泛化。
