# ablation-03 性能优化实施计划

## 通俗解释：究竟要解决什么

当前 ablation-03 像是让一名工人反复填写 800 张表：每张表内部可以用 8 个 CPU 线程，但 800 张表基本按顺序处理，RTX 4090 也没有参与核心模型训练。因此程序结果是可用的，等待时间却很长。

本计划只改变“怎样执行”，不改变“算什么”：训练数据、表示构造、抽样、折叠、lambda 候选、随机种子、指标、结果字段和统计解释都保持不变。优化后的程序首先应得到相同结果，再讨论速度提升。

## 专业判断：当前瓶颈和边界

- 主要瓶颈是 `R/ablation.R` 中的 learning-curve 三重循环（fraction × repeat × representation）以及每次 readout 的 3 lambda × 3 inner folds；当前配置约产生 800 次 `gblinear` XGBoost 拟合。
- 单个拟合通过 `nthread` 使用 CPU，当前 ablation-03 默认上限为 8；外层任务没有调度器，不能把独立拟合分配到其余 CPU 核心。
- 当前模型固定为 `gblinear + coord_descent`。XGBoost 的 GPU 加速主要针对树模型路径，不能把 `device="cuda"` 当作对现有模型的无风险替换；GPU 后端必须通过数值等价性和实际 GPU 利用率验收。
- Windows 下并行 worker 会复制大矩阵；并行度不能简单设为 32，必须同时约束 worker 数、每个 worker 的线程数和内存峰值。

## 目标与非目标

### 目标

- 保持当前公开结果契约和科学计算逻辑不变。
- 让 CPU 并行覆盖独立的 learning-curve、scaling 拟合任务，并避免线程过度订阅。
- 对 GPU 加速建立可验证的可选后端；只有在结果和性能都达标时才启用。
- 为高耗时节点增加可恢复的细粒度 checkpoint 和可观察进度。
- 在真实 ablation-03 输入上形成可重复的基线、加速后耗时和资源利用率报告。

### 非目标

- 不改用不同的分类器、不同的损失函数或 `gbtree` 来换取速度。
- 不改变样本过滤、cohort 抽样、inner-fold 划分、lambda 候选、seed、指标和报告解释。
- 不删除统计重复、null control、scaling 或 decoder；若要减少分析规模，另行审批参数变更。

## 改进方向

### 建立可重复的性能基线

先增加只读 benchmark 和 profiling 入口，记录每个阶段、每个拟合任务的 wall time、CPU time、峰值内存、线程数、cache hit/miss 和 GPU 利用率。基线必须覆盖：native geometry、retrieval、readout、learning-curve、cohort scaling、decoder。

用固定 seed 和当前缓存输入运行一个小型真实子集与一次完整配置，保存结果摘要和硬件信息。性能比较必须与同一输入、同一 cache 状态和同一 R/xgboost 版本对齐。

### 统一运行时配置，而不是把硬件策略散落在函数中

在表示分析配置中增加一个运行时执行上下文，至少包含：后端（`cpu`/`gpu`/`auto`）、总 worker 数、每 worker 线程数、GPU device、内存上限和 checkpoint 粒度。由一个运行时解析函数根据本机核数、可用内存和 GPU 能力生成最终计划。

默认值应保持当前 CPU 语义；ablation-03 的高性能配置通过显式环境变量或配置传入，不能依赖隐藏的全局选项。缓存键必须包含后端、线程/worker 配置和算法实现版本，避免把不同执行后端的结果混用。

### 优先做安全的 CPU 优化

- 将 learning-curve 和 scaling 的独立任务整理成带有 `fraction_index`、`repeat_id`、`representation` 或 `sequence_id` 的确定性 job table。
- 使用有界并行执行 job table；每个 worker 固定自己的 seed 和线程数，最终按 job key 排序合并，保证结果顺序、随机性和配对关系不变。
- 避免嵌套过度订阅：例如 4 个 worker 时每个 worker 分配 6 个 XGBoost 线程，而不是 4 个 worker 各自再开 32 个线程。实际组合通过 benchmark 选择，不预设“线程越多越快”。
- 在单个 readout 过程中复用已经转换好的矩阵、metadata、test `DMatrix` 和固定的 fold 索引，减少重复 `as.matrix`、切片和对象序列化；不得改变训练行集合和折叠边界。
- Windows worker 初始化时各自加载一次只读输入，避免每个任务反复传输约 1 GB 的表示对象；同时设置内存门槛，内存不足时自动退回较少 worker。

### 对 GPU 加速设置严格的可行性门槛

先用当前 xgboost 构建验证 CUDA 可用性，并确认 `gblinear/coord_descent` 是否真正产生 GPU kernel 和 GPU 利用率。若该路径不能 GPU 加速，不得仅添加 `device="cuda"` 就宣称完成。

在不改变线性目标、类别编码、class weight、lambda、nrounds、seed 和预测契约的前提下，评估两个可选后端：

- GPU 仅承担矩阵变换、decoder 的 PCA/回归等可证明等价的线性代数操作；分类拟合仍走 CPU。
- 实现一个与现有 `gblinear/coord_descent` 目标等价的 GPU 试验后端，作为独立实验，不直接替换默认路径。只有当预测、概率和指标在预设容差内一致，且真实任务 GPU 利用率和 wall time 均改善时才进入生产配置。

不得用 GPU `gbtree/hist` 直接替代当前线性模型；那会改变算法而不是性能优化。

### 细化 checkpoint 和可观测性

将 learning-curve 至少拆成 fraction/repeat/representation 级 checkpoint，scaling 拆成 sequence/module-count 级 checkpoint。每个 checkpoint 写入输入 hash、配置 hash、seed、后端、线程计划、状态和耗时；写入继续采用临时文件后原子替换。

日志应能回答“当前运行到哪个 job、已完成多少、最近一个 job 用时多久、预计剩余多久”，但不改变正式结果文件的 schema。中断后只重跑缺失或失效 job。

## 实施顺序

先完成基线和 profiling，再实现运行时配置与 CPU 任务调度；确认 CPU 版本结果等价后，再做矩阵复用和细粒度 checkpoint。GPU 后端最后单独验证，不让尚未证明的 CUDA 路径阻塞稳定 CPU 版本。

预期涉及的正式文件范围：

- `R/ablation.R`：运行时上下文、任务拆分、矩阵复用、checkpoint、可选后端接口。
- `test/ablation-03/05.00.00. 表示分析_functions.R`：性能配置入口和 benchmark 参数，不改变科学参数默认值。
- `test/ablation-03/tools/`：固定输入的性能基准与资源采样脚本。
- `test/ablation-03/tests/`：执行计划稳定性、checkpoint 恢复、结果等价和缓存键回归测试。
- `test/ablation-03/README.md`、`CHANGELOG.md`：记录运行参数、后端选择、兼容性和 `[Unreleased]` 性能改进。

## 验收标准

- 业务等价：同一输入、seed、配置和算法版本下，样本顺序、预测类别、概率、指标、抽样集合和正式结果字段一致；GPU 候选允许的浮点误差必须事先明确并逐项报告。
- 恢复正确：在任意 job 后中断，重新运行只补齐缺失 checkpoint，最终结果与一次完整运行一致。
- 并行正确：不同 worker 数和线程组合不会改变随机抽样、fold、lambda 选择或配对比较；无数据竞争、锁冲突和 cache 污染。
- 资源有效：完整运行期间能观察到多个 CPU worker 的实际负载；GPU 后端若启用，`nvidia-smi` 能看到 R 进程的计算使用，而不是只有桌面占用。
- 性能目标：以当前完整运行作为基线，优先目标是端到端至少 3 倍加速；若某一后端低于 1.5 倍或内存峰值不可接受，则不纳入默认配置。
- 质量门禁：定向 ablation-03 测试、`R CMD check`、结果契约检查和 BAC 验证均通过。

## 预期时间缩短幅度

以下预测以当前 i9-14900K（24 核 / 32 线程）、64 GB 内存、RTX 4090，以及当前 ablation-03 完整参数为基准。现有阶段时间来自本轮缓存时间戳；尚未完成的阶段按任务数量和已完成 readout 的单次拟合成本估算。正式实施前应通过 benchmark 校准，因此这些数字是决策区间，不是完成承诺。

| 阶段 | 当前预计耗时 | CPU 优化后 | GPU 条件达标后 | 主要依据 |
|---|---:|---:|---:|---|
| native geometry + retrieval | 约 3–5 分钟 | 约 2–4 分钟 | 约 1–3 分钟 | 当前已经较快，索引和距离计算只有有限优化空间 |
| readout | 实测约 39 分钟 | 约 12–20 分钟 | 约 8–15 分钟 | 两个表示可并发，矩阵和 fold 可复用 |
| learning-curve | 约 10–15 小时 | 约 2.5–4.5 小时 | 约 1.5–3 小时 | 800 次拟合可按 job 并行；受内存带宽和 Windows worker 复制限制，不按核心数线性缩短 |
| cohort scaling | 约 1–3 小时 | 约 20–60 分钟 | 约 15–45 分钟 | 约 36 次独立拟合可并行，且已有拟合缓存 |
| decoder | 约 5–30 分钟 | 约 3–15 分钟 | 约 1–8 分钟 | PCA、矩阵乘法和回归可能受益于并行 BLAS 或可选 GPU 线性代数 |
| 完整 05 表示分析 | 约 12–19 小时 | 约 3–6 小时 | 约 2–4 小时 | 包含阶段调度、checkpoint 和写盘开销 |

预期的第一阶段收益来自 CPU 调度和数据复用，目标是端到端缩短约 **60%–80%**，即约 **3–5 倍加速**。这是本计划的主要、可信度较高的性能目标。

GPU 的额外收益属于条件目标。若保持现有 `gblinear + coord_descent` 时不能证明 GPU 真正承担计算，则 GPU 不计入正式收益；此时交付以约 3–6 小时的 CPU 优化版本为目标。若等价 GPU 后端通过数值与资源验收，完整 05 阶段可进一步争取缩短到约 2–4 小时，相对当前约 **4–8 倍加速**。

性能验收采用三级标准：

- **最低可接受：** 完整阶段不超过基线的 50%，即至少 2 倍加速，同时结果等价、内存稳定。
- **计划目标：** 完整阶段约 3–6 小时，即 3–5 倍加速；learning-curve 不超过约 4.5 小时。
- **进取目标：** 等价 GPU 后端成立时约 2–4 小时，即 4–8 倍加速。

若某项优化使单阶段快于基线但导致完整运行内存峰值超过安全门槛、数值结果超出容差或恢复行为失效，则不计为有效加速并回退该项实现。

## 风险与回滚

- GPU 线性后端可能无法在数值和速度上同时达标；保留 CPU 后端作为默认和回滚路径。
- Windows 多进程复制大矩阵可能造成内存峰值；运行时必须限制 worker，并在超出门槛时自动降级。
- 并行结果的浮点归约顺序可能造成极小差异；合并顺序、容差和报告方式必须固定，不能把不可重复差异隐藏在缓存中。
- 新 checkpoint 与旧缓存 schema 不兼容时，使用新的 schema/算法 revision 让旧缓存安全失效，不覆盖旧正式结果。
- 任何正式实现前，先停止或等待当前 ablation-03 运行完成，避免新旧代码同时写入同一缓存根目录。
