# ablation-03 targets/renv/Rmd/并行重建计划

## 通俗解释：究竟发生了什么

当前 `R/ablation.R` 已经是一条能工作的消融分析流水线，但它同时承担了科学计算、流程编排、缓存、断点恢复和部分并行。`test/ablation-03` 又在脚本层维护了一套阶段状态和产品目录。这样同一件事由多层代码重复管理，长期修改时容易出现缓存、状态和实际计算不一致。

这次不做兼容迁移，而是基于 Git 中已保存的旧版本，重建一条新分析线：`targets` 成为唯一的依赖、缓存和恢复管理者；`renv` 作为 ablation-03 分析项目的环境锁定层，在串行骨架稳定后按正式复现需要启用；现有 Rmd 继续生成报告；crew 或 future 只承担一个明确的并行层。

改造后的直接变化是：改报告只重跑报告，改 retrieval 只影响 retrieval 及下游，learning curve 可以按独立任务恢复；不再维护旧 runner、旧 stage receipt 和新 targets store 的双重真相。

## 专业判断：问题在哪里

### 当前事实

- `R/ablation.R` 已包含 Direct/d1 准备、native geometry、retrieval、readout、learning curve、cohort scaling、decoder、结果汇总及自定义 checkpoint。
- `test/ablation-03` 已有 01–07 阶段脚本、Rmd、analysis plan、产品和报告资产。
- 当前自定义缓存和 workflow 不是错误，但职责与 targets 重叠；如果只在外层套 targets，会形成两套缓存、两套失效规则和潜在嵌套并行。
- learning curve 已经有可拆分的 `fraction × repeat × representation` 任务结构，适合改造成 targets dynamic branching。

### 本次决策

- 旧编号入口、旧产品、旧缓存和旧 workflow 不再作为运行时兼容目标；它们由 Git 和历史文件承担回退作用。
- 新流程使用新的目录和新的 targets store，不读取旧 workflow 的完成标记，也不把旧缓存提升为新缓存。
- CCS 包本身不嵌入 `renv`；包开发与质量门禁继续由 `DESCRIPTION`、`NAMESPACE`、`R CMD check` 和项目既有 R 版本约定负责。
- `renv` 只服务 `test/ablation-03` 的分析运行环境，且不是第一阶段 targets 串行骨架的硬前置；正式复现、跨机器运行或交付前再生成/提交 `renv.lock` 并执行 clean-library restore。
- 改造目标是“重用科学计算逻辑，重建运行边界”，不是重写统计方法。

## 要达到什么目标

- `R/ablation.R` 提供纯计算、可组合、可测试的分析节点和单任务函数。
- `test/ablation-03` 由 `_targets.R` 声明完整依赖图；不再依赖自定义编号 runner 执行主流程。
- ablation-03 在需要正式复现或跨机器交付时由项目级 `renv.lock` 固定 R 包环境；开发早期可以先不启用 renv。
- 现有 Rmd 保留为报告层，由 targets 的 `format = "file"` target 调用 `rmarkdown::render()`。
- learning curve 先保证串行恢复，再使用 crew dynamic branching；若环境或任务结构更适合函数内并行，使用 future，二者不默认嵌套。
- 新流程完成后，可以从一个全新的 targets store 独立完成数据准备、表示分析、生物分析、结构分析和报告渲染。

## 第一部分：如何改造 `R/ablation.R`

### 建立新的职责边界

把现有代码分成三层：

1. **纯计算层**：只接收 R 对象、参数和 seed，返回 R 对象；不写 workflow state、不创建 SUCCESS、不决定全局缓存路径。
2. **任务层**：封装一个可独立运行的科学任务，例如一个 learning-curve job 或一个 scaling design cell；保证 seed、输入行和输出字段确定。
3. **组合层**：由 targets 调用的组合函数，负责把节点输出拼接成最终结果；不再在内部模拟 targets 的节点失效。

### 暴露 representation 节点

将现有内部逻辑整理为以下可调用边界，具体命名可在实施时确定：

- `prepare_representation_inputs()`：Direct、既有 d1、metadata、feature/module manifest、query eligibility 和 selected blocks；
- `run_native_geometry()`：Direct/d1 native geometry 及必要诊断；
- `run_retrieval()`：癌种 retrieval、anchor retrieval、technical controls 和资格审计；
- `run_readout()`：两种 representation 的 supervised readout、配对结果和预测；
- `make_learning_curve_jobs()`：生成稳定排序的 fraction/repeat/representation job table；
- `run_learning_curve_job()`：只运行一个 job，返回固定 schema 的一行或一个小对象；
- `combine_learning_curve_jobs()`：按稳定 job key 合并并生成 paired curve；
- `make_scaling_jobs()` / `run_scaling_job()` / `combine_scaling_jobs()`：将 cohort scaling 的设计单元与汇总分开；
- `run_decoder()`：decoder 和 feature-type reconstruction；
- `assemble_representation_result()`：组合节点结果和审计摘要。

### 移除或旁路内部 workflow 责任

targets 接管后，以下责任从 `R/ablation.R` 移出：

- `.ablation_cached_node()` 及 retrieval/readout/learning-curve/decoder 的父级节点缓存；
- native geometry、Direct feature 和 scaling fit cache 的项目级路径管理；
- 以 `output.dir` 为中心的“目录非空/完成标记/恢复”判断；
- learning curve 内部的 PSOCK 调度。

科学函数仍可保留原子级临时对象，但不得将其作为 targets 节点完成状态。必要时先保留内部函数作为过渡实现，再删除重复缓存代码；新 targets 入口不能同时依赖旧 checkpoint 才能正确运行。

### 保持不变的科学契约

- Direct-GSClassifier 与 Cohort-d1 的定义、列顺序和 module block 合同；
- reference/query 样本资格、endpoint eligibility、tissue/cancer anchor 规则；
- native geometry、retrieval、readout、learning curve、scaling、decoder 的指标定义；
- 抽样、inner folds、lambda 候选、nrounds、seed 生成规则；
- 现有审计字段、paired comparison 和统计解释边界。

### `R/ablation.R` 的验证门槛

- 每个新节点都能在合成输入上独立执行，不需要 workflow lock、stage receipt 或外部产品目录。
- 单个 learning-curve job 的串行重复运行结果相同；job 合并顺序不依赖 worker 返回顺序。
- 新节点输出与旧实现进行固定 seed 的字段、维度、sample/feature hash 和数值容差比较。
- 不把 targets、crew 或 renv 的对象传入科学函数内部，避免基础设施和领域逻辑再次耦合。

## 第二部分：如何改造 `test/ablation-03`

### 建立新的项目布局

在 `test/ablation-03` 建立全新的运行入口和目录，不复用旧产品与旧外部缓存：

```text
test/ablation-03/
├── _targets.R
├── renv.lock                 # 正式复现/交付阶段生成；不是第一阶段硬前置
├── targets/
│   ├── functions.R
│   ├── config.R
│   └── report_helpers.R
├── Rmd/ 或根目录现有 Rmd
├── reports/
└── _targets/              # targets store；正式运行可通过配置放到外部磁盘
```

具体文件位置可遵循现有目录习惯，但必须只有一个 `_targets.R` 和一个 targets store 入口。旧 `products/`、旧 workflow helper 和旧 stage receipt 不进入新依赖图。

### 用 `_targets.R` 声明主依赖图

建议的 target 层级为：

```text
runtime_config
  ↓
data_inputs
  ↓
representation_inputs ──┬── native_geometry
                        ├── retrieval
                        ├── readout
                        ├── learning_jobs → learning_job_result → learning_curve
                        ├── scaling_jobs  → scaling_job_result  → cohort_scaling
                        └── decoder

biology_inputs → biology_result

representation_inputs + biology_inputs + native_geometry
  → structural_result

上述结果 → data/representation/biology/structural Rmd reports
```

目标图只拆独立失败边界、昂贵计算和多下游复用对象；不把每个临时 mutate、短生命周期矩阵或纯展示变量做成 target。

### 为 ablation-03 固定分析环境（可选分阶段启用）

`renv` 只在 `test/ablation-03` 项目初始化，不在 CCS 包源码或 `R/ablation.R` 内调用。第一阶段先让 targets 串行骨架和包接口通过小规模验证；确认依赖集合稳定后，再按实际调用锁定 targets、rmarkdown、knitr、xgboost、RcppAnnoy、irlba、CCS、GSClassifier、luckyBase、图表和报告依赖，提交 `renv.lock`，并在正式运行前用 clean library 做一次 restore 验证。

`renv` 只负责分析项目的包环境，不负责 target 失效；Git commit、参数文件、随机种子和（启用后）`renv.lock` 共同记录一次正式运行边界。CCS 包仍以固定 Git commit 或本地构建版本安装到分析环境中。

### 保留并改造现有 Rmd

04、05、06、07 Rmd 不迁移为 Quarto。每份报告变成一个 file target：

- 输入只列出相应的 target/product；
- Rmd 负责展示阈值、Top N、图表、表格和解读；
- 报告不重新训练模型、不自行判断缓存、不读取未声明临时文件；
- HTML、PDF、正式 figures/tables 的路径保持清晰且可追踪。

修改 Rmd 的颜色、文字或展示参数时，只失效报告 target，不重算上游分析。

### learning curve 并行化

第一步把 `make_learning_curve_jobs()` 和 `run_learning_curve_job()` 接入 targets，但先串行运行，确认 branch 输出和恢复语义。

第二步使用 targets dynamic branching：

```text
learning_jobs → learning_job_result pattern = map(learning_jobs)
```

默认优先 crew 作为 targets controller；每个 branch 内的 XGBoost 使用受控线程数。若选择 future，则只在一个 learning-curve target 内使用 `future.apply`，不再同时配置 crew。并行总预算、内存上限和随机数流必须通过真实小规模 benchmark 选择。

### scaling 与下游阶段

learning curve 稳定后，再判断 cohort scaling 是否值得拆成 dynamic branches；优先拆有昂贵拟合或明确恢复价值的 design cell，不为低成本 summary 增加 target。

biology 和 structural 分析作为独立 target 子图，直接消费声明的 biology/representation inputs 和结果，不读取旧 `ablation-experiment` 目录中的半成品。

## 共同实施顺序

1. **基线冻结**：保存当前测试、关键结果摘要、固定输入/seed、输出字段和允许数值容差；Git 作为唯一旧版本回退。
2. **改造 `R/ablation.R`**：先完成纯计算节点和单 job 接口；在 targets 之外用合成数据和小输入验证。
3. **建立 targets 串行骨架**：先接入 data、representation、native geometry、retrieval、readout、learning curve、scaling、decoder、biology、structural 和报告 target，不启用并行。
4. **决定是否启用项目级 renv**：在 targets 串行骨架验证依赖稳定后生成/更新 `renv.lock`；若只是当前机器的开发迭代，可暂缓，不阻塞接口改造。
5. **报告接入**：验证 Rmd 局部失效和报告产物契约。
6. **引入并行**：先串行 dynamic branch，再启用 crew 或 future 之一，完成结果等价和资源 benchmark。
7. **全流程验收**：新 store 从空目录运行完整 ablation-03，检查结果、恢复、报告、资源和日志；正式复现/交付场景还需通过 `renv::restore()` 后的同一验收。
8. **版本与安装门禁**：仅在新版 `R/ablation.R` 通过包检查、科学等价性验证、targets 串行全流程和必要的 renv restore 后，将 `DESCRIPTION` 的 patch 位递增 1；随后用 `C:\R\R-4.3.1` 对应 R 环境构建、检查并安装该版本 CCS 包，再记录安装库和版本证据。版本门禁前不得升版或安装“新版本”包。

## 如何确认完成

### `R/ablation.R`

- representation 的主要科学节点可独立调用；不需要旧 workflow 文件才能运行。
- 固定 seed 下新旧科学结果的样本顺序、feature 合同、抽样集合、预测和指标在预设容差内一致。
- learning-curve job 与 scaling job 可单独重跑并稳定合并。
- 科学函数不再负责 targets store、SUCCESS、stage receipt 或外部缓存路径。

### `test/ablation-03`

- `targets::tar_manifest()` 能显示完整目标图，`tar_outdated()` 的失效范围符合预期。
- 首次运行、相同输入重跑、单节点修改、报告修改、branch 中断和 branch 恢复均有可重复测试。
- 若本次启用了项目级环境锁定，`renv::restore()` 后可以执行 targets、分析节点和所有 Rmd 报告；未启用时，必须记录使用的 R 版本、CCS Git commit 和依赖版本快照。
- 只有所有前置验收通过后，`DESCRIPTION` 的 patch 位才从当前值递增 1；用 `C:\R\R-4.3.1` 对应 R 环境安装后，`packageVersion("CCS")` 与 `DESCRIPTION` 一致。
- crew/future 串行与并行结果一致；不会出现未预算的嵌套 worker、XGBoost 线程或内存复制。
- 新 store 可从空目录完成全流程，不读取旧产品和旧缓存。
- 现有 20–34 号科学/回归测试迁移到新入口或明确替换，不以旧 workflow 的通过结果代替新流程验收。

## 风险与待确认事项

- `R/ablation.R` 中原子函数与 workflow 编排边界较深，第一阶段必须避免边拆边改统计逻辑。
- targets dynamic branching 的对象序列化和 Windows 内存复制需要真实 benchmark；不能仅按逻辑 CPU 数设置 worker。
- crew 与当前 R/xgboost/Windows 环境的安装和兼容性需要实施前确认；如果不可用，先用串行 targets 或 future 过渡。
- `renv` 可能需要安装或编译包；初始化和 restore 需要明确的环境变更授权。
- patch 升级和安装属于发布门禁，不是开发循环；若任一科学等价性、包检查、targets 串行验收或 restore 验证失败，保持原版本并停止安装新版本。
- 新架构首次完整运行不应覆盖任何旧目录；只有新流程独立产物通过验收后，才可清理旧缓存。
- 计划重点是运行基础设施重建；科学结果若在双跑中出现差异，必须先停止迁移并定位差异来源。
