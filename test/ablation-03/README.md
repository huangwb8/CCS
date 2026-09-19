# ablation-03：表示性能、生物学锚点与结构可复现性

本目录按 `bensz-rmd-rules` 的编号分析单元组织。计算层只写可恢复数据产品，Rmd 只读取完成批次并负责展示、正式图表与证据锚定解读。

## 目录边界

- 工作流清单：根目录的 `analysis-plan.yaml` 是编号、依赖、输入与交付物的机器可读事实来源。
- 外部输入与科学配置：`raw/` 只保存输入边界说明和版本化的只读配置；真实基因组数据仍位于仓库外。
- 大体积缓存：`D:/cache/ccs/_ablation-03/`，可用 `CCS_ABLATION_CACHE_ROOT` 显式覆盖。
- 轻量产品清单：`products/main/<分析单元>/`，包含 `main.rds`、`summary.md`、`metadata.yaml` 与 `SUCCESS`。
- 正式图表与表格：`reports/figures/`、`reports/tables/`。
- 公共模板与辅助入口：`templates/` 保存工作流/Rmd helper，`scripts/` 保存可单独运行的维护和基准脚本。
- 质量保障：`tests/` 保存不依赖交互输入的合同测试。
- Rmd 与同名 HTML：本目录根部。
- 原有 `tmp/`：仅保留历史结果，不再作为新流程的读写入口。

## 分析单元

| 单元 | 职责 | 外部缓存阶段 |
|---|---|---|
| `01.00.00. 数据准备` | 输入、组织映射、CCS 对象与数据画像 | `01-data` |
| `02.00.00. 表示输入准备` | Direct/d1、样本契约与端点资格 | `01-representations` |
| `03.00.00. 生物输入准备` | 连续锚点与结构表达缓存 | `01-biology` |
| `04.00.00. 数据概览` | 输入审计报告 | 无 |
| `05.00.00. 表示分析` | 几何、retrieval、readout、learning curve、scaling、decoder | `ablation-experiment` |
| `06.00.00. 生物锚点分析` | 连续锚点效用与队列级推断 | `ablation-biology` |
| `07.00.00. 结构复现分析` | 双向结构复现与 matched-bank 敏感性 | `ablation-structural-reproducibility` |

## 运行

### targets 迁移入口（当前串行骨架）

新的运行入口是根目录的 `_targets.R`。它只加载目标 R 库中的 `CCS`，不
`source()` 仓库 `R/ablation.R`，正式 targets store 固定在
`D:/cache/ccs/_ablation-03/targets`（可用 `CCS_ABLATION_CACHE_ROOT`
覆盖）。首次运行需要准备一个 RDS 输入清单，并通过环境变量提供路径：

### renv 项目环境

`ablation-03` 的正式运行必须使用本目录的项目级 `renv` 环境。`renv` 锁定 R 包版本，`targets` 仍负责分析依赖、缓存和恢复；CCS 必须是已安装的 `0.8.3` 包，不能通过 `load_all()` 或 `source()` 代替。

首次初始化（只需在当前分析项目执行一次）：

```powershell
$r = 'C:/R/R-4.3.1/bin/Rscript.exe'
Push-Location 'test/ablation-03'
& $r --vanilla 'scripts/renv-ablation03.R' --command init
& $r --vanilla 'scripts/renv-ablation03.R' --command install-core
& $r --vanilla 'scripts/renv-ablation03.R' --command install-ccs
& $r --vanilla 'scripts/renv-ablation03.R' --command snapshot
Pop-Location
```

在干净 R 库或另一台机器复现：

```powershell
Push-Location 'test/ablation-03'
& $r --vanilla 'scripts/renv-ablation03.R' --command restore
& $r --vanilla 'scripts/renv-ablation03.R' --command check
& .\scripts\run-targets-renv.ps1 -Action manifest
& .\scripts\run-targets-renv.ps1 -Action make
Pop-Location
```

`run-targets-renv.ps1` 会先显式激活并检查 renv、`CCS 0.8.3`、lockfile 和核心依赖，再调用 targets。缺少 renv 或 CCS 版本不符时会提前失败。

```powershell
$env:CCS_ABLATION_INPUT_RDS = 'D:/cache/ccs/inputs/ablation-03-inputs.rds'
$env:CCS_ABLATION_RUN_ID = 'review-01'
$r = 'C:/R/R-4.3.1/bin/Rscript.exe'
Push-Location 'test/ablation-03'
& $r --vanilla -e "source('renv/activate.R'); targets::tar_manifest(script = '_targets.R')"
& $r --vanilla -e "source('renv/activate.R'); targets::tar_make(script = '_targets.R')"
Pop-Location
```

输入 RDS 至少包含 `object`（已安装 CCS 0.8.3 生成的 `CCS` 对象）、`data`
和 `metadata`；为兼容 `01.00.00. 数据准备.R` 的既有产物，
`resCCS_ablation`、`data_all`、`ablation_metadata` 也会分别作为这三个字段的
别名读取。运行入口会拒绝其它 CCS 版本，避免阶段对象与包实现不一致。
`biology_inputs` 与 `structural_inputs` 在相应 target 启用前也必须显式提供。
当前迁移阶段只建立串行依赖图，不启用 crew/future，
也不运行本目录原有的编号脚本。待串行结果经审阅并完成科学等价性验证后，
再引入动态分支和受控并行。

`targets/functions.R` 只负责路径、输入契约和运行元数据；计算通过
`CCS::ablation(step = ...)` 的 context → plan → run → result 阶段完成。开发期可以在 CCS 包项目中使用 `load_all()` 做
包测试，但不得把它作为 ablation-03 的正式运行方式。

以下编号脚本和历史产品暂时保留用于回溯，不属于新的 targets 依赖图：

在仓库根目录使用 Windows R 4.3.1：

```powershell
$r = 'C:/R/R-4.3.1/bin/Rscript.exe'
& $r --vanilla 'test/ablation-03/01.00.00. 数据准备.R'
& $r --vanilla 'test/ablation-03/02.00.00. 表示输入准备.R'
& $r --vanilla 'test/ablation-03/03.00.00. 生物输入准备.R'
& $r --vanilla 'test/ablation-03/05.00.00. 表示分析.R'
& $r --vanilla 'test/ablation-03/06.00.00. 生物锚点分析.R'
& $r --vanilla 'test/ablation-03/07.00.00. 结构复现分析.R'
```

任一阶段启动时先写 `run-state.rds(status = "running")` 并使旧完成标记失效；正式结果、来源凭据和项目内 `SUCCESS` 全部复验后才提交 `complete`。报告拒绝读取非 complete 批次。

流程支持 `BENSZ_FORCE_STEP=AA.BB.CC` 强制重算单元，以及 `BENSZ_RESUME_FROM=AA.BB.CC` 从指定单元恢复；恢复点之前若有任一无效 checkpoint，会立即停止并要求先修复前序单元。示例：

```powershell
$env:BENSZ_FORCE_STEP = '05.00.00'
& $r --vanilla 'test/ablation-03/05.00.00. 表示分析.R'
Remove-Item Env:BENSZ_FORCE_STEP

$env:BENSZ_RESUME_FROM = '05.00.00'
& $r --vanilla 'test/ablation-03/01.00.00. 数据准备.R'
& $r --vanilla 'test/ablation-03/02.00.00. 表示输入准备.R'
& $r --vanilla 'test/ablation-03/03.00.00. 生物输入准备.R'
& $r --vanilla 'test/ablation-03/05.00.00. 表示分析.R'
Remove-Item Env:BENSZ_RESUME_FROM
```

同一外部缓存根目录禁止并发运行多个阶段。阶段持有 `D:/cache/ccs/_ablation-03/.workflow-lock/`；正常完成后自动释放。若 R 进程异常退出，确认对应进程已经停止后，才可手动删除该锁目录并恢复运行。

表示分析在 `ablation-experiment/checkpoints/` 下按内容键保留 retrieval（含 anchor retrieval 与 controls）、readout、learning curve 和 decoder 多版本缓存。learning curve 还按 fraction/repeat/representation 保存 job 级 checkpoint，中断后只补算缺失或损坏的 job。相同输入、参数、seed、线程/worker 配置、算法修订、节点实际代码依赖与运行库版本命中缓存；改变 learning curve 设计不会使 retrieval 或 readout 失效，无关节点代码变化也不会清空其它节点。每个缓存同时具有独立 state 和 value hash；`running`、截断 RDS、schema/key/hash 不一致均视为 miss 并安全重算。正式 `manifest.rds` 的 `node_cache` 会记录本轮 `hit`/`miss`、写入状态及具体失效原因。

表示分析的 learning curve 支持显式的有界 PSOCK 并行。默认 `CCS_ABLATION_WORKERS=1` 保持原有 CPU 语义；在内存允许时可设置为 `2` 或 `4`，脚本会把 `CCS_ABLATION_CORES` 平均分配给 worker，避免 XGBoost 线程过度订阅。可选的 `CCS_ABLATION_MEMORY_GB` 设置 worker 总内存预算；入口按四个只读表示矩阵体积的三倍估计单 worker 峰值，并自动下调 worker 数。并行 job 按固定的 fraction/repeat/representation 顺序合并，seed、抽样和结果字段不变。建议先用小规模输入比较 wall time 与内存峰值，再用于完整运行：

```powershell
$env:CCS_ABLATION_WORKERS = '2'
$env:CCS_ABLATION_MEMORY_GB = '48'
& $r --vanilla 'test/ablation-03/05.00.00. 表示分析.R'
Remove-Item Env:CCS_ABLATION_WORKERS
Remove-Item Env:CCS_ABLATION_MEMORY_GB
```

若要启动一条与历史结果完全隔离的新分析线，可直接运行
`scripts/run-fresh-analysis.ps1`。它默认使用
`D:/cache/ccs/_ablation-03-fresh`，设置 Windows R 的可用 UTF-8 locale，先完成
01–03 准备，再运行 learning-curve 基准和 05–07 分析；不会删除或覆盖旧的
`_ablation-03` 缓存。可用 `-Workers`、`-TotalCores`、`-MemoryGB` 和
`-SkipBenchmark` 调整执行策略。

可先运行 `scripts/benchmark-learning-curve.R` 比较相同持久输入、seed 和 XGBoost 版本下的串行/并行耗时与结果等价性。默认只取第一个 fraction 和一次 repeat；可用 `CCS_ABLATION_BENCHMARK_FRACTIONS`、`CCS_ABLATION_BENCHMARK_REPEATS` 和 `CCS_ABLATION_BENCHMARK_WORKERS` 扩大基准规模。结果写入外部缓存阶段的 `performance/learning-curve-benchmark.csv`，不覆盖正式分析结果；若正式阶段持有 workflow lock，benchmark 会拒绝启动。

## 报告渲染

使用 `knit-rmd-html` 渲染全部报告；同名 HTML 写回本目录，矢量 PDF 写入 `reports/figures/`。JPG 检查预览优先写入当前 `BENSZ_TASK_ROOT`，普通人工运行则写入外部缓存的 `plot-previews/`。

外部 query 资格、表示尺度、生物锚点评分和结构复现的科学边界保持不变；本轮只重构执行边界、恢复语义和产物组织。
