# ablation-03：表示性能、生物学锚点与结构可复现性

ablation-03 是独立的可复现分析项目。`_targets.R` 是唯一的分析编排入口，负责依赖图、缓存、恢复、运行状态和 `crew` 并行调度。六个编号 R 文件保留为科学计算实现，由 targets 节点在同一 R 进程中消费；它们不再拥有 stage lock、stage receipt、SUCCESS 或手工 cache hit/miss 协议。

formal 分析与轻量验收使用完全相同的 targets 图、CCS 版本、参数、随机种子和科学代码。两者唯一允许的差异是 `CCS_ABLATION_INPUT_RDS` 指向的输入数据：轻量验收可以使用正式数据子集或同 schema 的自制 fixture。

## 分析单元

| 单元 | 科学职责 | targets 产物目录 |
|---|---|---|
| `01.01.00` | 输入对象、表达数据和元数据准备 | `01-data` |
| `01.02.00` | Direct/d1 表示输入与样本契约 | `01-representations` |
| `01.03.00` | 生物锚点与结构分析输入 | `01-biology` |
| `02.01.00` | 几何、检索、readout、learning curve、scaling、decoder | `ablation-experiment` |
| `02.02.00` | 生物锚点效用与队列级推断 | `ablation-biology` |
| `02.03.00` | 双向结构复现与 matched-bank 敏感性 | `ablation-structural-reproducibility` |

R Markdown 文件只消费已完成 target 产物，用于报告渲染，不参与科学计算调度。

## 正式运行

必须使用 Windows R 4.3.1、项目 renv 和已安装的 CCS 0.8.3：

```powershell
$r = 'C:/R/R-4.3.1/bin/Rscript.exe'
$input = 'D:/data/ccs/ablation-formal-inputs.rds'
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action make -InputRds $input -CacheRoot 'D:/cache/ccs/_ablation-03'
```

查看依赖图或监测运行：

```powershell
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action manifest -InputRds $input -CacheRoot 'D:/cache/ccs/_ablation-03'
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action watch -InputRds $input -CacheRoot 'D:/cache/ccs/_ablation-03'
```

`targets` 使用 `crew::crew_controller_local()`。worker 日志和 CPU/RAM 采样写入 cache root 下的 `logs/targets-crew/`，由 `resource_metrics` 与 `worker_health` target 读取。

输入 RDS 至少应包含 `object`、`data`、`metadata`；兼容旧字段名 `resCCS_ablation`、`data_all`、`ablation_metadata`。如果输入没有预计算的 biology/structural 输入，需通过 `CCS_FULL_EXPRESSION_RDS`、`CCS_GENE_SIGNATURE_RDS` 等环境变量提供同一数据边界所需的原始资源。

## 轻量验收

轻量验收不再是另一种 profile，也不减少 repeats、bootstrap、scaling、null controls 或 decoder。只需把 `-InputRds` 换成正式数据子集或同 schema 的小型 fixture，并使用相同 launcher、相同参数和相同 targets 图：

```powershell
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action make -InputRds 'test/ablation-03/raw/fixture-inputs.rds' `
  -CacheRoot 'test/ablation-03/tmp/targets-fixture'
```

样本减少可能使某些端点变为 `not_estimable`，但不得改变分析方法、target 图或参数语义。

## 环境与包边界

ablation-03 只依赖已安装的 CCS 包，不从仓库 `source()` `R/ablation.R`，也不使用 `pkgload::load_all()`。修改 CCS 源码后，必须先用 `C:/R/R-4.3.1` 构建并覆盖安装版本 `0.8.3`，再运行 targets。

targets store、科学产物和 observability 均位于显式 cache root；不得把大型 RDS 或个人/专有基因组数据提交到仓库。
