# ablation-03：表示性能、生物学锚点与结构可复现性

本目录是独立于 CCS 包构建内容的可复现分析项目。正式分析由 `_targets.R` 声明完整
依赖图，并由 `targets/` wiring、正式 targets store 和 `tar_make()` 统一调度、缓存与
恢复；不得使用逐脚本 runner 生成正式结果。六个编号 R 脚本是被 targets 节点消费的
科学实现，`run-ablation-03.R` 仅保留为历史/轻量诊断入口，不属于正式分析入口。

## 目录边界

- `analysis-plan.yaml`：编号、依赖、输入和交付物的机器可读清单。
- `raw/`：输入边界说明与版本化配置；真实基因组数据位于仓库外。
- `products/main/`：仅由 `formal` profile 写入的可审计 checkpoint 清单。
- `reports/`：正式图表与表格。
- `scripts/helpers/`：编号脚本共用的缓存、凭据和运行边界 helper。
- `tests/`：不依赖交互输入的合同测试。
- `tmp/`：本地轻量验收和临时输出；不属于 CCS 包。
- 大体积缓存不写入 CCS 包根目录；正式 targets store 固定为
  `D:/cache/ccs/_ablation-03/targets`，由 `_targets.yaml` 和 wiring 统一声明。

项目级 `renv/` 与 `renv.lock` 只属于本分析目录，不在 CCS 包根目录创建环境。
正式编号脚本使用本机已安装的 CCS `0.8.3`；不得 `source()` 仓库
`R/ablation.R` 或用 `pkgload::load_all()` 替代安装包。

## 分析单元

| 单元 | 职责 | 外部缓存阶段 |
|---|---|---|
| `01.01.00. 数据准备` | 输入、组织映射、CCS 对象与数据画像 | `01-data` |
| `01.02.00. 表示输入准备` | Direct/d1、样本契约与端点资格 | `01-representations` |
| `01.03.00. 生物输入准备` | 连续锚点与结构表达缓存 | `01-biology` |
| `02.01.00. 表示分析` | 几何、retrieval、readout、learning curve、scaling、decoder | `ablation-experiment` |
| `02.02.00. 生物锚点分析` | 连续锚点效用与队列级推断 | `ablation-biology` |
| `02.03.00. 结构复现分析` | 双向结构复现与 matched-bank 敏感性 | `ablation-structural-reproducibility` |

`01.04.00. 数据概览.Rmd` 及三个 `02.*.Rmd` 只消费已完成结果，不参与计算调度。

## 正式运行（targets）

在 CCS 仓库根目录使用 Windows R 4.3.1，通过 targets launcher 执行：

```powershell
$r = 'C:/R/R-4.3.1/bin/Rscript.exe'
& $r --vanilla 'test/ablation-03/scripts/renv-ablation03.R' --command check
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' -Action make
```

targets 使用 `crew::crew_controller_local()` 调度可并行的独立 target。默认最多启动 2 个
worker，可通过 `CCS_ABLATION_TARGET_WORKERS` 显式调整。每个 worker 的结构化日志和 CPU/RAM
采样写入正式 cache root 下的 `logs/targets-crew/workers/`，targets 主进程日志写入
`logs/targets-crew/main-process.log`。运行结束后，`resource_metrics` 和 `worker_health` target
读取这些记录；可用 `autometric::log_read()` 检查数据，并用
`autometric::log_plot(metrics, metric = "resident")` 绘制单个 worker 的内存曲线。

查看依赖图或监测正式运行：

```powershell
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' -Action manifest
& $r --vanilla -e "source('test/ablation-03/renv/activate.R'); targets::tar_watch(config='_targets.yaml', project='main', browse=TRUE)"
```

或使用 launcher 的 dashboard 操作：

```powershell
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' -Action watch
```

`tar_watch()` 必须读取 `D:/cache/ccs/_ablation-03/targets`；不得用旧 store 或自定义
learning-curve checkpoint 伪装成 targets 进度。

启动时写入 `.ablation03-root.rds`、`.ablation-entry-lock/owner.rds` 和
`runs/<run_id>/runtime.rds`。同一缓存根目录绑定一种 profile；活跃进程持有的缓存拒绝
并发，失主的本机锁先归档至 `stale-locks/` 再恢复。阶段异常会把对应 state 写为
`failed`，陈旧节点会标记为 `stale` 后重算。

完整数据位置仍可通过既有环境变量显式覆盖：`CCS_DATA_ROOT`、`CCS_SYNC_ROOT`、
`CCS_ABLATION_MODEL_ROOT`、`CCS_FULL_RESCCS_RDS`、`CCS_FULL_EXPRESSION_RDS` 和
`CCS_GENE_SIGNATURE_RDS`。入口本身不猜测缓存路径。

## 轻量真实样本验收

`lightweight` profile 仍执行完全相同的六个编号脚本，只通过已记录的运行参数限制
reference/query 样本数、learning-curve 重复和树轮数，降低结构 bootstrap 次数和
状态最小样本数，并关闭耗时的 scaling、null controls 与 decoder。样本上限仍保留
足够的 cohort 内重复，避免把结构端点拆成每个 cohort 仅一两个样本。测试缓存应放在
`test/ablation-03/tmp/` 的独立子目录：

轻量运行的 checkpoint 清单也写入显式 `--cache-root` 下的 `products/main/`，不会创建或
改写项目内的正式 `products/main/`。

```powershell
& $r --vanilla 'test/ablation-03/run-ablation-03.R' `
  --cache-root 'test/ablation-03/tmp/lightweight-real-20260920/cache' `
  --profile lightweight --cores 2 --workers 1 --memory-gb 16
```

如需改变轻量默认值，可在启动前显式设置
`CCS_ABLATION_MAX_REFERENCE_SAMPLES`、`CCS_ABLATION_MAX_QUERY_SAMPLES`、
`CCS_ABLATION_LEARNING_FRACTIONS`、`CCS_ABLATION_REPEATS`、
`CCS_ABLATION_NROUNDS`、`CCS_ABLATION_STRUCTURAL_MIN_ENTITY_N`、
`CCS_ABLATION_STRUCTURAL_BOOTSTRAP` 或 `CCS_ABLATION_MATCHED_REPEATS`。这些值
进入参数哈希，不会与正式缓存混用。

## 可选 benchmark

benchmark 是唯一额外运行脚本，只读一个已经完成的 formal/lightweight 缓存，并要求
显式、位于缓存根目录之外的输出目录：

```powershell
& $r --vanilla 'test/ablation-03/scripts/benchmark-learning-curve.R' `
  --cache-root 'test/ablation-03/tmp/lightweight-real-20260920/cache' `
  --output-dir 'test/ablation-03/tmp/lightweight-real-20260920/benchmark' `
  --workers 2
```

它比较相同输入、seed 和 XGBoost 版本下的串行/有界 PSOCK 结果，要求指标在容差内
一致，并写出 `learning-curve-benchmark.csv`。正式入口或阶段锁存在时 benchmark
拒绝运行，也不会向科学缓存写入结果。

## 报告与恢复

阶段完成后才提交 `SUCCESS` 和 stage receipt；报告拒绝读取非 `complete` 批次。
representation 节点和 learning-curve job 使用内容键、多版本 value/state 与 hash 校验。
改变 learning-curve 参数不会无故清空 retrieval/readout，失败或损坏的 job 只补算自身。

使用 `knit-rmd-html` 渲染报告；同名 HTML 写回本目录，矢量 PDF 写入
`reports/figures/`。外部 query 资格、表示尺度、生物锚点评分和结构复现的科学定义
不因入口收敛而改变。
