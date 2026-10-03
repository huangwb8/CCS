# ablation-03：表示性能、生物学锚点与结构可复现性

ablation-03 是独立的可复现分析项目。`_targets.R` 是唯一的分析编排入口，负责依赖图、缓存、恢复、运行状态和 `crew` 并行调度。六个编号 R 文件保留为科学计算实现，由 targets 节点在同一 R 进程中消费；它们不再拥有 stage lock、stage receipt、SUCCESS 或手工 cache hit/miss 协议。

formal 分析与轻量验收使用相同的 targets 图、CCS 版本、参数、随机种子和科学代码。测试子集重新计算生物输入；正式输入还包含预计算的生物和结构输入，小测不验证该分支。

## 分析单元

| 单元 | 科学职责 | targets 产物目录 |
|---|---|---|
| `01.01.00` | 输入对象、表达数据和元数据准备 | `01-data` |
| `01.02.00` | Direct/d1 表示输入与样本契约 | `01-representations` |
| `01.03.00` | 生物锚点与结构分析输入 | `01-biology` |
| `02.01.00` | 几何、检索、readout、learning curve、scaling、decoder | `ablation-experiment` |
| `02.02.00` | 生物锚点效用与队列级推断 | `ablation-biology` |
| `02.03.00` | 双向结构复现与 matched-bank 敏感性 | `ablation-structural-reproducibility` |
| `02.04.00` | decoder、局部 query、共享节点及固定训练设计的补充推断 | `statistical-inference` |

R Markdown 文件只消费已完成 target 产物，用于报告渲染，不参与科学计算调度。
表示报告另消费 `references/` 中 Figure S3 的历史 cohort Gain 汇总与来源记录，用于解释 d3 坐标重建 importance 与原始 d1 方差份额的区别。该输入不含患者级矩阵、不重拟合历史模型，纳入 `representation_report_sources` 文件依赖；当前方差份额仍来自正式 `representation_analysis` 产物。
四份 HTML 均由正式依赖图中的 file target 直出；不得在 `tar_make()` 之外用独立
runner 生成正式报告。
报告共用 `templates/` 下的 Liquid Glass 样式与交互模板；模板文件纳入四个报告的
targets 文件依赖。图表代码输出到 HTML 并默认折叠，可点击各代码块的 Code 按钮展开，
或用报告顶部的 Code → Show All Code 一次展开全部代码；setup 初始化代码不输出。
生物锚点报告附带配对效应、cohort 异质性和基因覆盖三张 PDF，并在
`reports/tables/` 导出队列级配对差值，供核对图中每个格子的样本数与方向。
配对效应图采用紧凑森林图与独立统计列，逐行显示差值及 95% CI、有效 cohort 数、
配对 query 人数和 BH 校正 P 值；统计列不参与效应横轴缩放，锚点顺序沿用冻结配置。
补充推断由 `statistical_inference`、`learning_query_inference` 和
`geometry_sensitivity` 三个 target 生成；报告分别依赖这些 target。
decoder 的 cohort 等权区间与原合并 query 逐特征分数分列；结构均值检验先检查共享节点
和模拟覆盖，稀疏网络保留 NA；100% 学习曲线的新区间只针对固定训练设计下
的 query cohort，原设计层 CI/P 值仍为 NA。几何诊断以 reference cohort
重采样给出 CKA 和距离排序的条件性 95% 区间；kNN Jaccard 固定完整近邻图后按 query cohort 重采样，P 值因没有预设零假设保持 NA；另报告
留一 reference cohort 敏感性，不把其范围标作 95% CI。

## 七组织深度设计

表示分析的 Figure 8 右列固定 ACC、BRCA、CRC、KIRC、PAAD、PRAD、STAD，每种 tissue 按 1–8 个冻结 cohort module 逐层增加，总模块数为 7、14、21、28、35、42、49、56。七种 tissue 与八层终点在看新结果前已固定；模块容量或 tissue 名称不满足配置时，正式设计直接失败。左列 breadth 继续逐次加入一种 tissue 的一个模块，Figure 9 继续使用原有同模块预算的 matched-size 配对。图中的区间以五次 bank 设计重复为单位，只描述冻结模块组合的探索性敏感性。

`ccs_code_files` 与 `ccs_code_identity` targets 跟踪项目 `R/ablation.R`、已安装 CCS 包的代码及版本描述，并在表示分析前核验关键函数体一致。`representation_sources` 同时跟踪科学参数和阶段脚本；源码或已安装包变化会触发 targets 失效，源码与已安装代码不一致时先停止正式计算。

## 正式运行

必须使用 Windows R 4.3.1、项目 renv 和与仓库 `DESCRIPTION` 版本一致的已安装 CCS 包。
正式输入快照来自已有的 01-data/inputs.rds，位于同一 cache root 但不会被数据准备阶段覆盖；切勿重新将 01-data/inputs.rds 作为输入。
当前 renv 锁文件中的 CCS、GSClassifier、luckyBase 缺少可恢复来源，现有机器已安装的包可运行，但跨机器恢复尚需补齐来源。

```powershell
$r = 'C:/R/R-4.3.1/bin/Rscript.exe'
$input = 'D:/cache/ccs/_ablation-03/formal-inputs.rds'
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action make -InputRds $input -CacheRoot 'D:/cache/ccs/_ablation-03'
```

launcher 会在调用 `tar_make()` 或 `tar_watch()` 前，先将 `_targets.yaml` 的 store 指向本次 `-CacheRoot/targets`。切换正式与隔离运行时须通过 launcher 设置 cache root；不要在另一运行仍活动时切换同一项目的 store 配置。正式报告默认写入分析项目根目录和 `reports/`，隔离运行须显式设置 `-OutputRoot`，把 HTML 与图表写入隔离目录。

查看依赖图或监测运行：

```powershell
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action manifest -InputRds $input -CacheRoot 'D:/cache/ccs/_ablation-03'
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action watch -InputRds $input -CacheRoot 'D:/cache/ccs/_ablation-03'
```

`targets` 使用 `crew::crew_controller_local()`。worker 日志和 CPU/RAM 采样写入 cache root 下的 `logs/targets-crew/`，由 `resource_metrics` 与 `worker_health` target 读取。

## 运行监控与卡点判读

以下命令只读取正式 targets store 和 observability 日志，不会启动第二套分析、修改缓存或干扰正在运行的 worker。先设置与正式运行一致的 cache root：

```powershell
$cacheRoot = 'D:/cache/ccs/_ablation-03'
$store = Join-Path $cacheRoot 'targets'
```

查看正式 target 状态：

```powershell
Get-Content (Join-Path $store 'meta/progress')
```

持续运行的 target 通常只显示为 `dispatched`。`targets` 不会把 target 内部每一步写入 store，因此可同时查看 store 中最近变化的文件：

```powershell
Get-ChildItem $store -Recurse -File |
  Sort-Object LastWriteTime -Descending |
  Select-Object -First 20 LastWriteTime, Length, FullName
```

target 完成、失败或提交对象时，`meta/` 与 `objects/` 才会出现新的变化。运行期间 store 暂时没有新时间戳，不等于 worker 已经卡死。

查看最近活跃的 crew worker 日志，并持续跟踪每秒一次的 CPU/RAM 心跳：

```powershell
$workerLog = Get-ChildItem (Join-Path $cacheRoot 'logs/targets-crew/workers/crew_log_*.log') |
  Sort-Object LastWriteTime -Descending |
  Select-Object -First 1
$workerLog.FullName
Get-Content $workerLog.FullName -Tail 5 -Wait
```

`__AUTOMETRIC__` 行中的 target 名应与 `meta/progress` 一致；target 名称之前的数值依次对应 `autometric::log_read()` 的 `core`、`cpu`、`resident` 和 `virtual` 字段。日志每秒新增且 `core` 长期接近 `100`，表示一个 CPU 核心仍在持续计算；日志停止更新、worker 进程消失或出现 error 时，才需要按失败或中断继续排查。主进程在等待 crew worker 时 `core` 接近 `0` 属于正常现象。

查看主进程心跳：

```powershell
Get-Content (Join-Path $cacheRoot 'logs/targets-crew/main-process.log') -Tail 5 -Wait
```

结束持续跟踪请按 `Ctrl+C`；这只退出日志查看，不会停止 `tar_make()` 或 crew worker。

输入 RDS 至少应包含 `object`、`data`、`metadata`；兼容旧字段名 `resCCS_ablation`、`data_all`、`ablation_metadata`。如果输入没有预计算的 biology/structural 输入，需通过 `CCS_FULL_EXPRESSION_RDS`、`CCS_GENE_SIGNATURE_RDS` 等环境变量提供同一数据边界所需的原始资源。

## 轻量验收

轻量验收不再是另一种 profile，也不减少 repeats、bootstrap、scaling、null controls 或 decoder。只需把 `-InputRds` 换成正式数据子集或同 schema 的小型 fixture，并使用相同 launcher、相同参数和相同 targets 图：

```powershell
$subsetInput = 'D:/path/to/existing-input-subset.rds'  # 换成已存在的同 schema 子集
$testRoot = 'tmp/tests/depth-seven-tissues'  # 相对于 test/ablation-03
& powershell -ExecutionPolicy Bypass -File 'test/ablation-03/scripts/run-targets-renv.ps1' `
  -Action make -InputRds $subsetInput `
  -CacheRoot $testRoot -OutputRoot $testRoot
```

样本减少可能使某些端点变为 `not_estimable`，但不得改变分析方法、target 图或参数语义。

## 环境与包边界

ablation-03 只依赖已安装的 CCS 包，不从仓库 `source()` `R/ablation.R`，也不使用 `pkgload::load_all()`。修改 CCS 源码后，必须先按项目版本门禁用 `C:/R/R-4.3.1` 构建、检查并安装 `DESCRIPTION` 指定的版本，再运行 targets。

targets store、科学产物和 observability 均位于显式 cache root；不得把大型 RDS 或个人/专有基因组数据提交到仓库。
报告渲染产生的 JPG 预览位于 cache root 的 `logs/report-previews/`，正式 HTML
与 PDF/CSV 图表材料仍分别位于分析项目根目录和 `reports/`。
