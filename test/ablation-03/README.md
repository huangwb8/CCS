# ablation-03：表示性能、gene signature 评估与结构可复现性

ablation-03 是独立的可复现分析项目。`_targets.R` 是唯一的分析编排入口，负责依赖图、缓存、恢复、运行状态和 `crew` 并行调度。六个编号 R 文件保留为科学计算实现，由 targets 节点在同一 R 进程中消费；它们不再拥有 stage lock、stage receipt、SUCCESS 或手工 cache hit/miss 协议。

formal 分析与轻量验收使用相同的 targets 图、CCS 版本、参数、随机种子和科学代码。测试子集重新计算生物输入；正式输入还包含预计算的生物和结构输入，小测不验证该分支。

## 分析单元

| 单元 | 科学职责 | targets 产物目录 |
|---|---|---|
| `01.01.00` | 输入对象、表达数据和元数据准备 | `01-data` |
| `01.02.00` | Direct/d1 表示输入与样本契约 | `01-representations` |
| `01.03.00` | gene signature 与结构分析输入 | `01-biology` |
| `02.01.00` | 几何、检索、readout、learning curve、scaling、decoder | `ablation-experiment` |
| `02.02.00` | gene signature 效用与队列级推断 | `ablation-biology` |
| `02.03.00` | 双向结构复现与 matched-bank 敏感性 | `ablation-structural-reproducibility` |
| `02.04.00` | decoder、局部 query、共享节点及固定训练设计的补充推断 | `statistical-inference` |

R Markdown 文件只消费已完成 target 产物，用于报告渲染，不参与科学计算调度。
表示报告另消费 `references/` 中 Figure S3 的历史 cohort Gain 汇总与来源记录，用于解释 d3 坐标重建 importance 与原始 d1 方差份额的区别。该输入不含患者级矩阵、不重拟合历史模型，纳入 `representation_report_sources` 文件依赖；当前方差份额仍来自正式 `representation_analysis` 产物。
四份 HTML 均由正式依赖图中的 file target 直出；不得在 `tar_make()` 之外用独立
runner 生成正式报告。
报告共用 `templates/` 下的 Liquid Glass 样式与交互模板；模板文件纳入四个报告的
targets 文件依赖。图表代码输出到 HTML 并默认折叠，可点击各代码块的 Code 按钮展开，
或用报告顶部的 Code → Show All Code 一次展开全部代码；setup 初始化代码不输出。
发表用图仅保留坐标、图例、分面名称及必要的数值标注；绘图符号解释、统计口径、
分析限制与解读写在正文或图注中，不嵌入图片的 caption 或 subtitle。
窄分面数值横轴采用稀疏刻度与适当小数位，并保证共同标尺不变；高覆盖队列图
加宽画布，避免五列刻度挤在一起。同预算模块比较保留全部整数设计点，横轴标签
旋转 45°。图表更新须核对正式 PDF 与 HTML，避免刻度或分面边界处的标签重叠。
gene signature 报告附带配对效应、cohort 异质性、基因覆盖与队列内区间四张 PDF，并在
`reports/tables/` 导出队列级配对差值，供核对图中每个格子的样本数与方向。
报告开篇、效应图解读与讨论明确核心结论：d1 改变了局部生物表达关系，但平均匹配
代价较小；对以跨队列表征为目标的 CCS，这是可以接受的方法取舍。效应量与队列
异质性共同限定这一判断，接近度变化不等同于分类准确率损失或正式非劣效结论。
报告“数据概览”从正式样本契约、表达缓存及 top-15 邻居名单动态列出 reference/query
的评分、检索和推断角色，展示来源与检测类型、全部队列及逐 gene signature 有效 query 数；
区分共享外部候选与癌种 readout 合格子集、严格诊断和低覆盖探索的资格。
同源汇总保存为 `02.02.00. cohort-sources.csv`、`02.02.00. cohort-usage.csv` 和
`02.02.00. query-usage.csv`，仅含队列元数据与计数；实际选中的 reference 数按不同
样本计数，不等同于邻居对次数，报告输入由 `biology_report_sources` 跟踪。
配对效应图采用紧凑森林图与独立统计列，逐行显示差值及 95% CI、有效 cohort 数、
配对 query 人数和 BH 校正 P 值；统计列不参与效应横轴缩放，gene signature 顺序沿用冻结配置。
数据概览与 gene signature 报告的五张正式 gene signature 图共用 `scripts/helpers/anchor_plot_labels.R`，
英文图统一显示 Proliferation、Immune TME、Stromal TME、IFNγ、IL6–JAK–STAT3；内部键为 proliferation、immune_tme、stromal_tme、ifn、il6。
队列内区间图直接消费 `statistical_inference$biology_local`，按 gene signature 分为五列，展示各
cohort 配对均值与点式 95% query bootstrap 区间，并导出同源 CSV；区间图与明细表
位于“队列间效应分布”的热图之后，连续核对同一格子的均值与不确定性。该区间未作同时
覆盖校正，不替代 cohort 间主效应区间，也不提供单格显著性结论。
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

高覆盖 gene signature 验证接入同一 DAG；正文以广队列基线、高覆盖验证及必要机制线索为主线，原 43-cohort 基线与各诊断产品继续保留。
设计冻结在 `config/biology-high-coverage.yml`，分析函数放在本分析的
`R/biology_high_coverage.R`，I/O 由 `targets/biology_high_coverage.R` 负责；直接复用
已安装 CCS 的检索、reference-only 读取器和 cohort 配对推断，不修改公共包 API 或版本。
`biology_high_coverage_inventory` 只提取测量信息，frontier/contracts 不依赖效应或读取器结果。
配置另放 `config/`，保持已有 `raw/` 只读。

每个 gene signature 在 60/70/80/90% 档独立搜索可兼容队列，冻结原 gene signature 共同基因及背景。
确定性搜索用基因可用频率排序生成全体、单队列及两队列交集种子，依次最大化 query
队列数、已知来源组数、query 数、共同 gene signature 数和 reference 队列数；这是有界搜索，
不宣称找到全局最优集合。主档优先选满足至少 8 个 query 队列、3 个来源组的最高覆盖档，
否则按预设 5 队列及 3 队列边界报告有限验证或敏感性证据。所有档完整保留。
每个 query 队列至少 20 人，同癌种 candidate pool 至少两个 reference 队列、15 个完整样本。
秩背景至少 3,000 个共同基因且样本背景完整；门槛不足保持不可估计。
实际排名使用完整的共同背景交集，不截断为 3,000 个基因。

高覆盖产品写入正式 cache root 的 `biology-high-coverage/`；患者级矩阵、邻居和读取器模型
留在该受控缓存。报告导出 `reports/tables/02.02.00. high-coverage-*.csv` 的聚合表和
`high-coverage-*.pdf` 图，包括测量 frontier、主效应、覆盖率稳定性、效用／读取器机制、
队列异质性，以及留一队列、来源组成／删除、主比较与共同队列对照、检索规则差中差
及有效人群支持五张敏感性图。原始明细表保留；空心点只表示描述性估计，缺少合格
配对的条件不在效应图中补零，其队列／患者数与资格完整显示在支持图。主效应与覆盖率效果图仅展示效应及区间均可估计的结果，省略无结果的行、
分面和人数标签；报告正文说明省略原因，完整状态与资格仍保留在聚合表中。
队列异质性图共享队列纵轴，仅在左侧显示名称，各 gene signature 按相同队列顺序对齐。
报告将“高覆盖 gene signature 验证”作为独立章节，紧接广队列基线；正文用一段说明全队列交集资格不足为何转向高覆盖队列，明确两种测量合同的区别。邻居重排与模块距离贡献留在正文“患者相似性重排的机制线索”中。
五个选定 gene signature 固定构成主检验族，缺失 gene signature 仍计入 BH 的五项分母；coverage
档的其它比较按完整五项族校正。完整表以 `BH_family` 逐行标明口径：主比较行沿用各 gene signature 选定档位组成的五项主族，其余行使用对应覆盖档、评分及候选池的五项族；同一覆盖档显示的 q 值可能属于不同族。95% cohort bootstrap 条件于固定 bank、合同及 atlas；没有独立
冻结的非劣效界值，不据 CI 跨零宣称无损。距离及技术限制使用配对患者差中差，留一及
C_core 只作预设敏感性估计。

`scripts/tests/biology-high-coverage-targets.R` 在隔离 `tmp/tests/high-coverage-*` 中用合成
数据执行同一 DAG 与报告，并核对效应列变化不改变合同、缺失／未知癌种、来源重复、
reference／秩背景不足、原正式文件不变及下游失效后的上游复用。该验收默认创建唯一
run root，并通过 `TAR_CONFIG` 使用隔离 YAML，不切换正式 `tar_watch()` 的 store。
验收的 `renv` sandbox 同样放入该 run root，避免并发 R 进程等待共享 sandbox 锁；
项目包库与锁文件保持原有配置。
第二个参数可指定已有隔离 run root，供同一 fixture 的增量恢复。
从仓库根目录执行验收，并使用 `--vanilla` 避免 `.Rprofile` 在路径隔离前激活环境：

```powershell
& "C:/R/R-4.3.1/bin/Rscript.exe" --vanilla test/ablation-03/scripts/tests/biology-high-coverage-targets.R
```

连续 gene signature 根因诊断由 `config/biology-diagnostics.yml` 冻结设计，计算层为
`R/ablation_biology.R`，同一 `_targets.R` 的 `biology_diagnostic_*` targets 调用
已安装 CCS。诊断结果写入 cache root 的 `biology-diagnostics/`，以
1e-10 容差核对广队列基线格子的效应和人数；冻结基因清单为独立 file target。
诊断节点直接依赖科学代码及共享标准化函数的身份；表示输入、gene signature 缓存和
基线文件均以 file target 跟踪实际内容，避免审计结果相同时截断必要重算。

诊断采用共同原尺度／样本内秩评分与全癌种／同癌种候选池四格，配套 reference
cohort 五折连续读取器及两个单因素 d1 距离对照，主评分为秩。严格秩背景至少
3,000 基因，高覆盖验证采用相同门槛；原尺度无背景门槛。直接输入不重叠、技术来源和对称留一
为严格补充。低覆盖探索产品独立保存，使用固定残余清单和有效 reference 标尺，
要求样本 gene signature 完整、秩背景完整及两臂 top-15 完整。多基因仅发布满足人数门槛的
cohort 条件区间，不发布 p/q；基质单基因只描述，不训练低覆盖读取器。
gene signature 保留原 50% 且不少于 8 个基因；同癌种候选池至少两个 reference cohort
和 15 个样本。Annoy 严格
ID 召回门禁为总体 0.95／每 cohort 0.90，不达标升级预算或精确检索；检测到
边界等距时按精确距离和 ID 决定名单。未知平台不作为真实平台。

同一比较固定共同患者／cohort，至少 20 个 query 才参与跨 cohort 推断。
2,000 次 bootstrap 至少 1,900 次有效、至少三个 cohort 才给 95% 区间；
相同有效 cohort 集合共享抽样索引，不同指标分别保存有效集合。MAE 差为
d1−Direct，正值较差；utility 差正值较好。XGBoost 只由预设 reference CV
规则触发，两臂同时运行，不使用外部 query 调参。输入不重叠审计采用完整
冻结 geneSet 的保守并集，不能据此宣称独立临床验证。

全队列共同测量资格、不可估计比较、读取器状态与检索核验完整保存在“附录：全队列共同测量审计”；固定残余基因的图表与解读另放“附录：固定残余基因的低覆盖探索”。这些状态属于全队列交集合同，不作为高覆盖验证失败的证据。全部原 R chunk、图表、科学参数与 targets 依赖继续保留，调整只涉及报告位置与解释。

生物报告在全队列共同交集效应可估计时输出评分候选池与读取器两张 PDF；低覆盖探索有有效
点估计时另输出 diagnostic-low-coverage-score-pool.pdf，并将全部预定评分／候选池与距离
切换差中差分别绘制为 diagnostic-low-coverage-score-pool-changes.pdf 和
diagnostic-low-coverage-distance-changes.pdf，保留实际人数、条件区间与单基因仅点估计。
测量门槛、邻居重叠与模块平方距离累计份额分别输出 diagnostic-measurement-gates.pdf、
diagnostic-neighbor-overlap.pdf 和 diagnostic-module-distance-concentration.pdf；这些描述图
不以严格效应可估计为前提，不把支持不足补为零效应。模块份额先在 cohort 内归一化、
再 cohort 等权汇总，各规则重新选择邻居，不作为固定邻居的纯权重分解。上述图经正式
biology_report target 更新 HTML，沿用原诊断产品与门槛。另补充候选池人数、严格比较支持矩阵、
读取器分支状态、残余代理可评分人数及检索严格召回五张诊断图，保留原表。
不可估计与未运行以状态展示；支持计数的零不当作零效应，人数图不替代推断资格。
检索图分别汇总 query 平均与最低队列平均，实际门禁仍按各候选组执行。导出
`reports/tables/02.02.00. diagnostic-*.csv` 聚合表，低覆盖推断、cohort、标尺和
有效人数分别保存在 diagnostic-low_coverage_*.csv，
包括分支状态、测量／候选池资格、推断、模块相关和单臂距离。患者级矩阵、
名单、读取器及共同 query 清单仅留在受控缓存。未知预处理保持未知。
邻居重叠图上下排列全癌种与同癌种实验，各组只保留实际有结果的条件，共用 [0,1]
横轴。行标签注明有效队列数，小点保留全部队列均值，菱形及短横线显示队列中位数
与四分位范围，右侧给出相应数值；该范围为描述性分布，不是置信区间。图中统一按 d1 与 Direct
标注；Jaccard 是对称的名单重叠度，不是有方向的效应差。解读区分局部关系重排、
生物效用与信息读取，说明不同候选池／人群间的重叠中位数不能直接作因果比较，
并提出共同患者配对及候选池大小参照的后续验证。
模块距离集中度图配套面向初学者的坐标／参照线说明，并从同源份额展示前 10% 模块
的贡献及累计一半距离所需模块数，解释组织等权为什么不等于实际模块贡献均分；
保留重新检索、不同有效队列和纯描述性结果的边界。
`tests/44-test-biology-diagnostics.R` 验证科学不变量；
`scripts/tests/biology-diagnostics-targets.R` 使用合成数据执行同一 DAG 和报告，
并在独立 `tmp/tests/biology-*` 中核对恢复、正式路径未变及 store 配置还原。
该验收还分别验证科学代码身份和表示输入文件内容的失效传播，恢复后保留
原始字节；可用第二个参数指定已有隔离 run root，让 targets 增量恢复。

五项 gene signature 的共同覆盖、实际人群及效应由生物报告动态展示。报告直接说明当前配置、测量支持和结果；分析变更记录见 Git 历史与 CHANGELOG。

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

## Gene signature 定义与结果

分析使用增殖、免疫微环境、基质微环境、IFNγ 和 IL6-JAK-STAT3 五项 gene signature，定义见 `config/biological-anchors.yml`。IFNγ 与 IL6-JAK-STAT3 分别按对应通路基因集评分，天然共享基因在各通路中保留。广队列接近度、高覆盖接近度与读取器采用各自固定的五项 BH 校正族，低覆盖探索不发布 p/q。数据概览、生物与结构报告由同一 targets DAG 生成。各项效应、95% CI、校正结果和不可估计原因见生物报告。
