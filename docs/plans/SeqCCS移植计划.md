## 测试环境

### 数据

```R
project_version <- 'PADv20240911'
data_version <- 'PADv20240810'

# Data
data_all <- readRDS(paste0('E:/iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_',data_version,'.rds'))
PAD <- readRDS(system.file("extdata", "PAD.train_20220916.rds", package = "GSClassifier")) 
geneSet <- PAD$geneSet
geneAnnotation <- PAD$geneAnnotation

# Parameters
ccs_params = list(
    device = "cpu",
    nfold = 5, nrounds = 50, nthread = 2,
    eta = 0.3, subsample = 1, 
    gamma = 0.05, alpha = 0.05,
    max_depth = 6, colsample_bytree = 0.8, min_child_weight = 1,
    n = 50, sampSize = 0.8, ptail = 0.1, nround.mode = c("fixed", "polling")[1]
)
```

### 测试空间

所有的测试过程和结果必须在这个文件夹下进行：

```
E:\iProjects\RCheck\SeqCCS\test01
```

## 目标

### 充分理解SeqCCS

在`E:\PythonCloud\Package\OmicEncoder`里有SeqCCS最原始的代码；尽管仍有优化空间，但是基本代表了SeqCCS的核心思想。你要充分阅读、理解SeqCCS的原理和实现过程。 OmicEncoder 是一个旧名，在新的测试或实现里不应该再使用。

### 充分理解 CCS

在`E:\RCloud\RFactory\ccs`里托管着CCS包的源代码。你要充分阅读、理解CCS的原理和实现过程。

这是CCS的核心过程：

- `data_all` 的不同肿瘤下有不同的cohort的数据（定义为`data.i.j`）
- `data_i` 里会定义好expr和subtype，然后基于它们， ccsSubMode函数就可以构建多个关于`data.i.j`的子模型`model.i.j`
- 因为`model.i.j`子模型学习了`data.i.j`的数据模式，可以构建一个函数`f(model.i.j, new_data)`用于生成关于`new_data`（可以是1个或多个样本）的声明性向量`V.i.j`（每个样本的声明性向量是独立的），从而实现了`new_data`的富有生物学意义的标准化
- 后续汇聚多个类似的`model.i.j`的`V.i.j`，这样对于任意的`new_data`，都可以有一个巨大的标准化向量`V`。它可以用于后续的降维、聚类等操作

`E:\RCloud\RFactory\ccs\R\ccs_model.R` 里定义了：当给定`data_all`，如何构建多个`model.i.j`（当然还会有其它代码的辅助，但核心的逻辑就在这）。目前，CCS包产生`model.i.j`的方式只有：

```R
# 假代码，仅演示业务逻辑
model.i.j = ccsSubModel_GSClassifier(
    data.i.j,
    geneSet,
    params, params_xg,
    seeds,
    numCores
)
```

它是由GSClassifier包所支持的（源代码在：`E:\RCloud\RFactory\GSClassifier`）。

在 `E:\iProjects\RCheck\GSClassifier\routine01\04.01.PADv20240911_Model training.R` 中，我使用ccs函数生成了多个子模型。然后，在`E:\iProjects\RCheck\GSClassifier\routine01\04.02.01.PADv20240911_Select01.R`中，

###  目前CCS产生`model.i.j`的局限性

目前`model.i.j`要生成，要有两个条件：

- **需要提前构建好data.i.j的subtype**： expr容易拿到，下载原始基因表达矩阵简单改改就行。 但subtype需要依赖特定的领域知识才可以获取，不够通用。这是最大的局限性。
- **需要有一个可以将data.i.j转化为model.i.j的方法**：目前暂时只有GSClassifier的方法；但很显然，我目前的设计有为未来接入更多算法做准备，所以才定义了`ccsSubModel_GSClassifier`这样的中间函数。日后，只要新算法的输出规范和它一样，就可以无缝衔接到CCS的流程里。`GSClassifier`是我很久以前设计的方法。当geneSet变得很大的时候，由于计算缓慢，GSClassifier将难以工作。而且，GSClassifier将data.i.j转化为model.i.j的方法不见得在生物学上、机器学习上是最优的。

### SeqCCS是如何用更好的策略构建`V`

我想开发一个新的方法来生成标准化向量`V`。当然，基于队列的方式生成`V`保持不变；我只是想动其中的关键环节：生成`V.i.j`。我希望，可以在**不需要提前构建好data.i.j的subtype的情况下生成V.i.j**。我的基本想法是这样的：

- SeqCCS会基于基因表达量转化为的rank sequence进行建模。建模时会通过随机mask某条sequence的某些基因，模型优化的目标就是尽可能准确地预测所有这些masked基因的原始位置。由于gene rank蕴含了一些生物学上的客观规律，通过SeqCCS的方法，模型就可以学习到一些生物学知识。
- SeqCCS有很多好处：
  - 基于sequence建模，可以充分地利用深度学习的方法，它们很高效，可以充分利用我电脑的强大显卡
  - 纳入多个异质性样本。而且每个异质性样本提供的rank长度可以不一样。
  - 在训练阶段，geneSet可以变得更大一点，也不用怕计算机资源跟不上，也不用怕速度很慢
  - 在预测阶段，geneSet可以小一点（总之就是不完整），这样模型可以预测完整的rank。因为预测阶段geneSet可以变小，这其实意味着在预测阶段SeqCCS对缺失值不敏感，这是很重要的；因为现实世界里的转录组经常有可能有缺失值，而它们可能会对模型性能有很大的影响。
  - 整个过程不需要显式地构建data.i.j的subtype，避免对特定领域知识的依赖。因此，SeqCCS少了生成`model.i.j`，因为它就是一个很大的模型。

SeqCCS的模型构建好后，输入`new_data`时（会转化为rank sequence），它会生成一个完整的序列。此时，我只要提取某个隐藏层，这个隐藏层可以作为`V`。然后就可以跑CCS包的下游分析，比如降维、聚类等。

### 我想做的事

SeqCCS（原称OmicEncoder）目前是python写的，而目前CCS是一个R包。 我希望你基于R的luz包和torch包重写类似的逻辑以和CCS相兼容。

你可以先在测试空间先构建一些新的函数或者计算流程来实现这一点；先不要动CCS的源代码。等我日后觉得时机成熟了，再更新CCS的源代码。

整个过程里，要严格遵守 bensz-rmd-rules skill 的规范开进行基于R的开发。

## 注意

- 对CCS包、GSClassifier包、OmicEncoder的相关信息、所有原有的代码或者数据都是只读模式，千万不要改，否则有严重后果！