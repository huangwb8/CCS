# 日常

---

基于 E:\iProjects\RCheck\GSClassifier\routine02\docs\plans\2026-08-01-ccs可新增消融实验计划.md 、 docs\plans\2026-08-01-metaccs-end-to-end-ablation.md 这2个计划，优化ccs包的 R\ablation.R ，进一步完善消融实验。

---

我准备在CCS包里添加一个新的`ablation.R`脚本（暂时放在E:\RCloud\RFactory\ccs\test\ablation里，等我测试完毕后，我再手动放在 ccs 的正式R脚本文件夹内），它里面设计一个新的函数`ablation`，它基于`CCS`这个Class进行设计（类似ccs包里的很多其它关键函数），专门用来落实 E:\iProjects\Manuscripts\CCS\paper\docs\plans\2026-07-27-ccs框架分层消融实验计划.md 里的消融实验。我希望，当用户运行函数`ablation`时，对于给定的CCS对象（以及必要的其它数据），该函数可以自动、准确、可重复地完成消融实验。这是我的大致想法：

- 本轮任务的工作目录是： E:\RCloud\RFactory\ccs\test\ablation，你可以在它及其子文件夹里生成、编辑文件。在此之外的文件，都是只读以保证安全。
- 你要设计好主函数和子函数的关系；主函数`ablation`放在最前面，子函数放在后面。不同的子函数代表不同的消融实验；这样以后我想到新的消融实验的时候，直接加一个子函数模块就行，非常方便。子函数也可以用来托管一些较通用的功能。总之你自己把握好。
- 和目前ccs包的函数、定义相协调，在函数的参数上尽量不要重复地造轮子，能复用就复用；实在没有再搞新的
- 函数`ablation`的参数的数量恰到好处就行，不要冗余。
- 不要动ccs包已有的函数体系；而是主动兼容、补充到目前ccs包的代码体系（`ablation`应该相当于ccs包的一个“挂件”）。代码风格上尽量同我的接近，以便我后期可以比较简单审查你写的代码。
- 加载"E:\RCloud\RFactory\ccs\test\ablation\ablation-test-data.R" （只读）后可以获得一个真实可用的CCS对象和训练样本的RNA原始表达量，里面有一些实际的数据。你运行消融实验代码进行测试时，可以利用这些数据；但千万不能编辑、改动这些数据的原始文件，以保证安全。

---

目前，CCS包有一个专门的github仓库 https://github.com/huangwb8/ccs.principle 提供教学功能；因此， CCS仓库的README要：

- 引用这个仓库
- 重点突出的同时尽量简洁；重点参考： https://github.com/huangwb8/GSClassifier 和 https://github.com/huangwb8/GSClassifier.principle
- CCS的论文正在投稿中。推荐引用文献应该留空，提示`Under审稿`之类的。 

# p1


## 函数

请优化`R/phenotype.R`中的`subtypePerformance.R`函数。我希望增加一个新功能：元分析（meta-analysis）汇总各分型/亚组的响应率：报告合并效应与 95% CI，使用随机效应模型（DerSimonian-Laird 或 REML），并给出 I^2（异质性指标）与 p-heterogeneity。新功能由一个新的子函数完成，名为`metaSubtypeRate`。

在`# Forest plot`和`# RR/NRR - scatter plot/box plot`之间插入，备注为`# Meta analysis`。如果需要添加新的子函数，请加到`newSubtype`这个子函数的后面、`####%%%%%%%%%%%%% Assistant functions %%%%%%%%%%%%%%%%%%%%####`的前面。

请注意，不要改动原有代码的主要部分，那些代码经常长期测试和验证，没有问题。

## 参数

- `data`：data.frame。包括"Cohort", "nResponse", "response_rate", "size", "Subtype", "tumor_type"等列（相当于原程序里的`df3`）。其中：
  - Cohort: Character。队列名
  - Subtype : int。分型对应的数字。
  - tumor_type：Character。肿瘤类型
  - size： int。特定的Cohort、Subtype的总样本量。
  - nResponse： int。特定的Cohort、Subtype的里reponse=1的样本量。
  - response_rate：float。等于nResponse/size。
- `method`：默认是随机效应模型。用于控制进行meta分析的模型方法
- `...`：一些其它的你认为比较重要的参数；比如控制图像的某些属性。

## 分析

- 要有一个meta分析的图，不同的Subtype（列）在不同的Cohort里（行），有一个response rate相关值及95%置信区间。
- 报告合并效应与 95% CI，给出 I^2、p-heterogeneity或其它比较重要的指标。

## 输出

一个List。包含下列子List：

  - `Data`：data.frame。和作图直接相关的数据。
  - `Plot`：ggplot2（或者其它更加适合meta analysis的图）。一个meta analysis的图。要求美观、大方、整洁，可以在Nature/Science/Cell等杂志上发表。

在最后的`l`的部分里，原本是这样：

```
Plot = list(
      ScatterBarPlot = plot_r,
      ForestPlot = plot_f,
      UtilityPlot = plot_utility
),
Data = list(
      ROC = data_roc,
      Normalization = data_norm,
      ClinicUtility = data_utility
)
```

现在要添加新的图和数据：

```
Plot = list(
      ScatterBarPlot = plot_r,
      ForestPlot = plot_f,
      MetaPlot = plot_m,
      UtilityPlot = plot_utility
),
Data = list(
      ROC = data_roc,
      Normalization = data_norm,
      MetaAnalysis = data_meta,
      ClinicUtility = data_utility
)
```

