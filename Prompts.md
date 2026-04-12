# 日常

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

