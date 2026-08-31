# Common

- 提交commit

```
用 git-commit skill生成commit。所有更新内容包含在1个commit里，分点罗列清楚即可。
```

- 学习分析

```
我对 test/ablation-02 的分析有一些疑问，想请你解答。 你准备好了吗？
```

# 日常

---

这么约束就行：

- "E:/iProjects/RCheck/GSClassifier/routine01/ccs/PADv20240911/resCCS.rds" 包含了“训练 metaCCS 的样本”和“外部数据样本”的 d1；但是应该没有d2-d3-metaCCS-normCCS那一套东西。 
- E:/Sync/Project/PADv20250720/models/resCCS_5ff3a2de76e6cf902e765e8224f9cb66.rds 仅包含“训练 metaCCS 的样本”，但是包含了d2-d3-metaCCS-normCCS那一套东西。 
- d1的计算是个耗时步骤，因此在  PADv20240911/resCCS.rds  这个阶段，其实就已经先把d1搞好了。 就是为了后续可以更加方便地使用。

但这里也有一个需要注意的点： d1是一个依赖cohort的值。如果cohort的数量增加，d1就会变长。 所以，显然  PADv20240911/resCCS.rds  里的d1是比resCCS_5ff3a2de76e6cf902e765e8224f9cb66.rds里的d1更长。如果在ablation-03的研究里要研究 “外部数据样本”的 d1，其实应该是仅截取由 “训练 metaCCS 的样本” 的那部分cohorts所call的那一部分（当然，剩余的部分也许是有用的；您也可以研究，只要你觉得对本次分析有用）。我的想法对吗？

---

test\ablation-03分析：

- 基于  docs\plans\2026-08-30-ablation-03生物学锚点与表示性质优化计划.md 做分析
- 分析的过程会引用一些不在 E:\RCloud\RFactory\ccs 里的文件；它们都是只读，严禁修改。你修改文件的范围要严格局限于  E:\RCloud\RFactory\ccs 及其子文件夹里
- 分析的时候可能会有一些中间数据已经存在，因此生成中间数据的代码可能不用重复运行。如果你觉得中间数据可能已经变化，应该重新运行生成这些中间数据的代码。
- 在 bensz-rmd-rules 的规范下工作（R环境可以用： C:\R\R-4.3.1）。 要出最后的html报告

---

这里我补充几点：

- E:\Sync\@Analysis\PanCan_Data\Level 1\PanCan_CancerSample_DataListForCCS_GEO+cBioPortal+UCXCXenav20240809.rds 是完整的矩阵；你可以看目前的数据里对应什么样本；然后回去这里找，肯定可以找到。 因此， 增殖、免疫/TME、PAM50、CMS 等生物 anchor 是可以做的
- "E:\RCloud\database\Signature\report\GeneSignature-HWB.rds" 是一个我之前整理过的已经发表的一些Gene Signature；对于你研究生物 anchor是有好处的。 你可以搞。
- 生存、治疗反应、临床结局这个在之前的讨论里已经说过了。 到了临床表型这一步，中间的影响因素太多，没有办法真的很好地研究。因此，还是把目标局限一下，能回答到“cohort-based representative对生物anchor的表征影响”已经足够了。从另外一个角度来说，临床表型有很多，你总是没有办法一一研究的；因此，关注“cohort-based representative的生物学意义”就足够了。
- d2/d3/metaCCS 稳定性的研究以前早就说过了。它们都是d1的衍生物，所以并不是研究重点。我们重点还是研究透d1的性质。
- 其它的探索方向你看着哪些合适的，尽量做一下；如果不能做就不要勉强。就算证据不是很强，有一些探索性的分析也可以接受；但是能做好的还是尽量做

本次优化将开启 ablation-03 这个分析，它是在02的基础上进一步优化。因为02分析虽然还不太全面，但是它基本上说明了“d1与Direct-GSClassifier明显不同”，也算是一个阶段性成果，因此我不打算覆盖它。03的分析应该包含02，然后有更多更好、更全面的分析，意义也要有所推进。如果你需要改CCS包的源代码，也可以；但我个人的建议是仅改 R\ablation.R ， 其它的不要大改；而且尽量保证对 ablation-02 的兼容（这个很容易；03主要是有新的分析，你直接添加一些新的部件，然后旧的部分不改就行了；就是这种“模块化增强”就可以保证兼容性）。

结合上述意见，你总结一下，写个优化计划在docs/plans 里。

---

test\ablation-02分析优化：

- 基于 docs\plans\2026-08-29-ablation-02-statistical-inference-optimization.md 优化源代码
- 重出rmd、html报告

---

你还是有一个小节；专门搞一个表格来做这件事吧，不要在正文里零散地做。 表格的任务：

- 指标
- 取值范围

- 意义：包括趋势变化意味着什么
- 通俗解读： 就是小孩子都可以看懂的

---

test\ablation-02分析优化： 

- 在  2fc352670840424b3b8961b3d5f2a50f568c5524 这一commit的基础上，基于 docs\plans\2026-08-08-ablation-02-cohort-bank-scaling优化计划.md 继续优化。

注意事项：

- 除了  test\ablation-02 、R\ablation.R 和 docs\关于CCS框架的讨论.md 可以修改，其它位置的文件都是严格只读的。

---

test\ablation-02 优化： 

- 做分析的时候，应该要先解释为什么要做；然后再开始分析、结果解读。目前解释“为什么”这一方面很薄弱
- 请在不改变分析代码的前提下，加强为什么要这么做。不然客户看不懂。
- 很多观测指标都比较抽象；在使用指标的时候，要用幼儿园小孩都能看懂的方式来解释这些指标的价值，深入浅出。

---

另外，我认为你在分析cohort bank效应的时候，大基础你没有把握住。cancer-type balanced accuracy、d1–GSClassifier TSP的变化（每翻倍斜率）、macro-AUROC变化（每翻倍斜率）等指标当然可以计算；但是，因为原始的d1本来就是基于150个训练集cohort module训练出来的，比较的时候，肯定很容易出现在150个训练集cohort module附近（你分析的时候是125-150）就出现平台期。这几乎是一个肯定的事，没什么特别。我觉得你的关注点要集中在scaling —— cohort scaling 对于 d1 的某些优势特性（在 2fc352670840424b3b8961b3d5f2a50f568c5524 里的 test\ablation-02 已经有一些分析的方向 ）的 scaling 是否存在？而且，scaling其实也有2个维度——会增加tissue多样性的scaling和不增加tissue多样性的scaling（即在已经存在的组织里单纯地增加cohort）对d1有什么影响？ ——也没有分析。你怎么看？

---

而且， cancer-type balanced accuracy 这个已经说过很多次；CCS框架并不是为了预测cancer-type，所以，这个balanced accuracy较GSClassifier-TSP低也不能说明CCS框架做了不好的事情。

---

test\ablation-02分析优化： 最近的cohort bank的分析，用的还是DIRECT TSP而不是 Direct-GSClassifier TSP，这个口径和之前的分析不一样。 我觉得还是要统一口径，采用之前用的基于GSClassifier计算的TSP。你的看法是？

---

test\ablation-02分析优化： 最近的cohort bank的分析，用的还是DIRECT TSP而不是 Direct-GSClassifier TSP，这个口径和之前的分析不一样。 我觉得还是要统一口径，采用之前用的基于GSClassifier计算的TSP。

---

test\ablation-02分析优化： 在目前测试数据里的32 gene signature下，评估随着cohort scaling，d1 vs.  Direct-GSClassifier TSP  的差异的变化趋势如何？这个问题还是很有意思的，因为d1是一个随着cohort会scalling的向量；但scalling后会加强还是削弱d1的优势暂时是未知的。我和reviewers们都对相关的趋势很感兴趣。你要优化分析代码以完成上述评估。请你：

- 必要时可以优化 R\ablation.R 以达成目标（目前它有 "scaling" 的模式，但可能有缺陷或未激活，不一定可以完成任务）。
- 最后记得在 docs\关于CCS框架的讨论.md 里记录你关于这方面的想法，特别是你的推测：如果gene signature不是32而是更大，比如2000之类的，你觉得相关的趋势会如何演化？

注意事项：

- 除了  test\ablation-02 、R\ablation.R 和 docs\关于CCS框架的讨论.md 可以修改，其它位置的文件都是严格只读的。

---

test\ablation-02分析优化： 我之前说， 代码重构一下，增强可读性。 可以有序号（比如 01 02），必要时可以有多个html/Rmd。这里的序号是指文件名要有序号。这样方便我知道代码的先后顺序。

---

test\ablation-02 分析优化： 增强代码的可视化

- data_all里，有一个组织类型是 undefined 。其实，里面有一些数据集是有真实的tissue信息的；你可以以 test\pre-train-info 的数据为准，调整一下。 
- 目前的代码可能能跑，但是可读性比较差。 我觉得在 bensz-rmd-rules 的规范下，你的代码要像“讲故事”一样：
  - 原始数据有哪些？它们有什么特点？可以做什么分析。
  - 原始数据的初始准备、分析后，会有一些重要的结果数据。它们具体是怎么算的？这些数据有哪些内涵需要向用户特别汇报？
  - 有了数据，就需要做可视化。 类似 test\ablation-02\ablation-experiment.Rmd 里的功能；最后出html。
  - 其它你觉得比较优秀的展示
- 代码重构一下，增强可读性。 可以有序号（比如 01 02），必要时可以有多个html/Rmd。
- 注意事项：
  - 除了  test\ablation-02 里的文件和 R\ablation.R可以修改，其它位置的文件都是严格只读的。

最终目标：让人类专家更加容易审核分析过程；让结果更适合直接用于发表；让分析更加详尽、通俗。

---

test\pre-train-info 分析优化：

- 结果里，有一些 tissue=undefined 的。其实，这些是历史遗留问题，因为当时由于某些原因，没有特别指定它们的组织类型；但它们其实也是有组织类型的。现在，在本分析中的xlsx结果文件中，我希望你给它们分配真实的tissue类型。

---

癌种结构同时被削弱，我感觉是预期内的，并不是真的是坏事。一般来说，CCS框架在设计的时候，其目标并不是为了还原癌种结构；毕竟，还原癌种没什么意思，也没什么意义；而且，我觉得就算可以还原癌种信息，也不是什么了不得的成就。我打个比方：一个样本的GSClassifier-TSP，肯定是会包含一点癌种特异性的信息在里面的；但是，我觉得CCS框架有点像是利用多癌种的信息构建一个“原始肿瘤”的虚拟元肿瘤，然后各种肿瘤是它的亚型； E:\iProjects\Manuscripts\NSFC_Young_2024 里有描述过这种思想。这个想法，其实并不以还原旧的癌种为目标，对吧？而且CCS的框架，其实是特征信息 + 队列信息的混合，所以，它不完全是TSP的信息，肯定是有队列的信息在里面，所以癌种结构被削弱是预期内的，不是一个坏事。 这也一定程度上可以和 E:\iProjects\Manuscripts\CCS\paper 里的口径统一，CCS框架就是想把feature和队列信息融合在一起，以方便构建全新的分型；这也是为什么“少量feature也可以scaling”的基本逻辑。你再调研一下，帮我分析一下，评估一下目前测试结果的解读要不要变。 当然，如果你有什么新的想法，也可以写。 最后，解读要更新到 test\ablation-02 的Rmd、html报告里。

---

根据  docs\plans\2026-08-05-ccs-d1-vs-tsp消融源码优化计划.md  优化 R\ablation.R。你可以在 test\ablation-02 里进行测试，保证 R\ablation.R 符合预期且可正常运行。你优化下 test\ablation-02\ablation-test-data.R ，引用  test\pre-train-info\DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810_BatchInfo.xlsx 这个数据，它对本次任务有重要作用。基于新的 test\ablation-02\ablation-test-data.R，根据新的 R\ablation.R（你直接在实验代码里 source 这个脚本；而不是构建新的CCS包；因为目前主要 是测试），完成消融实验，要求图文并茂、讲解清楚， 符合 bensz-rmd-rules 的规范。注意：

- resCCS 里应该有一些filtered cohorts，它们正好可以作为独立的外验证数据，你要好好利用。本次任务的最高目标： 基于  ablation-test-data.R 里的实际数据把预期的消融实验做到极致。
- 除了  test\ablation-02 里的文件和 R\ablation.R可以修改，其它位置的文件都是严格只读的。

---

根据 docs\关于CCS框架的讨论.md 里，其实最重要的问题，就是 d1 vs. TSP有什么区别。但是，如果观测指标拓展到d1的下游—— d2/d3/metaCCS/normCCS/ICI response，其实意义不大，因为这些部分涉及的超参数特别多。 比如，2个分组——d1 vs. TSP，其它参数基本一样，然后最后ICI reponse这里 d1 比 TSP好，也无法说明d1真的好，因为我可以说——TSP也许不是最佳的超参数组合，所以比不上基于cohort的d1。就算你设计再多的实验，也很难比较严格地证明这一点。所以，我觉得最重要的消融实验，应该把重点放在 d1 vs. TSP有什么区别 、d1改善了什么、d1妥协了什么，往这个方向去思考；当然，你有更好的意见。 你综合 docs\关于CCS框架的讨论.md、R\ablation.R 的设计和E:\iProjects\Manuscripts\CCS\paper\docs\plans\2026-07-27-ccs框架分层消融实验计划.md （这些文件有好的、也有不好的方面），写一个 R\ablation.R 的源代码优化计划，告诉我要怎么改良它。我们的目标就是：最符合CCS框架这个科研场景的最优化消融实验方案。 

---

预训练时我用的数据包是： E:/iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810.rds （只读） 。E:/iProjects/CCS_Data （只读） 里有生成这个数据的代码及对应的原始数据。 我希望你按图索骥，通过联网确认、适当编程或其它可靠、安全的方法，最终获得 DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810.rds 里的数据集的 docs/关于CCS框架的讨论.md 里 `## 证据层预训练语料的最小数据要求` 小节所描述的数据，以方便我进行后续分析。我的想法大致如下：
- 工作目录是： test/pre-train-info 。在本目录以外的所有文件，都是只读，严禁修改。
- 确认 `## 证据层预训练语料的最小数据要求` 小节所描述的数据的策略是多样的；有时候原始数据里可以直接获得；有时候可以需要联网找一下（比如，GEO数据集你联网找一下一般是有相关的数据的）；有时候可能需要通过某些API。一般来说，数据来源于GEO、cBioPortal和UCXC Xenav。你要灵活运用多种方法。
- 要求：信息保证准确；如果不确定，就用unknow代替。不能伪造数据。
- 一个一个数据集串行地搜索、确认，不要并行。每找一个数据集，关键中间数据要留痕；这样就算任务偶然中断，也可以继续下去。不管数据的来源如何，最终关键的中间数据要使用你自定义的规范方法来暂存。最后，直接用你开发的代码它们提取、merge，直接生成我需要的矩阵。
- 本次任务的出品：
  - ./test/pre-train-info/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810_BatchInfo.xlsx ：包含一个矩阵，记录我想要的上述数据信息。每一行是一个样本；每一列是一个信息。每个样本的ID、样本所属于的cohort名、tissue名，要和E:/iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810.rds的标签对齐；不要创造新的标签。
  - ./test/pre-train-info/human-should-check.md：系统地记录有哪些数据是不太确定的、需要人类进一步审核。
  - 其它重要的R脚本：工作时可以不遵守 bensz-rmd-rules 的规范，因为这不是一个标准的分析流程。可以按你自己的意思来。

---

设计CCS框架最初的灵感大致是：
- 在肿瘤学里，基于RNA的预测模型在临床转化过程中面临一些问题：（1）我们需要一个好的泛癌个体化转录组模型，因为它有希望解决一些小众癌种的特定药物（比如ICI）的疗效预测问题。但是，构建一个有临床应用价值的泛癌个体化转录组模型，一个重要前提是，这个模型在跨平台、跨临床中心的预测具有稳健性；只有这个“稳健性”得到了保证，这个地基打牢了，一个好的泛癌个体化转录组模型才是可能的。
- 不同的RNA队列会因为很多技术原因（比如属于不同的batch；不同测量平台[比如RNA-Seq或microarray; RNA-Seq也有可能导致]）；不同的标本类型导致RNA的绝对表达量无法被轻易地merge起来做事情（比如构建预测模型）。类似TSP这样的简洁技术很早就出现了——基于相对表达量。一般来说，很多论文、研究，包括我自己的实验，均表明TSP具有较强的跨平台、抗批次效应的能力。虽然绝对表达量到相对表达量的转换过程中，有一部分信息丢失了；但是却有希望把更多的RNA表达队列merge起来做事情。总之，一般认为TSP有这些优势：不依赖表达量绝对尺度、对单调变换不敏感、跨测序平台和标准化流程更容易迁移、每个特征都可以解释为具体的基因排序关系、单个新样本可以独立计算，不要求与训练队列共同标准化。
- 纯 TSP 存在三个根本问题：（1）全局聚类容易被组织来源支配：如果把所有癌种直接放在一起，最先聚出的通常可能是组织来源、细胞组成或测序差异，而不是期望的跨癌种机制 （2）TSP特征高度相关：共享基因的 TSP、互补 TSP 和满足传递关系的 TSP 会重复表达相似信息。
- 因此，我们想：利用大量异质转录组队列，在TSP的基础上构建一种可对单个新样本独立计算、能够跨平台迁移，并可服务不同临床任务的泛癌转录组证据层。

所以，我觉得你需要补充一个实验
- 对于本次任务的真实数据来说，2个TSP很相似的样本，在各自的队列里面，就PAD分型而言，是不是很不相似？

---

我忽然有一个灵光一闪： 目前，我在实际构建CCS框架的时候（E:\iProjects\Manuscripts\CCS [只读]），有一步是先构建metaCCS，然后再构建normCCS。现在，在你的提示下，我觉得metaCCS这一步既不必须的、也不重要。因为，我有d2向量（它比d1还是好一些，因此它基于tissue分别聚类了一下，避免样本量很大的tissue遮住样本量很少的tissue）就够了；因为metaCCS本来就是d3推导出来；但是求metaCCS的时候，我还要调各种参数，自讨苦吃。我就假设，我只要到d3就足够了。然后，我准备一些ICI队列准备构建normCCS。基于已训练的队列子模型和CCS框架，我直接就获取了这些ICI队列的样本的d2矩阵。然后，我已经有了ICI response这个表型标签。 如果我想专做ICI response，我就直接 ICI rsponse ~ d3矩阵构建模型就行了。我似乎不再需要normCCS。当然，这样做也有一些缺陷：
- 在叙事上面，这种直接的、功利的转录组模型的叙事可能不如metaCCS/normCCS，生物学上的可解释性较差些
- 在CCS的框架下， 对于已经完成训练的队列子模型来说，d1显然是稳定的；但d2/d3这2步的参数的限制，我不知道d3矩阵是不是稳定的；如果d3矩阵不是稳定的话，那 ICI rsponse ~ d3矩阵就行不通，我还是得走metaCCS~normCCS的老路；因为metaCCS是d1 ~ 自聚类的解码模型得到的，它是稳定的；metaCCS是通过确定的规则获得normCCS
结合你之前的分析、以及我这里的想法，你有什么意见、建议？

---

./test/tsp-based-clustering/tsp-based-clustering-test-data.R 里有一些数据，你看注释就知道它们的用途。 我希望做一些简单的分析：
- 基于训练集（这是个泛癌数据集），利用训练集里包含的30+基因构建一个纯TSP模型（TSP的计算最好用GSClassifier的，具体你可以参考下 ./R/ablation.R 是怎么做的）。先通过一些技术方法给纯TSP的表达矩阵进行聚类（这里的聚类结果就类似于metaCCS），然后通过 f: 聚类 ~ TSP矩阵构建机器学习模型。
- 使用 Validation Data - Used for building ICI-guided normCCS-like clustering 里对应的数据，构建一个类似4分型normCCS的二次分型系统 S
- 在 Validation Data - Used for validating normCCS-like clustering in ICI prediction 里对应的数据里验证二次分型系统 S对于ICI预测能力如何
- 所有的分析结果都在 ./test/tsp-based-clustering 这个文件夹里；在这个文件夹以外的所有文件都是只读的，严禁修改。 分析要基于 bensz-rmd-rules skill 的规范，结果的解读要足够详细。

---

基于 E:\iProjects\RCheck\GSClassifier\routine02\docs\plans\2026-08-01-ccs可新增消融实验计划.md 、 docs\plans\2026-08-01-metaccs-end-to-end-ablation.md 这2个计划，优化ccs包的 R\ablation.R ，进一步完善针对CCS架构的消融实验。这是我的大致想法：
- 本轮任务的工作目录是： E:\RCloud\RFactory\ccs\test\ablation，你可以在它及其子文件夹里生成、编辑文件。在此之外的文件，都是只读以保证安全。
- 和目前ccs包的函数、定义相协调，在函数的参数上尽量不要重复地造轮子，能复用就复用；实在没有再搞新的
- 函数`ablation`的参数的数量恰到好处就行，不要冗余。
- 在 test\ablation\ablation.md 里要更新本轮加上的消融实验。
- 不要动ccs包已有的函数体系；而是主动兼容、补充到目前ccs包的代码体系（`ablation`应该相当于ccs包的一个“挂件”）。代码风格上尽量同我的接近，以便我后期可以比较简单审查你写的代码。
- 加载"E:\RCloud\RFactory\ccs\test\ablation\ablation-test-data.R" （只读）后可以获得一个真实可用的CCS对象和训练样本的RNA原始表达量，里面有一些实际的数据。你运行消融实验代码进行测试时，可以利用这些数据；但千万不能编辑、改动这些数据的原始文件，以保证安全。

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

