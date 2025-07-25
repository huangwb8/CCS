
#' @description CSS Public Parameter Delivering
#' @param object a \code{\link{CCS-class}} object
#' @param size The plot size
#' @param width The plot width
#' @param height The plot height
#' @param geom Some of \code{'cancer_type'} and \code{'CCS'}.
#' @param hide.legend Some of \code{'cancer_type'} and \code{'CCS'}. Which type of data should be hided in the plot legend.
#' @param method Clustering methods. One of "\code{dbscan}", "\code{ward.D}", "\code{ward.D2}", "\code{single}", "\code{complete}", "\code{average}" (= UPGMA), "\code{mcquitty}" (= WPGMA), "\code{median}" (= WPGMC) or "\code{centroid}" (= UPGMC).
#' @param rm.zero Whether to remove zero value of CCS subtypes. In the strategy of \code{dbscan}, 0 means uncategorized, so you should set \code{rm.zero = TRUE}(default).
#' @param cover Whether to cover the existed result
#' @inheritParams ccs
CCSPublicParams <- function(
    object,
    size = 15,
    width = 12, height = 9,
    model.dir,
    geom,
    hide.legend,
    rm.zero = TRUE,
    method,
    verbose,
    numCores
  ){
  return(NULL)
}


#' @description softmax function
#' @param x a numeric vector
#' @return a normalized numeric vector
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @details
#' 尽管softmax是一个非线性函数，但softmax回归的输出仍然由输入特征的仿射变换决定。因此，softmax回归是一个线性模型（linear model）。
#' @examples
#' a <- softmax(c(1,2,3,4))
#' sum(a)
softmax <- function(x){
  return(exp(x)/sum(exp(x)))
}


#' @rdname CCS-method.createCCSData
#' @title createCCSData
#' @description Create training data for CCS package based on expression matrix and geneSet
#' @param dataList A nested list object. The first layer is the type of tumor; the second layer is the gene expression matrix.
#' @param ... parameters for \code{\link[GSClassifier]{subtypeVector}} function
#' @inheritParams GSClassifier::fitEnsembleModel
#' @inheritParams GSClassifier::subtypeVector
#' @importFrom GSClassifier subtypeVector
#' @import luckyBase
#' @return training data for \code{\link{ccs}} function
#' @seealso \code{\link{ccs}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @export
createCCSData <- function(
    dataList,
    geneSet,
    geneAnnotation = NULL,
    geneid = "ensembl",
    na.fill.method = c("quantile", "rpart", NULL)[1],
    na.fill.seed = 2024,
    verbose = TRUE,
    ...
){

  # Test
  if(F){
    library(luckyBase)
    Plus.library(c('GSClassifier'))
    dataList <- readRDS('E:/Sync/@Analysis/PanCan_Data/Level 1/PanCan_CancerSample_DataListForCCS_v20240614_GEO.rds')
    geneSet <- readRDS('E:/RCloud/database/Signature/report/Signature_PanCanWGCNA-Top10.rds')
    geneAnnotation <- NULL
    na.fill.method = c("quantile", "rpart", NULL)[1]
    na.fill.seed = 2024
    verbose = TRUE
  }

  # Environments
  if(is.null(geneAnnotation)){
    geneAnnotation <- common.annot[match(as.character(unique(unlist(geneSet))), common.annot$ENSEMBL),]
  }
  set.seed(na.fill.seed); seeds <- sample(1:100000, countListElement(dataList), replace = F)

  # Cancer type
  cancer_type <- names(dataList)

  # Calling subtype vector
  if(verbose) LuckyVerbose('createCCSData: Calling subtype vector...')
  l <- list(); seed_order <- 1
  for(i in 1:length(cancer_type)){ # i=1
    expr_type <- names(dataList[[cancer_type[i]]])
    for(j in 1:length(expr_type)){ # j=1
      expr.j <- dataList[[i]][[j]]
      expr.j.res <- GSClassifier::subtypeVector(
        expr = expr.j,
        geneSet = geneSet,
        verbose = verbose,
        ...
      )
      expr.j.res2 <- GSClassifier::geneMatch(
          X = expr.j.res$Data$expr,
          geneAnnotation = geneAnnotation,
          geneid = geneid,
          matchmode = c("fix", "free")[1]
        )
      # if(verbose) GSClassifier:::reportError(expr.j.res2)

      # Missing value imputation
      expr.j.res$Data$expr <- GSClassifier::na_fill(
        expr.j.res2$Subset,
        method = na.fill.method,
        seed = seeds[seed_order],
        verbose = verbose)
      seed_order <-  seed_order + 1

      # Other results
      expr.j.res$Data$matchError <- expr.j.res2$matchError
      expr.j.res$Data$missGenes <- expr.j.res2$missGenes
      l[[cancer_type[i]]][[expr_type[j]]] <- expr.j.res$Data
    }
  }

  # filtering out ineligable datasets
  data_filter <- NULL
  for(i in 1:length(l)){
    l_i <- l[[i]]
    for(j in 1:length(l_i)){
      l_j <- l_i[[j]]
      if(length(unique(l_j$subtype)) == 1){
        data_filter <- c(data_filter, names(l_i)[i])
        l[[i]][[j]] <- NULL
      }
    }
  }
  if(!is.null(data_filter)){
    if(verbose) LuckyVerbose('createCCSData: Datasets ',paste0(data_filter, collapse = ', '),' would be ignored due to lack of subtype diversity...')
  }

  # Output
  if(verbose) LuckyVerbose('createCCSData: Done!')
  return(l)
}


#' @rdname CCS-method.selectCCSData
#' @title selectCCSData
#' @description Advanced selection of CCS data list
#' @param dataList A nested list object. The first layer is the type of tumor; the second layer is the gene expression matrix.
#' @param coreTissueType Core tissue type you want for the CCS model establishment.
#' @param minSampleSize The lower margin of sample size of eligible cohorts.
#' @param maxMissingRate The upper margin of missing rate of eligible cohorts.
#' @inheritParams createCCSData
#' @importFrom GSClassifier geneMatch
#' @importFrom plyr laply
#' @importFrom purrr flatten
#' @import luckyBase
#' @return dataList for \code{\link{createCCSData}}.
#' @seealso \code{\link{ccs}}; \code{\link{createCCSData}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @export
selectCCSData <- function(
    dataList,
    geneSet,
    geneAnnotation = NULL,
    geneid = "ensembl",
    coreTissueType = c("ACC","AEG","BRCA","CRC","KIRC","PAAD","PRAD","STAD","Pediatric","Blood"),
    minSampleSize = 20,
    maxMissingRate = 0.3,
    verbose = TRUE
) {

  # Test
  if(F){

    library(luckyBase)
    Plus.library(c('GSClassifier', 'plyr', 'tidyr', 'purrr'))
    dataList = readRDS("E:/Sync/@Analysis/PanCan_Data/Level 1/PanCan_CancerSample_DataListForCCS_v20240809_GEO+cBioPortal+UCXCXena.rds")
    PAD <- readRDS(system.file("extdata", "PAD.train_20220916.rds", package = "GSClassifier"))
    geneSet <- PAD$geneSet
    geneAnnotation <- PAD$geneAnnotation
    geneid = "ensembl"
    coreTissueType = c("ACC","AEG","BRCA","CRC","KIRC","PAAD","PRAD","STAD","Pediatric","Blood")
    minSampleSize = 20
    maxMissingRate = 0.3
  }


  # Environments
  if(is.null(geneAnnotation)){
    geneAnnotation <- common.annot[match(as.character(unique(unlist(geneSet))), common.annot$ENSEMBL),]
  }


  # Filter small cohorts
  dataList_2 <- list()
  for(i in 1:length(dataList)){ # i=1
    cancer_type <- names(dataList)[i]
    dataList_i <- dataList[[i]]
    for(j in 1:length(dataList_i)){ # j=1
      dataList_i_j <- dataList_i[[j]]
      cohort_name <- names(dataList_i)[j]
      expr_i_j <- dataList[[cancer_type]][[cohort_name]]
      if(ncol(expr_i_j) >= minSampleSize){
        res_j <- GSClassifier::geneMatch(
          X = expr_i_j,
          geneAnnotation = geneAnnotation,
          geneid = geneid,
          matchmode = c("fix", "free")[1]
        )
        if(res_j$matchError <= maxMissingRate){
          dataList_2[[cancer_type]][[cohort_name]] <- res_j[["Subset"]]
        }
      } else {
        if(verbose) LuckyVerbose('selectCCSData: The sample size of ', names(dataList[[i]])[j], ' is lower than ', minSampleSize,'. Ignored!')
      }
    }
  }
  dataList_2 <- dataList_2[laply(dataList_2, length)>0]


  # Core Tissue Type
  dataList_3 <- list()
  for(i in 1:length(dataList_2)){ # i=1
    cancer_type_i <- names(dataList_2)[i]
    if(cancer_type_i %in% coreTissueType){
      if(!cancer_type_i %in% names(dataList_3)){
        dataList_3[[cancer_type_i]] <- dataList_2[[cancer_type_i]]
      } else {
        dataList_3[[cancer_type_i]] <- c(dataList_3[[cancer_type_i]], dataList_2[[cancer_type_i]])
      }
    } else {
      if(!'Undefined' %in% names(dataList_3)){
        dataList_3[['Undefined']] <- dataList_2[[cancer_type_i]]
      } else {
        dataList_3[['Undefined']] <- c(dataList_3[['Undefined']], dataList_2[[cancer_type_i]])
      }
    }
  }


  # Report
  raw_nSample <- laply(flatten(dataList), ncol) %>% sum(., na.rm = TRUE)
  new_nSample <- laply(flatten(dataList_3), ncol) %>% sum(., na.rm = TRUE)
  raw_nCohort <- length(flatten(dataList))
  new_nCohort <- length(flatten(dataList_3))
  diff_cohort_name <- setdiff(names(flatten(dataList)), names(flatten(dataList_3)))
  if(verbose) LuckyVerbose(raw_nSample, ' samples in the raw CCSDataList and ', new_nSample, ' (',round(new_nSample/raw_nSample,2)*100,"%)" ,' samples are kept. Filtered cohorts includes:    ', paste0(diff_cohort_name, collapse = ', '))


  # Output
  l <- list(
    Repeat = list(
      geneSet = geneSet,
      geneAnnotation = geneAnnotation,
      geneid = geneid,
      coreTissueType = coreTissueType,
      minSampleSize = minSampleSize,
      maxMissingRate = maxMissingRate
    ),
    Data = dataList_3,
    Report = list(
      raw_nSample = raw_nSample,
      raw_nCohort = raw_nCohort,
      new_nSample = new_nSample,
      new_nCohort = new_nCohort,
      diff_cohort_name = diff_cohort_name
    )
  )
  return(l)
}


#' @title Calinski-Harabasz index
#' @description Calinski-Harabasz index from fpc::cluster.stats
#' @inheritParams fpc::cluster.stats
#' @return Numeric. Calinski-Harabasz index.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @seealso \code{\link[fpc]{cluster.stats}}; \code{\link{CCS-method.optimizeCluster}}; \code{\link{CCS-method.exploreCluster}}
#' @export
CHindex <- function (
    d = NULL,
    clustering,
    noisecluster = TRUE
) {

  # Test
  if(F){

    library(luckyBase); library(CCS)

    # Data
    path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
    object <- readRDS(path_resCCS)
    object <- cluster(object, method = "dbscan", eps = 0.4, minPts = 5)
    # Parameters
    d <- dist(scale(object@Data$Probability$d3), method = "euclidean")
    clustering <- object@Data$CCS; table(clustering)
    noisecluster = TRUE
    verbose = T

    # Other parameters
    alt.clustering = NULL
    nndist = TRUE
    nnk = 2
    standardisation = "max"
    maxk = 10

    # ch = 831.1789
  }

  # Cluster message
  if (!is.null(d)){
    d <- as.dist(d)
  }
  cn <- max(clustering)
  clusteringf <- as.factor(clustering)
  clusteringl <- levels(clusteringf)
  cnn <- length(clusteringl)
  if (cn != cnn) {
    warning("clustering renumbered because maximum != number of clusters")
    for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
    cn <- cnn
  }
  n <- length(clustering)
  noisen <- 0
  cwn <- cn
  if (noisecluster) {
    noisen <- sum(clustering == cn)
    cwn <- cn - 1
  }

  # Data
  dmat <- as.matrix(d)
  nonnoise.ss <- sum(d^2)/n
  if (noisecluster){
    nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn, clustering <= cwn])^2)/sum(clustering <= cwn)
  }

  di <- list(); cluster.size <- within.cluster.ss <- 0
  for (i in 1:cn) {
    # if(verbose) LuckyVerbose(i)
    cluster.size[i] <- sum(clustering == i)
    di <- as.dist(dmat[clustering == i, clustering == i])
    if (i <= cwn) {
      within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
    }
  }
  between.cluster.ss <- nonnoise.ss - within.cluster.ss
  ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss * (cwn - 1))
  return(ch)
}


#' @description convert an expression data for t-SNE
#' @param res2 The output of \code{\link{oneCCSProbability}}
#' @param rm.dup.col Whether remove duplicated records across the colume (feature) level.
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom digest digest
#' @return a list with raw/cleaned data or md5 sum.
#' @details
#' DR algorithm like t-SNE doesn't allow the exact same rows (or samples), so \code{data_for_dr} is necessary.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
data_for_dr <- function(res2, rm.dup.col = F, verbose = F){

  res2_md5 <- apply(res2, 1, function(x)digest(x, algo="md5"))
  s <- !duplicated(res2_md5)
  if(rm.dup.col){
    res2_md5_col <- apply(res2, 2, function(x)digest(x, algo="md5"))
    s2 <- !duplicated(res2_md5_col)
  } else {
    s2 <- rep(TRUE, ncol(res2))
  }

  if(verbose){
    if(sum(!s)>=1) LuckyVerbose('No. of Duplicati sample = ',sum(!s),'. Remove them! No. of Kept sample = ', length(s) - sum(!s), ', Percentage=',round((length(s) - sum(!s))/length(s)*100, 2), '%.')
    if(sum(!s2)>=1) LuckyVerbose('No. of Duplicati feature = ',sum(!s2),'. Remove them! No. of Kept sample = ', length(s2) - sum(!s2), ', Percentage=',round((length(s2) - sum(!s2))/length(s2)*100, 2), '%.')
  }
  return(list(
    raw = list(md5=res2_md5, data=res2),
    cleaned = list(md5=res2_md5[s], data=res2[s,s2])
  ))
}


#' @description convert an expression data for t-SNE
#' @param expr Expression matrix based on \code{restSNE[['cleaned']][['data']]} with sample rows and feature cols.
#' @param restSNE Result from \code{data_for_tSNE}
#' @importFrom digest digest
#' @return a list with raw/cleaned data or md5 sum.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
repairCCS <- function(expr, restSNE){

  # Data
  md5_cleaned <- restSNE$cleaned$md5
  data_cleaned <- expr[names(md5_cleaned),]
  data_cleaned <- data_cleaned[names(md5_cleaned),]
  md5_raw <- restSNE$raw$md5
  md5_rest <- md5_raw[!names(md5_raw) %in% names(md5_cleaned)]

  # Index
  if(length(md5_rest)>0){
    index <- match(md5_rest, md5_cleaned)
    df <- NULL
    for(i in 1:length(index)){ # i=1
      df <- rbind(df, data_cleaned[index[i],])
      rownames(df)[i] <- names(md5_rest)[i]
    }
    df2 <- rbind(data_cleaned, df)
    df2 <- df2[names(md5_raw),]
  } else {
    df2 <- data_cleaned
  }


  # Output
  return(df2)

}


#' @description Dimensionality reduction for CCS probability matrix
#' @param SampleIDs Sample IDs like "GSM2411085". One or more IDs are approved.
#' @inheritParams ccs
#' @return Character. The cancer type of sample IDs.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
cancerType <- function(data,SampleIDs=c("GSM2411085","GSM2411084")){
  cancer_type <- lapply(data, function(x) as.character(unlist(lapply(x, function(y) colnames(y$expr)))))
  get_cancer_type <- function(SampleIDs,cancer_type){
    sapply(SampleIDs, function(x){
      for(i in 1:length(cancer_type)){
        if(x %in% cancer_type[[i]]){
          return(names(cancer_type)[i])
        }
      }
    })
  }
  return(get_cancer_type(SampleIDs, cancer_type))
}


#' @description Report a parameters group
#' @return Character
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
reportParams <- function(params.i){
  params.i.names <- names(params.i)
  report <- ''
  for(j in 1:length(params.i)){
    report <- paste0(report, params.i.names[j],"=", params.i[j], ifelse(j == length(params.i),"","; "))
  }
  return(report)
}


#' @description List parameters
#' @return List
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
listParams <- function(params.i){
  params.i.list <- list()
  for(j in 1:length(params.i)){
    params.i.list[[names(params.i)[j]]] <- params.i[1, j]
  }
  return(params.i.list)
}


#' @description Comparing real and pred
#' @param real Character. Real subtypes
#' @param pred Character. Predicted subtypes
#' @importFrom pROC multiclass.roc roc
#' @importFrom irr kappa2
#' @return data.frame
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
compareRealPred <- function(real, pred, cluster_translator=NULL){

  # Multi-ROC
  res.multi_auc <- multiclass.roc(real, pred, quiet = TRUE, direction = "auto")

  # Binary ROC
  res.binary <- list();
  response_subtype <- as.character(unique(real))
  for(i in 1:length(response_subtype)){ # i=1
    subtype.i <-  response_subtype[i]
    binary_i <- as.numeric(real %in% subtype.i)
    res.binary.i <- roc(binary_i, pred, quiet = TRUE)
    res.binary[[as.character(subtype.i)]] <- as.numeric(res.binary.i$auc)
  }
  if(!is.null(cluster_translator)){
    names(res.binary) <- convert(names(res.binary), 'adjust', 'raw', cluster_translator)
  }
  names(res.binary) <- paste('binary_auc_',names(res.binary), sep = '')

  # Kappa value
  res.kappa <- kappa2(data.frame(real, pred), "unweighted")
  res.all <- cbind(
    multi_auc = as.numeric(res.multi_auc$auc),
    accuracy = mean(pred == real),
    kappa = res.kappa$value,
    as.data.frame(res.binary)
  )

  # Output
  return(res.all)
}


#' @title Call ROC-AUC for Binary/multi-class response & multi-class predictor
#' @description Call ROC-AUC for Binary/multi-class response & multi-class predictor
#' @importFrom pROC multiclass.roc roc
#' @importFrom luckyBase convert LuckyVerbose
#' @importFrom dplyr arrange
#' @import tidyr
#' @return a result from \link[pROC]{multiclass.roc}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @export
multiXClass_roc <- function(response, predictor, verbose = T){

  # Test
  if(F){
    library(luckyBase)
    library(pROC)
    Plus.library(c('pROC', 'tidyr', 'dplyr'))
    # set.seed(123); response <- sample(c(0,1), 1000, replace = T)
    # set.seed(456); predictor <- sample(1:20, 1000, replace = T)
    data("myeloma", package = 'survminer')
    response <- myeloma$event
    predictor <- classify(
      myeloma$molecular_group,
      list(
        '1' = 'Cyclin D-1',
        '2' = 'Cyclin D-2',
        '3' = "Hyperdiploid",
        '4' = "Low bone disease",
        '5' = "MAF",
        '6' = "MMSET",
        '7' = "Proliferation",
        '8' = c(NA)
      ),
      useNA = F,
      cover = T
    ) %>% as.integer()
    verbose = T
  }

  # Data
  predictor_unique <- unique(predictor)

  # single ROC
  dat_auc <- NULL
  for(i in 1:length(predictor_unique)){ # i=1
    predictor_unique.i <- predictor_unique[i]
    predictor.i <- ifelse(predictor %in% predictor_unique.i, 1, 0)
    res.i <- roc(response, predictor.i, quiet = T)
    res.i.auc <- as.numeric(res.i$auc)
    dat_auc <- rbind(
      dat_auc,
      data.frame(
        class = predictor_unique.i,
        auc = res.i.auc,
        stringsAsFactors = F
      )
    )
  }
  dat_auc <- arrange(dat_auc, auc)
  class_translator <- data.frame(
    raw = dat_auc$class,
    adjust = 1:length(predictor_unique)
  )
  if(verbose){
    LuckyVerbose('multiXClass_roc: the relationship is as following:')
    print(class_translator)
  }
  predictor_adjust <- convert(predictor, 'raw', 'adjust', class_translator) %>% as.integer()

  # Data comparision
  res.adjust <- multiclass.roc(response, predictor_adjust, quiet = T)
  res.raw <- multiclass.roc(response, predictor, quiet = T)
  if(verbose){
    LuckyVerbose('multiXClass_roc: adjusted AUC=', round(res.adjust$auc, 5), '; raw AUC=', round(res.raw$auc, 5), '.')
  }

  # Output
  return(res.adjust)

}


#' @description Comparing real and pred - 2
#' @param real Integer containing 0 or 1 (Positive class should be 1)
#' @param pred Integer containing 0 or 1 (Positive class should be 1)
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc
#' @return data.frame
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
compareRealPred2 <- function(real, pred){

  #  Test
  if(F){
    library(luckyBase)
    Plus.library(c('irr','pROC','caret'))
    real = c(1, 0, 1, 1, 0, 1, 0, 0, 1, 0)
    pred = c(1, 0, 1, 0, 0, 1, 1, 0, 1, 0)
  }

  # Binary ROC
  res_roc<- roc(factor(real), pred, quiet = TRUE, ci=TRUE)

  # confusion matrix
  conf_matrix <- confusionMatrix(factor(pred), factor(real), positive = "1")

  # Summary
  res.all <- cbind(
    ROCAUC = as.numeric(res_roc$auc),
    ROCAUC_lower = as.numeric(res_roc$ci)[1],
    ROCAUC_upper = as.numeric(res_roc$ci)[3],
    accuracy = conf_matrix[["overall"]][["Accuracy"]],
    accuracy_lower = conf_matrix[["overall"]][["AccuracyLower"]],
    accuracy_upper = conf_matrix[["overall"]][["AccuracyUpper"]],
    as.data.frame(t(conf_matrix[["byClass"]]))
  )

  # Output
  return(res.all)
}


#' @description Comparing real and pred - 3
#' @param response Integer containing 0 or 1 (Positive class should be 1)
#' @param predictor Integer containing 0 or 1 (Positive class should be 1)
#' @inheritParams subtypePerformance
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc
#' @return data.frame
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
get_perform_markers <- function(
    response,
    predictor,
    n_norm_subtype,
    norm_roc_method
){
  if(n_norm_subtype == 2){
    CCS:::compareRealPred2(real = response, pred = predictor)
  } else if(n_norm_subtype > 2 & norm_roc_method == 'OptimizeRank'){
    data.frame(ROCAUC = CCS::multiXClass_roc(response, predictor, verbose = FALSE)[['auc']])
  } else if(n_norm_subtype > 2 & norm_roc_method == 'Raw'){
    data.frame(ROCAUC = multiclass.roc(response, predictor, quiet = T)[['auc']])
  } else if(!norm_roc_method %in% c('Raw','OptimizeRank')){
    stop('Not available ROC method. Please use one of "Raw" and "OptimizeRank"!')
  } else {
    NULL
  }
}


#' @description Comparing real and pred - 4. Optimized for CCS::subtypeROC function.
#' @param response Integer containing 0 or 1 (Positive class should be 1)
#' @param predictor Integer containing 0 or 1 (Positive class should be 1)
#' @inheritParams subtypePerformance
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc
#' @return data.frame
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
get_perform_markers_align <- function(
    response,
    predictor,
    n_norm_subtype,
    norm_roc_method
){
  if((n_norm_subtype == 2) |(n_norm_subtype > 2 & norm_roc_method == 'OptimizeRank')){
    data.frame(
      "ROCAUC" = CCS::multiXClass_roc(response, predictor, verbose = FALSE)[['auc']]
    )
  } else if(n_norm_subtype > 2 & norm_roc_method == 'Raw'){
    data.frame(
      "ROCAUC" = multiclass.roc(response, predictor, quiet = T)[['auc']]
    )
  } else if(!norm_roc_method %in% c('Raw','OptimizeRank')){
    stop('Not available ROC method. Please use one of "Raw" and "OptimizeRank"!')
  } else {
    NULL
  }
}

#' @description Generate scaller path for Deep Learning model
#' @importFrom digest digest
generate_scaller_path <- function(model.dir) {

  # Get time string with millisecond precision
  time_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S.%OS3")

  # Calculate and return the MD5 hash
  time_md5 <- digest::digest(time_str, algo = "md5", serialize = FALSE)

  scaller_path <- paste0(model.dir, '/.scaller/', time_md5)

  return(scaller_path)

}


