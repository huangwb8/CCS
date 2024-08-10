


#' @description Cohort-based probability
#' @import luckyBase
#' @import tidyr
#' @importFrom purrr flatten
#' @param parallel.method The strategy to to do parallel calculation via \code{\link{parCallEnsemble}}. If \code{parallel.method = 'ensemble'}(default), all dataset in \code{data} would be merge into one large dataset. If \code{parallel.method = 'discrete'}(legacy), every dataset would be calling respectively.
#' @inheritParams oneCCSProbability
#' @inheritParams ccsSubModel
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @description Cohort-based probability in an ensemble way
ccsProb <- function(
    data,
    model.dir,
    method,
    geneAnnotation, geneSet, geneid,
    parallel.method = c('ensemble','discrete')[1],
    numCores,verbose
){

  # Test
  if(F){
    library(luckyBase)
    np <- Plus.library(c('tidyr','dplyr','GSClassifier','purrr'))

    # data <- readRDS('E:/iProjects/CCS_Data/product/PanCan_CancerSample_DataListForCCS-Signature-PanCanWGCNA-Top10_v20240623_GEO.rds')[c(1,2)]
    # model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/Testv20240623'
    # geneSet <- readRDS('E:/RCloud/database/Signature/report/Signature_PanCanWGCNA-Top10_v20240623.rds')
    # geneAnnotation <- common.annot[match(as.character(unlist(geneSet)), common.annot$ENSEMBL),]

    data <- readRDS('E:/iProjects/CCS_Data/report/PanCan_CancerSample_DataListForCCS_v20240809_GEO+cBioPortal+UCXCXena_PAD_v20240810.rds')
    model.dir = 'E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810'
    PAD <- readRDS(system.file("extdata", "PAD.train_20220916.rds", package = "GSClassifier"))
    geneSet <- PAD$geneSet
    geneAnnotation <- PAD$geneAnnotation

    method = 'GSClassifier'
    numCores = 16
    verbose = T
    geneid = 'ensembl'
    parallel.method = 'ensemble'
  }

  if(verbose) LuckyVerbose('ccsProb: Calculate cohort-based probability...')
  path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', full.names = T, recursive = T) %>% .[grepl(paste0(names(flatten(data)), collapse = '|'), .)]
  path_resCBP <- paste0(model.dir, '/Cohort-based probability.rds')

  # All object
  if(!file.exists(path_resCBP)){

    res <- data.frame()
    for(i in 1:length(path_models)){ # i=1
      path_model1 <- path_models[i]
      name_model1 <- rev(Fastextra(path_model1, '[/]'))
      cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
      path_child <- paste0(model.dir,'/probability/',cancertype_model1, '/' , cohort_model1)
      dir.create(path_child, recursive = T, showWarnings = F)

      # Get & save every dataset-dataset probability as oneCCSProbabilityResult_*.rds
      if(parallel.method == 'ensemble'){
        CCS:::ccsProbEnsemble(
          method, data, model.dir, path_model1,
          geneAnnotation, geneSet, geneid,
          numCores, verbose
        )
      } else if(parallel.method == 'discrete'){
        CCS:::ccsProbDiscrete(
          method, data, model.dir, path_model1,
          geneAnnotation, geneSet, geneid,
          numCores, verbose
        )
      } else {
        stop('ccsProb: Wrong parallel method! Please use one of "ensemble" and "discrete"!')
      }

      # Merge Data
      # path_child = "E:/iProjects/RCheck/GSClassifier/test01/ccs/Testv20240623/probability/ACC/GSE143383"
      # path_child = "E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810/probability/ACC/GSE143383"; cancertype_model1 ='ACC'; cohort_model1 = 'GSE143383'
      a2 <- list.files(path_child, 'oneCCSProbabilityResult_', full.names = T, recursive = T) %>% lapply(., readRDS) %>% do.call("rbind", .)
      a2 <- a2[!(duplicated(a2$SampleIDs) | is.na(a2$SampleIDs)),] # Remove duplicated data
      colnames(a2)[2:ncol(a2)] <- paste(cancertype_model1, cohort_model1,  colnames(a2)[2:ncol(a2)], sep = '|')
      rownames(a2) <- as.character(a2$SampleIDs)

      # merge
      if(i==1){
        res <- a2
      } else {
        res <- cbind(res, a2[,-1])
      }
    }
    saveRDS(res, path_resCBP)

  } else {
    if(verbose) LuckyVerbose('ccsProb: The result of Cohort-based probability exists. Use it!')
    res <- readRDS(path_resCBP)
  }
  res2 <- res[,-1]; res2 <- matrix(as.numeric(as.matrix(res2)), nrow = nrow(res2), byrow = F, dimnames = list(as.character(res$SampleIDs), colnames(res2)))
  return(res2)
}


#### Other Assistant functions ####

#' @description Get target data
#' @importFrom dplyr full_join
#' @importFrom luckyBase Fastextra
getResData <- function(data, pred_i_res){
  data2 <- NULL
  for(n in 1:length(pred_i_res)){ # n=1
    text <- Fastextra(pred_i_res[n],'-')
    cancer_type.n <- text[1]; cohort.n <- gsub(paste0('^',cancer_type.n,'-',collapse = ''),'',pred_i_res[n])
    expr.n <- data[[cancer_type.n]][[cohort.n]][['expr']]
    colnames(expr.n) <- paste(cancer_type.n, cohort.n, colnames(expr.n), sep = '|') # Keep the sample name unique
    if(is.null(data2)){
      data2 <- data.frame(gene=rownames(expr.n), expr.n, check.names = FALSE)
    } else {
      data2 <- full_join(
        data2,
        data.frame(gene=rownames(expr.n), expr.n, check.names = FALSE),
        by = 'gene'
      )
    }
  }
  genes <- as.character(data2[,'gene'])
  data2 <- as.matrix(data2[-match('gene',colnames(data2))])
  rownames(data2) <- genes
  # colnames(data2) <- Fastextra(colnames(data2),'|',3)
  return(data2)
}


#' @description Split Probability data
#' @importFrom luckyBase Fastextra
splitProbData <- function(data2_prob, data, pred_i_res, path_child, cancertype_model1, cohort_model1){

  # Test
  if(F){
    test <- readRDS('E:/iProjects/RCheck/GSClassifier/test01/ccs/Testv20240623/probability/CRC/GSE27854/oneCCSProbabilityResult_Model-CRC-GSE27854_Data-AEG-GSE96667.rds')
    data2_prob <- readRDS('test/data2_prob.rds')
    data2_prob_withoutNA <- na.omit(data2_prob)
  }

  for(n in 1:length(pred_i_res)){ # n=1
    text <- Fastextra(pred_i_res[n],'-')
    cancer_type.n <- text[1]; cohort.n <- gsub(paste0('^',cancer_type.n,'-',collapse = ''),'',pred_i_res[n])
    sample.n.r <- colnames(data[[cancer_type.n]][[cohort.n]][['expr']])
    sample.n <- paste(cancer_type.n, cohort.n, sample.n.r, sep='|')
    path_a_tmp <- paste0(path_child, '/oneCCSProbabilityResult_Model-',cancertype_model1,'-',cohort_model1,'_Data-',cancer_type.n, '-',cohort.n,'.rds')
    index <- match(sample.n, as.character(data2_prob$SampleIDs)) # cannot use `match` because there might be duplicati samples!
    data2_prob_2 <- data2_prob[index,]
    data2_prob_2$SampleIDs <- as.character(Fastextra(data2_prob_2$SampleIDs,'[|]',3))
    saveRDS(data2_prob_2, path_a_tmp)
  }
}

#' @description Calculate probability in ensemble mode
#' @importFrom tidyr `%>%`
#' @import luckyBase
ccsProbEnsemble <- function(
    method, data, model.dir, path_model1,
    geneAnnotation, geneSet, geneid,
    numCores, verbose
){
  name_model1 <- rev(Fastextra(path_model1, '[/]'))
  cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
  path_child <- paste0(model.dir,'/probability/',cancertype_model1, '/' , cohort_model1)

  # All datasets to be predicted
  pred_all <- NULL
  for(i in 1:length(data)){
    for(j in 1:length(data[[i]])){
      pred_all <- c(pred_all, paste0(names(data)[i],'-',names(data[[i]])[j]))
    }
  }

  # check undo prediction by this submodel
  pred_i_done <- list.files(path_child, 'oneCCSProbabilityResult_', full.names = T, recursive = T) %>% Fastextra(split = '_Data-',2) %>% gsub('.rds$','',.)
  pred_i_res <- setdiff(pred_all, pred_i_done)

  # Call probability
  if(length(pred_i_res) > 0){

    data2 <- CCS:::getResData(data, pred_i_res)
    data2_prob <- CCS:::oneCCSProbability(
      method, list(expr=data2), path_model1,
      geneAnnotation, geneSet, geneid,
      numCores, dataName = paste(pred_i_res, collapse = ', '),
      verbose = T
    )
    data2_prob$SampleIDs <- gsub('[.]','|',data2_prob$SampleIDs)
    # saveRDS(data2_prob, paste0(path_child,'/','data2_prob.rds'))

    splitProbData(
      data2_prob, data, pred_i_res,
      path_child, cancertype_model1, cohort_model1
    )

  } else {
    if(verbose) LuckyVerbose('ccsProb: The result of Model ',paste0(cancertype_model1, ' - ' , cohort_model1),' exists. Use it!')
  }
}

#' @description Calculate probability in discrete mode
#' @import luckyBase
ccsProbDiscrete <- function(
    method, data, model.dir, path_model1,
    geneAnnotation, geneSet, geneid,
    numCores, verbose
){
  name_model1 <- rev(Fastextra(path_model1, '[/]'))
  cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
  path_child <- paste0(model.dir,'/probability/',cancertype_model1, '/' , cohort_model1)

  a <- NULL
  for(j in 1:length(data)){
    data_cancer <- data[[j]]; cancer_name <- names(data)[j]
    for(k in 1:length(data_cancer)){
      data_cohort <- data_cancer[[k]]; cohort_name <- names(data_cancer)[k]
      path_a_tmp <- paste0(path_child, '/oneCCSProbabilityResult_Model-',cancertype_model1,'-',cohort_model1,'_Data-',cancer_name, '-',cohort_name,'.rds')
      # Core function
      # 这里略作修改可以节省内存。有空再优化！
      if(!file.exists(path_a_tmp)){
        if(method == 'GSClassifier'){
          a_tmp <- CCS:::oneCCSProbability(
            method, data_cohort, path_model1,
            geneAnnotation, geneSet, geneid,
            numCores, dataName = paste0(cancer_name, ' - ',cohort_name), verbose = T )
        }
        saveRDS(a_tmp, path_a_tmp)
      } else {
        if(verbose) LuckyVerbose('ccsProb: The result of ', path_a_tmp,' exists. Use it!')
      }
    }
  }

}

#### Probability Core Functions ####


#' @description OneModel-OneData CSS process
#' @param data1 a list containing an expression matrix and its subtype vector
#' @param path_model1 the path of a (GSClassifier) model
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom luckyBase LuckyVerbose Fastextra
#' @return a data frame with sample IDS and the softmax probability.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
oneCCSProbability <- function(
    method = 'GSClassifier',
    data1,
    path_model1,
    geneAnnotation,
    geneSet,
    geneid,
    numCores,
    dataName = 'GSEXXX',
    verbose = T
){

  # Test
  if(F){
    library(luckyBase)
    Plus.library(c('GSClassifier'))
    data1 = data[[2]][[3]]
    path_model1 = "./ccs/test01/cohort.1.1/modelFit.rds"
    dataName = 'GSEXXX'
    verbose = T
    l <- list(A=list(A1=1,A2=2), B=list(B1=1,B2=2))
    l_name <- lapply(l, function(x) lapply(x, function(y) names(x)))

    # 20240302
    # Error in exp(x) : non-numeric argument to mathematical function
    data_all <- readRDS('E:/Sync/@Analysis/PanCan_Data/Level 1/PanCan_CancerSample_DataListForCCS-PAD_Train+Valid_v20231224.rds')
    PADi <- readRDS(system.file("extdata", paste0('PAD.train_20220916.rds'), package = "GSClassifier"))
    path_model1 = "E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240225/model/ACC/GSE33371/modelFit.rds"
    data <- readRDS('E:/Sync/@Analysis/PanCan_Data/Level 1/PanCan_CancerSample_DataListForCCS-PAD_Train+Valid_v20231224.rds')
    data1 = data[['PAAD']][['GSE21501']]
    # data1 = data[['AEG']][['GSE74553']]
    dataName = 'GSE21501'
    PADi <- readRDS(system.file("extdata", paste0('PAD.train_20220916.rds'), package = "GSClassifier"))
    geneAnnotation = PADi$geneAnnotation
    geneSet = PADi$geneSet
    geneid = "ensembl"
    numCores= 6
    verbose = T
    method = 'GSClassifier'
  }

  # Project
  x <- rev(Fastextra(path_model1, '[/]'))
  cancer <- x[3]; cohort <- x[2]
  project <- paste0(cancer, ' - ' ,cohort)
  if(verbose) LuckyVerbose('oneCCSProbability: Model ', project, '; Data ', dataName, ' ...')

  # Call probability score
  modelFit <- readRDS(path_model1)
  if(method == 'GSClassifier'){
    res <- oneCCSProbability_GSClassifier(
      data1,
      modelFit,
      geneSet,
      geneid,
      geneAnnotation,
      numCores
    )
  } else {
    if(verbose) LuckyVerbose('oneCCSProbability: Wrong model establishment methods. Stop CCS process! Σ( ° △ °|||)︴')
  }

  # Output
  # res <-
  #   cbind(
  #     SampleIDs = as.character(res[,'SampleIDs']),
  #     as.data.frame(t(apply(res[-c(1,2,3)], 1, softmax)))
  #   )

  # if(data1.name == cohort){
  #   sink(paste0(path_log,"/","submode_self prediction.txt"), append = TRUE)
  #   r <- paste0(round(mean(data1$subtype == res$BestCall_Max)*100, 2), '%')
  #   LuckyVerbose('Self Prediction Accuracy - ',cohort,': ', r, type = 'cat')
  #   sink()
  # }

  res <-
    cbind(
      SampleIDs = as.character(res[,'SampleIDs']),
      as.data.frame(res[-c(1,2,3)])
    )

  if(verbose) LuckyVerbose('oneCCSProbability: Done!')
  return(res)

}



#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @importFrom GSClassifier parCallEnsemble
oneCCSProbability_GSClassifier <- function(
    data1,
    modelFit,
    geneSet,
    geneid,
    geneAnnotation,
    numCores
){
  res <- parCallEnsemble(
    X = data1$expr,
    ens = modelFit$Model,
    geneAnnotation = geneAnnotation,
    geneSet = geneSet,
    geneid = geneid,
    scaller = NULL,
    subtype = NULL,
    numCores = numCores
  )
  return(res)
}


#### End ####


