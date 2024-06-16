


#' @description Cohort-based probability
#' @importFrom luckyBase Fastextra
#' @inheritParams oneCCSProbability
#' @inheritParams ccsSubModel
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
ccsProb <- function(
    data,
    model.dir,
    path_tmp,
    method,
    geneAnnotation, geneSet, geneid,
    numCores,verbose
  ){

  if(verbose) LuckyVerbose('ccsProb: Calculate cohort-based probability...')
  path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', full.names = T, recursive = T)
  path_resCBP <- paste0(model.dir, '/Cohort-based probability.rds')
  if(!file.exists(path_resCBP)){
    res <- data.frame()
    for(i in 1:length(path_models)){ # i=1
      path_model1 <- path_models[i]
      name_model1 <- rev(Fastextra(path_model1, '[/]'))
      cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
      path_child <- paste0(model.dir,'/probability/',cancertype_model1, '/' , cohort_model1)
      dir.create(path_child, recursive = T, showWarnings = F)
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
              a[[cancer_name]][[cohort_name]] <- a_tmp <- CCS:::oneCCSProbability(method, data_cohort, path_model1, geneAnnotation, geneSet, geneid, numCores, dataName = paste0(cancer_name, ' - ',cohort_name), verbose = T )
            }
            saveRDS(a_tmp, path_a_tmp)
          } else {
            if(verbose) LuckyVerbose('The result of ', path_a_tmp,' exists. Use it!')
            a[[cancer_name]][[cohort_name]] <- readRDS(path_a_tmp)
          }
        }
      }

      a2 <- NULL
      for(z in 1:length(a)){ # i=1
        a2 <- rbind(a2, do.call("rbind", a[[z]]))
      }
      a2 <- a2[!duplicated(a2$SampleIDs),] # Remove duplicated data
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
    if(verbose) LuckyVerbose('The result of Cohort-based probability exists. Use it!')
    res <- readRDS(path_resCBP)
  }
  res2 <- res[,-1]; res2 <- matrix(as.numeric(as.matrix(res2)), nrow = nrow(res2), byrow = F, dimnames = list(as.character(res$SampleIDs), colnames(res2)))
  return(res2)
}


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
  res <- as.data.frame(
    cbind(
      SampleIDs = as.character(res[,'SampleIDs']),
      t(apply(res[-c(1,2,3)], 1, softmax))
    )
  )
  if(verbose) LuckyVerbose('oneCCSProbability: Model ', project, '; Data ', dataName, ' Done!')
  return(res)

}



#### GSClassifier ####
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


