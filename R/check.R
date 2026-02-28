

#' @title Check whether CCS data list or submodels are well prepared
#' @description Check whether CCS data list or submodels are well prepared
#' @param mode One of \code{'dataset'} or \code{'submodel'}
#' @param nTest How many submodels you want to test
#' @param minAccuracy An acceptable accuracy
#' @inheritParams ccs
#' @importFrom purrr flatten
#' @importFrom luckyBase LuckyVerbose
#' @importFrom GSClassifier fitEnsembleModel
#' @importFrom tidyr `%>%`
#' @importFrom stringi stri_reverse
#' @return data frame
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @seealso \code{\link{ccs}}.
#' @export
ccsCheck <- function(
    data,
    geneSet,
    geneAnnotation,
    geneid = "ensembl",
    params = list(
      device = "cpu",
      nfold = 5, nrounds = 100,
      nthread = 2, eta = 0.3, gamma = 0, max_depth = 14, colsample_bytree = 1, min_child_weight = 1, subsample = 0.7,
      n = 30, sampSize = 0.7, ptail = 0.1, nround.mode = c("fixed", "polling")[1]
    ),
    train.proc = 0.8,
    seed = 145,
    mode = c('dataset','submodel')[1],
    nTest = 3,
    minAccuracy = 0.85,
    model.dir = "./ccs/PADv20240810",
    numCores = 16,
    verbose = TRUE
){

  # Test
  if(F){

    library(luckyBase)
    Plus.library(c('CCS', 'GSClassifier', 'plotly','cowplot','tidyr','ggplot2','purrr','stringi','digest', 'pROC','ComplexHeatmap','scales','plyr','dplyr','forestplot','ggrepel','writexl','readxl','patchwork','gtable','grid'))

    mode = c('dataset','submodel')[1]

    project <- 'GibbsPanCanv20240909'
    data <- readRDS(paste0('E:/iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_',project,'.rds'))
    ImmuneSubtype <- readRDS(system.file("extdata", "ImmuneSubtype.rds", package = "GSClassifier"))
    geneSet <- ImmuneSubtype$geneSet %>% llply(.,function(x)convert(x, fromtype = 'SYMBOL',totype = 'ENSEMBL'))
    geneAnnotation <- common.annot[match(as.character(unique(unlist(geneSet))), common.annot$ENSEMBL),]

    geneid = "ensembl"
    params = list(
      device = "cpu",
      nfold = 5, nrounds = 100,
      nthread = 2, eta = 0.3, gamma = 0, max_depth = 14, colsample_bytree = 1, min_child_weight = 1, subsample = 0.7,
      n = 30, sampSize = 0.7, ptail = 0.1, nround.mode = c("fixed", "polling")[1]
    )
    seed = 145;
    train.proc = 0.8; nTest = 3; minAccuracy = 0.85
    model.dir = paste("./ccs/",project, sep = '')
    verbose = TRUE
    numCores = 16
  }

  # Parameters
  # params_2 <- params[match(intersect(c("device", "eta", "gamma", "max_depth", "min_child_weight", "nfold", "nrounds", "nthread"), names(params)), names(params))]
  params_2 <- params[setdiff(names(params), c('n', 'sampSize', 'ptail', 'nround.mode'))]


  if(mode == 'dataset'){

    set.seed(seed); seeds <- sample(1:20000, 2, replace = T)
    set.seed(seeds[1]); data_test <- flatten(data) %>% .[sample(1:length(.), nTest, replace = F)]
    data_test_res <- get_train_valid(data_test, train.proc, seed = seeds[2])
    data_test_train <- data_test_res[['data_test_train']]
    data_test_valid <- data_test_res[['data_test_valid']]

    # submodel <- readRDS("E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810/model/ACC/GSE143383/modelFit.rds")
    submodels <- list()
    for(i in 1:length(data_test_train)){
      submodel_data <- data_test_train[[i]]
      cohortName_i <- names(data_test_train)[i]
      if(verbose) LuckyVerbose('ccsCheck: build submodels based on cohort - ', cohortName_i, type = 'cat')
      submodels[[cohortName_i]] <- fitEnsembleModel(
        Xs = submodel_data$expr,
        Ys = submodel_data$subtype,
        geneSet = geneSet,
        na.fill.method = c('quantile','rpart',NULL)[1],
        na.fill.seed=2022,
        n = params$n,
        sampSize = params$sampSize,
        sampSeed = 2020,
        breakVec = c(0, 0.25, 0.5, 0.75, 1.0),
        params = params_2,
        nround.mode = params$nround.mode,
        xgboost.seed = 105,
        ptail = params$ptail,
        verbose = TRUE,
        numCores = numCores
      )
    }

    accuracy_1 <- ccsCheck_accuracy(data_test_valid, submodels, geneSet, geneAnnotation, geneid, numCores, verbose)
    accuracy_2 <- ccsCheck_accuracy(data_test_train, submodels, geneSet, geneAnnotation, geneid, numCores, verbose)

  } else if(mode == 'submodel'){

    set.seed(seed); path_submodel <- list.files(paste0(model.dir,'/model', collapse = ''), pattern = 'modelFit.rds$', recursive = T, full.names = T) %>% .[sample(1:length(.), nTest, replace = F)]

    data_test <- list(); submodels <- list()
    for(i in 1:length(path_submodel)){ # i=1
      path_submodel_i <- path_submodel[i]
      cancerType_i <- stri_reverse(Fastextra(stri_reverse(path_submodel_i),'[/]',3))
      cohortName_i <- stri_reverse(Fastextra(stri_reverse(path_submodel_i),'[/]',2))
      if(verbose) LuckyVerbose('ccsCheck: load the submoel of cohort - ', cohortName_i, type = 'cat')
      data_test[[cohortName_i]] <- data[[cancerType_i]][[cohortName_i]]
      submodels[[cohortName_i]] <- readRDS(path_submodel_i)
    }

    accuracy_1 <- ccsCheck_accuracy(data_test, submodels, geneSet, geneAnnotation, geneid, numCores, verbose)
    accuracy_2 <- NULL

  } else {
    stop('ccsCheck: wrong mode. Please select one of "dataset" or "submodel"!')
  }

  if(verbose) LuckyVerbose('ccsCheck-valid: Accuracy=', paste0(accuracy_1, collapse = ', '),', medianAccuracy=', median(accuracy_1, na.rm = T), type = 'cat')
  if(!is.null(accuracy_2)){
    if(verbose) LuckyVerbose('ccsCheck-train: Accuracy=', paste0(accuracy_2, collapse = ', '),', medianAccuracy=', median(accuracy_2, na.rm = T), type = 'cat')
  } else {
    if(verbose) LuckyVerbose('ccsCheck-train: Unavailable!')
  }

  if(median(accuracy_1, na.rm = T) > minAccuracy){
    if(verbose) LuckyVerbose('ccsCheck: pass!', type = 'cat')
  } else {
    LuckyVerbose('ccsCheck: Attention! Median accuracy is too low!')
  }

  return(data.frame(cohort = names(data_test),trainAccuracy = accuracy_2, validAccuracy = accuracy_1, stringsAsFactors = F))

}


####%%%%%%%%%%%%%% Assistant functions %%%%%%%%%%%%%%%%%%####


#' @importFrom GSClassifier modelData
#' @importFrom plyr ldply llply
#' @importFrom tidyr `%>%`
get_train_valid <- function(data_test, train.proc = 0.8, seed = seeds[2]){

  data_test_df <- ldply(data_test, function(x){
    data.frame(
      SampleIDs = colnames(x$expr),
      Subtypes = x$subtype,
      stringsAsFactors = F)
  }, .id = "Cohort")

  data_test_df_modelRes <- modelData(
    design = data_test_df,
    id.col = "SampleIDs",
    variable = c("Cohort", "Subtypes"),
    Prop = train.proc,
    seed = seed
  )

  trainSamples <- data_test_df_modelRes$Data$Train$SampleIDs
  validSamples <- data_test_df_modelRes$Data$Valid$SampleIDs

  data_test_train <- llply(data_test, function(x){
    index <- colnames(x[['expr']]) %in% trainSamples
    x[['expr']] <- as.matrix(x[['expr']]) %>% .[,index,drop=FALSE]
    x[['subtype']] <- x[['subtype']] %>% .[index]
    return(x)
  })

  data_test_valid <- llply(data_test, function(x){
    index <- colnames(x[['expr']]) %in% validSamples
    x[['expr']] <- as.matrix(x[['expr']]) %>% .[,index,drop=FALSE]
    x[['subtype']] <- x[['subtype']] %>% .[index]
    return(x)
  })

  return(list(
    data_test_train = data_test_train,
    data_test_valid = data_test_valid
  ))

  # Legacy
  if(F){
    data_test_train <- list(); data_test_valid <- list()
    for(i in 1:length(data_test)){ # i=1
      cohort_i <- names(data_test)[i]
      expr_i <- data_test[[cohort_i]][['expr']]
      subtype_i <- data_test[[cohort_i]][['subtype']]
      set.seed(seeds[i+1]); data_test_train[[cohort_i]][['expr']] <- expr_i %>% .[,sample(1:ncol(.), round(ncol(.)*train.proc), replace = F),drop=FALSE]
      data_test_valid[[cohort_i]][['expr']] <- expr_i %>% .[,setdiff(colnames(.), colnames(data_test_train[[cohort_i]][['expr']]))]
      data_test_train[[cohort_i]][['subtype']] <- subtype_i[colnames(data_test_train[[cohort_i]][['expr']])]
      data_test_valid[[cohort_i]][['subtype']] <- subtype_i[colnames(data_test_valid[[cohort_i]][['expr']])]
    }
  }

}


#' @importFrom GSClassifier callEnsemble parCallEnsemble
#' @importFrom luckyBase LuckyVerbose
ccsCheck_accuracy <- function(data_test, submodels, geneSet, geneAnnotation, geneid, numCores, verbose){

  accuracy <- NULL

  for(i in 1:length(submodels)){ # i=1

    cohortName_i <- names(data_test)[i]

    submodel_data <- data_test[[i]]

    submodel <- submodels[[i]]

    if(verbose) LuckyVerbose('ccsCheck: subtype calling of cohort - ', cohortName_i, type = 'cat')

    if(ncol(submodel_data$expr) > numCores){
      submodel_res <- parCallEnsemble(
        X = submodel_data$expr,
        ens = submodel$Model,
        geneAnnotation = geneAnnotation,
        geneSet = geneSet,
        geneid = geneid,
        scaller = NULL,
        subtype = NULL,
        verbose = verbose,
        numCores = numCores
      )
    } else {
      submodel_res <- callEnsemble(
        X = submodel_data$expr,
        ens = submodel$Model,
        geneAnnotation = geneAnnotation,
        geneSet = geneSet,
        geneid = geneid,
        scaller = NULL,
        subtype = NULL,
        verbose = verbose,
      )
    }

    submodel_real <- submodel_data$subtype
    submodel_pred <- submodel_res$BestCall_Max

    accuracy <- c(accuracy, mean(submodel_real == submodel_pred))
  }
  return(accuracy)
}
