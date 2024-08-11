

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
#' @return NULL
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
    seed = 145,
    mode = c('dataset','submodel')[1],
    nTest = 3,
    minAccuracy = 0.9,
    model.dir = "./ccs/PADv20240810",
    verbose = TRUE
){

  # Test
  if(F){
    mode = c('dataset','submodel')[1]
    data = data
    geneSet = geneSet
    geneAnnotation = geneAnnotation
    geneid = "ensembl"
    params = list(
      device = "cpu",
      nfold = 5, nrounds = 100,
      nthread = 2, eta = 0.3, gamma = 0, max_depth = 14, colsample_bytree = 1, min_child_weight = 1, subsample = 0.7,
      n = 30, sampSize = 0.7, ptail = 0.1, nround.mode = c("fixed", "polling")[1]
    )
    seed = 145; nTest = 3
    model.dir = "./ccs/PADv20240810"
    verbose = TRUE
  }

  # Parameters
  params_2 <- params[match(intersect(c("device", "eta", "gamma", "max_depth", "min_child_weight", "nfold", "nrounds", "nthread"), names(params)), names(params))]

  if(mode == 'dataset'){

    set.seed(seed); data_test <- flatten(data) %>% .[sample(1:length(.), nTest, replace = F)]

    # submodel <- readRDS("E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810/model/ACC/GSE143383/modelFit.rds")
    submodels <- list()
    for(i in 1:length(data_test)){
      submodel_data <- data_test[[i]]
      cohortName_i <- names(data_test)[i]
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
        numCores = 16
      )
    }

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

  } else {
    stop('ccsCheck: wrong mode. Please select one of "dataset" or "submodel"!')
  }

  accuracy <- ccsCheck_accuracy(data_test, submodels, geneSet, geneAnnotation, geneid, verbose)

  if(verbose) LuckyVerbose('ccsCheck: Accuracy=', paste0(accuracy, collapse = ', '),', medianAccuracy=', median(accuracy, na.rm = T), type = 'cat')

  if(median(accuracy, na.rm = T) > minAccuracy){
    if(verbose) LuckyVerbose('ccsCheck: pass!', type = 'cat')
  } else {
    stop('ccsCheck: median accuracy is too low! Stopped')
  }

}


####%%%%%%%%%%%%%% Assistant functions %%%%%%%%%%%%%%%%%%####

#' @importFrom GSClassifier callEnsemble
#' @importFrom luckyBase LuckyVerbose
ccsCheck_accuracy <- function(data_test, submodels, geneSet, geneAnnotation, geneid, verbose){

  accuracy <- NULL

  for(i in 1:length(submodels)){ # i=1

    cohortName_i <- names(data_test)[i]

    submodel_data <- data_test[[i]]

    submodel <- submodels[[i]]

    if(verbose) LuckyVerbose('ccsCheck: subtype calling of cohort - ', cohortName_i, type = 'cat')

    submodel_res <- callEnsemble(
      X = submodel_data$expr,
      ens = submodel$Model,
      geneAnnotation = geneAnnotation,
      geneSet = geneSet,
      geneid = geneid,
      scaller = NULL,
      subtype = NULL
    )

    submodel_real <- submodel_data$subtype
    submodel_pred <- submodel_res$BestCall_Max

    accuracy <- c(accuracy, mean(submodel_real == submodel_pred))
  }
  return(accuracy)
}
