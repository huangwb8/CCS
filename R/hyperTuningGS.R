

#' @title Hyperparameter turning for GSClassifier models
#' @description Hyperparameter turning for GSClassifier models
#' @param Xs_train Expression matrix of the training cohort
#' @param Ys_train Subtype vector of the training cohort
#' @param Xs_valid Expression matrix of the internal validation cohort
#' @param Ys_valid Subtype vector of the internal validation cohort
#' @param para.grid a data frame of alternative parameters
#' @param model.dir the data of pr
#' @inheritParams GSClassifier::parCallEnsemble
#' @import GSClassifier
#' @importFrom digest digest
#' @importFrom xgboost xgb.DMatrix xgb.cv xgboost
#' @details
#' abc \cr
#' cde \cr
#' @return NULL. Data in \code{model.dir}
#' @seealso \code{\link[GSClassifier]{parCallEnsemble}}; \code{\link[GSClassifier]{fitEnsembleModel}};
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' @export
hyperTuningGS <- function(
    Xs_train, Ys_train,
    Xs_valid, Ys_valid,
    geneSet = NULL,
    geneAnnotation = NULL,
    geneid = "ensembl",
    para.grid = NULL,
    seed = 489,
    model.dir = './hyperTuningGS/test01',
    verbose = TRUE,
    numCores = 2
){

  # Test
  if(F){

    seed = 897; numCores = 24
    library(luckyBase)
    nPac <- c('GSClassifier', 'xgboost', 'pROC', 'digest', 'irr'); Plus.library(nPac)
    testData <- readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))
    expr <- testData$PanSTAD_expr_part
    design <- testData$PanSTAD_phenotype_part
    modelInfo <- modelData(
      design,
      id.col = "ID",
      variable = c("platform", "PAD_subtype"),
      Prop = 0.05,
      seed = 145
    )
    modelInfo <- modelData(
      modelInfo$Data$Train,
      id.col = "ID",
      variable = c("platform", "PAD_subtype"),
      Prop = 0.7,
      seed = 258
    )

    Xs_train <- expr[,modelInfo$Data$Train$ID]
    y <- modelInfo$Data$Train
    y <- y[colnames(Xs_train),]
    Ys_train <- ifelse(y$PAD_subtype == 'PAD-I',1,ifelse(y$PAD_subtype == 'PAD-II',2,ifelse(y$PAD_subtype == 'PAD-III',3,ifelse(y$PAD_subtype == 'PAD-IV',4,NA)))); table(Ys_train)/length(Ys_train)

    Xs_valid <- expr[,modelInfo$Data$Valid$ID]
    y <- modelInfo$Data$Valid
    y <- y[colnames(Xs_valid),]
    Ys_valid <- ifelse(y$PAD_subtype == 'PAD-I',1,ifelse(y$PAD_subtype == 'PAD-II',2,ifelse(y$PAD_subtype == 'PAD-III',3,ifelse(y$PAD_subtype == 'PAD-IV',4,NA)))); table(Ys_valid)/length(Ys_valid)

    PADi <- readRDS(system.file("extdata", paste0('PAD.train_20220916.rds'), package = "GSClassifier"))
    geneSet <- PADi$geneSet
    geneAnnotation <- PADi$geneAnnotation

    para.grid  = NULL
    geneid = "ensembl"
    verbose = T

    model.dir = './hyperTuningGS/test01'

  }

  # Parameters grid
  if(is.null(para.grid)){

    # https://xgboost.readthedocs.io/en/latest/parameter.html

    # https://medium.com/@rithpansanga/optimizing-xgboost-a-guide-to-hyperparameter-tuning-77b6e48e289d

    para.grid <- expand.grid(

      ## No. of cross validation subcohorts
      nfold = 5,

      ## Max number of boosting iterations. GSClassifier have optimized nrounds, so here I set a large value.
      nrounds = 100,

      ## No. of CPU cores used
      nthread = 10,

      ## Step size shrinkage used in update to prevent overfitting.  [0,1]
      eta = c(0.5, 0.3, 0.1),

      ## The larger gamma is, the more conservative the algorithm will be. [0,∞]
      gamma = c(0, 0.1, 0.2, 0.5),

      ## maximum depth of a tree. Deeper trees can capture more complex patterns in the data, but may also lead to overfitting. [0,∞]
      max_depth = c(6, 10, 14),

      ## percentage of columns used for each tree construction. Lowering this value can prevent overfitting by training on a subset of the features. (0, 1]
      colsample_bytree = 1,

      ## The larger min_child_weight is, the more conservative the algorithm will be. [0,∞]
      min_child_weight = 1,

      ## for preventing overfitting
      subsample = c(0.7, 1),

      # fitEnsembleModel
      n = c(100, 500),
      sampSize = 0.7,
      ptail = c(0.4, 0.2, 0.3, 0.5)
    )


    set.seed(2023); s <- sample(1:nrow(para.grid), nrow(para.grid), replace = F)

    para.grid <- para.grid[s,]

  }

  # Grobal seeds
  set.seed(seed); seeds <- sample(1:10000, 20, replace = T)
  dir.create(model.dir, showWarnings = F, recursive = T)

  # Test parameters
  for(i in 1:nrow(para.grid)){ # i=6

    # para
    params <- list()
    for(j in 1:ncol(para.grid)){
      params[[colnames(para.grid)[j]]] <- para.grid[i,j]
    }
    params_xg <- params[-match(c('n','sampSize','ptail'), names(params))]
    params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg))]

    # dir
    project <- digest(params, algo="md5")
    LuckyVerbose('New project: ', project)
    path_child <- paste0(model.dir,'/', project)
    dir.create(path_child, recursive = T, showWarnings = F)
    saveRDS(params, paste0(path_child, '/params.rds'))

    # fit
    path_fit <- paste0(path_child, '/modelFit.rds')
    if(!file.exists(path_fit)){
      modelFit <- fitEnsembleModel(
        Xs = Xs_train,
        Ys = Ys_train,
        geneSet = geneSet,
        na.fill.method = "quantile",
        na.fill.seed = seeds[1],
        n = params$n,
        sampSize = params$sampSize,
        sampSeed = seeds[2],
        breakVec = c(0, 0.25, 0.5, 0.75, 1),
        params = params_xg,
        xgboost.seed = seeds[3],
        caret.grid = NULL,
        ptail = params$ptail,
        verbose = T,
        numCores = numCores
      )
      saveRDS(modelFit, path_fit)
    } else {
      LuckyVerbose(project, ': modelFit exists. Ignored!')
      modelFit <- readRDS(path_fit)
    }

    # pred_train
    path_scaller <- paste0(path_child, '/scaller.rds')
    if(!file.exists(path_scaller)){

      resTrain <- parCallEnsemble(
        X = Xs_train,
        ens = modelFit$Model,
        geneAnnotation = geneAnnotation,
        geneSet = geneSet,
        geneid = geneid,
        scaller = NULL,
        subtype = NULL,
        numCores = numCores
      )

      # xgboost via best interation
      nSubtype <- length(unique(Ys_train))
      dtrain <- xgb.DMatrix(as.matrix(resTrain[4:(3 + nSubtype)]), label = Ys_train-1)
      set.seed(seeds[4])
      cvRes <- xgb.cv(
        params = params_xg2,
        data = dtrain,
        nrounds=params_xg$nrounds,
        nfold=params_xg$nfold,
        num_class = nSubtype,
        early_stopping_rounds=10,
        objective = "multi:softmax",
        verbose = 0
      )

      best_iteration <- cvRes$best_iteration

      # xgboost via best interation
      set.seed(seeds[4])
      scaller <- xgboost(
        params = params_xg2,
        data = dtrain,
        nrounds = best_iteration,
        num_class = nSubtype,
        objective = "multi:softmax",
        verbose = 0
      )

      # cvRes <- xgb.cv(data = dtrain,
      #                 nrounds=100,
      #                 nthread=10,
      #                 nfold=5,
      #                 max_depth=5,
      #                 eta=0.5,
      #                 early_stopping_rounds=100,
      #                 num_class = 4,
      #                 objective = "multi:softmax",
      #                 verbose = 1)

      # scaller <- xgboost(
      #   data = dtrain,
      #   max_depth=5,
      #   eta=0.5,
      #   nrounds = cvRes$best_iteration,
      #   nthread=10,
      #   num_class = 4,
      #   objective = "multi:softmax"
      # )

      # output
      saveRDS(scaller, path_scaller)
    } else {
      LuckyVerbose(project, ': scaller exists. Ignored!')
      scaller <- readRDS(path_scaller)
    }


    # fitness for training cohort
    path_fitness_train <- paste0(path_child, '/fitness_train.rds')
    if(!file.exists(path_fitness_train)){
      # multi-roc; binary-roc; accuracy; kappa;
      Ys_train_pred <- predict(scaller, as.matrix(resTrain[4:7])) + 1
      res.all <- compareRealPred(Ys_train, Ys_train_pred)
      names(res.all) <- paste('train_', names(res.all), sep = '')
      saveRDS(res.all, path_fitness_train)
    } else {
      LuckyVerbose(project, ': fitness for training cohort exists. Ignored!')
    }

    # fitness for internal validation cohort
    path_fitness_valid <- paste0(path_child, '/fitness_valid.rds')
    if(!file.exists(path_fitness_valid)){
      # multi-roc; binary-roc; accuracy; kappa;
      resValid <- parCallEnsemble(
        X = Xs_valid,
        ens = modelFit$Model,
        geneAnnotation = geneAnnotation,
        geneSet = geneSet,
        geneid = geneid,
        scaller = scaller,
        subtype = NULL,
        numCores = numCores
      )
      Ys_valid_pred <- resValid$BestCall
      res.all <- compareRealPred(Ys_valid, Ys_valid_pred)
      names(res.all) <- paste('valid_', names(res.all), sep = '')
      saveRDS(res.all, path_fitness_valid)
    } else {
      LuckyVerbose(project, ': fitness for internal validation cohort exists. Ignored!')
    }

  }

  # Output
  return(NULL)

}









