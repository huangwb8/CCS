

setGeneric("scaller", function(object, ...) {
  standardGeneric("scaller")
})


#' @rdname CCS-method.scaller
#' @title CCS method: scaller
#' @description \code{scaller} method for \code{CCS} class
#' @param ... Parameters for \code{\link[xgboost]{xgboost}}.
#' @param best_iteration if \code{best_iteration = NULL}(default), cross validation would be used to decide the best iteration round.
#' @inheritParams CCSPublicParams
#' @inheritParams GSClassifier::modelData
#' @import xgboost
#' @importFrom GSClassifier modelData
#' @importFrom pROC multiclass.roc roc
#' @importFrom irr kappa2
#' @import tidyr
#' @details
#' The CCS object must be modified by \code{\link{CCS-method.cluster}} or \code{\link{CCS-method.optimizeCluster}}.
#' @seealso \code{\link{ccs}};\code{\link[xgboost]{xgboost}};
#' @examples
#' resCCS <- scaller(resCCS)
#' @exportMethod scaller
setMethod(
  "scaller",
  signature(object='CCS'),
  function(
    object,
    model.dir = NULL,
    Prop = 0.8,
    best_iteration = NULL,
    seed = 456,
    cover = FALSE,
    rm.zero = TRUE,
    numCores = 6,
    verbose = TRUE,
    ...
  ){

    # Test
    if(F){
      library(luckyBase)
      Plus.library(c('xgboost', 'GSClassifier','pROC','irr','tidyr'))
      model.dir <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240106'
      path_resCCS <- paste0(model.dir,'/resCCS.rds')
      object <- readRDS(path_resCCS)
      # Renew CCS manually
      classifier <- list(
        "13" = c("13","15","20"),
        "6" = c("48","36","38","3"),
        "45" = c("45","23"),
        "75" = c("75","73"),
        "26" = c("26","25","58","65","27","24"),
        "34" = c("34", "31"),
        "40" = c("56", "43","40"),
        "55" = c("55", "57"),
        "16" = c("16", "66"),
        "82" = c("80", "82"),
        "44" = c("44", "37"),
        "76" = c("76","77","42"),
        "4" = c("4","10"),
        "8" = c("8","17"),
        "18" = c("18","9"),
        "7" = c("7","1","60","5","61"),
        "35" = c("35", "67", "47","41","83","64","22","14","11","12", "69"),
        "70" = c("70", "59")
      )
      # object  <- CCS::manualCluster(object , classifier)
      numCores = 6; seed = 456; Prop = 0.8
      verbose = TRUE; cover = TRUE; rm.zero = TRUE
    }

    # CCS normalization
    set.seed(seed); seeds <- sample(1:10000, 10, replace = T)
    y2 <- object@Data$CCS
    if(rm.zero){
      target_sample <- !y2 %in% 0
      y2 <- y2[target_sample]
      if(verbose) LuckyVerbose('scaller: Remove ', sum(!target_sample), ' zero value in CCS subtypes...')
    }
    y2_uniqueSubtype <- unique(y2)
    cluster_translator <- data.frame(
      raw = y2_uniqueSubtype,
      adjust = 0:(length(y2_uniqueSubtype) - 1),
      stringsAsFactors = F
    )
    y2_name <- names(y2)
    y2 <- convert(y2, 'raw', 'adjust', cluster_translator) %>% as.integer()
    names(y2) <- y2_name

    # Data
    dat <- data.frame(ID=y2_name, subtype = y2, stringsAsFactors = F)
    dat_modelData <- GSClassifier::modelData(
      design = dat,
      id.col = "ID",
      variable = c("subtype"),
      Prop = Prop,
      seed = seeds[1])
    id_train <- as.character(dat_modelData$Data$Train[,'ID'])
    id_valid <- as.character(dat_modelData$Data$Valid[,'ID'])
    res2 <- object@Data$Probability$d1[y2_name,]
    res2_train <- res2[id_train,] # dim(res2_train)
    y2_train <- y2[id_train]
    res2_valid <- res2[id_valid,] # dim(res2_valid)
    y2_valid <- y2[id_valid]

    # Self-defined parameters
    if(verbose) LuckyVerbose('scaller: adjust parameters...')
    params <- object@Repeat$params
    params_name <- names(params)
    args <- list(...)
    args_name <- names(args)
    for(i in 1:length(args_name)){ # i=1
      args_name_i <- args_name[i]
      if(args_name_i %in% params_name){
        if(args[[args_name_i]] != params[[args_name_i]]){
          if(verbose) LuckyVerbose('Change `', args_name_i, '` from ', params[[args_name_i]], ' to ', args[[args_name_i]], '...', levels = 3)
        }
      } else {
        if(verbose) LuckyVerbose('Add new parameter `', args_name_i, '` = ', args[[args_name_i]], levels = 3)
      }
      params[[args_name_i]] <- args[[args_name_i]]
    }


    # Model parameters
    if(is.null(model.dir)){
      model.dir = object@Repeat$model.dir
    }
    path_scaller <- paste0(model.dir, '/scaller.rds')

    params_xg <- params[-match(c('n','sampSize','ptail','nround.mode'), names(params))]
    params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg))]


    # Check old scaller
    if(!is.null(object@Data[['scaller']]) & !cover){
      stop('scaller: With scaller object in the CCS object. Ignore! You can set `cover = TRUE` to make `scaller` run and cover the old result.')
    }
    if(verbose) LuckyVerbose('scaller: Build subtype caller ...')


    # Subtype Caller
    if((!file.exists(path_scaller)) | cover){

      nSubtype <- length(unique(y2_train))
      dtrain <- xgb.DMatrix(res2_train, label = y2_train)

      # Parameters
      params_xg3 <- params_xg2
      params_xg3[['nthread']] <- numCores

      # find best interation
      if(is.null(best_iteration)){
        # xgboost cv for best interation exploration
        if(verbose) LuckyVerbose('scaller: xgboost cv for best interation exploration...')
        set.seed(seeds[2])
        cvRes <- xgb.cv(
          params = params_xg3,
          data = dtrain,
          nrounds=params_xg$nrounds,
          nfold=params_xg$nfold,
          num_class = nSubtype,
          early_stopping_rounds=10,
          objective = "multi:softmax",
          verbose = ifelse(verbose, 1, 0)
        )
        best_iteration <- cvRes$best_iteration
        # mymusic()
      } else {
        if(verbose) LuckyVerbose('scaller: use self-defined best interaction round...')
      }

      # xgboost via best interation
      if(verbose) LuckyVerbose('scaller: xgboost model based on the best interation...')
      set.seed(seeds[2])
      scaller <- xgboost(
        params = params_xg3,
        data = dtrain,
        nrounds = best_iteration,
        num_class = nSubtype,
        objective = "multi:softmax",
        verbose = ifelse(verbose, 1, 0)
      )
      saveRDS(scaller, path_scaller)
    } else {
      if(verbose) LuckyVerbose('scaller: The subtype caller exists. Use it!')
      scaller <- readRDS(path_scaller)
      params_xg3 <- object@Data[['scaller.parameters']][['Params']]
    }
    object@Data[['scaller']] <- scaller


    # Validation
    if(verbose) LuckyVerbose('scaller: validate `scaller` in the internal cohort...')
    y2_valid_pred <- predict(scaller, res2_valid) # adjusted value
    # label_valid_pred <- CCS:::adjustXGBoostSubtype(object, label_valid_pred, verbose) # Do not use this step
    res.all <- CCS:::compareRealPred(y2_valid, y2_valid_pred, cluster_translator)
    if(verbose) LuckyVerbose('scaller: Performance -- multi_auc=', round(res.all$multi_auc[1], 4), '; accuracy=', round(res.all$accuracy[1], 4), '; kappa=', round(res.all$kappa[1], 4))
    object@Data[['scaller.performance']] <- res.all


    # Output
    params_xg4 <- params_xg3
    params_xg4[['best_iteration']] <- best_iteration
    object@Data[['scaller.parameters']] <- list(
      Prop = Prop,
      Seed = seed,
      cluster_translator = cluster_translator,
      Params = params_xg4
    )
    if(verbose) LuckyVerbose('scaller: All done!')
    return(object)

  }
)
