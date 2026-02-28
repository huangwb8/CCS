
#' @description Core function of sub-model establishment in CCS
#' @param method method for the establishment of sub-models. One of \code{'GSClassifier'}
#' @inheritParams ccsSubModel_GSClassifier
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
ccsSubModel <- function(
    data,
    model.dir,
    method = 'GSClassifier',
    geneSet,
    params,
    seeds,
    numCores,
    parallel.method = c('ensemble','discrete')[1],
    verbose
  ){

  # Information
  if(verbose) LuckyVerbose('ccsSubModel: Build ',method,' models...')

  # Parameters
  params_xg <- params[-match(c('n','sampSize','ptail'), names(params))]
  # params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg))]

  # Process
  for(i in 1:length(data)){ # i=1

    data.i <- data[[i]]; cancer_type.i <- names(data)[i]

    for(j in 1:length(data.i)){ # j=1

      data.i.j <- data.i[[j]];

      cohort.j <- names(data.i)[j]

      project <- paste0(cancer_type.i, ' - ' ,cohort.j)

      if(verbose) LuckyVerbose('ccsSubModel: New project - ', project)

      path_child <- paste0(model.dir,'/model/',cancer_type.i, '/' ,cohort.j)
      dir.create(path_child, showWarnings = F, recursive = T)

      # fit
      path_fit <- paste0(path_child, '/modelFit.rds')
      if(!file.exists(path_fit)){

        if(method == 'GSClassifier'){
          modelFit <- CCS:::ccsSubModel_GSClassifier(
            data.i.j,
            geneSet,
            params, params_xg,
            seeds,
            numCores
          )
          saveRDS(modelFit, path_fit)
        } else {
          if(verbose) LuckyVerbose('ccsSubModel: Wrong model establishment methods. Stop CCS process! Σ( ° △ °|||)︴')
        }

      } else {
        if(verbose) LuckyVerbose('ccsSubModel: ',project, ' - modelFit exists. Ignored!')
      }

    }

  }

}


#### GSClassifier ####
#' @param data a list with components (expression matrix + subtype vector)
#' @param model.dir a character. the path of model series.
#' @param params parameters of \code{\link[xgboost]{xgb.cv}}, \code{\link[xgboost]{xgboost}}, and \code{\link[GSClassifier]{fitEnsembleModel}}. Some important options are:
#' \itemize{
#'   \item \code{nfold} No. of cross validation subcohorts.
#'   \item \code{nrounds} Max number of boosting iterations. GSClassifier have optimized nrounds, so here I set a large value.
#'   \item \code{nthread} No. of CPU cores used in \code{\link[xgboost]{xgboost}}.
#'   \item \code{eta} Step size shrinkage used in update to prevent overfitting. Range = [0,1]
#'   \item \code{gamma} The larger gamma is, the more conservative the algorithm will be. Range = [0,∞]
#'   \item \code{max_depth} Maximum depth of a tree. Deeper trees can capture more complex patterns in the data, but may also lead to overfitting. Range = [0,∞]
#'   \item \code{colsample_bytree} Percentage of columns used for each tree construction. Lowering this value can prevent overfitting by training on a subset of the features. Range = (0, 1]
#'   \item \code{min_child_weight} The larger min_child_weight is, the more conservative the algorithm will be. Range = [0,∞]
#'   \item \code{subsample} Preventing overfitting
#'   \item \code{n} \code{\link[GSClassifier]{fitEnsembleModel}} parameter
#'   \item \code{sampSize} \code{\link[GSClassifier]{fitEnsembleModel}} parameter. Range = (0,1]
#'   \item \code{ptail} \code{\link[GSClassifier]{fitEnsembleModel}} parameter. Range = (0,0.5]
#' }
#' @inheritParams GSClassifier::fitEnsembleModel
#' @importFrom GSClassifier fitEnsembleModel
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
ccsSubModel_GSClassifier <- function(
    data.i.j,
    geneSet,
    params, params_xg,
    seeds,
    numCores
){

  nround.mode = {
    if(is.null(params$nround.mode)){
      nround.mode = 'polling'
    } else {
      nround.mode = params$nround.mode
    }
    nround.mode
  }

  modelFit <- GSClassifier::fitEnsembleModel(
    Xs = data.i.j$expr,
    Ys = data.i.j$subtype,
    geneSet = geneSet,
    na.fill.method = "quantile",
    na.fill.seed = seeds[1],
    n = params$n,
    sampSize = params$sampSize,
    sampSeed = seeds[2],
    breakVec = c(0, 0.25, 0.5, 0.75, 1),
    params = params_xg,
    nround.mode = nround.mode,
    xgboost.seed = seeds[3],
    caret.grid = NULL,
    ptail = params$ptail,
    verbose = T,
    numCores = numCores
  )

  # model size control
  # modelFit[['Repeat']] <- NA

  # Output
  return(modelFit)
}
