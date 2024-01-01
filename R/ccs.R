
#' @name CCS-class
#' @title CCS-class
#' @docType class
#' @description CCS class
#' @slot Repeat The data for repeatability
#' @slot Model a ggplot plot for CCS subtype visualization
#' @slot Data Data about CCS probability and CCS subtypes
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @exportClass CCS
#' @keywords classes
setClass(
  "CCS",
  slots = c(
    Repeat="list",
    Model="list",
    Data="list"
  ),
  prototype = list(
    Repeat = list(
      geneSet = list(),
      geneAnnotation = data.frame(),
      geneid = character(),
      params = list(),
      seed = numeric(),
      params.NbClust = list(),
      model.dir = character()
    ),
    Model = list(),
    Data = list(
      Probability = list(),
      CCS = integer()
    )
  )
)


#' @title Cohort Congress System
#' @description Cohort Congress System
#' @param data a list with components (expression matrix + subtype vector)
#' @param model.dir a character. the path of model series.
#' @param params parameters of \code{\link[xgboost]{xgb.cv}}, \code{\link[xgboost]{xgboost}}, and \code{\link[GSClassifier]{fitEnsembleModel}}.
#' @param dimension the hyperparameters for 2 level t-SNE.
#' @param random.cluster Whether to use random samples to decide the best subtype number. Default is \code{TRUE}, which is faster but more vague comparing with \code{FALSE} when \code{data} contains lots of samples (>10000, for example).
#' @param params.NbClust Parameters of \code{CCS:::NbClust2}. Possible options are:
#' \itemize{
#'   \item \code{distance} the distance measure to be used to compute the dissimilarity matrix. This must be one of: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski" or "NULL". By default, distance="euclidean". If the distance is "NULL", the dissimilarity matrix (diss) should be given by the user. If distance is not "NULL", the dissimilarity matrix should be "NULL".
#'   \item \code{method} the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
#'   \item \code{index} the index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP, Gamma, Gplus and Tau), "alllong" (all indices with Gap, Gamma, Gplus and Tau included). Cory's recommendation: one of "ch", "duda", "cindex", "gamma", "beale", "gap", and "all".
#'  \item \code{min.nc} minimal number of clusters, between 1 and (number of objects - 1)
#' }
#' @inheritParams GSClassifier::parCallEnsemble
#' @inheritParams drCCSProbability
#' @importFrom luckyBase LuckyVerbose is.one.true Fastextra
#' @importFrom Rtsne Rtsne
#' @import GSClassifier
#' @import xgboost
#' @return a CCS object
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @seealso \code{\link[xgboost]{xgb.cv}}; \code{\link[xgboost]{xgboost}}; \code{\link[NbClust]{NbClust}};
#' @references
#' \href{https://cran.r-project.org/web/packages/NbClust/index.html}{1. NbCluster package} \cr
#' \href{http://www.sthda.com/english/wiki/wiki.php?id_contents=7940}{2. A tutorial - DBSCAN: density-based clustering for discovering clusters in large datasets with noise} \cr
#' \href{https://satijalab.org/seurat/articles/get_started_v5_new}{Getting Started with Seurat • Seurat}
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' @export
ccs <- function(
    data,
    geneSet,
    geneAnnotation,
    geneid = "ensembl",
    params = list(
      # No. of cross validation subcohorts
      nfold = 5,
      # Max number of boosting iterations. GSClassifier have optimized nrounds, so here I set a large value.
      nrounds = 100,
      # No. of CPU cores used
      nthread = 2,
      # Step size shrinkage used in update to prevent overfitting.  [0,1]
      eta = 0.5,
      # The larger gamma is, the more conservative the algorithm will be. [0,∞]
      gamma = 0,
      # maximum depth of a tree. Deeper trees can capture more complex patterns in the data, but may also lead to overfitting. [0,∞]
      max_depth = 10,
      # percentage of columns used for each tree construction. Lowering this value can prevent overfitting by training on a subset of the features. (0, 1]
      colsample_bytree = 1,
      # The larger min_child_weight is, the more conservative the algorithm will be. [0,∞]
      min_child_weight = 1,
      # for preventing overfitting
      subsample = 0.7,
      # fitEnsembleModel
      n = 100,
      sampSize = 0.7,
      ptail = 0.2
    ),
    seed = 489,
    dimension = c(2,2),
    perplexity = 30,
    theta = 0.3,
    random.cluster = TRUE,
    params.NbClust = list(
      distance = "euclidean",
      min.nc = 2,
      max.nc= 10,
      method = "ward.D2",
      index = "all"
    ),
    model.dir = './ccs/project_01',
    verbose = TRUE,
    numCores = 4
){

  # Test
  if(F){

    library(luckyBase)
    nPac <- c('GSClassifier', 'xgboost', 'pROC', 'digest', 'irr','Rtsne','ggplot2'); Plus.library(nPac)
    testData <- readRDS(system.file("extdata", "testData.rds", package = "GSClassifier"))
    expr <- testData$PanSTAD_expr_part
    design <- testData$PanSTAD_phenotype_part

    # Simulated data
    cancer_ntype <- 4
    cancer_nCohort <- 10
    set.seed(123); cohort_nSample <- sample(100:500, cancer_nCohort*cancer_ntype, replace = T)
    set.seed(123); cohort_seeds <- sample(1:10000, cancer_nCohort*cancer_ntype, replace = T)
    data <- list()
    for(i in 1:cancer_ntype){ # i=1
      for(j in 1:cancer_nCohort){ # j=1
        set.seed(cohort_seeds[i*j])
        select_sample <- sample(1:ncol(expr), cohort_nSample[i*j],replace = F)
        data[[paste('cancer',i, sep = '_')]][[paste('cohort',i,j, sep = '.')]][['expr']] <- expr[,select_sample]
        y <- design[select_sample,]
        data[[paste('cancer',i, sep = '_')]][[paste('cohort',i,j, sep = '.')]][['subtype']] <- ifelse(y$PAD_subtype == 'PAD-I',1,ifelse(y$PAD_subtype == 'PAD-II',2,ifelse(y$PAD_subtype == 'PAD-III',3,ifelse(y$PAD_subtype == 'PAD-IV',4,NA))))
      }
    }

    data_test <- list(
      cancer_1 = list(
        cohort.1.1 = data[[1]][[1]],
        cohort.1.2 = data[[1]][[2]]
      ),
      cancer_2 = list(
        cohort.2.1 = data[[2]][[1]],
        cohort.2.2 = data[[2]][[2]]
      )
    )

    # Other parameters
    PADi <- readRDS(system.file("extdata", paste0('PAD.train_20220916.rds'), package = "GSClassifier"))
    geneSet <- PADi$geneSet
    geneAnnotation <- PADi$geneAnnotation
    geneid <- "ensembl"
    model.dir = './test/ccs/project_01'
    numCores = 6
    seed = 489
    verbose = T
    random.cluster = T
    min.nc = 2
    max.nc= 10
    dimension = c(2,2)
    perplexity = 30
    theta = 0.3

    # Model parameters
    params <- list(
      ## No. of cross validation subcohorts
      nfold = 5,

      ## Max number of boosting iterations. GSClassifier have optimized nrounds, so here I set a large value.
      nrounds = 100,

      ## No. of CPU cores used
      nthread = 2,

      ## Step size shrinkage used in update to prevent overfitting.  [0,1]
      eta = 0.5,

      ## The larger gamma is, the more conservative the algorithm will be. [0,∞]
      gamma = 0,

      ## maximum depth of a tree. Deeper trees can capture more complex patterns in the data, but may also lead to overfitting. [0,∞]
      max_depth = 10,

      ## percentage of columns used for each tree construction. Lowering this value can prevent overfitting by training on a subset of the features. (0, 1]
      colsample_bytree = 1,

      ## The larger min_child_weight is, the more conservative the algorithm will be. [0,∞]
      min_child_weight = 1,

      ## for preventing overfitting
      subsample = 0.7,

      # fitEnsembleModel
      n = 100,
      sampSize = 0.7,
      ptail = 0.2
    )

  }

  # Data name check
  data_allName <- c(names(data), as.character(unlist(lapply(data, names))))
  check <- is.one.true(grepl('[|]', data_allName))
  if(check){
    stop("Character '|' is used in your data names, which is not allowed. Stop ccs! Σ(°△ °|||)︴")
  }

  # Grobal seeds
  set.seed(seed); seeds <- sample(1:10000, 20, replace = T)
  dir.create(model.dir, showWarnings = F, recursive = T)

  # Parameters
  params_xg <- params[-match(c('n','sampSize','ptail'), names(params))]
  params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg))]
  path_ccs <- paste0(model.dir, "/resCCS.rds")

  # Model
  if(verbose) LuckyVerbose('Build GSClassifier models...')
  for(i in 1:length(data)){ # i=1

    data.i <- data[[i]]; cancer_type.i <- names(data)[i]

    for(j in 1:length(data.i)){ # j=1

      data.i.j <- data.i[[j]];

      cohort.j <- names(data.i)[j]

      project <- paste0(cancer_type.i, ' - ' ,cohort.j)

      if(verbose) LuckyVerbose('New project: ', project)

      path_child <- paste0(model.dir,'/model/',cancer_type.i, '/' ,cohort.j)
      dir.create(path_child, showWarnings = F, recursive = T)

      # fit
      path_fit <- paste0(path_child, '/modelFit.rds')
      if(!file.exists(path_fit)){
        modelFit <- fitEnsembleModel(
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
          xgboost.seed = seeds[3],
          caret.grid = NULL,
          ptail = params$ptail,
          verbose = T,
          numCores = numCores
        )
        saveRDS(modelFit, path_fit)
      } else {
        if(verbose) LuckyVerbose(project, ': modelFit exists. Ignored!')
      }

    }

  }

  # Cohort-based probability
  if(verbose) LuckyVerbose('Calculate cohort-based probability...')
  if(T){
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
        path_prob <- paste0(path_child,'/ccsProb.rds')
        if(!file.exists(path_prob)){
          a <- lapply(data, function(x) lapply(x, function(y) oneCCSProbability(y, path_model1, geneAnnotation, geneSet, geneid, numCores)))
          a2 <- NULL
          for(z in 1:length(a)){ # i=1
            a2 <- rbind(a2, do.call("rbind", a[[z]]))
          }
          a2 <- a2[!duplicated(a2$SampleIDs),] # Remove duplicated data
          colnames(a2)[2:ncol(a2)] <- paste(cancertype_model1, cohort_model1,  colnames(a2)[2:ncol(a2)], sep = '|')
          rownames(a2) <- as.character(a2$SampleIDs)
          saveRDS(a2, path_prob)
        } else {
          if(verbose) LuckyVerbose(cancer_type.i, '-' ,cohort.j, ': The result of ccsProb exists. Use it!')
          a2 <- readRDS(path_prob)
        }

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
  }

  # Dimensionality reduction
  if(verbose) LuckyVerbose('Dimensionality reduction via t-SNE...')
  if(T){
    path_dr <- paste0(model.dir, '/drMatrix.rds')
    if(!file.exists(path_dr)){
      # Dimensionality reduction - Level 1
      d2 <- drCCSProbability(
        res2,
        reference = Fastextra(colnames(res2), '[|]', 1),
        dims = dimension[1],
        perplexity = perplexity,
        theta = theta,
        seed = seeds[4],
        numCores = numCores,
        verbose = verbose
      ) # ; mymusic()

      # Dimensionality reduction - Level 2
      d3 <- drCCSProbability(
        d2,
        reference = NULL,
        dims = dimension[2],
        perplexity = perplexity,
        theta = theta,
        seed = seeds[5],
        numCores = numCores,
        verbose = verbose
      )
      saveRDS(list(d2=d2,d3=d3), path_dr)
    } else {
      if(verbose) LuckyVerbose('The result of dimensionality reduction exists. Use it!')
      dat_dr <- readRDS(path_dr)
      d2 <- dat_dr[['d2']]
      d3 <- dat_dr[['d3']]
    }
  }

  # CCS subtypes
  if(verbose) LuckyVerbose('Call CCS subtypes...')
  if(T){
    if(nrow(d3) > 1000 & random.cluster){
      set.seed(seeds[6])
      dat_cluster <- scale(d3[sample(1:nrow(d3), 1000),], center = T, scale = T)
      numComplete <- NbClust2(
        data = dat_cluster,
        distance = params.NbClust[['distance']],
        min.nc = params.NbClust[['min.nc']],
        max.nc = params.NbClust[['max.nc']],
        method = params.NbClust[['method']],
        index = params.NbClust[['index']],
        verbose = verbose
      )
      dis <- dist(scale(d3, center = T, scale = T), method = params.NbClust[['distance']])
      hc <- hclust(dis, method =  params.NbClust[['method']])
      y2 <- cutree(hc, length(unique(numComplete$Best.partition)))
    } else {
      dat_cluster <- scale(d3, center = T, scale = T)
      set.seed(seeds[6]);
      numComplete <- NbClust2(
        data = dat_cluster,
        distance = params.NbClust[['distance']],
        min.nc = params.NbClust[['min.nc']],
        max.nc = params.NbClust[['max.nc']],
        method = params.NbClust[['method']],
        index = params.NbClust[['index']],
        verbose = verbose
      )
      y2 <- numComplete$Best.partition
    }
  }

  # Subtype Caller
  if(verbose) LuckyVerbose('Build subtype caller ...')
  if(T){
    path_scaller <- paste0(model.dir, '/scaller.rds')
    if(!file.exists(path_scaller)){

      nSubtype <- length(unique(y2))
      dtrain <- xgb.DMatrix(res2, label = y2-1)

      # Parameters
      params_xg3 <- params_xg2
      params_xg3[['nthread']] <- numCores

      # xgboost cv for best interation exploration
      set.seed(seeds[7])
      cvRes <- xgb.cv(
        params = params_xg3,
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
      set.seed(seeds[7])
      scaller <- xgboost(
        params = params_xg3,
        data = dtrain,
        nrounds = best_iteration,
        num_class = nSubtype,
        objective = "multi:softmax",
        verbose = 0
      )
      saveRDS(scaller, path_scaller)
    } else {
      if(verbose) LuckyVerbose('The subtype caller exists. Use it!')
      scaller <- readRDS(path_scaller)
    }
  }

  # Output
  l <- new(
    'CCS',
    Repeat = list(
      geneSet = geneSet,
      geneAnnotation = geneAnnotation,
      geneid = geneid,
      params = params,
      seed = seed,
      params.NbClust = params.NbClust,
      model.dir = model.dir
    ),
    Model = list('light mode'),
    Data = list(
      Probability = list(
        d1 = res2,
        d2 = d2,
        d3 = d3
      ),
      CCS = y2,
      CancerType = cancerType(data,names(y2))
    )
  )
  saveRDS(l, path_ccs)
  if(verbose) LuckyVerbose('CCS: All done!')
  return(l)
}


#' @rdname CCS-method.predict
#' @title CCS method: predict
#' @description \code{predict} method for \code{CCS} class
#' @param object a CCS object
#' @inheritParams ccs
#' @inheritParams GSClassifier::parCallEnsemble
#' @import GSClassifier
#' @import xgboost
#' @importFrom luckyBase Fastextra
#' @return predict: A list object with CCS subtype prediction.
#' @seealso \code{\link{ccs}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ccs_pred <- predict(resCCS, model.dir = "./ccs/project_01")
#' @export
predict.CCS <- function(
    object, X,
    model.dir = NULL,
    verbose = T,
    numCores = 4){

  # Test
  if(F){
    library(luckyBase); library(GSClassifier)
    work.space = 'E:/RCloud/RFactory/CCS/test/ccs/project_01'
    object = readRDS(paste0(work.space, '/', 'resCCS.rds'))
    X = data[["cancer_1"]][["cohort.1.1"]][["expr"]][,1:5]
    model.dir = './test/ccs/project_01'
    numCores = 4
    verbose = F
  }

  # Model parameters
  if(is.null(model.dir)){
    model.dir = object@Repeat$model.dir
  }
  geneAnnotation = object@Repeat$geneAnnotation
  geneSet = object@Repeat$geneSet
  geneid = object@Repeat$geneid
  scaller = readRDS(paste0(model.dir,'/scaller.rds'))
  models = object@Model

  # Expression matrix
  X <- GSClassifier:::rightX(X)
  nSample <- ncol(X)

  # Call CCS probability
  X_CCSprobability <- NULL
  if(identical(list('light mode'), models)){
    if(verbose) LuckyVerbose('Light mode CCS. Use external model...')
    path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', recursive = T, full.names = T)
    for(i in 1:length(path_models)){ # i=1

      # A GSClassifier model
      path_models.i <- path_models[i]
      modelFit <- readRDS(path_models.i)
      name_model1 <- rev(Fastextra(path_models.i, '[/]'))
      cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
      if(verbose) LuckyVerbose('Load model of ', cancertype_model1, ' - ', cohort_model1 ," cohort : ")

      # Call CCS probability
      if(nSample > 100){
        X_CCSprobability.i <- parCallEnsemble(
          X = X,
          ens = modelFit$Model,
          geneAnnotation = geneAnnotation,
          geneSet = geneSet,
          geneid = geneid,
          scaller = NULL,
          subtype = NULL,
          verbose = verbose,
          numCores = numCores
        )
      } else {
        X_CCSprobability.i <- callEnsemble(
          X = X,
          ens = modelFit$Model,
          geneAnnotation = geneAnnotation,
          geneSet = geneSet,
          geneid = geneid,
          scaller = NULL,
          subtype = NULL,
          verbose = verbose
        )
      }

      colnames(X_CCSprobability.i)[3:ncol(X_CCSprobability.i)] <- paste(cancertype_model1, cohort_model1,  colnames(X_CCSprobability.i)[3:ncol(X_CCSprobability.i)], sep = '|')
      res.i <- t(apply(X_CCSprobability.i[-c(1,2,3)], 1, softmax))

      # Merge results
      if(i == 1){
        X_CCSprobability <- res.i
      } else {
        X_CCSprobability <- cbind(X_CCSprobability, res.i)
      }
    }
  } else {
    if(verbose) LuckyVerbose('Complete mode CCS. Use internal model...')
    for(i in 1:length(models)){ # i=1
      cancertype_model1 <- names(models)[i]
      model.i <- models[[i]]
      for(j in 1:length(model.i)){
        cohort_model1 <- names(model.i)[[j]]
        modelFit <- model.i[[j]]

        # Call CCS probability
        if(nSample > 100){
          X_CCSprobability.i <- parCallEnsemble(
            X = X,
            ens = modelFit$Model,
            geneAnnotation = geneAnnotation,
            geneSet = geneSet,
            geneid = geneid,
            scaller = NULL,
            subtype = NULL,
            verbose = verbose,
            numCores = numCores
          )
        } else {
          X_CCSprobability.i <- callEnsemble(
            X = X,
            ens = modelFit$Model,
            geneAnnotation = geneAnnotation,
            geneSet = geneSet,
            geneid = geneid,
            scaller = NULL,
            subtype = NULL,
            verbose = verbose
          )
        }

        colnames(X_CCSprobability.i)[3:ncol(X_CCSprobability.i)] <- paste(cancertype_model1, cohort_model1,  colnames(X_CCSprobability.i)[3:ncol(X_CCSprobability.i)], sep = '|')
        res.i <- t(apply(X_CCSprobability.i[-c(1,2,3)], 1, softmax))

        # Merge results
        if(i == 1 & j == 1){
          X_CCSprobability <- res.i
        } else {
          X_CCSprobability <- cbind(X_CCSprobability, res.i)
        }
      }
    }
  }

  X_CCS_Pred <- predict(scaller, X_CCSprobability) + 1
  names(X_CCS_Pred) <- colnames(X)


  # Output
  l <- list(
    X = X,
    model.dir = model.dir,
    CCS = list(
      Probability = cbind(SampleIDs = colnames(X), as.data.frame(X_CCSprobability)),
      Prediction = X_CCS_Pred
    )
  )
  if(verbose) LuckyVerbose('All done!')
  return(l)
}


#' @rdname CCS-method.plot
#' @title CCS method: plot
#' @description \code{plot} method for \code{CCS} class
#' @param object a CCS object
#' @param size the size of ggplot
#' @param CCS Character. A vector of samples' CCS subtype.
#' @param geom Some of \code{'cancer_type'} and \code{'CCS'}.
#' @import ggplot2
#' @importFrom luckyBase mycolor
#' @return plot: A ggplot2 object
#' @seealso \code{\link{ccs}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ccs_plot <- plot(resCCS)
#' @export
plot.CCS <- function(
    object,
    CCS = NULL,
    geom = c('cancer_type','CCS')[2],
    size = 15){

  # Test
  if(F){
    library(luckyBase)
    Plus.library(c('ggplot2'))
    model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225'
    object = readRDS(paste0(model.dir, '/resCCS.rds'))
    CCS = NULL
    geom = c('cancer_type','CCS')
    size = 15
  }

  # Data
  dat_plot <- as.data.frame(object@Data[["Probability"]][["d3"]])
  if(is.null(CCS)){
    y2 <- object@Data[["CCS"]]
  } else {
    y2 <- CCS
  }

  # plot head
  if('CCS' %in% geom){
    dat_plot <- cbind(dat_plot, CCS = paste('CCS',y2,sep = ''))
    gghead <-
      ggplot(dat_plot, aes(x = `all|D1`, y = `all|D2`, group = CCS)) +
      geom_point(aes(color=CCS), size = size/15*3, shape = 1, stroke = size/15*1.5) +
      labs(title = "t-SNE Visualization", color = 'CCS')
  }
  if('cancer_type' %in% geom){
    dat_plot <- cbind(dat_plot, cancer_type = object@Data[["CancerType"]])
    unique_cancertype <- unique(dat_plot$cancer_type)
    unique_csstype <- unique(dat_plot$CCS)
    # http://www.sthda.com/english/wiki/ggplot2-point-shapes
    default_shape <- c(1,2,4); default_shape <- c(default_shape, setdiff(1:25, default_shape))
    gghead <-
      ggplot(dat_plot, aes(x = `all|D1`, y = `all|D2`, group = CCS)) +
      geom_point(aes(color=CCS, shape = cancer_type), size = size/15*3) +
      scale_color_manual(values = mycolor[1:length(unique_csstype)], breaks = unique_csstype) +
      scale_shape_manual(values = default_shape[1:length(unique_cancertype)], breaks = unique_cancertype) +
      labs(title = "t-SNE Visualization", color = 'CCS', shape = 'Cancer')
  }

  # Plot complete
  if(T){
    p <- gghead +
      theme_bw() +
      theme(
        axis.text = element_text(size = size/15*12,colour = "black",face = "bold"),
        plot.title = element_text(size = size,colour = "black",face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = size,colour = "black",face = "bold"),
        axis.title.y = element_text(size = size,colour = "black",face = "bold"),
        legend.text = element_text(size = size/15*12,colour = "black",face = "bold"),
        legend.title = element_text(size = size/15*12,colour = "black",face = "bold"),
        legend.position='right',
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = size/15*12,colour = "black",face = "bold")); # win.graph(10,10); print(p)
  }

  # Output
  return(p)
}



setGeneric("completeModel", function(object, ...) {
  standardGeneric("completeModel")
})

#' @rdname CCS-method.completeModel
#' @title CCS method: completeModel
#' @description \code{completeModel} method for \code{CCS} class
#' @inheritParams ccs
#' @importFrom luckyBase Fastextra
#' @return completeModel: a complete CCS class.
#' @exportMethod completeModel
setMethod(
  "completeModel",
  signature(object='CCS'),
  function(object,
           model.dir = NULL){

    # Test
    if(F){
      model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225'
      object = readRDS(paste0(model.dir, '/resCCS.rds'))
    }

    # Model parameters
    if(is.null(model.dir)){
      model.dir = object@Repeat$model.dir
    }

    # Load model
    path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', recursive = T, full.names = T)
    object@Model <- list()
    for(i in 1:length(path_models)){ # i=1
      path_model1 <- path_models[i]
      name_model1 <- rev(Fastextra(path_model1, '[/]'))
      cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
      object@Model[[cancertype_model1]][[cohort_model1]] <- readRDS(path_model1)
    }

    # Output
    return(object)

  }
)



