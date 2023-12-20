

#' @description OneModel-OneData CSS process
#' @param data1 a list containing an expression matrix and its subtype vector
#' @param path_model1 the path of a (GSClassifier) model
#' @importFrom GSClassifier parCallEnsemble
#' @return a data frame with sample IDS and the softmax probability.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
oneCCSProbability <- function(data1, path_model1){

  # Test
  if(F){
    data1 = data[[1]][[1]]
    path_model1 = "./ccs/test01/cohort.1.1/modelFit.rds"
  }

  # Call probability score
  modelFit <- readRDS(path_model1)
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

  # Output
  res <- as.data.frame(
    cbind(
      SampleIDs = as.character(res[,'SampleIDs']),
      t(apply(res[-c(1,2,3)], 1, softmax))
    )
  )

  # return(res[,-match(c("BestCall","BestCall_Max"), colnames(res))])
  return(res)

}


#' @description convert an expression data for t-SNE
#' @param res2 The output of \code{\link{oneCCSProbability}}
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom digest digest
#' @return a list with raw/cleaned data or md5 sum.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
data_for_tSNE <- function(res2,verbose){

  res2_md5 <- apply(res2, 1, function(x)digest(x, algo="md5"))
  s <- !duplicated(res2_md5)
  if(verbose){
    LuckyVerbose('No. of Duplicati sample = ',sum(!s),'. Remove them!')
  }
  return(list(
    raw = list(md5=res2_md5, data=res2),
    cleaned = list(md5=res2_md5[s], data=res2[s,])
  ))
}


#' @title Cohort Congress System
#' @description Cohort Congress System
#' @param data a list with components (expression matrix + subtype vector)
#' @param model.dir a character. the path of model series.
#' @param params parameters like X function
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom luckyBase LuckyVerbose
#' @importFrom Rtsne Rtsne
#' @import ggplot2
#' @import GSClassifier
#' @return a CCS object
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
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
    min.nc = 2,
    max.nc= 10,
    model.dir = './ccs/project_01',
    verbose = T,
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
    model.dir = './ccs/test01'
    numCores = 6
    seed = 489
    verbose = T
    min.nc = 2
    max.nc= 10

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

  # Grobal seeds
  set.seed(seed); seeds <- sample(1:10000, 20, replace = T)
  dir.create(model.dir, showWarnings = F, recursive = T)

  # Parameters
  params_xg <- params[-match(c('n','sampSize','ptail'), names(params))]
  params_xg2 <- params_xg[-match(c('nfold','nrounds'), names(params_xg))]

  # Model
  for(i in 1:length(data)){ # i=1

    data.i <- data[[i]]

    for(j in 1:length(data.i)){ # j=1

      data.i.j <- data.i[[j]]

      project <- names(data.i)[j]

      LuckyVerbose('New project: ', project)

      path_child <- paste0(model.dir,'/', project)
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
        LuckyVerbose(project, ': modelFit exists. Ignored!')
        modelFit <- readRDS(path_fit)
      }

    }

  }

  # Cohort-based probability
  path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', full.names = T, recursive = T)
  res <- data.frame()
  for(i in 1:length(path_models)){ # i=1
    path_model1 <- path_models[i]
    name_model1 <- rev(Fastextra(path_model1, '[/]'))[2]
    # a <- lapply(data_test, function(x) lapply(x, function(y) css_one(y, path_model1)))
    a <- lapply(data, function(x) lapply(x, function(y) oneCCSProbability(y, path_model1)))
    a2 <- do.call("rbind", do.call("rbind", a))
    colnames(a2)[2:ncol(a2)] <- paste(name_model1, colnames(a2)[2:ncol(a2)], sep = '_')
    if(i==1){
      res <- a2
    } else {
      res <- cbind(res, a2[,-1])
    }

  }
  saveRDS(res, paste0(model.dir, '/Cohort-based probability.rds'))
  # res <- readRDS(paste0(model.dir, '/Cohort-based probability.rds'))
  res2 <- res[,-1]; res2 <- matrix(as.numeric(as.matrix(res2)), nrow = nrow(res2), byrow = F, dimnames = list(rownames(res2), colnames(res2)))


  # Dimensionality reduction
  d1 <- data_for_tSNE(res2, verbose)
  set.seed(seeds[4])
  tsne_result <- Rtsne(
    # scale(d1$cleaned$data, center = T, scale = T),
    d1$cleaned$data,
    dims = 2, perplexity = 30, theta = 0,
    verbose = verbose
  )
  tsne_data <- as.data.frame(tsne_result$Y)
  colnames(tsne_data) <- c("Dimension 1", "Dimension 2")
  p <- ggplot(tsne_data, aes(x = `Dimension 1`, y = `Dimension 2`)) +
    geom_point() +
    ggtitle("t-SNE Visualization")
  # win.graph(); print(p)


  # CCS subtypes
  if(nrow(tsne_data) > 500){
    set.seed(178);
    numComplete <- NbClust2(
      data = tsne_data[sample(1:nrow(tsne_data), 500),],
      distance = "euclidean",
      min.nc = min.nc,
      max.nc= max.nc,
      method = "ward.D2",
      index = "all",
      verbose = verbose
    )
    dis <- dist(tsne_data, method = "euclidean")
    hc <- hclust(dis, method =  "ward.D2")
    y2 <- cutree(hc, length(unique(numComplete$Best.partition)))
  } else {
    numComplete <- NbClust2(
      data = tsne_data,
      distance = "euclidean",
      min.nc = min.nc,
      max.nc= max.nc,
      method = "ward.D2",
      index = "all",
      verbose = verbose
    )
    y2 <- numComplete$Best.partition
  }
  index <- match(d1$raw$md5, d1$cleaned$md5)
  y3 <- NULL
  for(i in 1:length(index)){
    y3 <- c(y3, y2[index[i]])
  }


  # Output
  l <- list(
    Repeat = list(
      data = data,
      geneSet = geneSet,
      geneAnnotation = geneAnnotation,
      geneid = geneid,
      params = params,
      seed = seed,
      min.nc = min.nc,
      max.nc= max.nc,
      model.dir = model.dir
    ),
    Data = list(
      Probability = res,
      CCS = y3
    ),
    Plot = p
  )
  return(l)
}
