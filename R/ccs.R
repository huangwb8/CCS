
#' @import ggplot2
#' @author Weibin Huang<\email{654751191@@qq.com}>
setOldClass(c("gg","ggplot"))

#' @name CCS-class
#' @title CCS-class
#' @docType class
#' @description CCS class
#' @slot Repeat The data for repeatability
#' @slot Data Data about CCS probability and CCS subtypes
#' @slot Plot a ggplot plot for CCS subtype visualization
#' @import ggplot2
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @exportClass CCS
#' @keywords classes
setClass("CCS",
         slots = c(
           Repeat="list",
           Data="list",
           Plot="gg"
         ),
         prototype = list(
           Repeat = list(
             data = list(),
             geneSet = list(),
             geneAnnotation = data.frame(),
             geneid = character(),
             params = list(),
             seed = numeric(),
             min.nc = numeric(),
             max.nc= numeric(),
             model.dir = character()
           ),
           Data = list(
             Probability = list(),
             CCS = integer(),
             scaller = list()
           ),
           Plot = ggplot()
         )
)


#' @title Cohort Congress System
#' @description Cohort Congress System
#' @param data a list with components (expression matrix + subtype vector)
#' @param model.dir a character. the path of model series.
#' @param params parameters like X function
#' @param dimension the hyperparameters for 2 level t-SNE.
#'
#' @inheritParams GSClassifier::parCallEnsemble
#' @inheritParams drCCSProbability
#' @inheritParams NbClust2
#' @importFrom luckyBase LuckyVerbose is.one.true Fastextra
#' @importFrom Rtsne Rtsne
#' @import ggplot2
#' @import GSClassifier
#' @import xgboost
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
    dimension = c(2,2),
    perplexity = 30,
    theta = 0.3,
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
    model.dir = './test/ccs/project_01'
    numCores = 6
    seed = 489
    verbose = T
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

  # Model
  for(i in 1:length(data)){ # i=1

    data.i <- data[[i]]; cancer_type.i <- names(data)[i]

    for(j in 1:length(data.i)){ # j=1

      data.i.j <- data.i[[j]];

      cohort.j <- names(data.i)[j]

      project <- paste0(cancer_type.i, ' - ' ,cohort.j)

      if(verbose) LuckyVerbose('New project: ', project)

      path_child <- paste0(model.dir,'/',cancer_type.i, '/' ,cohort.j)
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
        modelFit <- readRDS(path_fit)
      }

    }

  }

  # Cohort-based probability
  path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', full.names = T, recursive = T)
  path_resCBP <- paste0(model.dir, '/Cohort-based probability.rds')
  if(!file.exists(path_resCBP)){
    res <- data.frame()
    for(i in 1:length(path_models)){ # i=1
      path_model1 <- path_models[i]
      name_model1 <- rev(Fastextra(path_model1, '[/]'))
      cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
      a <- lapply(data, function(x) lapply(x, function(y) oneCCSProbability(y, path_model1)))
      a2 <- do.call("rbind", do.call("rbind", a))
      # a2 <- res[c(1:5)]; colnames(a2)[2:5] <- 1:4
      colnames(a2)[2:ncol(a2)] <- paste(cancertype_model1, cohort_model1,  colnames(a2)[2:ncol(a2)], sep = '|')
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
  res2 <- res[,-1]; res2 <- matrix(as.numeric(as.matrix(res2)), nrow = nrow(res2), byrow = F, dimnames = list(rownames(res2), colnames(res2)))
  d1 <- data_for_tSNE(res2, verbose);


  # Dimensionality reduction
  if(T){
    # Dimensionality reduction - Level 1
    d2 <- d1$cleaned$data
    reference <- Fastextra(colnames(d2), '[|]', 1)
    d3 <- drCCSProbability(
      d2,
      reference = reference,
      dims = dimension[1],
      perplexity = perplexity,
      theta = theta,
      seed = seeds[4],
      verbose = verbose
    )

    # Dimensionality reduction - Level 2
    d4 <- drCCSProbability(
      d3,
      reference = NULL,
      dims = dimension[2],
      perplexity = perplexity,
      theta = theta,
      seed = seeds[5],
      verbose = verbose
    )
  }


  # CCS subtypes
  if(T){
    if(nrow(d4) > 500){
      set.seed(seeds[6]);
      numComplete <- NbClust2(
        data = d4[sample(1:nrow(d4), 500),],
        distance = "euclidean",
        min.nc = min.nc,
        max.nc= max.nc,
        method = "ward.D2",
        index = "all",
        verbose = verbose
      )
      dis <- dist(d4, method = "euclidean")
      hc <- hclust(dis, method =  "ward.D2")
      y2 <- cutree(hc, length(unique(numComplete$Best.partition)))
    } else {
      set.seed(seeds[6]);
      numComplete <- NbClust2(
        data = d4,
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
    names(y3) <- as.character(res$SampleIDs)
  }


  # Plot
  if(verbose) LuckyVerbose('Dimensionality reduction visulization via ggplot2...')
  if(T){
    dat_plot <- cbind(d4, CCS = paste('CCS',y2,sep = ''))
    size = 15 # plot size
    p <- ggplot(dat_plot, aes(x = `all|D1`, y = `all|D2`, color = CCS)) +
      geom_point() +
      labs(title = "t-SNE Visualization") +
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


  # Subtype Caller
  if(verbose) LuckyVerbose('Build subtype caller ...')
  path_scaller <- paste0(model.dir, '/scaller.rds')
  if(!file.exists(path_scaller)){

    nSubtype <- length(unique(y2))
    dtrain <- xgb.DMatrix(d2, label = y2-1)

    # xgboost cv for best interation exploration
    set.seed(seeds[7])
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
    set.seed(seeds[7])
    scaller <- xgboost(
      params = params_xg2,
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

  # Output
  # l <- list(
  #   Repeat = list(
  #     data = data,
  #     geneSet = geneSet,
  #     geneAnnotation = geneAnnotation,
  #     geneid = geneid,
  #     params = params,
  #     seed = seed,
  #     min.nc = min.nc,
  #     max.nc= max.nc,
  #     model.dir = model.dir
  #   ),
  #   Data = list(
  #     Probability = list(
  #       raw = res,
  #       d2 = d3,
  #       d3 = d4
  #     ),
  #     CCS = y3
  #   ),
  #   Plot = p
  # )

  # Make a new Animal object
  l <- new(
    'CCS',
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
      Probability = list(
        raw = res,
        d2 = d3,
        d3 = d4
      ),
      CCS = y3,
      scaller = scaller
    ),
    Plot = p
  )
  if(verbose) LuckyVerbose('All done!')
  return(l)
}


#' @title CCS-method
#' @name CCS-method
#' @description \code{predict} method for \code{CCS} class
#' @param object a CCS object
#' @inheritParams ccs
#' @inheritParams GSClassifier::parCallEnsemble
#' @import GSClassifier
#' @import xgboost
#' @importFrom luckyBase Fastextra
#' @return A list object with CCS subtype prediction.
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
  scaller = object@Data$scaller

  # Expression matrix
  X <- GSClassifier:::rightX(X)
  nSample <- ncol(X)

  # Call CCS probability
  path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', recursive = T, full.names = T)
  X_CCSprobability <- NULL
  for(i in 1:length(path_models)){ # i=1

    # A GSClassifier model
    path_models.i <- path_models[i]
    modelFit <- readRDS(path_models.i)
    name_model1 <- rev(Fastextra(path_models.i, '[/]'))
    cohort_model1 <- name_model1[2]; cancertype_model1 <- name_model1[3]
    LuckyVerbose('Load model of ', cancertype_model1, ' - ', cohort_model1 ," cohort : ")

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
        verbose = F,
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
        verbose = F
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
  return(l)
}




