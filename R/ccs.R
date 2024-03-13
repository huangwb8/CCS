
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
#' @inheritParams GSClassifier::parCallEnsemble
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
      nfold = 5,
      nrounds = 100,
      nthread = 2,
      eta = 0.5,
      gamma = 0,
      max_depth = 10,
      colsample_bytree = 1,
      min_child_weight = 1,
      subsample = 0.7,
      n = 100,
      sampSize = 0.7,
      ptail = 0.2
    ),
    seed = 489,
    model.dir = './ccs/project_01',
    verbose = TRUE,
    numCores = 4
){

  # Data name check
  data_allName <- c(names(data), as.character(unlist(lapply(data, names))))
  check <- is.one.true(grepl('[|]', data_allName))
  if(check){
    stop("Character '|' is used in your data names, which is not allowed. Stop ccs! Σ(°△ °|||)︴")
  }

  # Grobal seeds
  set.seed(seed); seeds <- sample(1:10000, 20, replace = T)
  dir.create(model.dir, showWarnings = F, recursive = T)
  path_tmp <- paste0(model.dir,'/tmp')
  dir.create(path_tmp, showWarnings = F, recursive = T)

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
          # a <- lapply(data, function(x) lapply(x, function(y) oneCCSProbability(y, path_model1, geneAnnotation, geneSet, geneid, numCores, verbose = T))) # No detail. Work, but not I want.
          a <- NULL
          for(j in 1:length(data)){
            data_cancer <- data[[j]]; cancer_name <- names(data)[j]
            for(k in 1:length(data_cancer)){
              data_cohort <- data_cancer[[k]]; cohort_name <- names(data_cancer)[k]
              path_a_tmp <- paste0(path_tmp, '/oneCCSProbabilityResult_Model-',cancertype_model1,'-',cohort_model1,'_Data-',cancer_name, ' - ',cohort_name,'.rds')
              if(!file.exists(path_a_tmp)){
                a[[cancer_name]][[cohort_name]] <- a_tmp <- oneCCSProbability(data_cohort, path_model1, geneAnnotation, geneSet, geneid, numCores, dataName = paste0(cancer_name, ' - ',cohort_name),verbose = T)
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

  # Output
  l <- new(
    'CCS',
    Repeat = list(
      geneSet = geneSet,
      geneAnnotation = geneAnnotation,
      geneid = geneid,
      params = params,
      seed = seed,
      model.dir = model.dir
    ),
    Model = list(NA),
    Data = list(
      Probability = list(
        d1 = res2,
        d2 = NA,
        d3 = NA
      ),
      CCS = NA,
      CancerType = cancerType(data, rownames(res2))
    )
  )
  saveRDS(l, path_ccs)
  if(verbose) LuckyVerbose('CCS: All done!')
  return(l)
}



#' @rdname CCS-method.predict
#' @title CCS method: predict
#' @description \code{predict} method for \code{CCS} class
#' @inheritParams CCSPublicParams
#' @inheritParams GSClassifier::parCallEnsemble
#' @import GSClassifier
#' @import xgboost
#' @import tidyr
#' @importFrom luckyBase Fastextra convert
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
  scaller = object@Data$scaller
  # scaller = readRDS(paste0(model.dir,'/scaller.rds'))
  models = object@Model
  cluster_translator = object@Data[['scaller.parameters']][['cluster_translator']]

  # Check integrity
  if(is.null(scaller)){
    stop('Lack subtype caller. Stop!')
  }

  # Expression matrix
  X <- GSClassifier:::rightX(X)
  nSample <- ncol(X)

  # Call CCS probability
  X_CCSprobability <- NULL
  if(identical(models, list(NA))){
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

  # Prediction
  X_CCS_Pred <- predict(scaller, X_CCSprobability)
  # X_CCS_Pred <- CCS:::adjustXGBoostSubtype(object, X_CCS_Pred, verbose)
  X_CCS_Pred <- convert(X_CCS_Pred, 'adjust', 'raw', cluster_translator) %>% as.integer()
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
#' @description \code{plot} method for \code{CCS} class. Plot dimension scatter plot for CCS results.
#' @param CCS Character. A vector of samples' CCS subtype.
#' @inheritParams CCSPublicParams
#' @import ggplot2
#' @import luckyBase
#' @return plot: A ggplot2 object
#' @seealso \code{\link{ccs}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ccs_plot <- plot(resCCS)
#' @export
plot.CCS <- function(
    object,
    CCS = NULL,
    geom = c('cancer_type','CCS'),
    hide.legend = c('cancer_type','CCS')[2],
    rm.zero = TRUE,
    size = 15,
    verbose = TRUE){

  # Test
  if(F){
    library(luckyBase)
    Plus.library(c('ggplot2','tidyr','dplyr','circlize','ComplexHeatmap'))
    model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225'
    object = readRDS(paste0(model.dir, '/resCCS.rds'))
    CCS = NULL
    geom = c('cancer_type','CCS')[2]
    hide.legend = c('cancer_type','CCS')
    size = 15
    rm.zero = TRUE
    verbose = TRUE
  }

  # Data
  dat_plot <- as.data.frame(object@Data[["Probability"]][["d3"]])
  if(is.null(CCS)){
    y2 <- object@Data[["CCS"]]
  } else {
    y2 <- CCS
  }
  cancer_type <- object@Data[["CancerType"]]

  # Remove zero value
  if(rm.zero){
    target_sample <- !y2 %in% 0
    dat_plot <- dat_plot[target_sample,]
    y2 <- y2[target_sample]
    cancer_type <- cancer_type[target_sample]
    if(verbose) LuckyVerbose('plot.CCS: Remove ', sum(!target_sample), ' zero value in CCS subtypes...')
  }

  # plot head
  dat_plot <- cbind(dat_plot, CCS = paste('CCS',y2,sep = ''))
  unique_csstype <- unique(dat_plot$CCS)
  n_ccstype <- length(unique_csstype)
  if(n_ccstype > length(mycolor)){
    default_color <- c(mycolor, setdiff(scales::hue_pal()(n_ccstype), mycolor))[1:n_ccstype]
  } else {
    default_color <- mycolor[1:n_ccstype]
  }
  if('CCS' %in% geom){
    gghead <-
      ggplot(dat_plot, aes(x = `all|D1`, y = `all|D2`, group = CCS)) +
      geom_point(aes(color=CCS), size = size/15*3, shape = 16, stroke = size/15*1.5) + # shape = 1
      scale_color_manual(values = default_color, breaks = unique_csstype) +
      labs(title = "", color = 'CCS')
  }
  if('cancer_type' %in% geom){
    dat_plot <- cbind(dat_plot, cancer_type = cancer_type)
    unique_cancertype <- unique(dat_plot$cancer_type)
    # http://www.sthda.com/english/wiki/ggplot2-point-shapes
    default_shape <- c(1,2,4); default_shape <- c(default_shape, setdiff(1:25, default_shape))
    gghead <-
      ggplot(dat_plot, aes(x = `all|D1`, y = `all|D2`, group = CCS)) +
      geom_point(aes(color=CCS, shape = cancer_type), size = size/15*3) +
      scale_color_manual(values = default_color, breaks = unique_csstype) +
      scale_shape_manual(values = default_shape[1:length(unique_cancertype)], breaks = unique_cancertype) +
      labs(title = "", color = "CCS", shape = 'Cancer')
  }


  # Legend
  if('CCS' %in% hide.legend){
    gghead <- gghead + guides(color = "none", shape = "legend")
  }
  if('cancer_type' %in% hide.legend){
    gghead <- gghead + guides(color = "legend", shape = "none")
  }
  if(all(c('CCS','cancer_type') %in% hide.legend)){
    gghead <- gghead + guides(color = "none", shape = "none")
  }

  # Plot complete
  if(T){
    p1 <- gghead +
      labs(x = 'Dimension 1', y = 'Dimension 2') +
      theme_bw() +
      theme(
        axis.text = element_text(size = size/15*12,colour = "black",face = "bold"),
        plot.title = element_text(size = size,colour = "black",face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = size,colour = "black",face = "bold"),
        axis.title.y = element_text(size = size,colour = "black",face = "bold"),
        legend.text = element_text(size = size/15*12,colour = "black",face = "bold"),
        legend.title = element_text(size = size/15*12,colour = "black",face = "bold"),
        legend.position='right',
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = size/15*2),
        axis.ticks = element_line(colour = "black",linewidth = size/15*1.4,linetype = 1, lineend = 'square'),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = size/15*12,colour = "black",face = "bold")); # win.graph(10,7); print(p1)
  }

  # Output
  if(verbose) LuckyVerbose('plot.CCS: All done!')
  return(p1)
}



setGeneric("plotCancerSubtype", function(object, ...) {
  standardGeneric("plotCancerSubtype")
})

#' @rdname CCS-method.plotCancerSubtype
#' @title CCS method: plotCancerSubtype
#' @description \code{plotCancerSubtype} method for \code{CCS} class. Plot CCS subtype across cancer types.
#' @inheritParams CCSPublicParams
#' @import ComplexHeatmap
#' @importFrom luckyBase LuckyVerbose
#' @importFrom grid gpar grid.text
#' @return plotCancerSubtype: a complete CCS class.
#' @exportMethod plotCancerSubtype
setMethod(
  "plotCancerSubtype",
  signature(object='CCS'),
  function(object,
           CCS = NULL,
           rm.zero = TRUE,
           size = 15,
           verbose = TRUE){

    # Test
    if(F){
      library(luckyBase)
      Plus.library(c('tidyr','dplyr','circlize','ComplexHeatmap'))
      model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240106'
      object = readRDS(paste0(model.dir, '/resCCS.rds'))
      CCS = NULL
      rm.zero = TRUE
      size = 15
      verbose = TRUE
    }

    # Data
    if(is.null(CCS)){
      y2 <- object@Data[["CCS"]]
    } else {
      y2 <- CCS
    }
    cancer_type <- object@Data[["CancerType"]]
    # y2 <- y2[names(cancer_type)]

    # Remove zero value
    if(rm.zero){
      target_sample <- !y2 %in% 0
      y2 <- y2[target_sample]
      cancer_type <- cancer_type[target_sample]
      if(verbose) LuckyVerbose('plotCancerSubtype: Remove ', sum(!target_sample), ' zero value in CCS subtypes...')
    }

    # Data
    hm.dat <- table(y2, cancer_type) %>% as.data.frame() %>% spread(cancer_type, Freq)
    hm_x <- as.matrix(hm.dat[,-1]); rownames(hm_x) <- as.character(hm.dat[,1]); hm_x[hm_x == 0] <- NA;
    # hm_x2 <- as.matrix(as.vector(hm_x)/sum(hm_x, na.rm = T)) %>% scale() %>% as.vector()
    # hm_x3 <- matrix(hm_x2, nrow = nrow(hm_x), byrow = F, dimnames = list(rownames(hm_x), colnames(hm_x)))
    hm_x3 <- scale(hm_x)

    # Cell function
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(ifelse(is.na(hm_x[i, j]), '', hm_x[i, j]), x, y, gp = gpar(fontsize = size/15*5))
    }

    # Heatmap
    # col = colorRamp2(c(0,0.01,0.03,0.05,0.06,0.07,0.08,0.11), c("#3300FF","#0066FF","#00FFFF","#00FF66","#33FF00","#CCFF00","#FF9900","#FF0000"))
    p2 <- Heatmap(
      hm_x3,
      # col = col,
      name = 'Percentage',
      na_col = "grey",

      cluster_columns = F,
      show_column_names = T,
      show_column_dend = F,
      column_names_rot = 45,
      # column_names_max_height = unit(6, "cm"),
      column_names_gp = gpar(fontsize = size/15*8, fontface = "bold"),

      cluster_rows = F,
      show_row_names = F,
      row_names_side = 'left',
      show_row_dend = F,
      row_names_gp = gpar(fontsize = size/15*8, fontface = "bold"),

      show_heatmap_legend = F,

      cell_fun = cell_fun
    ) # ;print(p2)

    # Output
    if(verbose) LuckyVerbose('plotCancerSubtype: All done!')
    return(p2)

  }
)


setGeneric("completeModel", function(object, ...) {
  standardGeneric("completeModel")
})

#' @rdname CCS-method.completeModel
#' @title CCS method: completeModel
#' @description \code{completeModel} method for \code{CCS} class
#' @inheritParams CCSPublicParams
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



