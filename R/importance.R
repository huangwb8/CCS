

setGeneric("importance", function(object, ...) {
  standardGeneric("importance")
})

#' @rdname CCS-method.importance
#' @title CCS method: importance
#' @description \code{importance} method for \code{CCS} class
#' @param ... parameters for \code{\link[xgboost]{xgb.cv}} and \code{\link[xgboost]{xgboost}}.
#' @param seed random seed for \code{xgboost}
#' @inheritParams CCSPublicParams
#' @inheritParams ccs
#' @inheritParams dr
#' @importFrom plyr llply
#' @seealso \code{\link{ccs}};\code{\link{plotImportance}}.
#' @examples
#' To be continued!
#' @exportMethod importance
setMethod(
  "importance",
  signature(object='CCS'),
  function(
    object,
    seed = 2024,
    verbose = TRUE,
    numCores = 16,
    ...
  ){

    # Test
    if(F){
      library(luckyBase)
      np <- c('xgboost','tidyr','plyr','dplyr','purrr','furrr'); Plus.library(np)
      path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240526/resCCS.rds'
      object <- readRDS(path_resCCS)
      object@Repeat$method <- 'GSClassifier'
      object@Repeat$model.dir <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240526'
      numCores = 32
      seed = 2024
      numCores = 16
      verbose = TRUE
    }

    # Environment
    method <- object@Repeat[["method"]]
    if(!method %in% c('GSClassifier')){
      stop('importance: Error method Σ(°△ °|||)︴ Please choose one of "GSClassifier"!')
    }

    # Self-defined parameters
    if(T){
      if(verbose) LuckyVerbose('importance: adjust parameters...')
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
      params_xg <- params[-match(c('n','sampSize','ptail', 'nround.mode'), names(params))]
      params_xg[['objective']] <- "reg:squarederror"
      params_xg[['nthread']] <- numCores # numCores
      # params_xg[['device']] <- 'cpu'
      # params_xg[['eta']] <- 0.2
    }

    # Cohort & tissue importance
    if(verbose) LuckyVerbose('importance: Calculate cohort & tissue importance...')
    t1 <- system.time(
      l <- ctImportance(object, params_xg, seed, verbose)
    )
    if(verbose) LuckyVerbose('importance: ',round(sum(t1, na.rm = T), 2),'s consuming.')

    # Feature importance
    if(method == 'GSClassifier'){
      if(verbose) LuckyVerbose('importance: Calculate feature importance...')
      t2 <- system.time(
        l[['importance_feature']] <- featureImportance_GSClassifier(object, l, numCores)
      )
      if(verbose) LuckyVerbose('importance: ',round(sum(t2, na.rm = T), 2),'s consuming.')
    } else {
      stop('Wrong model method. Please set one of "GSClassifier"!')
    }

    # Get merged importance
    l.merge <- llply(l, mergeImportance)

    # Output
    object@Data[['importance']] <- l.merge
    return(object)
  }
)



setGeneric("plotImportance", function(object, ...) {
  standardGeneric("plotImportance")
})

#' @rdname CCS-method.plotImportance
#' @title CCS method: plotImportance
#' @description \code{plotImportance} method for \code{CCS} class
#' @param type One of \code{"tissue"}, \code{"cohort"}, or \code{"feature"}.
#' @param convertGene Whether to convert names like \code{"ENSG00000185811:ENSG00000187189"} to names like \code{"IKZF1:TSPYL4"}.
#' @param ... parameters to \code{\link[ggplot2]{geom_bar}}.
#' @inheritParams CCSPublicParams
#' @inheritParams ccs
#' @inheritParams dr
#' @importFrom plyr llply ldply
#' @importFrom dplyr slice_max
#' @importFrom tidyr `%>%`
#' @importFrom ggpubr rotate_x_text
#' @import ggplot2
#' @seealso \code{\link{ccs}}; \code{\link{importance}}.
#' @examples
#' To be continued!
#' @exportMethod plotImportance
setMethod(
  "plotImportance",
  signature(object='CCS'),
  function(
    object,
    type = c('tissue','cohort','feature')[1],
    convertGene = F,
    nTop = 10,
    size = 10,
    ...
  ){

    # Test
    if(F){
      library(luckyBase)
      np <- c('plyr','dplyr','ggplot2', 'ggpubr'); Plus.library(np)
      path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20240526/resCCS.rds'
      object <- readRDS(path_resCCS)
      type = c('tissue','cohort','feature')[3]
      nTop = 10
      size = 10
      convertGene = T
    }

    # Data
    df <- object@Data[['importance']][[paste0('importance_',type)]] %>%
      llply(., function(x) {
      if(nTop <= 1){
        x <- x %>% slice_max(Gain, prop = nTop)
      } else {
        x <- x %>% slice_max(Gain, n = nTop)
      }

      # Gene Symbol
      if(type %in% 'feature' & convertGene){
        x$Feature <- sapply(x$Feature, function(x){
          if(grepl('^ENSG',x)){
            # Name like "ENSG00000185811:ENSG00000187189"
            return(Fastextra(x, '[:]') %>% convert(., fromtype = 'ENSEMBL',totype = 'SYMBOL') %>% paste0(.,collapse = ':'))
          } else {
            # Keep it as raw
            return(x)
          }
        }, USE.NAMES = F)
      }

      # Re-factor
      x$Feature <- factor(x$Feature, level = rev(as.character(x$Feature)))
      return(x)
    }) %>%
      ldply(.id = "Target")

    # Plot
    y.name <- ifelse(type=='tissue','Tissue', ifelse(type=='cohort','Cohort','Feature'))
    p <- ggplot(df, aes(x = Gain, y = Feature)) +
      # geom_bar(stat = "identity", color = "black") +
      geom_bar(stat = "identity", color = "black", ...) +
      theme_bw() +
      labs(title = NULL,
           x = "Importance",
           y = y.name) +
      facet_grid(. ~ Target) +
      theme(
        axis.text = element_text(size = size, colour = "black",face = "bold"),
        axis.title = element_text(size = size*1.1, colour = "black",face = "bold"),
        strip.text = element_text(size = size*1.1, colour = "black",face = "bold"),
        strip.background = element_rect(fill="white")
      ) +
      rotate_x_text(angle = 45)

    # Output result
    return(p)
  }
)


####=================Assistant functions===============####

#' @description Cohort & Tissue Importance estimation
#' @importFrom xgboost xgb.DMatrix xgb.cv xgboost xgb.importance
#' @importFrom luckyBase Fastextra
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
ctImportance <- function(object, params_xg, seed, verbose){

  d1 <- object@Data[["Probability"]][["d1"]]
  d3 <- object@Data[["Probability"]][["d3"]]

  l <- list()
  for(i in 1:ncol(d3)){ # i=1

    featureName.i <- colnames(d3)[i]

    # Data
    x <- d1[rownames(d3),]
    y <- as.numeric(d3[,i])

    # xgboost
    dtrain <- xgb.DMatrix(data = x, label = y)

    set.seed(seed)
    cvRes <- xgb.cv(
      params = {
        p <- params_xg;
        p[c('nfold','nrounds')] <- NULL;
        p
      },
      data = dtrain,
      nrounds=params_xg$nrounds,
      nfold=params_xg$nfold,
      early_stopping_rounds=10,
      verbose = 0
    )

    set.seed(seed)
    bst_model <- xgboost(
      params = {
        p <- params_xg;
        p[c('nfold','nrounds')] <- NULL;
        p['device'] <- 'gpu';
        p
      },
      data = dtrain,
      nrounds = 105,
      verbose = 0
    )

    # 提取变量重要性
    importance_matrix <- xgb.importance(
      feature_names = colnames(x),
      model = bst_model
    )
    x2 <- cbind(N=1, as.matrix(importance_matrix[,-1]))

    # Cohort Importance
    importance_matrix_cohort <-
      mergeMatrixDup(
        x2,
        mergeCol = F,
        mergeRow = T,
        fun_row = function(x)sum(x, na.rm = T),
        refRow = Fastextra(importance_matrix$Feature, '[|]', 2),
        verbose = F
      )

    # Tissue Importance
    importance_matrix_tissue <-
      mergeMatrixDup(
        x2,
        mergeCol = F,
        mergeRow = T,
        fun_row = function(x)sum(x, na.rm = T),
        refRow = Fastextra(importance_matrix$Feature, '[|]', 1),
        verbose = F
      )

    # Output
    l[['importance_cohort']][[featureName.i]] <- cbind(Feature=rownames(importance_matrix_cohort), as.data.frame(importance_matrix_cohort)) %>% arrange(desc(Gain))
    l[['importance_tissue']][[featureName.i]] <- cbind(Feature = rownames(importance_matrix_tissue), as.data.frame(importance_matrix_tissue)) %>% arrange(desc(Gain))
  }

  return(l)

}


#' @description Feature(genes) Importance estimation
#' @importFrom xgboost xgb.DMatrix xgb.cv xgboost xgb.importance
#' @importFrom luckyBase Fastextra
#' @importFrom purrr map_dfr imap_dfr
#' @importFrom dplyr full_join arrange summarize
#' @importFrom plyr ddply
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
featureImportance_GSClassifier <- function(object, resctImportance, numCores){

  # Test
  if(F){
    resctImportance <- l
  }

  path_submodels <- list.files(
    object@Repeat[["model.dir"]],
    pattern = 'modelFit.rds',
    full.names = T,
    recursive = T
  )

  if(is.na(object@Model[[1]])){

    process_model <- function(model, subtype, tissue, dataset) {
      if (is.null(model)) return(NULL)

      importance_matrix <- xgb.importance(
        feature_names = model$genes,
        model = model$bst
      )

      if (nrow(importance_matrix) == 0) return(NULL)

      cbind(
        tissue = tissue,
        dataset = dataset,
        subtype = subtype,
        importance_matrix
      )
    }

    df <- path_submodels %>%
      map_dfr(~{
        model_i <- readRDS(.x)[['Model']]
        tissue <- rev(Fastextra(.x, '[/]'))[3]
        dataset <- rev(Fastextra(.x, '[/]'))[2]

        model_i %>%
          imap_dfr(~{
            imap_dfr(.x, ~process_model(.x, .y, tissue, dataset))
          })
      })

    # df <- readRDS('./test/Demo-ImportanceResults.rds')[["importance_feature"]]

    # Summary
    df2 <- ddply(
      df,
      c('dataset','Feature'),
      dplyr::summarize,
      Gain = median(Gain,na.rm = T),
      Cover = median(Cover,na.rm = T),
      Frequency = median(Frequency, na.rm = T)
    )

    # Get median feature importance
    importance_cohort <- resctImportance[['importance_cohort']]
    importance_feature <- list()
    for(i in 1:length(importance_cohort)){ # i=1
      importance_cohort_i <- importance_cohort[[i]]
      importance_cohort_i <- cbind(dataset = rownames(importance_cohort_i), as.data.frame(importance_cohort_i[, c('Gain','Cover','Frequency')]))
      df.i <- full_join(df2, importance_cohort_i, by='dataset')
      df.i.2 <- df.i[,c("Gain.x","Cover.x","Frequency.x")] * df.i[,c("Gain.y","Cover.y","Frequency.y")]
      colnames(df.i.2) <- c(c("Gain","Cover","Frequency"))
      df.i.2 <- cbind(df.i[,c("dataset","Feature")],df.i.2)
      df.i.3 <- ddply(
        df.i.2,
        "Feature",
        dplyr::summarize,
        Gain = median(Gain, na.rm = T),
        Cover = median(Cover, na.rm = T),
        Frequency = median(Frequency, na.rm = T)
      )
      for(j in 2:ncol(df.i.3)){
        df.i.3[,j] <- df.i.3[,j]/sum(df.i.3[,j], na.rm = T)
      }
      df.i.3 <- arrange(df.i.3, desc(Gain))
      importance_feature[[names(importance_cohort)[i]]] <- df.i.3
    }
  }

  # Output the result
  return(importance_feature)

  # Trash
  if(F){

    # Legacy version
    df <- data.frame()
    for(i in 1:length(path_submodels)){ # i=1
      model_i <- readRDS(path_submodels[i])[['Model']]
      for(j in 1:length(model_i)){ # j = 1
        model_j <- model_i[[j]]
        for(k in 1:length(model_j)){ # k=1
          model_k <- model_j[[k]]
          # LuckyVerbose('i=',i,'; j=',j,'; k=',k)
          if(!is.null(model_k)){
            importance_matrix <- xgb.importance(
              feature_names = model_k$genes,
              model = model_k$bst
            );
            if(nrow(importance_matrix) > 0){
              df.k <- cbind(
                tissue = rev(Fastextra(path_submodels[i],'[/]'))[3],
                dataset = rev(Fastextra(path_submodels[i],'[/]'))[2],
                subtype = names(model_j)[k],
                importance_matrix
              )
              df <- rbind(df, df.k)
            }
          }
        }
      }
    }


    # parallel version: 看似没问题，但实际运行时有报错： list contains fewer than 2 elements。 可能与.options = furrr_options(seed = TRUE)有关。
    if(T){
      plan(multisession, workers = numCores)

      process_model <- function(model, subtype, tissue, dataset) {
        if (is.null(model)) return(NULL)

        importance_matrix <- xgb.importance(
          feature_names = model$genes,
          model = model$bst
        )

        if (nrow(importance_matrix) == 0) return(NULL)

        cbind(
          tissue = tissue,
          dataset = dataset,
          subtype = subtype,
          importance_matrix
        )
      }

      df <- future_map_dfr(path_submodels, ~{
        model_i <- readRDS(.x)[['Model']]
        tissue <- rev(Fastextra(.x, '[/]'))[3]
        dataset <- rev(Fastextra(.x, '[/]'))[2]

        future_imap_dfr(model_i, ~{
          future_imap_dfr(.x, ~process_model(.x, .y, tissue, dataset))
        }, .options = furrr_options(seed = TRUE))
      }, .options = furrr_options(seed = TRUE)
      )

      # End parallel
      plan(sequential)
    }
  }

}


#' @description Merge all|D1 and all|D2
#' @importFrom dplyr full_join arrange
#' @importFrom luckyBase mergeMatrixDup Fastextra
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
mergeImportance <- function(lst){


  # Test
  if(F){
    lst <- l[[1]]
  }

  idCol <- c('Feature')
  targetCol <- c("Gain","Cover","Frequency")
  for(i in 1:length(lst)){ # i=1
    if(i == 1){
      m <- lst[[i]][c(idCol,targetCol)]
    } else {
      m.i <- lst[[i]][c(idCol,targetCol)]
      m.i2 <- full_join(m, m.i,by = 'Feature')
      m.i3 <- m.i2[-match(idCol, colnames(m.i2))]
      m.i3 <- mergeMatrixDup(
        x = m.i3,
        mergeCol = T,
        fun_col = function(x)sum(x, na.rm = T),
        refCol = Fastextra(colnames(m.i3),'[.]',1),
        mergeRow = F,
        verbose = F
      )
      m <- cbind(m.i2[,idCol],as.data.frame(m.i3)); colnames(m)[1:length(idCol)] <- idCol
    }
  }
  m[,targetCol] <- m[,targetCol] / length(lst)
  m <- arrange(m,desc(Gain)) %>% as.data.frame()
  lst <- c(merge=list(m[c(idCol,targetCol)]),lst)
  return(lst)

}


#### End ####

