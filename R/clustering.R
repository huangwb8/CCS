

setGeneric("cluster", function(object, ...) {
  standardGeneric("cluster")
})

#' @rdname CCS-method.cluster
#' @title CCS method: cluster
#' @description \code{cluster} method for \code{CCS} class
#' @param distance Distance method. Not work when \code{method="dbscan"}
#' @param k An integer scalar or vector with the desired number of groups.
#' @param ... Other parameters for \code{\link[dbscan]{dbscan}} or \code{\link[stats]{hclust}}.
#' @inheritParams CCSPublicParams
#' @importFrom dbscan dbscan
#' @seealso \code{\link{ccs}}.
#' @examples
#' resCCS <- cluster(resCCS)
#' @exportMethod cluster
setMethod(
  "cluster",
  signature(object='CCS'),
  function(
    object,
    method = 'dbscan',
    distance = NULL,
    k = NULL,
    ...
  ){

    # Test
    if(F){
      library(dbscan)
      path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
      object <- readRDS(path_resCCS)
    }

    # Data
    d3 <- object@Data$Probability$d3
    dat_cluster <- scale(d3, center = T, scale = T)
    hclust_method <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

    # Clustering
    if(method == 'dbscan'){
      # d3_dbscan <- dbscan(d3, eps = 1, minPts = 5)
      d3_dbscan <- dbscan(dat_cluster, ...)
      d3_CCS <- d3_dbscan$cluster; names(d3_CCS) <- rownames(d3)
    } else if(method %in% hclust_method){
      if(is.null(distance))
        stop('Please set distance. One of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".')
      dis <- dist(dat_cluster, method = distance)
      hc <- hclust(dis, method = method)
      if(is.null(k))
        stop('Please set k,
an integer scalar or vector with the desired number of groups.')
      d3_CCS <- cutree(hc, k=k)
    } else {
      stop('Please set method. One of "dbscan", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".')
    }

    # Results
    object@Data$CCS <- d3_CCS
    if('scaller' %in% names(object@Data)){
      object@Data[['scaller']] <- NULL
    }
    return(object)
  }
)


setGeneric("optimizeCluster", function(object, ...) {
  standardGeneric("optimizeCluster")
})

#' @rdname CCS-method.optimizeCluster
#' @title CCS method: optimizeCluster
#' @description \code{optimizeCluster} method for \code{CCS} class
#' @param method Clustering methods. One of "\code{dbscan}", "\code{ward.D}", "\code{ward.D2}", "\code{single}", "\code{complete}", "\code{average}" (= UPGMA), "\code{mcquitty}" (= WPGMA), "\code{median}" (= WPGMC) or "\code{centroid}" (= UPGMC).
#' @param params a data frame containing parameters for \code{CCS::cluster}.
#' @param cover Boolean. Whether to cocer the existed result of \code{optimizeCluster}. Default if \code{FALSE}.
#' @inheritParams CCSPublicParams
#' @importFrom dbscan dbscan
#' @importFrom dplyr arrange
#' @seealso \code{\link{ccs}}.
#' @examples
#' resCCS <- optimizeCluster(resCCS)
#' @exportMethod optimizeCluster
setMethod(
  "optimizeCluster",
  signature(object='CCS'),
  function(
    object,
    method = 'dbscan',
    params = expand.grid(
      eps = c(seq(0.02,0.1,0.01),seq(0.2,1,0.1)),
      minPts = c(5,6,7,8,9,10,15,20)
    ),
    cover = FALSE,
    verbose = TRUE
  ){

    # Test
    if(F){
      library(luckyBase)
      Plus.library(c('dbscan', 'fpc', 'clue', 'dplyr'))
      path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
      object <- readRDS(path_resCCS)
      method = 'dbscan'
      params = expand.grid(
        eps = c(0.02,0.05,0.1,seq(0.2,1,0.1)),
        minPts = c(5,10,15,20,25,30)
      )
      cover = FALSE
      verbose = TRUE

    }

    if(!cover & !is.null(object@Data[['optimizeCluster']])){
      if(verbose) stop('optimizeCluster: The result exist. Ignored! You can set `cover = TRUE` to make `optimizeCluster` run and cover the old result.', call. = FALSE)
    }

    if(verbose) LuckyVerbose('optimizeCluster: Start...')

    # Distance
    d3_dist <- dist(scale(object@Data$Probability$d3), method = "euclidean")

    # Calinski-Harabasz index
    df <- NULL
    for(i in 1:nrow(params)){ # i=1
      params.i <- params[i,]

      # Verbose management
      if(verbose) LuckyVerbose(reportParams(params.i))

      # Params list
      params.i.list <- listParams(params.i)
      params.i.list[['method']] <- method
      params.i.list[['object']] <- object

      # CCS subtypes based on the params
      object <- do.call(CCS::cluster, params.i.list)

      # Calculate Calinski-Harabasz index
      ch.i <- CCS::CHindex(d = d3_dist,
                           clustering = object@Data$CCS,
                           noisecluster = TRUE)

      # Output
      df <- rbind(
        df,
        cbind(
          CHindex = ch.i,
          params.i
        )
      )
    }
    df <- arrange(df, desc(CHindex))
    if(verbose) LuckyVerbose('The best parameter group is: ', reportParams(df[1,]))
    # mymusic()

    # Renew object
    params.i.list <- listParams(df[1,])[-1] # Remove CHindex
    params.i.list[['method']] <- method
    params.i.list[['object']] <- object
    object <- do.call(CCS::cluster, params.i.list)

    # Output
    object@Data[['optimizeCluster']] <- list(
      method = method,
      CHindex = df
    )
    if(verbose) LuckyVerbose('optimizeCluster: All done!')
    return(object)
  }
)

setGeneric("exploreCluster", function(object, ...) {
  standardGeneric("exploreCluster")
})

#' @rdname CCS-method.exploreCluster
#' @title CCS method: exploreCluster
#' @description \code{exploreCluster} method for \code{CCS} class
#' @inheritParams CCSPublicParams
#' @param nTop Integer. The count of top parameters you want to plot.
#' @seealso \code{\link{ccs}}.
#' @exportMethod exploreCluster
setMethod(
  "exploreCluster",
  signature(object='CCS'),
  function(object,
           nTop = 10,
           geom = c('cancer_type','CCS'),
           hide.legend = c('cancer_type','CCS')[2],
           model.dir = NULL,
           size = 15){

    # Test
    if(F){
      model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225'
      object = readRDS(paste0(model.dir, '/resCCS.rds'))
      nTop = 10
      geom = c('cancer_type','CCS')
      size = 15
    }

    # Model parameters
    if(is.null(model.dir)){
      model.dir = object@Repeat$model.dir
    }
    path_child <- paste0(model.dir, '/clustering')
    dir.create(path_child, showWarnings = F, recursive = T)
    params <- object@Data$optimizeCluster$CHindex[1:nTop,]

    # Plot
    l <- list()
    for(i in 1:nrow(params)){ # i=1

      # Parameters
      params.i <- params[i,]
      params.i.list <- CCS:::listParams(params.i)
      params.i.list[['method']] <- object@Data$optimizeCluster$method
      params.i.list[['object']] <- object
      params.i.list <- params.i.list[-match('CHindex', names(params.i.list))]

      # New CCS subtype
      object.i <- do.call(CCS::cluster, params.i.list)

      # Plot
      name.i <- CCS:::reportParams(params.i)
      cairo_pdf(paste0(path_child, "/DimPlot_",name.i,'.pdf'), width = size/15*12, height = size/15*10)
      l[[name.i]] <- plot(object.i,
                          CCS = object.i@Data$CCS,
                          geom = geom,
                          size = size)
      print(l[[name.i]])
      dev.off()

    }

    # Output
    # saveRDS(l, paste0(path_child, '/Results-of-exploreCluster.rds')) # Very large in the local disk. Ignored.
    return(l)
  }
)
