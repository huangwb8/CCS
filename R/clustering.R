

setGeneric("cluster", function(object, ...) {
  standardGeneric("cluster")
})

#' @rdname CCS-method.cluster
#' @title CCS method: cluster
#' @description \code{cluster} method for \code{CCS} class
#' @param object a \code{\link{CCS-class}} object
#' @param method Clustering methods. One of "\code{dbscan}", "\code{ward.D}", "\code{ward.D2}", "\code{single}", "\code{complete}", "\code{average}" (= UPGMA), "\code{mcquitty}" (= WPGMA), "\code{median}" (= WPGMC) or "\code{centroid}" (= UPGMC).
#' @param distance Distance method. Not work when \code{method="dbscan"}
#' @param k An integer scalar or vector with the desired number of groups.
#' @param ... Other parameters for \code{\link[dbscan]{dbscan}} or \code{\link[stats]{hclust}}.
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
    return(object)
  }
)


setGeneric("optimizeCluster", function(object, ...) {
  standardGeneric("optimizeCluster")
})

#' @rdname CCS-method.optimizeCluster
#' @title CCS method: optimizeCluster
#' @description \code{optimizeCluster} method for \code{CCS} class
#' @param object a \code{\link{CCS-class}} object
#' @param ... Other parameters for
#' @importFrom dbscan dbscan
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
      eps = c(0.02,0.05,0.1,seq(0.2,1,0.1)),
      minPts = c(5,10,15,20,25,30)
    ),
    ...
  ){

    # Test
    if(F){
      library(dbscan); library(fpc); library(clues)
      path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
      object <- readRDS(path_resCCS)
      params = expand.grid(
        eps = c(0.02,0.05,0.1,seq(0.2,1,0.1)),
        minPts = c(5,10,15,20,25,30)
      )
    }

    # Calinski-Harabasz index
    for(i in 1:nrow(params)){ # i=1
      params.i <- params[i,]
      object <- cluster(
        object,
        method = "dbscan",
        eps = params.i[,'eps'],
        minPts = params.i[,'minPts']
      )

      # Extract cluster assignments and noise points
      cluster_assignments <- object@Data$CCS

      # Calculate Calinski-Harabasz index
      d3_dist <- dist(scale(object@Data$Probability$d3))
      res <- CHindex(d = d3_dist, clustering = cluster_assignments)
      # Very slow Σ( ° △ °|||)︴
    }
  }
)

setGeneric("exploreCluster", function(object, ...) {
  standardGeneric("exploreCluster")
})

#' @rdname CCS-method.exploreCluster
#' @title CCS method: exploreCluster
#' @description \code{exploreCluster} method for \code{CCS} class
#' @seealso \code{\link{ccs}}.
#' @exportMethod exploreCluster
setMethod(
  "exploreCluster",
  signature(object='CCS'),
  function(object,
           nCluster = 2:20,
           model.dir = NULL,
           size = 15){

    # Test
    if(F){
      model.dir = 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225'
      object = readRDS(paste0(model.dir, '/resCCS.rds'))
      nCluster = 4:5
    }

    # Model parameters
    if(is.null(model.dir)){
      model.dir = object@Repeat$model.dir
    }
    distance  = object@Repeat$params.NbClust$distance
    method = object@Repeat$params.NbClust$method
    path_child <- paste0(model.dir, '/clustering')
    dir.create(path_child, showWarnings = F, recursive = T)

    # Data
    d3 <- object@Data$Probability$d3
    dat_cluster <- scale(d3, center = T, scale = T)
    dis <- dist(dat_cluster, method = distance)
    hc <- hclust(dis, method = method)

    # Plot
    l <- list()
    for(i in nCluster){ # i=15
      y2 <- cutree(hc, i)
      y2 <- y2[rownames(d3)]
      name.i <- paste0("distance=",distance,'_method=',method,'_nCluster=',i)
      cairo_pdf(paste0(path_child, "/DimPlot_",name.i,'.pdf'), width = size/15*12, height = size/15*10)
      l[[name.i]] <- plot(object,
                          CCS = y2,
                          geom = c('cancer_type','CCS'),
                          size = size)
      print(l[[name.i]])
      dev.off()
    }

    # Output
    # saveRDS(l, paste0(path_child, '/Results-of-exploreCluster.rds')) # Very large in the local disk. Ignored.
    return(l)
  }
)
