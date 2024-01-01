

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
