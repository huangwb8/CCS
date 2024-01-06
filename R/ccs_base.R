
#' @description CSS Public Parameter Delivering
#' @param object a \code{\link{CCS-class}} object
#' @param size The plot size
#' @param geom Some of \code{'cancer_type'} and \code{'CCS'}.
#' @param hide.legend Some of \code{'cancer_type'} and \code{'CCS'}. Which type of data should be hided in the plot legend.
#' @param method Clustering methods. One of "\code{dbscan}", "\code{ward.D}", "\code{ward.D2}", "\code{single}", "\code{complete}", "\code{average}" (= UPGMA), "\code{mcquitty}" (= WPGMA), "\code{median}" (= WPGMC) or "\code{centroid}" (= UPGMC).
#' @inheritParams ccs
CCSPublicParams <- function(
    object,
    size = 15,
    model.dir,
    geom,
    verbose,
    numCores){
  return(NULL)
}


#' @description softmax function
#' @param x a numeric vector
#' @return a normalized numeric vector
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @details
#' 尽管softmax是一个非线性函数，但softmax回归的输出仍然由输入特征的仿射变换决定。因此，softmax回归是一个线性模型（linear model）。
#' @examples
#' a <- softmax(c(1,2,3,4))
#' sum(a)
softmax <- function(x){
  return(exp(x)/sum(exp(x)))
}

#' @description OneModel-OneData CSS process
#' @param data1 a list containing an expression matrix and its subtype vector
#' @param path_model1 the path of a (GSClassifier) model
#' @importFrom GSClassifier parCallEnsemble
#' @return a data frame with sample IDS and the softmax probability.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
oneCCSProbability <- function(
    data1,
    path_model1,
    geneAnnotation,
    geneSet,
    geneid,
    numCores
){

  # Test
  if(F){
    data1 = data[[2]][[3]]
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
  return(res)

}


#' @description convert an expression data for t-SNE
#' @param res2 The output of \code{\link{oneCCSProbability}}
#' @param rm.dup.col Whether remove duplicated records across the colume (feature) level.
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom digest digest
#' @return a list with raw/cleaned data or md5 sum.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
data_for_tSNE <- function(res2, rm.dup.col = F, verbose = F){

  res2_md5 <- apply(res2, 1, function(x)digest(x, algo="md5"))
  s <- !duplicated(res2_md5)
  if(rm.dup.col){
    res2_md5_col <- apply(res2, 2, function(x)digest(x, algo="md5"))
    s2 <- !duplicated(res2_md5_col)
  } else {
    s2 <- rep(TRUE, ncol(res2))
  }

  if(verbose){
    if(sum(!s)>=1) LuckyVerbose('No. of Duplicati sample = ',sum(!s),'. Remove them!')
    if(sum(!s2)>=1) LuckyVerbose('No. of Duplicati feature = ',sum(!s2),'. Remove them!')
  }
  return(list(
    raw = list(md5=res2_md5, data=res2),
    cleaned = list(md5=res2_md5[s], data=res2[s,s2])
  ))
}


#' @description convert an expression data for t-SNE
#' @param expr Expression matrix based on \code{restSNE[['cleaned']][['data']]} with sample rows and feature cols.
#' @param restSNE Result from \code{data_for_tSNE}
#' @importFrom digest digest
#' @return a list with raw/cleaned data or md5 sum.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
repairCCS <- function(expr, restSNE){

  # Data
  md5_cleaned <- restSNE$cleaned$md5
  data_cleaned <- expr[names(md5_cleaned),]
  data_cleaned <- data_cleaned[names(md5_cleaned),]
  md5_raw <- restSNE$raw$md5
  md5_rest <- md5_raw[!names(md5_raw) %in% names(md5_cleaned)]

  # Index
  if(length(md5_rest)>0){
    index <- match(md5_rest, md5_cleaned)
    df <- NULL
    for(i in 1:length(index)){ # i=1
      df <- rbind(df, data_cleaned[index[i],])
      rownames(df)[i] <- names(md5_rest)[i]
    }
    df2 <- rbind(data_cleaned, df)
    df2 <- df2[names(md5_raw),]
  } else {
    df2 <- data_cleaned
  }


  # Output
  return(df2)

}


#' @description Dimensionality reduction for CCS probability matrix
#' @param SampleIDs Sample IDs like "GSM2411085". One or more IDs are approved.
#' @inheritParams ccs
#' @return Character. The cancer type of sample IDs.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
cancerType <- function(data,SampleIDs=c("GSM2411085","GSM2411084")){
  cancer_type <- lapply(data, function(x) as.character(unlist(lapply(x, function(y) colnames(y$expr)))))
  get_cancer_type <- function(SampleIDs,cancer_type){
    sapply(SampleIDs, function(x){
      for(i in 1:length(cancer_type)){
        if(x %in% cancer_type[[i]]){
          return(names(cancer_type)[i])
        }
      }
    })
  }
  return(get_cancer_type(SampleIDs, cancer_type))
}

#' @title Calinski-Harabasz index
#' @description Calinski-Harabasz index from fpc::cluster.stats
#' @inheritParams fpc::cluster.stats
#' @return Numeric. Calinski-Harabasz index.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @seealso \code{\link[fpc]{cluster.stats}}; \code{\link{CCS-method.optimizeCluster}}; \code{\link{CCS-method.exploreCluster}}
#' @export
CHindex <- function (
    d = NULL,
    clustering,
    noisecluster = TRUE
) {

  # Test
  if(F){

    library(luckyBase); library(CCS)

    # Data
    path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
    object <- readRDS(path_resCCS)
    object <- cluster(object, method = "dbscan", eps = 0.4, minPts = 5)
    # Parameters
    d <- dist(scale(object@Data$Probability$d3), method = "euclidean")
    clustering <- object@Data$CCS; table(clustering)
    noisecluster = TRUE
    verbose = T

    # Other parameters
    alt.clustering = NULL
    nndist = TRUE
    nnk = 2
    standardisation = "max"
    maxk = 10

    # ch = 831.1789
  }

  # Cluster message
  if (!is.null(d)){
    d <- as.dist(d)
  }
  cn <- max(clustering)
  clusteringf <- as.factor(clustering)
  clusteringl <- levels(clusteringf)
  cnn <- length(clusteringl)
  if (cn != cnn) {
    warning("clustering renumbered because maximum != number of clusters")
    for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
    cn <- cnn
  }
  n <- length(clustering)
  noisen <- 0
  cwn <- cn
  if (noisecluster) {
    noisen <- sum(clustering == cn)
    cwn <- cn - 1
  }

  # Data
  dmat <- as.matrix(d)
  nonnoise.ss <- sum(d^2)/n
  if (noisecluster){
    nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn, clustering <= cwn])^2)/sum(clustering <= cwn)
  }

  di <- list(); cluster.size <- within.cluster.ss <- 0
  for (i in 1:cn) {
    # if(verbose) LuckyVerbose(i)
    cluster.size[i] <- sum(clustering == i)
    di <- as.dist(dmat[clustering == i, clustering == i])
    if (i <= cwn) {
      within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
    }
  }
  between.cluster.ss <- nonnoise.ss - within.cluster.ss
  ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss * (cwn - 1))
  return(ch)
}


#' @description Report a parameters group
#' @return Character
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
reportParams <- function(params.i){
  params.i.names <- names(params.i)
  report <- ''
  for(j in 1:length(params.i)){
    report <- paste0(report, params.i.names[j],"=", params.i[j], ifelse(j == length(params.i),"","; "))
  }
  return(report)
}

#' @description List parameters
#' @return List
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
listParams <- function(params.i){
  params.i.list <- list()
  for(j in 1:length(params.i)){
    params.i.list[[names(params.i)[j]]] <- params.i[1, j]
  }
  return(params.i.list)
}






