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


#' @description Calinski-Harabasz index from fpc::cluster.stats
CHindex <- function (d = NULL, clustering, alt.clustering = NULL, noisecluster = TRUE, nndist = TRUE, nnk = 2, standardisation = "max", maxk = 10, cvstan = sqrt(length(clustering)))
{
  lweight <- function(x, md) (x < md) * (-x/md + 1)
  if (!is.null(d))
    d <- as.dist(d)
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
  parsimony <- cn/maxk
  diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
  for (i in 1:cn) cluster.size[i] <- sum(clustering == i)
  if (is.numeric(standardisation))
    stan <- standardisation
  else stan <- switch(standardisation, max = max(d), ave = mean(d),
                      q90 = quantile(d, 0.9), 1)
  pk1 <- cluster.size/n
  pk10 <- pk1[pk1 > 0]
  h1 <- -sum(pk10 * log(pk10))
  if (!(standardisation == "none")) {
    pkmax <- rep(1, cn)/cn
    h1 <- -h1/sum(pkmax * log(pkmax))
  }
  corrected.rand <- vi <- NULL
  if (!is.null(alt.clustering)) {
    choose2 <- function(v) {
      out <- numeric(0)
      for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2,
                                              choose(v[i], 2), 0)
      out
    }
    cn2 <- max(alt.clustering)
    clusteringf <- as.factor(alt.clustering)
    clusteringl <- levels(clusteringf)
    cnn2 <- length(clusteringl)
    if (cn2 != cnn2) {
      warning("alt.clustering renumbered because maximum != number of clusters")
      for (i in 1:cnn2) alt.clustering[clusteringf == clusteringl[i]] <- i
      cn2 <- cnn2
    }
    nij <- table(clustering, alt.clustering)
    dsum <- sum(choose2(nij))
    cs2 <- numeric(0)
    for (i in 1:cn2) cs2[i] <- sum(alt.clustering == i)
    sum1 <- sum(choose2(cluster.size))
    sum2 <- sum(choose2(cs2))
    pk2 <- cs2/n
    pk12 <- nij/n
    corrected.rand <- (dsum - sum1 * sum2/choose2(n))/((sum1 +
                                                          sum2)/2 - sum1 * sum2/choose2(n))
    pk20 <- pk2[pk2 > 0]
    h2 <- -sum(pk20 * log(pk20))
    icc <- 0
    for (i in 1:cn) for (j in 1:cn2) if (pk12[i, j] > 0)
      icc <- icc + pk12[i, j] * log(pk12[i, j]/(pk1[i] *
                                                  pk2[j]))
    vi <- h1 + h2 - 2 * icc
  }

  dmat <- as.matrix(d)
  within.cluster.ss <- 0
  overall.ss <- nonnoise.ss <- sum(d^2)/n
  if (noisecluster)
    nonnoise.ss <- sum(as.dist(dmat[clustering <= cwn,
                                    clustering <= cwn])^2)/sum(clustering <= cwn)
  ave.between.matrix <- separation.matrix <- matrix(0,
                                                    ncol = cn, nrow = cn)
  nnd <- numeric(0)
  cvnndc <- rep(NA, cn)
  mnnd <- cvnnd <- NULL
  di <- list()
  for (i in 1:cn) {
    cluster.size[i] <- sum(clustering == i)
    di <- as.dist(dmat[clustering == i, clustering ==
                         i])
    if (i <= cwn) {
      within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
      within.dist <- c(within.dist, di)
    }
    if (length(di) > 0) {
      diameter[i] <- max(di)
      average.distance[i] <- mean(di)
      median.distance[i] <- median(di)
    }
    else diameter[i] <- average.distance[i] <- median.distance[i] <- NA
    bv <- numeric(0)
    for (j in 1:cn) {
      if (j != i) {
        sij <- dmat[clustering == i, clustering ==
                      j]
        bv <- c(bv, sij)
        if (i < j) {
          separation.matrix[i, j] <- separation.matrix[j,
                                                       i] <- min(sij)
          ave.between.matrix[i, j] <- ave.between.matrix[j,
                                                         i] <- mean(sij)
          if (i <= cwn & j <= cwn)
            between.dist <- c(between.dist, sij)
        }
      }
    }
    separation[i] <- min(bv)
    average.toother[i] <- mean(bv)
  }
  if (nndist) {
    kenough <- cluster.size > nnk
    for (i in (1:cn)[kenough]) {
      nndi <- apply(dmat[clustering == i, clustering ==
                           i], 1, sort, partial = nnk + 1)[nnk + 1, ]
      nnd <- c(nnd, nndi)
      cvnndc[i] <- sd(nndi)/mean(nndi)
    }
    cvnnd <- weighted.mean(cvnndc, pk1, na.rm = TRUE)
    mnnd <- mean(nnd, na.rm = TRUE)
    if (!standardisation == "none") {
      maxnnd <- max(apply(dmat, 1, sort, partial = nnk +
                            1)[nnk + 1, ])
      mnnd <- mnnd/maxnnd
      cvnnd <- cvnnd/cvstan
    }
  }
  average.between <- mean(between.dist)/stan
  average.within <- weighted.mean(average.distance, cluster.size,
                                  na.rm = TRUE)/stan
  nwithin <- length(within.dist)
  nbetween <- length(between.dist)
  between.cluster.ss <- nonnoise.ss - within.cluster.ss
  ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss *
                                                   (cwn - 1))
  return(ch)
}




