

setGeneric("dr", function(object, ...) {
  standardGeneric("dr")
})

#' @rdname CCS-method.dr
#' @title CCS method: dr
#' @description \code{dr} method for \code{CCS} class. Dimensionality reduction (DR) for CCS probability matrix.
#' @param object a \code{\link{CCS-class}} object
#' @param method One of \code{'UMAP'} and \code{'t-SNE'}.
#' @param dimension the hyperparameters for 2 level t-SNE.
#' @param ... Parameters of the core DR function.
#' @inheritParams ccs
#' @inheritParams drCCSProbability_UMAP
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom luckyBase LuckyVerbose
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @seealso \code{\link{ccs}}.
#' @exportMethod dr
setMethod(
  "dr",
  signature(object='CCS'),
  function(
    object,
    method = c('UMAP','t-SNE','PCA')[1],
    dimension = c(2,2),
    model.dir = NULL,
    seed = 456,
    numCores = 6,
    verbose = T,
    ...
  ){

  args <- list(...)

  # Data
  data <- object@Data$Probability$d1
  reference <- Fastextra(colnames(data), '[|]', 1)
  if(is.null(model.dir)){
    model.dir <- object@Repeat$model.dir
  }
  path_ccs <- paste0(model.dir, "/resCCS.rds")
  path_dr <- paste0(model.dir, '/drMatrix.rds')

  # Dimensionality reduction (DR)
  if(!file.exists(path_dr)){
    if (method == 'UMAP') {
      if(verbose) LuckyVerbose('dr: use UMAP methods...')
      d2 <- drCCSProbability_UMAP(
        data = data,
        reference = reference,
        dims = dimension[1],
        verbose = verbose,
        ...
      )
      d3 <- drCCSProbability_UMAP(
        data = d2,
        reference = NULL,
        dims = dimension[2],
        verbose = verbose,
        ...
      )
    } else if (method == 't-SNE') {
      if(verbose) LuckyVerbose('dr: use t-SNE methods...')
      d2 <- drCCSProbability_tSNE(
        data = data,
        reference = reference,
        dims = dimension[1],
        seed = seed,
        numCores = numCores,
        verbose = verbose,
        ...
      )
      d3 <- drCCSProbability_tSNE(
        data = d2,
        reference = NULL,
        dims = dimension[2],
        seed = seed,
        numCores = numCores,
        verbose = verbose,
        ...
      )
    } else if (method == 'PCA') {
      if(verbose) LuckyVerbose('dr: use PCA methods...')
      d2 <- drCCSProbability_PCA(
        data = data,
        reference = reference,
        dims = dimension[1],
        ...
      )
      d3 <- drCCSProbability_PCA(
        data = d2,
        reference = NULL,
        dims = dimension[2],
        ...
      )
    } else {
      stop("Wrong method. One of 'UMAP', 't-SNE', and 'PCA'. Please retry!")
    }
    saveRDS(list(d2=d2,d3=d3), path_dr)
  } else {
    if(verbose) LuckyVerbose('The result of dimensionality reduction exists. Use it!')
    dat_dr <- readRDS(path_dr)
    d2 <- dat_dr[['d2']]
    d3 <- dat_dr[['d3']]
  }

  # Output
  object@Data$Probability$d2 <- d2
  object@Data$Probability$d3 <- d3
  if(verbose) LuckyVerbose('dr: All done!')
  return(object)
}
)



#' @description Dimensionality reduction for CCS probability matrix
#' @param data a matrix with feature columns and sample rows.
#' @param reference Character. Cancer type or other references.
#' @param dims The target number of dimensionality reduction. Default is 2.
#' @inheritParams ccs
#' @importFrom dplyr full_join
#' @importFrom luckyBase LuckyVerbose
#' @importFrom umap umap
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
drCCSProbability_UMAP <- function(
    data,
    reference,
    dims = 2,
    verbose = F,
    ...
){

  # Test
  if(F){
    library(digest); library(GSClassifier); library(umap); library(dplyr)
    model.dir <- "./test/ccs/project_01"
    path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
    object <- readRDS(path_resCCS)
    data <- object@Data$Probability$d2
    reference <- Fastextra(colnames(data), '[|]', 1)
    dims = 2
    verbose = T
  }

  # Data
  d2 <- data
  if(is.null(reference)){
    reference <- rep('all', ncol(d2))
  }
  reference_unique <- unique(reference)
  nCancerType <- length(reference_unique)

  # Dimensionality reduction
  d2_dr_data <- data.frame()
  for(i in 1:nCancerType){ # i=1
    subtype.i <- reference_unique[i]
    d2_i <- d2[, reference %in% subtype.i]
    if(verbose) LuckyVerbose("drCCSProbability: ",subtype.i)
    d2_dr <- data_for_dr(d2_i, rm.dup.col = F, verbose = verbose)
    d2_dr_cleanedData <- d2_dr[['cleaned']][['data']]
    # d2_i_dr <- umap(d2_dr_cleanedData, n_components = dims)
    d2_i_dr <- umap(
      d2_dr_cleanedData,
      n_components = dims,
      ...)
    d2_i_dr_data <- as.data.frame(d2_i_dr$layout)
    rownames(d2_i_dr_data) <- names(d2_dr[['cleaned']][['md5']])
    d2_i_dr_data <- repairCCS(expr = d2_i_dr_data, restSNE = d2_dr)
    colnames(d2_i_dr_data) <- paste(subtype.i, paste('D',1:dims,sep=''), sep='|')
    d2_i_dr_data <- cbind(SampleIDs = rownames(d2_i_dr_data), d2_i_dr_data)

    if(i == 1){
      d2_dr_data <- d2_i_dr_data
    } else {
      d2_dr_data <- full_join(d2_dr_data, d2_i_dr_data, by = "SampleIDs")
    }
  }

  # Output
  rownames(d2_dr_data) <- as.character(d2_dr_data$SampleIDs)
  if(verbose) LuckyVerbose("drCCSProbability: Done!")
  return(as.matrix(d2_dr_data[,-1]))
}


#' @description Dimensionality reduction for CCS probability matrix
#' @param data a matrix with feature columns and sample rows.
#' @param reference Character. Cancer type or other references.
#' @param dims The target number of dimensionality reduction. Default is 2.
#' @inheritParams ccs
#' @importFrom dplyr full_join
#' @importFrom luckyBase LuckyVerbose
#' @importFrom stats prcomp
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
drCCSProbability_PCA <- function(
    data,
    reference,
    dims = 2,
    verbose = F,
    ...
){

  # Test
  if(F){
    library(digest); library(GSClassifier); library(dplyr); library(stats)
    model.dir <- "./test/ccs/project_01"
    path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
    object <- readRDS(path_resCCS)
    data <- object@Data$Probability$d1
    reference <- Fastextra(colnames(data), '[|]', 1)
    dims = 2
    verbose = T
  }

  # Data
  d2 <- data
  if(is.null(reference)){
    reference <- rep('all', ncol(d2))
  }
  reference_unique <- unique(reference)
  nCancerType <- length(reference_unique)

  # Dimensionality reduction
  d2_dr_data <- data.frame()
  for(i in 1:nCancerType){ # i=1
    subtype.i <- reference_unique[i]
    d2_i <- d2[, reference %in% subtype.i]
    if(verbose) LuckyVerbose("drCCSProbability: ",subtype.i)
    d2_dr <- CCS:::data_for_dr(d2_i, rm.dup.col = F, verbose = verbose)
    d2_dr_cleanedData <- d2_dr[['cleaned']][['data']]
    d2_i_dr <- prcomp(
      d2_dr_cleanedData,
      center = TRUE, scale. = TRUE,
      rank. = dims,
      ...
    )
    d2_i_dr_data <- as.data.frame(d2_i_dr$x)
    rownames(d2_i_dr_data) <- names(d2_dr[['cleaned']][['md5']])
    d2_i_dr_data <- CCS:::repairCCS(expr = d2_i_dr_data, restSNE = d2_dr)
    colnames(d2_i_dr_data) <- paste(subtype.i, paste('D',1:dims,sep=''), sep='|')
    d2_i_dr_data <- cbind(SampleIDs = rownames(d2_i_dr_data), d2_i_dr_data)

    if(i == 1){
      d2_dr_data <- d2_i_dr_data
    } else {
      d2_dr_data <- full_join(d2_dr_data, d2_i_dr_data, by = "SampleIDs")
    }
  }

  # Output
  rownames(d2_dr_data) <- as.character(d2_dr_data$SampleIDs)
  if(verbose) LuckyVerbose("drCCSProbability: Done!")
  return(as.matrix(d2_dr_data[,-1]))
}




#' @description Dimensionality reduction for CCS probability matrix
#' @param seed a seed for random process management. Not needed in \code{drCCSProbability_UMAP}.
#' @inheritParams drCCSProbability_UMAP
#' @inheritParams ccs
#' @importFrom dplyr full_join
#' @importFrom luckyBase LuckyVerbose
#' @importFrom Rtsne Rtsne
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
drCCSProbability_tSNE <- function(
    data,
    reference,
    dims = 2,
    seed = 46,
    numCores = 6,
    verbose = F,
    ...
){

  # Test
  if(F){
    library(digest); library(GSClassifier); library(Rtsne)
    model.dir <- "./test/ccs/project_01"
    path_resCBP <- paste0(model.dir, '/Cohort-based probability.rds')
    res <- readRDS(path_resCBP)
    res2 <- res[,-1]; res2 <- matrix(as.numeric(as.matrix(res2)), nrow = nrow(res2), byrow = F, dimnames = list(rownames(res2), colnames(res2)))
    d1 <- data_for_dr(res2, verbose = T)
    d2 <- d1$cleaned$data
    reference <- Fastextra(colnames(res2), '[|]', 1)
    dims = 2
    perplexity = 30
    theta = 0.3
    verbose = T
    seed = 46
  }

  # Data
  d2 <- data
  if(is.null(reference)){
    reference <- rep('all', ncol(d2))
  }
  reference_unique <- unique(reference)
  nCancerType <- length(reference_unique)
  set.seed(seed); seeds_dr <- sample(1:10000, nCancerType, replace = F)

  # Dimensionality reduction
  d2_tsne_data <- data.frame()
  for(i in 1:nCancerType){ # i=1
    subtype.i <- reference_unique[i]
    d2_i <- d2[, reference %in% subtype.i]
    if(verbose) LuckyVerbose("drCCSProbability: ",subtype.i)
    d2_tSNE <- data_for_dr(d2_i, rm.dup.col = F, verbose = verbose)
    d2_tSNE_cleanedData <- d2_tSNE[['cleaned']][['data']]
    set.seed(seeds_dr[i])
    d2_i_tsne <- Rtsne(
      d2_tSNE_cleanedData,
      dims = dims,
      verbose = verbose,
      num_threads = numCores,
      ...
    )
    d2_i_tsne_data <- as.data.frame(d2_i_tsne$Y)
    rownames(d2_i_tsne_data) <- names(d2_tSNE[['cleaned']][['md5']])
    d2_i_tsne_data <- repairCCS(expr = d2_i_tsne_data, restSNE = d2_tSNE)
    colnames(d2_i_tsne_data) <- paste(subtype.i, paste('D',1:dims,sep=''), sep='|')
    d2_i_tsne_data <- cbind(SampleIDs = rownames(d2_i_tsne_data), d2_i_tsne_data)

    if(i == 1){
      d2_tsne_data <- d2_i_tsne_data
    } else {
      d2_tsne_data <- full_join(d2_tsne_data, d2_i_tsne_data, by = "SampleIDs")
    }
  }

  # Output
  rownames(d2_tsne_data) <- as.character(d2_tsne_data$SampleIDs)
  if(verbose) LuckyVerbose("drCCSProbability: Done!")
  return(as.matrix(d2_tsne_data[,-1]))
}

