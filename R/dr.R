

setGeneric("dr", function(object, ...) {
  standardGeneric("dr")
})

#' @rdname CCS-method.dr
#' @title CCS method: dr
#' @description \code{dr} method for \code{CCS} class. Dimensionality reduction (DR) for CCS probability matrix.
#' @param object a \code{\link{CCS-class}} object
#' @param dimension the hyperparameters for 2 level t-SNE.
#' @param cover Whether to create an exact new dr result
#' @param ... Parameters of DR functions.
#' @inheritParams ccs
#' @inheritParams drCCSProbability
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom luckyBase LuckyVerbose
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @details
#' ## Parameters
#'
#' ### Method == 't-SNE'
#'
#' - \code{seed}: random seed
#' - \code{numCores}: CPU threads
#'
#' @seealso \code{\link{ccs}}.
#' @exportMethod dr
setMethod(
  "dr",
  signature(object='CCS'),
  function(
    object,
    method = c('UWOT','UMAP','t-SNE','PCA')[1],
    dimension = c(2,2),
    cover = F,
    model.dir = NULL,
    seed = 2024,
    verbose = T,
    ...
  ){

  # Data
  data <- object@Data$Probability$d1
  reference <- Fastextra(colnames(data), '[|]', 1)
  if(is.null(model.dir)){
    model.dir <- object@Repeat$model.dir
  }
  path_ccs <- paste0(model.dir, "/resCCS.rds")
  path_dr <- paste0(model.dir, '/drMatrix.rds')

  # Dimensionality reduction (DR)
  if(!file.exists(path_dr) | cover){
    d2 <- drCCSProbability(
      data = data,
      method = method,
      reference = reference,
      dims = dimension[1],
      seed = seed,
      verbose = verbose,
      ...)
    d3 <- drCCSProbability(
      data = d2,
      method = method,
      reference = NULL,
      dims = dimension[2],
      seed = seed,
      verbose = verbose,
      ...)
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
#' @param seed random seed
#' @inheritParams ccs
#' @inheritParams CORE_DR
#' @importFrom dplyr full_join
#' @importFrom luckyBase LuckyVerbose
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
drCCSProbability <- function(
    data,
    method = c('UWOT','UMAP','t-SNE','PCA')[1],
    reference,
    dims = 2,
    verbose = F,
    seed = 2024,
    ...
){

  # Test
  if(F){
    library(digest); library(GSClassifier); library(umap); library(dplyr)
    # model.dir <- "./test/ccs/project_01"
    # path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test01/ccs/v20231225/resCCS.rds'
    model.dir <- "E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810"
    path_resCCS <- 'E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810/resCCS.rds'
    method = c('UWOT','UMAP','t-SNE','PCA')[1]
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

  # Seeds management
  set.seed(seed); seeds_dr <- sample(1:10000, nCancerType, replace = F)

  # Dimensionality reduction
  d2_dr_data <- data.frame()
  for(i in 1:nCancerType){ # i=1

    # Data Preparation
    subtype.i <- reference_unique[i]
    d2_i <- d2[, reference %in% subtype.i]
    if(verbose) LuckyVerbose("drCCSProbability - Method_",method,": ",subtype.i)
    d2_dr <- CCS:::data_for_dr(d2_i, rm.dup.col = F, verbose = verbose)
    d2_dr_cleanedData <- d2_dr[['cleaned']][['data']]

    # Core function for DR
    d2_i_dr_data <- CCS:::CORE_DR(
      method = method,
      data = d2_dr_cleanedData,
      dims = dims,
      seed = seeds_dr[i],
      ...
    )

    # Repair data
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
  if(verbose) LuckyVerbose("drCCSProbability - Method_",method,": Done!")
  return(as.matrix(d2_dr_data[,-1]))
}



#' @description Core Function for Dimensionality reduction in CCS
#' @param method One of \code{'UWOT'}, \code{'UMAP'},\code{'t-SNE'} and \code{'PCA'}
#' @param data cleaned data from the result of \code{CCS:::data_for_dr}
#' @param dims The target number of dimensionality reduction.
#' @importFrom Rtsne Rtsne
#' @importFrom stats prcomp
#' @importFrom umap umap.knn
#' @importFrom uwot load_uwot
#' @return data frame
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
CORE_DR <- function(method, data, dims, seed = seeds_dr[i],...){

  if(method == 'UMAP'){
    # [Frequently Asked Questions — umap 0.5 documentation](https://umap-learn.readthedocs.io/en/latest/faq.html)
    set.seed(seed)
    d2_i_dr <- umap::umap(
      data,
      n_components = dims,
      ...)
    d2_i_dr_data <- as.data.frame(d2_i_dr$layout)

  } else if(method == 'UWOT'){
    # [Frequently Asked Questions — umap 0.5 documentation](https://umap-learn.readthedocs.io/en/latest/faq.html)
    d2_i_dr <- uwot::umap(
      data,
      n_components = dims,
      seed = seed,
      ...)
    d2_i_dr_data <- as.data.frame(d2_i_dr)

  } else if(method == 'PCA'){
    set.seed(seed)
    d2_i_dr <- prcomp(
      data,
      center = TRUE, scale. = TRUE,
      rank. = dims,
      ...
    )
    d2_i_dr_data <- as.data.frame(d2_i_dr$x)

  } else if(method == 't-SNE'){
    set.seed(seed)
    d2_i_dr <- Rtsne(
      data,
      dims = dims,
      verbose = verbose,
      ...
    )
    d2_i_dr_data <- as.data.frame(d2_i_dr$Y)
  } else {
    stop("Wrong method. One of 'UMAP', 't-SNE', and 'PCA'. Please retry!")
  }

  # Output
  return(d2_i_dr_data)

}




