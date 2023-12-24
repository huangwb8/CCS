#' @description OneModel-OneData CSS process
#' @param data1 a list containing an expression matrix and its subtype vector
#' @param path_model1 the path of a (GSClassifier) model
#' @importFrom GSClassifier parCallEnsemble
#' @return a data frame with sample IDS and the softmax probability.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
oneCCSProbability <- function(data1, path_model1){

  # Test
  if(F){
    data1 = data[[1]][[1]]
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
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom digest digest
#' @return a list with raw/cleaned data or md5 sum.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
data_for_tSNE <- function(res2,verbose){

  res2_md5 <- apply(res2, 1, function(x)digest(x, algo="md5"))
  s <- !duplicated(res2_md5)
  if(verbose){
    LuckyVerbose('No. of Duplicati sample = ',sum(!s),'. Remove them!')
  }
  return(list(
    raw = list(md5=res2_md5, data=res2),
    cleaned = list(md5=res2_md5[s], data=res2[s,])
  ))
}


#' @description Dimensionality reduction for CCS probability matrix
#' @param d2 Data frame. Cleaned CCS probability data after \code{data_for_tSNE}.
#' @param reference Character. Cancer type or other references.
#' @param seed Seeds for \code{\link[Rtsne]{Rtsne}} function.
#' @inheritParams Rtsne::Rtsne
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom Rtsne Rtsne
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
drCCSProbability <- function(
    d2,
    reference,
    dims = 2,
    perplexity = 30,
    theta = 0.3,
    seed = 46,
    verbose = F
){

  # Test
  if(F){
    library(digest); library(GSClassifier); library(Rtsne)
    model.dir <- "./test/ccs/project_01"
    path_resCBP <- paste0(model.dir, '/Cohort-based probability.rds')
    res <- readRDS(path_resCBP)
    res2 <- res[,-1]; res2 <- matrix(as.numeric(as.matrix(res2)), nrow = nrow(res2), byrow = F, dimnames = list(rownames(res2), colnames(res2)))
    d1 <- data_for_tSNE(res2, verbose = T)
    d2 <- d1$cleaned$data
    reference <- Fastextra(colnames(res2), '[|]', 1)
    dims = 2
    perplexity = 30
    theta = 0.3
    verbose = T
    seed = 46
  }

  # Data
  if(is.null(reference)){
    reference <- rep('all', ncol(d2))
  }
  reference_unique <- unique(reference)
  nSubtype <- length(reference_unique)
  set.seed(seed); seeds_dr <- sample(1:10000, nSubtype, replace = F)

  # Dimensionality reduction
  d2_tsne_data <- data.frame()
  for(i in 1:nSubtype){ # i=1
    subtype.i <- reference_unique[i]
    d2_i <- d2[, reference %in% subtype.i]
    set.seed(seeds_dr[i])
    d2_i_tsne <- Rtsne(
      d2_i,
      dims = dims, perplexity = perplexity, theta = theta,
      verbose = verbose
    )
    d2_i_tsne_data <- as.data.frame(d2_i_tsne$Y)
    colnames(d2_i_tsne_data) <- paste(subtype.i, paste('D',1:dims,sep=''), sep='|')
    if(i == 1){
      d2_tsne_data <- d2_i_tsne_data
    } else {
      d2_tsne_data <- cbind(d2_tsne_data, d2_i_tsne_data)
    }
  }

  # Output
  return(d2_tsne_data)
}
