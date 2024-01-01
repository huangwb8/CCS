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
#' @param d2 Data frame. Cleaned CCS probability data after \code{data_for_tSNE}.
#' @param reference Character. Cancer type or other references.
#' @param seed Seeds for \code{\link[Rtsne]{Rtsne}} function.
#' @inheritParams Rtsne::Rtsne
#' @inheritParams GSClassifier::parCallEnsemble
#' @importFrom dplyr full_join
#' @importFrom Rtsne Rtsne
#' @importFrom luckyBase LuckyVerbose
#' @return data frame. Series of PC1, PC2 and etc after dimensionality reduction
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
drCCSProbability <- function(
    d2,
    reference,
    dims = 2,
    perplexity = 30,
    theta = 0.3,
    seed = 46,
    numCores = 6,
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
  nCancerType <- length(reference_unique)
  set.seed(seed); seeds_dr <- sample(1:10000, nCancerType, replace = F)

  # Dimensionality reduction
  d2_tsne_data <- data.frame()
  for(i in 1:nCancerType){ # i=1
    subtype.i <- reference_unique[i]
    d2_i <- d2[, reference %in% subtype.i]
    if(verbose) LuckyVerbose("drCCSProbability: ",subtype.i)
    d2_tSNE <- data_for_tSNE(d2_i, rm.dup.col = F, verbose = verbose)
    d2_tSNE_cleanedData <- d2_tSNE[['cleaned']][['data']]
    set.seed(seeds_dr[i])
    d2_i_tsne <- Rtsne(
      d2_tSNE_cleanedData,
      dims = dims, perplexity = perplexity, theta = theta,
      verbose = verbose,
      num_threads = numCores
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





