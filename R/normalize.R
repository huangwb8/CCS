

#' @title Normalization
#' @description Normalization for \code{CCS} class and \code{d1}
#' @param object \code{CCS} class or \code{d1}
#' @param .fun a function for normalization
#' @importFrom luckyBase Fastextra
#' @return an object after normalization
#' @seealso \code{\link{ccs}}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' @export
normalize <- function(object, .fun = CCS:::softmax){

  # Test
  if(F){
    library(luckyBase)
    object <- readRDS("E:/iProjects/RCheck/GSClassifier/test01/ccs/Testv20240712/resCCS_raw.rds")
    .fun = CCS:::softmax
  }

  # d1
  if('CCS' %in% class(object)){
    is.normalized <- object@Data$Probability$d1_normalized$.true
    if(!is.null(is.normalized)){
      if(is.normalized){
        stop("normilize: You're doing repeated normalization, which is dangerous!")
      }
    }
    d1 <- object@Data$Probability$d1
  } else {
    d1 <- object
  }

  # Reference
  ref <- unique(paste(Fastextra(colnames(d1),'[|]',1), Fastextra(colnames(d1),'[|]',2), sep = '|'))

  # Normalization
  d1_norm <- NULL
  for(i in 1:length(ref)){ # i=1
    ref_i <- ref[i]
    cancer_i <- Fastextra(ref_i,'[|]',1); cohort_i <- Fastextra(ref_i,'[|]',2)
    select <- grepl(cancer_i, colnames(d1)) & grepl(cohort_i,colnames(d1))
    d1_i <- d1[ ,select]
    d1_i_norm <- t(apply(d1_i,1,.fun))
    if(i == 1){
      d1_norm <- d1_i_norm
    } else {
      d1_norm <- cbind(d1_norm, d1_i_norm)
    }
  }

  # Output
  if('CCS' %in% class(object)){
    object@Data$Probability$d1 <- d1_norm
    object@Data$Probability$d1_normalized$.true <- TRUE
    object@Data$Probability$d1_normalized$.fun <- .fun
    object@Data$Probability$d2 <- NA
    object@Data$Probability$d3 <- NA
    return(object)
  } else {
    return(d1_norm)
  }
}
