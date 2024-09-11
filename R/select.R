

setGeneric("filterCCS", function(object, ...) {
  standardGeneric("filterCCS")
})


#' @rdname CCS-method.filterCCS
#' @title CCS method: filterCCS
#' @description \code{filterCCS} method for \code{CCS} class. filterCCS specified submodels
#' @param object a \code{\link{CCS-class}} object
#' @param pattern if \code{method='match'}, it should be characters with the format \code{tissue-cohort}(such as \code{STAD-GSE22377}).
#' @param method One of \code{'grep'} or \code{match}.
#' @inheritParams ccs
#' @importFrom luckyBase LuckyVerbose Fastextra Fastgrep
#' @importFrom stringi stri_reverse
#' @importFrom tidyr %>%
#' @return a new \code{\link{CCS-class}} object
#' @author Weibin Huang</email{hwb2012@@qq.com}>
#' @seealso \code{\link{ccs}}.
#' @exportMethod filterCCS
setMethod(
  "filterCCS",
  signature(object='CCS'),
  function(
    object,
    model.dir = NULL,
    pattern = 'blood',
    method = c('grep','match')[1]
  ){

    # Test ----
    if(F){
      library(luckyBase)
      Plus.library(c('stringi','tidyr'))
      object <- readRDS("E:/Sync/Project/GSClassifier_test2/PADv20240810/optimizeDR/636a52381159cbec9787a0171db7911b/resCCS.rds")
      model.dir <- 'E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810'
      patten = 'Blood/Broad-Blood-2022'
      method = c('grep','match')[1]
      object_raw <- readRDS("E:/iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810/resCCS.rds")

      # path_models_res <- list.files(path = model.dir, pattern = 'modelFit.rds$', recursive = T, full.names = T) %>% .[!grepl(paste0(models_filtered_name,collapse = '|'), .)]

    }

    # Path of submodels ----
    if(is.null(model.dir)){
      model.dir = object@Repeat$model.dir
    }
    path_models <- list.files(path = model.dir, pattern = 'modelFit.rds$', recursive = T, full.names = T)

    # Path of filtered submodels ----
    if(method == 'grep'){
      path_models_filtered <- path_models[grepl(paste0(pattern, collapse = '|'), path_models, ignore.case = T)]
      models_filtered_name <- sapply(path_models_filtered, function(x){
        tissue <- stri_reverse(x) %>% Fastextra(.,'/',3) %>% stri_reverse()
        cohort <- stri_reverse(x) %>% Fastextra(.,'/',2) %>% stri_reverse()
        return(paste0(tissue,'/',cohort))
      }) %>% unlist() %>% as.character()
    } else if(method == 'match'){
      models_filtered_name <- pattern
      path_models_filtered <- path_models[Fastgrep(Fastextra(models_filtered_name,'[/]',2), path_models)]
    } else {
      stop('select.CCS: Wrong method. Please use one of "grep" and "match"!')
    }

    # filtered samples ----
    samples_filtered <- sapply(path_models_filtered, function(x){
      submodel <- readRDS(x)
      return(colnames(submodel$Repeat$Xs))
    }) %>% unlist() %>% as.character()

    # d1 matrix ----
    d1 <- object@Data$Probability$d1
    index <- Fastgrep(gsub('[/]', '\\|', models_filtered_name), colnames(d1))
    d1_res <- d1[!rownames(d1) %in% samples_filtered, -index, drop = FALSE]

    # summary ----
    Data <- list(
      Probability = list(
        d1 = d1_res,
        d2 = NA,
        d3 = NA
      ),
      CCS = NA,
      CancerType = object@Data$CancerType %>% .[!names(.) %in% samples_filtered],
      filtered.cohort = models_filtered_name
    )
    object@Data <- Data

    # output ----
    return(object)
  }
)

