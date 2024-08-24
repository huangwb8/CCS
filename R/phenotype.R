

#' @title Tumor clinical parameters based on molecular subtypes
#' @description Tumor clinical parameters based on molecular subtypes
#' @param data Data frame
#' @param col_subtype Colname of subtype
#' @param col_cohort Colname of cohort
#' @param col_response Colname of effect (positive=1, negative=0)
#' @param col_id Colname of sample ID
#' @param col_tissue Colname of source tissue
#' @importFrom plyr ddply rbind.fill
#' @importFrom dplyr filter reframe arrange desc summarize
#' @return Data frame containing useful clinical parameters:
#' \itemize{
#'   \item \code{RRMean} Mean response rate
#'   \item \code{RRMean_lower} Lower boundary of response rate
#'   \item \code{RRMean_upper} Upper boundary of response rate
#'   \item \code{NRRMean} Mean normalized response rate. The larger (absolute) the better.
#'   \item \code{NRRMean_lower} Lower boundary of normalized response rate
#'   \item \code{NRRMean_upper} Upper boundary of normalized response rate
#'   \item \code{H_index} Over \code{x}% cohorts with response rate over \code{x}%. The larger (raw) the better.
#'   \item \code{tissue_consistence_*} Consistency across cohorts in specified tissues. The larger (raw) the better.
#' }
#' @author Weibin Huang</email{hwb2012@@qq.com}>
#' @export
subtypeEffect <- function(
    data,
    col_subtype = 'Subtype', col_cohort = 'Cohort',
    col_response = 'response', col_id = 'SampleIDs', col_tissue = 'tumor_type'
){

  # Test
  if(F){
    library(luckyBase)
    Plus.library(c('GSClassifier','CCS', 'pROC','ComplexHeatmap','scales','plotly','plyr','dplyr','ggplot2','forestplot','cowplot','ggrepel'))
    data <- readRDS('E:/iProjects/RCheck/GSClassifier/test01/test/Data-for-DecisionEstimateAnalysis.rds')
    col_subtype = 'CCS'
    col_cohort = 'Cohort'
    col_response = 'response'
    col_id = 'SampleIDs'
    col_tissue = 'tumor_type'

    data = df %>% filter(grepl('Immunotherapy', drug_category))
    col_subtype = "CCS"
    col_cohort = "Cohort"
    col_response = "response"
    col_id = "SampleIDs"
    col_tissue = "tissue"
  }

  # Colnames
  col_selected <- c(col_cohort, col_id, col_tissue, col_subtype, col_response)
  df <- data[, col_selected]
  colnames(df) <- c('Cohort','SampleIDs','tumor_type','Subtype','response')

  # Load data
  df2 <- ddply(df, .variables = c('Cohort', 'Subtype', 'response'), plyr::summarize, nSample = length(SampleIDs), tumor_type = unique(tumor_type))
  df2_type <- ddply(df2, .variables	= c('Cohort', 'Subtype'), dplyr::summarize)
  df3_x <- data.frame()
  for(i in 1:nrow(df2_type)){ # i=1
    df2.i <- dplyr::filter(df2, Cohort %in% df2_type[i,'Cohort'], Subtype %in% df2_type[i,'Subtype'])
    df3.i <- data.frame(
      tissue = unique(df2.i$tumor_type),
      size = sum(df2.i$nSample),
      nResponse = ifelse(1 %in% df2.i$response, df2.i$nSample[df2.i$response == 1], 0)
    )
    df3.i$response_rate = df3.i$nResponse/df3.i$size
    df3_x <- rbind(df3_x, df3.i)
  }
  df3 <- cbind(df2_type, df3_x)
  df4 <- ddply(df3, .variables = c('Cohort'), reframe, Subtype = 'All', tissue = unique(tissue),size = sum(size), nResponse = sum(nResponse), response_rate = nResponse/size)
  df6 <- NULL
  for(i in 1:nrow(df3)){
    df3.i <- df3[i,] # i=1
    df3.i_all_response_rate <-  df4[df4$Cohort %in% df3.i$Cohort & df4$Subtype == 'All','response_rate']
    df3.i$normalized_response_rate <- (df3.i$response_rate - df3.i_all_response_rate)/df3.i_all_response_rate
    df6 <- rbind(df6, df3.i)
  }
  df6 <- arrange(df6, Subtype)


  # Index
  all_Subtype <- unique(df6$Subtype)
  all_cohort <- unique(df6$Cohort)
  all_nSample <- sum(df6$size, na.rm = TRUE)
  df7 <- NULL
  for(i in 1:length(all_Subtype)){ # i=15

    Subtype_i <- all_Subtype[i]

    df6.i <- df6[df6$Subtype %in% Subtype_i,]

    # 受测试队列数量及百分比
    # 衡量该亚型的广度
    Subtype_i_AppearCohortNum <- length(unique(df6.i$Cohort))
    Subtype_i_AppearCohortPercent <- round(Subtype_i_AppearCohortNum*100/length(all_cohort), 2)

    # 亚型总病例数及百分比
    # 衡量该亚型的深度
    Subtype_i_PatientNum <- sum(df6.i$size, na.rm = TRUE)
    Subtype_i_PatientPercent <- round(Subtype_i_PatientNum*100/all_nSample, 2)


    # Response rate & 95CI
    # 衡量该亚型的绝对效应
    Subtype_i_vectorRR <- df6.i$response_rate
    # 计算均值
    Subtype_i_vectorRR_Mean <- mean(Subtype_i_vectorRR, na.rm = TRUE)
    if(length(Subtype_i_vectorRR) >= 3){
      n <- length(Subtype_i_vectorRR); alpha <- 1 - 0.95
      t_value <- qt(1 - alpha/2, df = n - 1) # 计算t分布的临界值
      standard_error <- sd(Subtype_i_vectorRR) / sqrt(n) # 计算标准误差
      margin_of_error <- t_value * standard_error # 计算误差边界
      Subtype_i_vectorRR_lower <- max(c(0, Subtype_i_vectorRR_Mean - margin_of_error))
      Subtype_i_vectorRR_upper <- min(c(1, Subtype_i_vectorRR_Mean + margin_of_error))
      # Subtype_i_vectorRR_int <- calculate_bounded_confidence_interval(Subtype_i_vectorRR)
      # Subtype_i_vectorRR_lower <- Subtype_i_vectorRR_int[1]
      # Subtype_i_vectorRR_upper <- Subtype_i_vectorRR_int[2]
    } else {
      Subtype_i_vectorRR_lower <- Subtype_i_vectorRR_upper <- NA
    }

    # Normalized response rate &95CI
    # 衡量该亚型的相对效应
    Subtype_i_vectorNRR <- df6.i$normalized_response_rate
    Subtype_i_vectorNRR_Mean <- mean(Subtype_i_vectorNRR, na.rm = TRUE)
    if(length(Subtype_i_vectorNRR) >= 3){
      n <- length(Subtype_i_vectorNRR); alpha <- 1 - 0.95
      t_value <- qt(1 - alpha/2, df = n - 1) # 计算t分布的临界值
      standard_error <- sd(Subtype_i_vectorNRR) / sqrt(n) # 计算标准误差
      margin_of_error <- t_value * standard_error # 计算误差边界
      Subtype_i_vectorNRR_lower <- Subtype_i_vectorNRR_Mean - margin_of_error
      Subtype_i_vectorNRR_upper <- Subtype_i_vectorNRR_Mean + margin_of_error
      # Subtype_i_vectorNRR_int <- calculate_bounded_confidence_interval(Subtype_i_vectorNRR)
      # Subtype_i_vectorNRR_lower <- Subtype_i_vectorNRR_int[1]
      # Subtype_i_vectorNRR_upper <- Subtype_i_vectorNRR_int[2]
    } else {
      Subtype_i_vectorNRR_lower <- Subtype_i_vectorNRR_upper <- NA
    }


    # Rate-based h-index
    Subtype_i_hindex <- NULL; explorer <- 0.0001; add_unit <- 0.0001
    for(j in 1:round(1/add_unit)){ # j=4117
      j.check <- mean(Subtype_i_vectorRR - explorer >= 0)
      if(j.check >= explorer){
        Subtype_i_hindex <- explorer * 100
      }
      explorer <- explorer + add_unit
    }
    if(is.null(Subtype_i_hindex)){
      Subtype_i_hindex = 0
    }


    # 同种tissue一致性评价
    df_consistence_i <- NULL
    tissue_i <- unique(df6.i$tissue)
    for(j in 1:length(tissue_i)){ # j=1
      tissue_j <- tissue_i[j]
      df6.j <- df6.i[df6.i$tissue %in% tissue_j,]
      if(nrow(df6.j) == 1){
        next
      }
      index_consistence_j <- (max(mean(df6.j$normalized_response_rate > 0), mean(df6.j$normalized_response_rate < 0)) - 0.5) * 2
      res.j <- data.frame(test=index_consistence_j); colnames(res.j) <- paste0('tissue_consistence_',gsub(' ','',tissue_j))
      if(is.null(df_consistence_i)){
        df_consistence_i <- res.j
      } else {
        df_consistence_i <- cbind(df_consistence_i, res.j)
      }
    }

    # Result
    df7_i <- data.frame(
      Subtype = Subtype_i,
      AppearCohortNum = Subtype_i_AppearCohortNum,
      AppearCohortPercent = Subtype_i_AppearCohortPercent,
      PatientNum = Subtype_i_PatientNum,
      PatientPercent = Subtype_i_PatientPercent,
      RRMean = Subtype_i_vectorRR_Mean,
      RRMean_lower = Subtype_i_vectorRR_lower,
      RRMean_upper = Subtype_i_vectorRR_upper,
      NRRMean = Subtype_i_vectorNRR_Mean,
      NRRMean_lower = Subtype_i_vectorNRR_lower,
      NRRMean_upper = Subtype_i_vectorNRR_upper,
      H_index = Subtype_i_hindex,
      stringsAsFactors = F
    )

    if(!is.null(df_consistence_i)){
      df7_i <- cbind(
        df7_i,
        df_consistence_i
      )
    }

    df7 <- plyr::rbind.fill(df7, df7_i)
  }
  df7 <- arrange(df7, desc(AppearCohortNum), desc(H_index))

  return(df7)
}


####%%%%%%%%%%%%% Assistant functions %%%%%%%%%%%%%%%%%%%%####

calculate_bounded_confidence_interval <- function(data, conf_level = 0.95, lower_limit = 0, upper_limit = 1) {
  # 计算样本均值
  sample_mean <- mean(data, na.rm = T)

  # 计算样本标准误差
  standard_error <- sd(data, na.rm = T) / sqrt(length(data), na.rm = T)

  # 计算t值（自由度为n-1）
  degrees_of_freedom <- length(data) - 1
  t_value <- qt((1 + conf_level) / 2, df = degrees_of_freedom)

  # 计算置信区间
  margin_of_error <- t_value * standard_error
  lower_bound <- max(lower_limit, sample_mean - margin_of_error)
  upper_bound <- min(upper_limit, sample_mean + margin_of_error)

  # 返回结果
  return(c(lower_bound, upper_bound))
}






















