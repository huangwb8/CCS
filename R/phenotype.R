


#' @title Downstream analysis: Molecular subtypes & Clinical features
#' @description Downstream analysis: Molecular subtypes & Clinical features
#' @inheritParams subtypeEffect
#' @inheritParams CCSPublicParams
#' @param cohort_level The level of cohort used in scater+bar plot
#' @param color_tissue a named vector like \code{c(Brain = "#984EA3",Oral = "#BF5B17")}
#' @param color_cohort Default is \code{NULL}, which means that all cohort names were colored as black. You can also use a named vector like \code{c('Kim2018_PD-1_gastric' = "#E31A1C")}.
#' @param color_subtype Default is \code{NULL}, which means that all raw subtype names were colored as \code{scales::hue_pal()(n)}. You can also use a named vector like \code{c('0' = "#E31A1C")}.
#' @param n_norm_subtype The number of normalized subtypes in discrete mode,
#'   or the maximum number in continuous mode (which also includes K = 1).
#' @param norm_mode Either \code{"discrete"} (default, search only the specified
#'   K) or \code{"continuous"} (search K = 1 through the specified maximum).
#'   Continuous search fixes the state order once, exhaustively evaluates K = 1,
#'   and uses constrained DP within each previous top parent's partitions for
#'   K >= 2. All legal single-split extensions are scored with the public
#'   objective before pruning. Every child preserves all parent cutpoints.
#'   This mode always uses inherited DP, irrespective of \code{norm_search_method}.
#' @param norm_top_n Positive integer, default 10. Maximum number of distinct
#'   candidates retained at each K in continuous mode. Only retained candidates
#'   appear in Raw/Summary; each is checked with the reference ROC evaluator.
#'   Summary combines all K and keeps the existing column schema; the number of
#'   grades can be read from Info or the per-layer SearchDiagnostics. K = 1 has
#'   constant predictor AUC 0.5. Finite-width pruning is approximate and neither
#'   guarantees a global optimum nor monotone improvement as K increases.
#' @param norm_search_method Search algorithm: \code{"bruteforce"} (the default,
#'   retaining all partitions) or \code{"dp"} (deterministic multi-start dynamic
#'   programming followed by local refinement of the selected objective).
#'   DP is exact only for its fixed-direction weighted mean Raw AUC surrogate,
#'   not for any public objective under automatic direction or \code{OptimizeRank}. It runs on one
#'   core, does not change the parallel backend, and returns only evaluated final
#'   candidates plus \code{SearchDiagnostics} under \code{Data$Normalization}.
#'   Oversized exhaustive searches are rejected before allocation. Every cohort
#'   must contain both response classes; required input columns must be complete.
#' @param norm_roc_method The method to call ROC AUC. One of \code{'OptimizeRank'} and \code{'Raw'}.
#' @param norm_objective Candidate scoring strategy (exact names only).
#'   \code{"median"} (default) retains the legacy median including pooled
#'   \code{All}, with missing values removed. \code{"cohort_median"} and
#'   \code{"mean"} use only real cohorts, with equal cohort weights.
#'   \code{"lower_tail"} averages the lowest 25\% of cohort probability mass,
#'   including a fractional boundary cohort. \code{"mean_lower_tail"} averages
#'   that lower tail and the cohort mean with weights 0.5 each. Non-default
#'   strategies require finite AUCs for every real cohort and exclude \code{All}.
#'   Both search methods use the same objective throughout selection/refinement.
#'   New strategies add \code{objective}, \code{objectiveScore},
#'   \code{cohortMedianROCAUC}, \code{cohortMeanROCAUC} and
#'   \code{cohortLowerTailROCAUC} to the normalization Summary; rows sort by
#'   score, cohort mean (descending), then Info (ascending). The existing
#'   \code{medianROCAUC} always retains its legacy definition.
#'   Non-default settings are recorded in \code{Repeat}; absence of
#'   \code{Repeat$norm_objective} means legacy \code{"median"}, including an
#'   explicit median call. With \code{n_norm_subtype = NULL}, the objective is
#'   validated and recorded but no normalization search or diagnostics are made.
#'   These are selection criteria; external predictive gains require held-out
#'   validation with training-derived ranks, cutpoints, direction and mapping.
#' @param norm_rank_marker The marker for normalization rank. Interpretation:
#' \itemize{
#'   \item \code{NRR} Ranked only by normalized response rate.
#'   \item \code{RR} Ranked only by raw response rate.
#'   \item \code{NRR+RR} Ranked by normalized response rate and then by raw response rate.
#'   \item \code{RR+NRR} Ranked by raw response rate and then by normalized response rate.
#'   \item \code{MeanNRR} Equal-cohort mean normalized response rate.
#'   \item \code{LowerTailNRR} Mean of the lowest 25\% cohort NRR mass,
#'     including a fractional boundary cohort.
#'   \item \code{ShrunkLogOR} Equal-cohort mean of Jeffreys-smoothed state
#'     log-odds minus Jeffreys-smoothed overall cohort log-odds.
#'   \item \code{RankFusion} Equal-weight average ranks of median NRR and RR.
#' }
#'   New rules break ties by median NRR, median RR, then state identifier.
#'   Missing states contribute no observation in that cohort; these are
#'   training-derived ordering heuristics, not claims of external superiority.
#' @param plot_layout_widths \code{layout_widths} for \code{\link[patchwork]{plot_layout}}
#' @param plot_ROCCutoff The cut-off value for visualization of ROC AUC
#' @param plot_silence_cohort Whether to silence cohort names in the plot
#' @param plot_silence_subtype Whether to silence subtype names in the plot
#' @param meta_method Random-effects estimator for the meta-analysis. One of \code{'DL'}, \code{'REML'}, \code{'ML'}, \code{'PM'} (Paule-Mandel), \code{'SJ'} (Sidik-Jonkman), \code{'HE'} (Hedges), \code{'EB'} (Empirical Bayes), or \code{'HK'} (Hartung-Knapp using REML + Knapp-Hartung adjustment).
#' @param meta_cohort_font_size Text size (in points) for cohort labels on the meta-analysis plot. Defaults to 15\% larger than the previous styling.
#' @param meta_subtitle_font_size Text size (in points) for the overall summary subtitle (I^2, p, confidence interval) in each meta-analysis panel.
#' @importFrom tidyr `%>%`
#' @importFrom dplyr summarize arrange desc filter
#' @importFrom plyr ddply
#' @importFrom ggtext element_markdown
#' @return A list
#' \itemize{
#'   \item \code{Repeat} Data or parameters for repeatability
#'   \item \code{Plot} 1. Scatter plot + bar plot for RR & NRR; 2. Forest plot for RR.
#'   \item \code{Data} 1. ROC data; 2. Normalized subtype results 3. Clinic utility from \code{\link{subtypeEffect}}.
#' }
#' @examples
#' \dontrun{
#' df <- expand.grid(Cohort = c("A", "B", "C"), Subtype = 1:8,
#'                   response = 0:1, replicate = 1:10)
#' df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
#' df$tumor_type <- "Lung"
#' fit <- subtypePerformance(df, col_subtype = "Subtype",
#'   col_tissue = "tumor_type", n_norm_subtype = 5,
#'   norm_search_method = "dp", norm_roc_method = "Raw", numCores = 1)
#' fit$Data$Normalization$SearchDiagnostics
#' }
#' @author Weibin Huang</email{hwb2012@@qq.com}>
#' @export
subtypePerformance <- function(
    data,
    col_subtype = "CCS",
    col_cohort = "Cohort",
    col_response = "response",
    col_id = "SampleIDs",
    col_tissue = "tissue",
    cohort_level = NULL,
    color_cohort = NULL,
    color_subtype = NULL,
    color_tissue = c(
      Brain = "#984EA3",
      Oral = "#BF5B17",
      Nasopharynx = "#8DD3C7",
      Lung = "#FF6868",
      Gastric = "#FF9843",
      Breast = "#A6D854",
      Colon = "#E31A1C",
      Kedney = "#86A7FC",
      Bladder = "#3468C0",
      Ovary = "#B3B3B3",
      Blood = "#F781BF",
      Skin = "#4DAF4A"
    ),
    n_norm_subtype = 2,
    norm_roc_method = c('OptimizeRank', 'Raw')[1],
    norm_rank_marker = c('NRR', 'RR', 'NRR+RR', 'RR+NRR',
                         'MeanNRR', 'LowerTailNRR', 'ShrunkLogOR', 'RankFusion')[1],
    norm_search_method = c("bruteforce", "dp"),
    norm_mode = c("discrete", "continuous")[1],
    plot_ROCCutoff = 0.8,
    plot_layout_widths = c(9, 1),
    plot_silence_cohort = F,
    plot_silence_subtype = F,
    meta_method = c('DL','REML','ML','PM','SJ','HE','EB','HK'),
    meta_cohort_font_size = 11.5,
    meta_subtitle_font_size = 10,
    numCores = NULL,
    verbose = TRUE,
    norm_objective = "median",
    norm_top_n = 10L
){

  norm_objective <- norm_validate_objective(norm_objective)
  norm_search_method <- match.arg(norm_search_method)
  norm_mode <- norm_validate_mode(norm_mode)
  norm_validate_top_n(norm_top_n)
  norm_validate_rank_marker(norm_rank_marker)

  # Test
  if(F){

    library(luckyBase)
    Plus.library(c('CCS', 'GSClassifier', 'plotly','cowplot','tidyr','ggplot2','ggpubr','purrr','furrr','stringi','digest', 'pROC','ComplexHeatmap','scales','plyr','dplyr','forestplot','ggrepel','writexl','readxl','patchwork','gtable','grid',"reshape2","circlize","parallel","foreach","doParallel","pROC","ggtext"))
    source('./R/ccs_base.R',encoding = 'utf-8')


    # data = read_xlsx('E:/iProjects/RCheck/GSClassifier/routine01/test/DataFrame_CCS+ClinicFeature_PanIMTv20240726+CDSDBv20240726.xlsx')
    data = readRDS('E:/iProjects/RCheck/GSClassifier/routine01/test/df2_v20250704.rds')
    col_subtype = "CCS"
    col_cohort = "Cohort"
    col_response = "response"
    col_id = "SampleIDs"
    col_tissue = "tissue"

    cohort_level = c('Carol2020_PD1_Melanoma','E-MTAB-4030_Nivolumab_Melanoma','Eliezer2015_CTLA-4_Melanoma','Gide2019_PD1_Melanoma','Gide2019_PD1+CTLA4_Melanoma','GSE35640_MAGE-A3_melanoma','GSE78220_PD-1_melanoma','GSE91061_Nivolumab_Melanoma','Hugo2016_PD1_Melanoma','Noam2018_CTLA-4+PD-1_Melanoma','Riaz2017_PD1_Melanoma','GSE162137_Pembrolizumab_Mature T-cell and NK-cell lymphoma','IMvigor210_PD-L1_metastatic urothelial cancer','Braun2020_PD-1_kedney clear cell renal','CheckMate025_PD1_kedney clear cell renal','E-MTAB-3218_Nivolumab_Clear cell renal cell carcinoma','Kim2018_PD-1_gastric', 'GSE154538_Nivolumab_Esophageal cancer','Cho2020_PD-1_Nonsmall-cell lung cancer','Hyunchul2019_PD-1_Nonsmall-cell lung cancer','GSE179730_Nivolumab_Oral squamous cell carcinoma','Zhao2019_PD-1_Glioblastoma')
    # extravalid_cohort = c('LSY01_PD1XCTLA4_LUNG','LSY01_PD1XCTLA4_NPC')

    color_tissue <- c(
      Brain = "#984EA3",
      Oral = "#BF5B17",
      Nasopharynx = "#8DD3C7",
      Lung = "#FF6868",
      Gastric = "#FF9843",
      Breast = "#A6D854",
      Colon = "#E31A1C",
      Kedney = "#86A7FC",
      Bladder = "#3468C0",
      Ovary = "#B3B3B3",
      Blood = "#F781BF",
      Skin = "#4DAF4A"
    )

    color_cohort = c(
      'Kim2018_PD-1_gastric' = "#E31A1C"
    )

    color_subtype = c(
      '9' = "#FB8072"
    )

    n_norm_subtype <- 3
    # n_norm_subtype <- 2
    norm_roc_method = c('OptimizeRank', 'Raw')[2]
    norm_rank_marker = c('NRR', 'RR')[2]
    plot_ROCCutoff = 0.8
    plot_silence_cohort = T
    plot_silence_subtype = T

    numCores = NULL
    verbose = T
    dodge.width = 0.8
    plot_layout_widths = c(9, 1)

  }

  # Data
  if(T){

    # Data: df
    col_selected <- c(col_cohort, col_id, col_tissue, col_subtype, col_response)
    df <- data[col_selected]
    colnames(df) <- c('Cohort','SampleIDs','tumor_type','Subtype','response')
    if(is.null(cohort_level)){
      cohort_level = as.character(unique(df$Cohort))
    } else {
      cohort_level = cohort_level[cohort_level %in% as.character(unique(df$Cohort))]
    }
    if(is.null(cohort_level) || length(cohort_level) == 0){
      cohort_level = as.character(unique(df$Cohort))
    }
    df <- df %>% filter(Cohort %in% cohort_level)

    # Data for performance analysis
    df2 <- ddply(df, c('Cohort', 'Subtype', 'response'), dplyr::summarize, nSample = length(SampleIDs), tumor_type = unique(tumor_type))
    df3 <- ddply(
      df2, c('Cohort', 'Subtype', 'tumor_type'), dplyr::summarize,
      size = sum(nSample, na.rm = TRUE),
      nResponse = ifelse(length(nSample[response %in% c(1,"1")]) == 0, 0, nSample[response %in% c(1,"1")]),
      response_rate = nResponse/size
    )
    df4 <- ddply(df3, c('Cohort'), dplyr::reframe, Subtype = 'All', size = sum(size, na.rm = TRUE), nResponse = sum(nResponse, na.rm = TRUE), response_rate = nResponse/size)
    df6 <- NULL
    for(i in 1:nrow(df3)){
      df3.i <- df3[i,] # i=1
      df3.i_all_response_rate <-  df4[df4$Cohort %in% df3.i$Cohort & df4$Subtype == 'All','response_rate']
      df3.i$normalized_response_rate <- (df3.i$response_rate - df3.i_all_response_rate)/df3.i_all_response_rate
      df6 <- rbind(df6, df3.i)
    }

    # Data: RR/NRR - scatter plot/box plot
    dat.plot <- arrange(df6, Subtype)
    dat.plot$normalized_response_rate <- ifelse(dat.plot$normalized_response_rate >1, 1, ifelse(dat.plot$normalized_response_rate < -1, -1, dat.plot$normalized_response_rate))
    dat.plot$Cohort <- factor(dat.plot$Cohort, levels = rev(cohort_level))
    dat.plot$Subtype <- factor(dat.plot$Subtype)
    nSubtype <- length(unique(dat.plot$Subtype))

  }

  # ROC analysis
  if(verbose) LuckyVerbose('subtypePerformance: ROC analysis...')
  data_roc <- subtypeROC(df, norm_roc_method)

  # Forest plot
  if(verbose) LuckyVerbose('subtypePerformance: Forest plot...')
  plot_f <- forestPlotSubtypeRate(df3)

  # Meta analysis
  if(verbose) LuckyVerbose('subtypePerformance: Meta analysis...')
  valid_meta_methods <- c('DL','REML','ML','PM','SJ','HE','EB','HK')
  if(length(meta_method) == 0 || is.null(meta_method)){
    meta_method <- valid_meta_methods[1]
  } else {
    meta_method <- toupper(trimws(meta_method[1]))
  }
  if(!meta_method %in% valid_meta_methods){
    stop(sprintf(
      "subtypePerformance: meta_method must be one of %s.",
      paste(valid_meta_methods, collapse = ", ")
    ))
  }
  meta_res <- metaSubtypeRate(
    data = df3,
    method = meta_method,
    cohort_font_size = meta_cohort_font_size,
    summary_font_size = meta_subtitle_font_size,
    cohort_levels = cohort_level
  )

  # RR/NRR - scatter plot/box plot
  if(verbose) LuckyVerbose('subtypePerformance: RR/NRR - scatter plot/box plot...')
  plot_r <- plotSubtypeRate(
    dat.plot, data_roc, cohort_level,
    color_tissue, color_cohort, color_subtype,
    dodge.width = 0.8, plot_ROCCutoff, plot_layout_widths,
    plot_silence_cohort, plot_silence_subtype
  ) # plot_r$`Normalized response rate`

  # Performance of normalized subtype
  if(verbose) LuckyVerbose('subtypePerformance: Normailized subtype...')
  if(!is.null(n_norm_subtype)){
    norm_function <- if (norm_mode == "continuous") subtypeNorm_continuous else subtypeNorm
    norm_args <- list(
      df,
      n_norm_subtype = n_norm_subtype,
      norm_roc_method = norm_roc_method,
      norm_rank_marker = norm_rank_marker,
      numCores = numCores,
      verbose = verbose,
      norm_search_method = norm_search_method,
      norm_objective = norm_objective
    )
    if (norm_mode == "continuous") norm_args$norm_top_n <- norm_top_n
    data_norm <- do.call(norm_function, norm_args)
  } else {
    data_norm <- NULL
  }

  # Clinical utility
  if(verbose) LuckyVerbose('subtypePerformance: Clinical utility...')
  data_utility <- subtypeEffect(
    data = df,
    col_subtype = 'Subtype', col_cohort = 'Cohort',
    col_response = 'response', col_id = 'SampleIDs', col_tissue = 'tumor_type'
  )
  plot_utility <- plotSubtypeEffect(data_utility)

  # Output
  l <- list(
    Repeat = list(
      data = data,
      col_subtype = col_subtype,
      col_cohort = col_cohort,
      col_response = col_response,
      col_id = col_id,
      col_tissue = col_tissue,
      cohort_level = cohort_level,
      n_norm_subtype = n_norm_subtype,
      norm_search_method = norm_search_method
    ),
    Plot = list(
      ScatterBarPlot = plot_r,
      ForestPlot = plot_f,
      MetaPlot = meta_res$Plot,
      UtilityPlot = plot_utility
    ),
    Data = list(
      ROC = data_roc,
      Normalization = data_norm,
      MetaAnalysis = meta_res$Data,
      ClinicUtility = data_utility
    )
  )
  if (norm_objective != "median") {
    l$Repeat <- c(l$Repeat, norm_objective_settings(norm_objective))
  }
  if (norm_mode == "continuous") {
    l$Repeat$norm_mode <- norm_mode
    l$Repeat$norm_top_n <- norm_top_n
    l$Repeat$norm_search_method <- "dp"
  }
  if (norm_mode == "continuous" ||
      !norm_rank_marker %in% c("NRR", "RR", "NRR+RR", "RR+NRR")) {
    l$Repeat$norm_rank_marker <- norm_rank_marker
    l$Repeat$norm_roc_method <- norm_roc_method
  }
  if(verbose) LuckyVerbose('subtypePerformance: All done!')
  return(l)

}



####%%%%%%%%%%%%% sub-module functions %%%%%%%%%%%%%%%%%%%%####

#' @importFrom dplyr summarize arrange desc
#' @importFrom plyr ddply rbind.fill
#' @importFrom tidyr `%>%`
subtypeROC <- function(df, norm_roc_method){

  df_roc_01 <- cbind(
    Cohort = 'All',
    nSample = nrow(df),
    get_perform_markers_align(
      df$response, df$Subtype,
      n_norm_subtype = length(unique(df$Subtype)),
      norm_roc_method = norm_roc_method),
    stringsAsFactors = FALSE
  )

  df_roc_02 <- df %>%
    ddply(
      .variables = c('Cohort'),
      .fun = function(sub_df) {
        nSample = nrow(sub_df)
        if (length(unique(sub_df$response)) < 2 | length(unique(sub_df$Subtype)) < 2) {
          sub_df_perf <- data.frame(
            "accuracy" = NA,
            "accuracy_lower"= NA,
            "accuracy_upper"= NA,
            "Balanced Accuracy"= NA,
            "Detection Prevalence"= NA,
            "Detection Rate"= NA,
            "F1" = NA,
            "Neg Pred Value"= NA,
            "Pos Pred Value"= NA,
            "Precision"= NA,
            "Prevalence"= NA,
            "Recall"= NA,
            "ROCAUC" =NA,
            "ROCAUC_lower"= NA,
            "ROCAUC_upper"= NA,
            "Sensitivity"= NA,
            "Specificity"= NA
          )
        } else {
          sub_df_perf <- get_perform_markers_align(
            response = sub_df$response,
            predictor = sub_df$Subtype,
            n_norm_subtype = length(unique(sub_df$Subtype)),
            norm_roc_method = norm_roc_method
          )
        }
        return(cbind(nSample = nSample, sub_df_perf))
      }
    ) %>%
    arrange(desc(ROCAUC))


  df_roc <- rbind.fill(df_roc_01, df_roc_02)
  return(df_roc)
}


#' @importFrom plyr ddply
#' @import ggplot2
#' @import patchwork
#' @importFrom ggpubr rotate_x_text
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @importFrom ggtext element_markdown
plotSubtypeRate <- function(
    dat.plot, data_roc, cohort_level,
    color_tissue, color_cohort, color_subtype,
    dodge.width = 0.8, plot_ROCCutoff = 0.8, plot_layout_widths = c(9, 1),
    plot_silence_cohort, plot_silence_subtype
){

  # Barplot: Tissue
  if(F){

    unique_tissue <- unique(dat.plot$tumor_type)

    dat.plot.bar <- ddply(dat.plot, c('tumor_type', 'Cohort'), dplyr::summarize)
    dat.plot.bar$Cohort <- factor(dat.plot.bar$Cohort, levels = rev(cohort_level))

    stacked_bar <- ggplot(dat.plot.bar, aes(x = 1, y = Cohort, fill = tumor_type)) +
      geom_tile() +
      geom_text(aes(label = tumor_type), color = "white") +
      scale_fill_manual(values = color_tissue[names(color_tissue) %in% unique_tissue]) +
      theme_void() +
      theme(legend.position = "none")
  }

  # Barplot: ROC AUC
  if(T){

    data_roc2 <- data_roc %>% filter(!Cohort %in% 'All') %>% inner_join(., ddply(dat.plot, c('tumor_type', 'Cohort'), dplyr::summarize), by = 'Cohort')
    data_roc2$Cohort <- factor(data_roc2$Cohort, levels = rev(cohort_level))
    unique_tissue <- unique(data_roc2$tumor_type)

    stacked_bar <- ggplot(data_roc2, aes(x = ROCAUC, y = Cohort, fill = tumor_type)) +
      geom_bar(stat="identity", color = "black", linewidth = 1) + # fill='transparent'
      scale_fill_manual(values = color_tissue[names(color_tissue) %in% unique_tissue]) +
      geom_vline(xintercept = plot_ROCCutoff, linetype = "dashed", color = "black", linewidth = 1) +
      scale_x_continuous(
        breaks = c(0,0.2,0.4,0.6,0.8,1),
        labels = c(0,0.2,0.4,0.6,0.8,1),
        limits = c(0, 1), expand = c(0, 0)
      ) +
      labs(title = 'ROCAUC') +
      # theme_minimal() +
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        legend.position = "none"
      )

  }

  # Config colors
  if(T){
    original_labels <- unique(as.character(dat.plot$Cohort))
    styled_labels <- sapply(original_labels, function(label) {
      if (label %in% names(color_cohort)) {
        return(paste0("<span style='color:", as.character(color_cohort[label]) ,";'>", label, "</span>"))
      } else {
        return(label)
      }
    }, USE.NAMES = FALSE)

    if(!is.null(color_subtype)){
      subtype_raw <- levels(dat.plot$Subtype)
      color_subtype_raw <- scales::pal_hue()(length(unique(as.character(dat.plot$Subtype)))); names(color_subtype_raw) <- levels(dat.plot$Subtype)
      color_subtype_update <- sapply(subtype_raw, function(label) {
        if (label %in% names(color_subtype)) {
          return(color_subtype[names(color_subtype) %in% label])
        } else {
          return(color_subtype_raw[label])
        }
      }, USE.NAMES = FALSE)
      boxplot_fill_color <- scale_fill_manual(
        values = as.character(color_subtype_update),
        breaks = names(color_subtype_update)
      )
    } else {
      boxplot_fill_color <- NULL
    }
  }

  # Dodge
  dodge <- position_dodge(width = dodge.width)

  # Scatterplot + boxplot: Response rate
  if(T){
    # Boxplot
    p1 <- ggplot(dat.plot,
                 aes(x=Subtype,
                     y=response_rate,
                     fill = Subtype)) +
      stat_boxplot(geom ='errorbar', width = 0.3,linewidth = 1, position = dodge) +
      geom_boxplot(outlier.shape = NA,width = 0.6,linewidth = 1, color='black', position = dodge) +
      boxplot_fill_color +
      geom_point(position = position_jitterdodge(jitter.width = 0.3,dodge.width = dodge.width), aes(group = Subtype), alpha = 0.35, size = 1.5) +
      labs(x = NULL,y = 'Response rate') +
      guides(color = "none", fill = "none") +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", linewidth=1.5),
        axis.line = element_line(colour = "black", linewidth=0, linetype = 1)
      )

    # Raw response rate - Scatter plot
    p2 <- ggplot(dat.plot, aes(x=Subtype, y=Cohort, colour = response_rate, size = size)) +
      geom_point(alpha = 1) +
      scale_y_discrete(
        breaks = original_labels,
        labels = styled_labels
      ) +
      geom_text_repel(aes(x=Subtype, y=Cohort, label = size), size = 4, color = 'black', box.padding = unit(0.5, "lines")) +
      labs(y = NULL, title = NULL, colour = 'Response rate') +
      scale_colour_gradientn(
        colours  = c("#5F8D4E","#D3D3D3","#D24545"),  #387ADF
        values = scales::rescale(c(0, 0.25, 1)),
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) +
      scale_size(
        range = c(3,12),
        limits = c(1, max(dat.plot$size, na.rm = T)),
        transform = c('reverse',"identity")[2]
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_markdown(),
        panel.border = element_rect(colour = "black", linewidth=1.5),
        axis.line = element_line(colour = "black", linewidth=0, linetype = 1),
        legend.position = 'bottom'
      ) +
      rotate_x_text(angle = 90)

    # Whether to hide labels of cohorts/subtypes
    if(plot_silence_cohort){
      p2 <- p2 + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
    }
    if(plot_silence_subtype){
      p2 <- p2 + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    }

    # Plot: patchwork grid
    p12 <- (
      p1 + theme(axis.title = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()) +
        plot_spacer() +
        p2 +
        stacked_bar
    ) +
      plot_layout(ncol = 2, nrow = 2,
                  widths = plot_layout_widths,
                  heights = c(2, 8)) &
      theme(plot.margin = margin(1, 1, 1, 1))

    # Rotation of p12
    if(T){
      p1_rotated <- p1 + coord_flip()
      p2.2 <- p2; p2.2$layers[[2]] <- NULL
      p2_rotated <- p2.2 +
        labs(x = NULL, y = 'Cohort') +
        coord_flip() +
        theme(
          legend.position = 'bottom',
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
          )
      stacked_bar_rotated <- stacked_bar +
        coord_flip() +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_text()
        )
      p12_rotated <- (
        plot_spacer() + stacked_bar_rotated +
          p1_rotated + theme(axis.title.x = element_blank()) +
          p2_rotated
      ) +
        plot_layout(
          ncol = 2,
          nrow = 2,
          widths = c(2, 8),
          heights = rev(plot_layout_widths)
        ) &
        theme(plot.margin = margin(1, 1, 1, 1))
      # print(p12_rotated)
    }

  }

  # Scatterplot + boxplot: Normalized response rate
  if(T){
    p3 <- ggplot(dat.plot,
                 aes(x=Subtype,
                     y=normalized_response_rate,
                     fill = Subtype)) +
      stat_boxplot(geom ='errorbar', width = 0.3,linewidth = 1, position = dodge) +
      boxplot_fill_color +
      geom_boxplot(outlier.shape = NA,width = 0.6,linewidth = 1, color='black', position = dodge) +
      geom_point(position = position_jitterdodge(jitter.width = 0.3,dodge.width = dodge.width), aes(group = Subtype), alpha = 0.35, size = 1.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = NULL,y = 'Normalized response rate') +
      guides(color = "none", fill = "none") +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", linewidth=1.5),
        axis.line = element_line(colour = "black", linewidth=0, linetype = 1)
      )

    # Normalized response rate - Scatter plot
    p4 <- ggplot(dat.plot, aes(x=Subtype, y=Cohort, colour = normalized_response_rate, size = size)) +
      geom_point(alpha = 1) +
      scale_y_discrete(
        breaks = original_labels,
        labels = styled_labels
      ) +
      geom_text_repel(aes(x=Subtype, y=Cohort, label = size), size = 4, color = 'black', box.padding = unit(0.5, "lines")) +
      labs(y = NULL, title = NULL, colour = 'Normalized response rate') +
      scale_colour_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        limits = c(-1, 1),  midpoint = 0,
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour",
        labels = c('≤-1.0','-0.5','0','0.5','≥1.0')
      ) +
      scale_size(
        range = c(3,12),
        limits = c(1,max(dat.plot$size, na.rm = T)),
        transform = c('reverse',"identity")[2]
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_markdown(),
        panel.border = element_rect(colour = "black", linewidth=1.5),
        axis.line = element_line(colour = "black", linewidth=0, linetype = 1),
        legend.position = 'bottom'
      ) + rotate_x_text(angle = 90)


    # Whether to hide labels of cohorts/subtypes
    if(plot_silence_cohort){
      p4 <- p4 + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank())
    }
    if(plot_silence_subtype){
      p4 <- p4 + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    }

    # Plot: patchwork grid
    p34 <- (
      p3 + theme(axis.title = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()) +
        plot_spacer() +
        p4 +
        stacked_bar
    ) +
      plot_layout(ncol = 2, nrow = 2,
                  widths = plot_layout_widths,
                  heights = c(2, 8)) &
      theme(plot.margin = margin(1, 1, 1, 1))

    # Rotation of p34
    if(T){
      p3_rotated <- p3 + coord_flip()
      p4.2 <- p4; p4.2$layers[[2]] <- NULL
      p4_rotated <- p4.2 +
        labs(x = NULL, y = 'Cohort') +
        coord_flip() +
        theme(
          legend.position = 'bottom',
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      stacked_bar_rotated <- stacked_bar +
        coord_flip() +
        theme(
          axis.text.x = element_blank(),
          axis.text.y = element_text()
        )
      p34_rotated <- (
        plot_spacer() + stacked_bar_rotated +
          p3_rotated + theme(axis.title.x = element_blank()) +
          p4_rotated
      ) +
        plot_layout(
          ncol = 2,
          nrow = 2,
          widths = c(2, 8),
          heights = rev(plot_layout_widths)
        ) &
        theme(plot.margin = margin(1, 1, 1, 1))
      # print(p34_rotated)
    }

  }

  # Merge
  if(T){
    p1234 <- (
      p1 + theme(axis.title = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()) +
      p3 + theme(axis.title = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()) +
        scale_y_continuous(position = "right") +
      plot_spacer() +
      p2 + guides(size = "none") +
      p4 + guides(size = "none") + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank()) +
      stacked_bar
    ) +
      plot_layout(
        ncol = 3, nrow = 2,
        widths = c(plot_layout_widths[1]/2,plot_layout_widths[1]/2,plot_layout_widths[2]),
        heights = c(2, 8)
      ) &
      theme(plot.margin = margin(1, 1, 1, 1))

    p1234_rotated <- (
      plot_spacer() +
      stacked_bar_rotated +
      stacked_bar_rotated + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      plot_spacer() +
      p1_rotated + theme(axis.title.x = element_blank()) +
      p2_rotated + guides(size = "none") +
      p4_rotated + guides(size = "none") +
      p3_rotated + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_blank())
    ) +
      plot_layout(
        ncol = 4, nrow = 2,
        widths = c(1,4,4,1),
        heights = c(2, 8)
      ) &
      theme(plot.margin = margin(1, 1, 1, 1)); # print(p1234_rotated)

  }

  # https://wilkelab.org/cowplot/reference/plot_grid.html
  # Not work well. Deprecated
  # p12 <- plot_grid(
  #   p1 + theme(axis.title = element_blank(),
  #              axis.text.x = element_blank(),
  #              axis.ticks.x = element_blank()),
  #   NULL,
  #   p2,
  #   stacked_bar,
  #   align = c("hv"),
  #   axis =  "tblr",
  #   ncol = 2, nrow = 2,
  #   rel_widths = c(0.7, 0.3),
  #   rel_heights = c(0.2, 0.8)
  # )

  return(list(
    'Response rate' = p12,
    'Normalized response rate' = p34,
    'Merge' = p1234,
    'Rotated response rate' = p12_rotated,
    'Rotated normalized response rate' = p34_rotated,
    'Rotated merge' = p1234_rotated
  ))

}


#' @importFrom plyr ddply
#' @importFrom tidyr `%>%`
#' @importFrom dplyr summarize arrange desc filter
#' @importFrom forestplot forestplot fpColors fpTxtGp
#' @importFrom grid unit
#' @importFrom stats binom.test fisher.test
forestPlotSubtypeRate <- function(df3){

  # Data: statistics
  if(T){

    # confidential intervals
    df3_conf <-  df3[c("Cohort", "Subtype", "size", "nResponse")] %>% ddply(., c('Subtype'), dplyr::summarize, size = sum(size, na.rm = T), nResponse=sum(nResponse, na.rm = T))
    No_allPatient <- sum(df3_conf$size, na.rm = T)
    No_allresponder <- sum(df3_conf$nResponse, na.rm = T)
    # binom.test(No_allresponder, No_allPatient)
    p_all <- No_allresponder / No_allPatient
    df3_conf2 <- data.frame()
    for(i in 1:nrow(df3_conf)){ # i=1
      df3_conf_i <- df3_conf[i,]
      # H0: ICI治疗对该亚组的患者没有积极影响。
      res.i <- binom.test(df3_conf_i$nResponse,df3_conf_i$size,conf.level= 0.95, p = p_all) # p_all
      df3_conf_i$rr <- res.i[["estimate"]]
      res.i.conf <- res.i[["conf.int"]]
      df3_conf_i$rr_lower <- res.i.conf[[1]]
      df3_conf_i$rr_upper <- res.i.conf[[2]]
      df3_conf_i$binom_P <- res.i$p.value

      # Fisher's Exact Test
      # 干预组-阳性，对照组-阳性， 干预组-阴性， 对照组-阴性
      fisher_i <- c(df3_conf_i$nResponse, No_allresponder, df3_conf_i$size - df3_conf_i$nResponse,  No_allPatient - No_allresponder)
      dim(fisher_i) <- c(2,2)
      fisher_i_res <- fisher.test(fisher_i)
      df3_conf_i$fisher_or <- fisher_i_res$estimate
      df3_conf_i$fisher_or_lower <- fisher_i_res$conf.int[1]
      df3_conf_i$fisher_or_upper <- fisher_i_res$conf.int[2]
      df3_conf_i$fisher_P <- fisher_i_res$p.value

      # Output
      df3_conf2 <- rbind(df3_conf2, df3_conf_i)
    }

  }

  # Plot
  if(T){

    # data
    # df3_conf3 <- df3_conf2 %>% filter(size >= 10) %>% arrange(binom_P, desc(size))
    df3_conf3 <- df3_conf2 %>% filter(size >= 10) %>% arrange(desc(rr), binom_P, desc(size))
    df3_conf3_colsum <- colSums(df3_conf3)
    tabletext <- cbind(
      c("Subtype","\n", df3_conf3$Subtype),
      # c(paste0("No. of patients \n n=",as.integer(df3_conf3_colsum['size'])),"\n", df3_conf3$size),
      # c(paste0("No. of responders \n n=", as.integer(df3_conf3_colsum['nResponse'])),"\n", df3_conf3$nResponse),
      c(paste0("No. of patients \n n=", No_allPatient),"\n", df3_conf3$size),
      c(paste0("No. of responders \n n=", No_allresponder),"\n", df3_conf3$nResponse),
      c(paste0("95% CI"),"\n", paste(sprintf("%.3f",df3_conf3$rr),' ','[', sprintf("%.3f",df3_conf3$rr_lower),'-',  sprintf("%.3f",df3_conf3$rr_upper), ']', sep = '')),
      c('P Value',"\n",sprintf("%.4f",df3_conf3$binom_P))
    )

    # Dynamic box size
    # dynamic_boxsize <- function(study_sizes){
    #
    #   # Define the desired minimum and maximum box sizes for visual appeal
    #   min_display_boxsize <- 0.15 # Smallest box size for the smallest study
    #   max_display_boxsize <- 0.4  # Largest box size for the largest study
    #
    #   # Handle cases:
    #   # 1. All studies have the same size.
    #   # 2. Studies have different sizes (requires scaling).
    #   # 3. Only one study (assign an average or max box size).
    #   if (length(study_sizes) == 1) {
    #     # If there's only one study/subgroup, assign a medium or max box size
    #     dynamic_box_sizes <- rep((min_display_boxsize + max_display_boxsize) / 2, length(study_sizes))
    #     # Or simply: dynamic_box_sizes <- max_display_boxsize
    #   } else if (min(study_sizes) == max(study_sizes)) {
    #     # If all studies have the same sample size, use an average box size
    #     dynamic_box_sizes <- rep((min_display_boxsize + max_display_boxsize) / 2, length(study_sizes))
    #   } else {
    #     # Linearly scale sample sizes to the range [min_display_boxsize, max_display_boxsize]
    #     # Formula: new_value = new_min + (value - old_min) * (new_max - new_min) / (old_max - old_min)
    #     dynamic_box_sizes <- min_display_boxsize +
    #       (study_sizes - min(study_sizes)) *
    #       (max_display_boxsize - min_display_boxsize) /
    #       (max(study_sizes) - min(study_sizes))
    #   }
    #   return(dynamic_box_sizes)
    # }

    # plot
    if(T){

      p_f <- forestplot(
        labeltext=tabletext, graph.pos=5,
        mean=c(NA,NA,df3_conf3$rr),
        lower=c(NA,NA,df3_conf3$rr_lower),
        upper=c(NA,NA,df3_conf3$rr_upper),
        xlog = F,
        xticks = seq(0.1,0.9,0.2), xticks.digits = 1,
        xlab="Response rate",
        col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
        txt_gp=fpTxtGp(label=gpar(cex = 1.1),
                       ticks=gpar(cex = 1.1),
                       xlab=gpar(cex = 1.4),
                       title=gpar(cex = 1.1)),
        zero = p_all,
        cex = 0.9, lineheight = "auto",
        colgap = unit(8,"mm"),
        lwd.ci=2,
        # boxsize=0.3*nrow(df3_conf3)/25,
        boxsize = 0.275,
        ci.vertices=TRUE, ci.vertices.height = 0.25
      )
      # α=0.0; cairo_pdf(paste0(resPath, '/forestplot_response rate across Subtype subtypes.pdf'),width = 10*(1-α), height = 6*(1-α)*nrow(tabletext)/15); print(p_f); dev.off()
    }

  }

  return(p_f)
}


#' @importFrom plyr ddply
#' @importFrom tidyr `%>%`
#' @importFrom dplyr summarize arrange desc filter group_by ungroup
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach `%dopar%`
#' @import luckyBase
subtypeNorm <- function(
    df,
    n_norm_subtype=2,
    norm_roc_method = 'Raw',
    norm_rank_marker = 'RR',
    numCores = NULL, verbose = TRUE,
    norm_search_method = c("bruteforce", "dp"),
    norm_objective = "median"
){

  norm_objective <- norm_validate_objective(norm_objective)
  # Test
  if (FALSE) {
    df <- expand.grid(Cohort = c("A", "B"), Subtype = 1:4, response = 0:1)
    df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
    df$tumor_type <- "Lung"
    subtypeNorm(df, 3, "Raw", "RR", 1, FALSE, "dp", "mean_lower_tail")
  }
  norm_search_method <- match.arg(norm_search_method)
  norm_validate_input(df, n_norm_subtype, norm_roc_method, numCores)
  norm_validate_rank_marker(norm_rank_marker)
  search_started <- proc.time()[["elapsed"]]

  # Data
  if(T){
    df2 <- ddply(df, c('Cohort', 'Subtype', 'response'), dplyr::summarize, nSample = length(SampleIDs), tumor_type = unique(tumor_type))
    df3 <- ddply(
      df2, c('Cohort', 'Subtype', 'tumor_type'), dplyr::summarize,
      size = sum(nSample, na.rm = TRUE),
      nResponse = ifelse(length(nSample[response %in% c(1,"1")]) == 0, 0, nSample[response %in% c(1,"1")]),
      response_rate = nResponse/size
    )
    df4 <- ddply(df3, c('Cohort'), dplyr::reframe, Subtype = 'All', size = sum(size, na.rm = TRUE), nResponse = sum(nResponse, na.rm = TRUE), response_rate = nResponse/size)
    df6 <- NULL
    for(i in 1:nrow(df3)){
      df3.i <- df3[i,] # i=1
      df3.i_all_response_rate <-  df4[df4$Cohort %in% df3.i$Cohort & df4$Subtype == 'All','response_rate']
      df3.i$normalized_response_rate <- (df3.i$response_rate - df3.i_all_response_rate)/df3.i_all_response_rate
      df6 <- rbind(df6, df3.i)
    }
    df_cutoff <- df6 %>% ddply(
      c('Subtype'), dplyr::summarize,
      medianNRR = median(normalized_response_rate, na.rm = T),
      medianRR = median(response_rate, na.rm = T)
    )

    # Rank
    df_cutoff <- norm_rank_subtypes(df, norm_rank_marker, df_cutoff)

    if (n_norm_subtype > nrow(df_cutoff)) {
      stop("n_norm_subtype must not exceed the number of observed subtypes.")
    }
    if (norm_search_method == "dp") {
      result <- subtype_norm_dp(df, df_cutoff, n_norm_subtype,
                                 norm_roc_method, search_started, verbose,
                                 objective = norm_objective)
      result$SearchDiagnostics$rank_marker <- norm_rank_marker
      return(result)
    }
    n_partitions <- choose(nrow(df_cutoff) - 1, n_norm_subtype - 1)
    estimated_result_bytes <- n_partitions * (length(unique(df$Cohort)) + 1) * 192
    if (!is.finite(n_partitions) || n_partitions > 1e6 || estimated_result_bytes > 2 * 1024^3) {
      stop("Brute-force search exceeds the partition/result-memory budget; use norm_search_method = 'dp'.")
    }
    combinations <- combn(seq_len(nrow(df_cutoff) - 1), n_norm_subtype - 1)
  }

  # Parallel processing: Find best combination
  if(T){

    # Parallel: Start
    if(is.null(numCores)){
      numCores <- detectCores() - 1  # Use all but one core
    }
    registerDoParallel(cores = numCores)
    on.exit(stopImplicitCluster(), add = TRUE)

    if(verbose) LuckyVerbose('subtypeNorm: Parallel processing - Find best combination...',levels = 2)

    # Splite data
    if(ncol(combinations) >= 10*numCores){
      XL <- cut_vector(1:ncol(combinations), nsplit = numCores)
    } else {
      XL <- cut_vector(1:ncol(combinations), nsplit = 4)
    }

    system.time(
      df_perform <- foreach(
        i = 1:length(XL), .combine = 'rbind',
        .packages = c("plyr", "dplyr", "dplyr", "luckyBase","pROC")
      ) %dopar% {
        batch_results <- lapply(XL[[i]], function(k){ # k=1

          # 基于某个特定的标准化方案（某个Info）计算

          x <- combinations[, k]
          medianNRR_category <- rep(NA, nrow(df_cutoff))

          # Calculate categories
          for (j in 0:length(x)) {
            lower_limit_j <- ifelse(j == 0, 1, (x[j] + 1))
            upper_limit_j <- ifelse(j == length(x), nrow(df_cutoff), x[j + 1])
            medianNRR_category[lower_limit_j:upper_limit_j] <- n_norm_subtype - j
          }
          df_cutoff$medianNRR_category <- medianNRR_category - 1

          # Summarize cutoff data
          df_cutoff_summary <- summarize(
            group_by(df_cutoff, medianNRR_category),
            annotation = paste0(Subtype, collapse = ',')
          )

          df$normSubtype <- convert(df$Subtype, "Subtype", "medianNRR_category", df_cutoff) %>% as.numeric()

          # Group relationship
          info_str <- paste(
            paste(
              df_cutoff_summary$medianNRR_category,
              df_cutoff_summary$annotation,
              sep = "="
            ),
            collapse = "; "
          )

          # Create summary for all data
          df_perform_i1 <- cbind(
            Cohort = 'All',
            Info = info_str,
            nSample = nrow(df),
            get_perform_markers(
              df$response, df$normSubtype,
              n_norm_subtype = n_norm_subtype,
              norm_roc_method = norm_roc_method),
            stringsAsFactors = FALSE
          )

          # Create summary by cohort
          df_perform_i2 <- df %>% ddply(
            ., c('Cohort'),
            dplyr::summarize,
            Info = info_str,
            nSample = length(Cohort),
            get_perform_markers(
              response, normSubtype,
              n_norm_subtype = n_norm_subtype,
              norm_roc_method = norm_roc_method)
          )

          rbind(df_perform_i1, as.data.frame(df_perform_i2))
        })
        do.call('rbind', batch_results)
      }
    )
    # 用户   系统   流逝
    # 0.86   0.52 249.47
    stopImplicitCluster()

  }

  # Summary
  # cl <- makeCluster(numCores)
  # registerDoParallel(cl)
  # df_perform <- df_perform %>% arrange(Cohort, desc(ROCAUC))
  if(verbose) LuckyVerbose('subtypeNorm: Summary...',levels = 2)
  system.time(
    df_perform_summary <- df_perform %>%
      group_by(Info) %>%
      dplyr::summarize(medianROCAUC = median(ROCAUC, na.rm = TRUE)) %>%
      ungroup() %>%
      arrange(desc(medianROCAUC))
  )
  if (norm_objective != "median") {
    df_perform_summary <- norm_objective_summary(df_perform, norm_objective)
  }
  # system.time({
  #   df_perform_summary <- ddply(
  #     df_perform,
  #     .variables = "Info",
  #     .fun = dplyr::summarize,
  #     medianROCAUC = median(ROCAUC, na.rm = TRUE)
  #   ) %>% arrange(desc(medianROCAUC))
  # })
  # stopCluster(cl)


  # Output
  if(verbose) LuckyVerbose('subtypeNorm: Done!',levels = 2)
  return(
    list(
      Raw = df_perform,
      Summary = df_perform_summary
    )
  )

}


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
  df2 <- ddply(df, .variables = c('Cohort', 'Subtype', 'response'), dplyr::summarize, nSample = length(SampleIDs), tumor_type = unique(tumor_type))
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
    # 0为趋势完全不一致，1为趋势完全一致。>0或<1为趋势部分一致。
    # 只要不是0都可以接受。
    df_consistence_i <- NULL
    tissue_i <- unique(df6.i$tissue)
    for(j in 1:length(tissue_i)){ # j=1
      tissue_j <- tissue_i[j]
      df6.j <- df6.i[df6.i$tissue %in% tissue_j,]
      if(nrow(df6.j) == 1){
        next
      }
      index_consistence_j <- (max(mean(df6.j$normalized_response_rate > 0, na.rm = T), mean(df6.j$normalized_response_rate < 0, na.rm = T)) - 0.5) * 2
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


#' @title Plot for subtypeEffect
#' @description Plot for subtypeEffect
#' @param data Data from \code{\link{subtypeEffect}} function
#' @importFrom dplyr %>% mutate select all_of distinct
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_remove
#' @importFrom ggplot2 ggplot aes geom_col geom_text scale_y_continuous labs coord_cartesian theme element_blank element_text element_line geom_hline geom_errorbar geom_point scale_fill_gradientn geom_tile margin
#' @importFrom patchwork plot_layout
#' @importFrom scales number rescale squish
#' @return Combined ggplot object
#' @author Weibin Huang</email{hwb2012@@qq.com}>
#' @export
plotSubtypeEffect <- function(data) {
  # Dependencies: ggplot2, dplyr, tidyr, stringr, patchwork, scales
  # suppressPackageStartupMessages({
  #   library(dplyr)
  #   library(tidyr)
  #   library(stringr)
  #   library(ggplot2)
  #   library(patchwork)
  #   library(scales)
  # })

  # Ensure Subtype is an ordered factor (numeric order if numeric-like, else natural order)
  df <- data %>%
    mutate(Subtype = {
      if (suppressWarnings(all(!is.na(as.numeric(as.character(Subtype)))))) {
        factor(Subtype, levels = sort(unique(as.numeric(as.character(Subtype)))))
      } else {
        factor(Subtype, levels = unique(Subtype))
      }
    })

  # Publication-friendly shared theme
  theme_shared <- theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 12),
      axis.title.y = element_text(margin = margin(r = 8)),
      # axis.text.x = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10)
    )

  # Helper: add a top expansion to make space for text labels above bars
  top_expand <- expansion(mult = c(0.02, 0.10))

  # 1) AppearCohortPercent bar + label AppearCohortNum
  p1 <- ggplot(df, aes(x = Subtype, y = AppearCohortPercent)) +
    geom_col(fill = "#2E86AB", width = 0.7) +
    geom_text(aes(label = AppearCohortNum), vjust = -0.35, size = 3.3) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = top_expand) +
    labs(title = "Cohort (%)", y = "Percentage", x = NULL) +
    coord_cartesian(clip = "off") +
    theme_shared +
    theme(plot.margin = margin(t = 8, r = 6, b = 4, l = 6))

  # 2) PatientPercent bar + label PatientNum
  p2 <- ggplot(df, aes(x = Subtype, y = PatientPercent)) +
    geom_col(fill = "#F39C12", width = 0.7) +
    geom_text(aes(label = PatientNum), vjust = -0.35, size = 3.3) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), expand = top_expand) +
    labs(title = "Patient (%)", y = "Percentage", x = NULL) +
    coord_cartesian(clip = "off") +
    theme_shared +
    theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 6))

  # 3) RRMean with CI as point + errorbar
  p3 <- ggplot(df, aes(x = Subtype, y = RRMean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_errorbar(aes(ymin = RRMean_lower, ymax = RRMean_upper),
                  width = 0.2, color = "#27AE60", linewidth = 0.7) +
    geom_point(size = 2.6, color = "#27AE60") +
    labs(title = "RR (Mean [95% CI])", y = "Value", x = NULL) +
    theme_shared +
    theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 6))

  # 4) NRRMean with CI as point + errorbar
  p4 <- ggplot(df, aes(x = Subtype, y = NRRMean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_errorbar(aes(ymin = NRRMean_lower, ymax = NRRMean_upper),
                  width = 0.2, color = "#E74C3C", linewidth = 0.7) +
    geom_point(size = 2.6, color = "#E74C3C") +
    labs(title = "NRR (Mean [95% CI])", y = "Value", x = NULL) +
    theme_shared +
    theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 6))

  # 5) H_index bar + label (NEW: add value labels above bars; keep consistent width)
  p5 <- ggplot(df, aes(x = Subtype, y = H_index)) +
    geom_col(fill = "#6C5CE7", width = 0.7) +
    geom_text(
      aes(label = scales::number(H_index, accuracy = 0.01)), # adjust accuracy as needed
      vjust = -0.35, size = 3.3, color = "black"
    ) +
    scale_y_continuous(expand = top_expand) +
    labs(title = "H_index", y = "Score", x = NULL) +
    coord_cartesian(clip = "off") +
    theme_shared +
    theme(plot.margin = margin(t = 4, r = 6, b = 4, l = 6))

  # 6) Heatmap for tissue_consistence_* columns (auto-detect)
  tissue_cols <- grep("^tissue_consistence_", names(df), value = TRUE)

  tissue_long <- df %>%
    select(Subtype, all_of(tissue_cols)) %>%
    pivot_longer(cols = -Subtype, names_to = "tissue", values_to = "value") %>%
    mutate(
      tissue = stringr::str_remove(tissue, "^tissue_consistence_"),
      tissue = factor(tissue, levels = sort(unique(tissue)))
    )

  # Compute a sensible height ratio for the heatmap
  n_tissues <- tissue_long %>% distinct(tissue) %>% nrow()
  heatmap_height <- max(2.8, n_tissues * 0.35)

  p6 <- ggplot(tissue_long, aes(x = Subtype, y = tissue, fill = value)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = c("#FFFFFF", "#FFE6DC" ,"#FFAD95", "#FF4D40", "#FE0000"),
      values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
      limits = c(0, 1),
      oob = squish,
      na.value = "grey",
      name = NULL
      # name = "Consistency"
    ) +
    labs(title = "Tissue consistency", x = NULL, y = "Tissue") +
    theme_shared +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(t = 4, r = 6, b = 6, l = 6),
      axis.text.x = element_text(size = 15, face = "bold")
    )

  # Stack all plots vertically using patchwork
  layout_heights <- c(2.2, 2.2, 1.8, 1.8, 1.8, heatmap_height)

  # Collect guides and place legend at bottom to keep panel widths consistent
  combined <- (p1 / p2 / p3 / p4 / p5 / p6) +
    plot_layout(heights = layout_heights, guides = "collect") &
    theme(legend.position = "bottom")

  return(combined)
}


#' @title Get new normalized subtypes based on a record from \code{subtypePerformance}
#' @description Get new normalized subtypes based on a record from \code{subtypePerformance}
#' @param x a character
#' @param record a record from \code{subtypePerformance}
#' @importFrom luckyBase Fastextra convert
#' @importFrom tidyr `%>%`
#' @return a character
#' @seealso \code{\link{subtypePerformance }}.
#' @author Weibin Huang<\email{hwb2012@@qq.com}>
#' @export
newSubtype <- function(x, record){
  record_sub <- Fastextra(record, '; ')

  record_annot <- lapply(record_sub, function(x){
    data.frame(
      normSubtype = Fastextra(x,'=',1),
      Subtype = Fastextra(x,'=',2) %>% Fastextra(., ','),
      stringsAsFactors = F
    )
  }) %>% do.call('rbind', .)

  x2 <- convert(x, 'Subtype','normSubtype', record_annot) %>% as.numeric()
  return(x2)
}


#' @title Meta-analysis of subtype response rates
#' @description Perform a random-effects meta-analysis across cohorts for every subtype and prepare publication-ready data and plots.
#' @param data A data.frame with columns \code{Cohort}, \code{Subtype}, \code{size}, \code{nResponse}, and optionally \code{response_rate} and \code{tumor_type}.
#' @param method Random-effects estimator code. One of \code{'DL'}, \code{'REML'}, \code{'ML'}, \code{'PM'} (Paule-Mandel), \code{'SJ'} (Sidik-Jonkman), \code{'HE'} (Hedges), \code{'EB'} (Empirical Bayes), or \code{'HK'} (Hartung-Knapp using REML with Knapp-Hartung adjustment).
#' @param conf.level Confidence level for study-wise and pooled intervals.
#' @param summary_label Label used for the pooled row in the forest-style plot.
#' @param min_studies Minimum number of cohorts required to report heterogeneity statistics.
#' @param cohort_font_size Text size (in points) for cohort labels on the plot.
#' @param summary_font_size Text size (in points) for the subtitle (overall CI/I^2/p info).
#' @param cohort_levels Optional character vector specifying the display order of cohorts (top to bottom, excluding the summary row).
#' @importFrom stats binom.test pchisq qnorm
#' @return A list with \code{Data} (data.frame) and \code{Plot} (ggplot).
metaSubtypeRate <- function(
    data,
    method = c('DL','REML','ML','PM','SJ','HE','EB','HK'),
    conf.level = 0.95,
    summary_label = "Overall",
    min_studies = 1,
    cohort_font_size = 11.5,
    summary_font_size = 10,
    cohort_levels = NULL
){

  valid_methods <- c('DL','REML','ML','PM','SJ','HE','EB','HK')
  if(length(method) == 0 || is.null(method)){
    method <- valid_methods[1]
  } else {
    method <- toupper(trimws(method[1]))
  }
  if(!method %in% valid_methods){
    stop(sprintf(
      "metaSubtypeRate: method must be one of %s.",
      paste(valid_methods, collapse = ", ")
    ))
  }
  if(!is.numeric(cohort_font_size) || length(cohort_font_size) == 0 || !is.finite(cohort_font_size)){
    cohort_font_size <- 11.5
  }
  if(!is.numeric(summary_font_size) || length(summary_font_size) == 0 || !is.finite(summary_font_size)){
    summary_font_size <- 10
  }
  if(!is.null(cohort_levels)){
    cohort_levels <- as.character(cohort_levels)
  }
  required_cols <- c("Cohort", "Subtype", "nResponse", "size")
  if(!all(required_cols %in% colnames(data))){
    stop("metaSubtypeRate: data must contain Cohort, Subtype, nResponse, and size columns.")
  }

  df <- data
  if(!"tumor_type" %in% colnames(df)){
    df$tumor_type <- NA_character_
  }
  if(!"response_rate" %in% colnames(df)){
    df$response_rate <- with(df, ifelse(size > 0, nResponse/size, NA))
  } else {
    idx_fill <- which(is.na(df$response_rate) & df$size > 0)
    df$response_rate[idx_fill] <- df$nResponse[idx_fill]/df$size[idx_fill]
  }

  df <- df[is.finite(df$size) & df$size > 0, , drop = FALSE]
  if(nrow(df) == 0){
    return(list(Data = NULL, Plot = NULL))
  }

  df$Subtype <- as.character(df$Subtype)
  df$Cohort <- as.character(df$Cohort)
  df$response_rate <- pmax(pmin(df$response_rate, 1 - 1e-6), 1e-6)

  inv_logit <- function(x){
    exp(x)/(1 + exp(x))
  }

  use_hk <- identical(method, "HK")
  rma_method <- if(use_hk) "REML" else method
  method_label <- if(use_hk) "REML+HK" else method

  estimate_tau2_dl <- function(yi, vi){
    k <- length(yi)
    if(k <= 1){
      return(0)
    }
    wi <- 1/vi
    mu_fixed <- sum(wi * yi)/sum(wi)
    Q <- sum(wi * (yi - mu_fixed)^2)
    C <- sum(wi) - sum(wi^2)/sum(wi)
    if(C <= 0){
      return(0)
    }
    max(0, (Q - (k - 1))/C)
  }

  format_p <- function(x){
    if(is.na(x)){
      return("NA")
    }
    if(x < 1e-3){
      return("<0.001")
    }
    sprintf("%.3f", x)
  }

  split_df <- split(df, df$Subtype)
  data_meta <- lapply(names(split_df), function(subtype_key){
    df_sub <- split_df[[subtype_key]]
    df_sub <- df_sub[order(df_sub$Cohort), , drop = FALSE]
    events <- df_sub$nResponse
    totals <- df_sub$size
    prop_adj <- (events + 0.5)/(totals + 1)
    yi <- log(prop_adj/(1 - prop_adj))
    vi <- 1/(totals * prop_adj * (1 - prop_adj))
    ci_bounds <- t(vapply(
      seq_len(nrow(df_sub)),
      function(i){
        bt <- stats::binom.test(events[i], totals[i], conf.level = conf.level)
        c(bt$conf.int[1], bt$conf.int[2])
      },
      FUN.VALUE = numeric(2)
    ))
    df_sub$ci_lower <- ci_bounds[, 1]
    df_sub$ci_upper <- ci_bounds[, 2]
    df_sub$meta_y <- yi
    df_sub$meta_v <- vi
    df_sub$method <- method_label

    k <- nrow(df_sub)
    hetero_valid <- k >= max(min_studies, 2)
    se_mu <- NA_real_
    mu_hat <- NA_real_
    ci_lower_mu <- NA_real_
    ci_upper_mu <- NA_real_
    tau2 <- 0
    I2 <- 0
    p_het <- NA_real_
    level_pct <- conf.level * 100

    if(k == 1){
      tau2 <- 0
      se_mu <- sqrt(vi)
      mu_hat <- yi
      z_value <- stats::qnorm(0.5 + conf.level/2)
      ci_lower_mu <- mu_hat - z_value * se_mu
      ci_upper_mu <- mu_hat + z_value * se_mu
    } else {
      meta_fit <- tryCatch(
        metafor::rma.uni(
          yi = yi,
          vi = vi,
          method = rma_method,
          test = if(use_hk) "knha" else "z",
          level = level_pct
        ),
        error = function(e) NULL
      )

      if(!is.null(meta_fit)){
        mu_hat <- as.numeric(meta_fit$b)
        se_mu <- as.numeric(meta_fit$se)
        ci_lower_mu <- meta_fit$ci.lb
        ci_upper_mu <- meta_fit$ci.ub
        tau2 <- ifelse(is.na(meta_fit$tau2), 0, meta_fit$tau2)
        if(hetero_valid){
          I2 <- ifelse(is.null(meta_fit$I2), NA_real_, meta_fit$I2)
          p_het <- meta_fit$QEp
        }
      } else {
        tau2 <- estimate_tau2_dl(yi, vi)
        w_star_tmp <- 1/(vi + tau2)
        mu_hat <- sum(w_star_tmp * yi)/sum(w_star_tmp)
        se_mu <- sqrt(1/sum(w_star_tmp))
        z_value <- stats::qnorm(0.5 + conf.level/2)
        ci_lower_mu <- mu_hat - z_value * se_mu
        ci_upper_mu <- mu_hat + z_value * se_mu
        if(hetero_valid){
          wi_fixed <- 1/vi
          mu_fixed <- sum(wi_fixed * yi)/sum(wi_fixed)
          Q <- sum(wi_fixed * (yi - mu_fixed)^2)
          if(Q > 0){
            I2 <- max(0, (Q - (k - 1))/Q) * 100
            p_het <- stats::pchisq(Q, df = k - 1, lower.tail = FALSE)
          } else {
            I2 <- 0
            p_het <- NA_real_
          }
        }
      }
    }

    tau2 <- max(tau2, 0)
    w_star <- 1/(vi + tau2)
    weights_norm <- w_star/sum(w_star)
    df_sub$weight <- weights_norm
    df_sub$effect_type <- 'Cohort'
    df_sub$tau2 <- NA_real_
    df_sub$I2 <- NA_real_
    df_sub$p_heterogeneity <- NA_real_
    df_sub$nStudy <- nrow(df_sub)

    summary_rate <- inv_logit(mu_hat)
    summary_lower <- inv_logit(ci_lower_mu)
    summary_upper <- inv_logit(ci_upper_mu)
    if(!hetero_valid){
      I2 <- 0
      p_het <- NA_real_
    }

    summary_row <- df_sub[1, , drop = FALSE]
    summary_row$Cohort <- summary_label
    summary_row$tumor_type <- NA_character_
    summary_row$size <- sum(totals, na.rm = TRUE)
    summary_row$nResponse <- sum(events, na.rm = TRUE)
    summary_row$response_rate <- summary_rate
    summary_row$ci_lower <- summary_lower
    summary_row$ci_upper <- summary_upper
    summary_row$meta_y <- mu_hat
    summary_row$meta_v <- se_mu^2
    summary_row$weight <- 1
    summary_row$effect_type <- 'Summary'
    summary_row$tau2 <- tau2
    summary_row$I2 <- I2
    summary_row$p_heterogeneity <- p_het
    summary_row$nStudy <- k

    rbind(df_sub, summary_row)
  })

  data_meta <- do.call(rbind, data_meta)
  if(is.null(data_meta) || nrow(data_meta) == 0){
    return(list(Data = NULL, Plot = NULL))
  }

  subtype_levels <- unique(data_meta$Subtype)
  suppressWarnings(subtype_numeric <- as.numeric(subtype_levels))
  if(all(!is.na(subtype_numeric))){
    subtype_levels <- subtype_levels[order(subtype_numeric)]
  }
  data_meta$Subtype <- factor(data_meta$Subtype, levels = subtype_levels)

  cohorts_present <- unique(df$Cohort)
  if(!is.null(cohort_levels)){
    cohort_order <- cohort_levels[cohort_levels %in% cohorts_present]
  } else {
    cohort_order <- cohorts_present
  }
  if(length(cohort_order) == 0){
    cohort_order <- cohorts_present
  }
  cohort_levels_base <- c(summary_label, cohort_order)
  cohort_levels_plot <- rev(cohort_levels_base)

  data_meta$effect_type <- factor(data_meta$effect_type, levels = c('Cohort','Summary'))
  data_meta$Cohort_display <- ifelse(
    data_meta$effect_type == "Summary",
    summary_label,
    data_meta$Cohort
  )
  data_meta$Cohort_display <- factor(data_meta$Cohort_display, levels = cohort_levels_plot)
  cohort_label_values <- setNames(rev(cohort_levels_base), cohort_levels_plot)

  color_palette <- c(Cohort = "#1F77B4", Summary = "#D62728")
  size_palette <- c(Cohort = 2.5, Summary = 3.5)
  shape_palette <- c(Cohort = 21, Summary = 23)

  axis_label <- function(x){
    paste0(round(x * 100), "%")
  }

  build_panel <- function(df_sub, subtype_label, show_legend = FALSE, show_axis_labels = FALSE){
    df_sub$Cohort_display <- factor(df_sub$Cohort_display, levels = cohort_levels_plot)
    df_summary <- subset(df_sub, effect_type == "Summary")
    summary_text <- NULL
    if(nrow(df_summary) > 0){
      df_summary$label <- sprintf(
        "%.1f%% [%.1f%%, %.1f%%]\nI^2 = %s\np = %s",
        df_summary$response_rate * 100,
        df_summary$ci_lower * 100,
        df_summary$ci_upper * 100,
        ifelse(is.na(df_summary$I2), "NA", sprintf("%.1f%%", df_summary$I2)),
        vapply(df_summary$p_heterogeneity, format_p, character(1))
      )
      summary_text <- df_summary$label[1]
    }

    axis_labels <- if(show_axis_labels) cohort_label_values else rep("", length(cohort_levels_plot))
    axis_text <- if(show_axis_labels) element_text(size = cohort_font_size, color = "#333333") else element_blank()
    axis_ticks <- if(show_axis_labels) element_line(color = "#D0D0D0", linewidth = 0.3) else element_blank()

    g <- ggplot(
      df_sub,
      aes(x = response_rate, y = Cohort_display)
    ) +
      geom_errorbarh(
        aes(
          xmin = ci_lower,
          xmax = ci_upper,
          color = effect_type
        ),
        height = 0.2,
        linewidth = 0.7,
        alpha = 0.9
      ) +
      geom_point(
        aes(
          color = effect_type,
          fill = effect_type,
          shape = effect_type,
          size = effect_type
        ),
        stroke = 0.4
      ) +
      scale_color_manual(values = color_palette) +
      scale_fill_manual(values = color_palette) +
      scale_shape_manual(values = shape_palette) +
      scale_size_manual(values = size_palette) +
      scale_x_continuous(
        labels = axis_label,
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.2))
      ) +
      scale_y_discrete(
        limits = cohort_levels_plot,
        labels = axis_labels
      ) +
      labs(
        title = subtype_label,
        subtitle = summary_text,
        x = NULL,
        y = NULL,
        color = NULL,
        shape = NULL,
        fill = NULL,
        size = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = axis_text,
        axis.ticks.y = axis_ticks,
        axis.title.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(face = "bold", hjust = 0.5, size = summary_font_size, margin = margin(b = 6)),
        plot.margin = margin(t = 25, r = 20, b = 25, l = if(show_axis_labels) 10 else 5)
      ) +
      coord_cartesian(clip = "off")

    if(!show_legend){
      g <- g +
        guides(
          color = "none",
          fill = "none",
          shape = "none",
          size = "none"
        )
    }

    return(g)
  }

  subtype_order <- levels(data_meta$Subtype)
  panel_list <- lapply(seq_along(subtype_order), function(i){
    subtype_i <- subtype_order[i]
    df_sub <- subset(data_meta, Subtype == subtype_i)
    build_panel(
      df_sub,
      subtype_label = subtype_i,
      show_legend = (i == 1),
      show_axis_labels = (i == 1)
    )
  })

  combined_panels <- patchwork::wrap_plots(panel_list, nrow = 1, guides = "collect")
  combined_panels <- combined_panels &
    theme(legend.position = "bottom")

  x_label_plot <- ggplot() +
    annotate(
      geom = "text",
      x = 0.5,
      y = 0.5,
      label = "Response rate",
      fontface = "bold",
      size = 6,
      color = "#333333"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_void() +
    theme(plot.margin = margin(t = -5, r = 0, b = 10, l = 0))

  plot_m <- combined_panels / x_label_plot
  plot_m <- plot_m +
    patchwork::plot_layout(heights = c(1, 0.08))
  plot_m <- plot_m +
    patchwork::plot_annotation(
      title = "Meta-analysis of subtype response rates",
      subtitle = paste0("Random-effects (", method_label, ")")
    )

  return(list(
    Data = data_meta,
    Plot = plot_m
  ))
}


####%%%%%%%%%%%%% Assistant functions %%%%%%%%%%%%%%%%%%%%####


# For subtypeEffect
calculate_bounded_confidence_interval <- function(data, conf_level = 0.95, lower_limit = 0, upper_limit = 1) {
  # 计算样本均值
  sample_mean <- mean(data, na.rm = T)

  # 计算样本标准误差
  standard_error <- sd(data, na.rm = T) / sqrt(length(data))

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




# Shared candidate x cohort scoring for both normCCS search algorithms.
norm_validate_objective <- function(objective) {
  choices <- c("median", "cohort_median", "mean", "lower_tail", "mean_lower_tail")
  if (!is.character(objective) || length(objective) != 1L ||
      is.na(objective) || !objective %in% choices) {
    stop("norm_objective must be one of: ", paste(choices, collapse = ", "), ".")
  }
  objective
}

norm_objective_score <- function(auc, cohorts, objective = "median") {
  # Test
  if (FALSE) {
    auc <- matrix(c(0.9, 0.5, 0.7, 0.8, 0.9, 1), nrow = 1L)
    norm_objective_score(auc, c("All", LETTERS[1:5]), "lower_tail")
  }
  objective <- norm_validate_objective(objective)
  if (!is.matrix(auc) || !is.numeric(auc) || ncol(auc) != length(cohorts) ||
      !is.character(cohorts) || anyNA(cohorts) || any(!nzchar(cohorts)) ||
      anyDuplicated(cohorts)) {
    stop("auc must be a numeric candidate-by-cohort matrix with unique, nonempty cohort identifiers.")
  }
  if (objective == "median") {
    return(apply(auc, 1L, stats::median, na.rm = TRUE))
  }
  real <- auc[, cohorts != "All", drop = FALSE]
  if (!ncol(real)) stop("norm_objective requires at least one real cohort (excluding 'All').")
  if (any(!is.finite(real))) {
    bad <- which(!is.finite(real), arr.ind = TRUE)[1L, ]
    stop("Non-finite AUC at candidate ", bad[1L], ", cohort '",
         cohorts[cohorts != "All"][bad[2L]], "'.")
  }
  if (objective == "mean") return(rowMeans(real))
  if (objective == "cohort_median") return(apply(real, 1L, stats::median))
  mass <- 0.25 * ncol(real)
  whole <- floor(mass)
  fraction <- mass - whole
  tail <- apply(real, 1L, function(x) {
    x <- sort(x)
    (sum(x[seq_len(whole)]) + if (fraction > 0) fraction * x[whole + 1L] else 0) / mass
  })
  if (objective == "lower_tail") return(tail)
  0.5 * rowMeans(real) + 0.5 * tail
}

norm_objective_settings <- function(objective) {
  settings <- list(norm_objective = objective)
  if (objective %in% c("lower_tail", "mean_lower_tail")) settings$tail_fraction <- 0.25
  if (objective == "mean_lower_tail") settings$mean_weight <- 0.5
  settings
}

# Preserve legacy columns; never relabel a new score as medianROCAUC.
norm_objective_summary <- function(raw, objective) {
  infos <- sort(unique(as.character(raw$Info)))
  cohorts <- sort(unique(as.character(raw$Cohort)))
  key <- paste(match(raw$Info, infos), match(raw$Cohort, cohorts), sep = ":")
  if (anyDuplicated(key) || nrow(raw) != length(infos) * length(cohorts)) {
    stop("Each candidate must have exactly one AUC for every cohort.")
  }
  auc <- matrix(NA_real_, length(infos), length(cohorts))
  auc[cbind(match(raw$Info, infos), match(raw$Cohort, cohorts))] <- raw$ROCAUC
  summary <- data.frame(Info = infos,
                        medianROCAUC = norm_objective_score(auc, cohorts, "median"),
                        stringsAsFactors = FALSE)
  summary$objective <- objective
  summary$objectiveScore <- norm_objective_score(auc, cohorts, objective)
  summary$cohortMedianROCAUC <- norm_objective_score(auc, cohorts, "cohort_median")
  summary$cohortMeanROCAUC <- norm_objective_score(auc, cohorts, "mean")
  summary$cohortLowerTailROCAUC <- norm_objective_score(auc, cohorts, "lower_tail")
  summary <- summary[order(-summary$objectiveScore, -summary$cohortMeanROCAUC, summary$Info), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}


# Internal normCCS search kernels. States are ordered from high to low RR/NRR.
# No RNG, persistent cache, native code, or parallel backend is used here.

norm_validate_input <- function(df, k, method, cores) {
  required <- c("Cohort", "Subtype", "response", "SampleIDs", "tumor_type")
  if (!is.data.frame(df) || !all(required %in% names(df)) || !nrow(df)) {
    stop("df must be a nonempty data frame with Cohort, Subtype, response, SampleIDs and tumor_type.")
  }
  if (anyNA(df[required]) || !all(as.character(df$response) %in% c("0", "1"))) {
    stop("Required columns must be complete; response must contain only 0 and 1.")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2 || k != floor(k)) {
    stop("n_norm_subtype must be a single integer >= 2.")
  }
  if (length(method) != 1L || !method %in% c("Raw", "OptimizeRank")) {
    stop("norm_roc_method must be 'Raw' or 'OptimizeRank'.")
  }
  if (!is.null(cores) && (!is.numeric(cores) || length(cores) != 1L ||
      !is.finite(cores) || cores < 1 || cores != floor(cores))) {
    stop("numCores must be NULL or a positive integer.")
  }
  if (any(vapply(split(df$response, df$Cohort, drop = TRUE),
                 function(x) length(unique(x)) != 2L, logical(1)))) {
    stop("Every cohort must contain both response classes (0 and 1).")
  }
  if (any(as.character(df$Cohort) == "All")) {
    stop("Cohort name 'All' is reserved for the pooled summary.")
  }
}

# Prefix counts and interval loss relative to keeping all M distinct scores.
# L[a,b] = U[a,b] - P[a,b]*N[a,b]/2; negative loss is legal.
norm_prepare_counts <- function(df, ranked) {
  m <- nrow(ranked)
  cohorts <- sort(unique(as.character(df$Cohort)))
  state <- match(as.character(df$Subtype), as.character(ranked$Subtype))
  rows <- c(list(seq_len(nrow(df))), lapply(cohorts, function(x) which(df$Cohort == x)))
  names(rows) <- c("All", cohorts)
  # pROC respects factor level order, even when '0' is the case level.
  roc_levels <- levels(droplevels(as.factor(df$response)))
  p <- n <- matrix(0, length(rows), m)
  loss <- first <- array(0, c(length(rows), m, m))
  median_states <- matrix(0L, length(rows), 4L)
  for (c in seq_along(rows)) {
    ids <- rows[[c]]
    positive <- as.character(df$response[ids]) == roc_levels[2L]
    p[c, ] <- tabulate(state[ids[positive]], nbins = m)
    n[c, ] <- tabulate(state[ids[!positive]], nbins = m)
    cp <- c(0, cumsum(p[c, ]))
    cn <- c(0, cumsum(n[c, ]))
    pair_prefix <- c(0, cumsum(p[c, ] * (cn[-1L] - n[c, ] / 2)))
    first_state <- match(seq_len(m), state[ids])
    first_state[is.na(first_state)] <- Inf
    for (a in seq_len(m)) {
      b <- a:m
      pg <- cp[b + 1L] - cp[a]
      ng <- cn[b + 1L] - cn[a]
      loss[c, a, b] <- (pg * cn[b + 1L] -
        (pair_prefix[b + 1L] - pair_prefix[a]) - pg * ng / 2) / (sum(p[c, ]) * sum(n[c, ]))
      first[c, a, b] <- cummin(first_state[b])
    }
    for (j in 1:2) {
      counts <- if (j == 1L) p[c, ] else n[c, ]
      total <- sum(counts)
      median_states[c, (2L * j - 1L):(2L * j)] <-
        vapply(c(ceiling(total / 2), floor(total / 2) + 1),
               function(x) which(cumsum(counts) >= x)[1L], integer(1))
    }
  }
  prefix <- function(x) t(apply(cbind(0, x), 1, cumsum))
  base_auc <- rowSums(p * (rowSums(n) - t(apply(n, 1, cumsum)) + n / 2)) /
    (rowSums(p) * rowSums(n))
  list(m = m, cohorts = names(rows), rows = rows, p = p, n = n,
       roc_levels = roc_levels,
       cp = prefix(p), cn = prefix(n), loss = loss, first = first,
       median_states = median_states, base_auc = base_auc)
}

# Exact O(K*M^2) solution for a weighted fixed-direction mean AUC.
norm_solve_dp <- function(cost, k) {
  m <- nrow(cost)
  value <- matrix(Inf, k, m)
  parent <- matrix(0L, k, m)
  value[1L, ] <- cost[1L, ]
  for (g in 2:k) {
    for (b in g:m) {
      ends <- (g - 1L):(b - 1L)
      choices <- value[g - 1L, ends] + cost[cbind(ends + 1L, b)]
      best <- which.min(choices)
      value[g, b] <- choices[best]
      parent[g, b] <- ends[best]
    }
  }
  cuts <- integer(k - 1L)
  b <- m
  for (g in k:2) {
    cuts[g - 1L] <- parent[g, b]
    b <- cuts[g - 1L]
  }
  list(cuts = cuts, loss = value[k, m])
}

norm_median_bin <- function(counts) {
  cumulative <- counts
  if (ncol(counts) > 1L) for (j in 2:ncol(counts)) {
    cumulative[, j] <- cumulative[, j - 1L] + counts[, j]
  }
  total <- rowSums(counts)
  (rowSums(cumulative < ceiling(total / 2)) +
     rowSums(cumulative < (floor(total / 2) + 1))) / 2 + 1
}

# Batched legacy AUC, including pROC's median-based automatic direction.
# OptimizeRank preserves first-observation tie order in one-vs-rest AUCs.
norm_auc_batch <- function(cuts, model, method = "Raw") {
  if (is.null(dim(cuts))) cuts <- matrix(cuts, nrow = 1L)
  nr <- nrow(cuts)
  k <- ncol(cuts) + 1L
  m <- model$m
  left <- cbind(1L, cuts + 1L)
  right <- cbind(cuts, m)
  interval <- left + (right - 1L) * m
  result <- matrix(0, nr, length(model$cohorts), dimnames = list(NULL, model$cohorts))
  for (c in seq_along(model$cohorts)) {
    if (method == "Raw" || k == 2L) {
      auc <- model$base_auc[c] - rowSums(matrix(model$loss[c, , ][as.vector(interval)], nr, k))
      median_group <- matrix(vapply(model$median_states[c, ],
                                    function(s) rowSums(cuts < s), numeric(nr)), nr, 4L)
      reverse <- rowSums(median_group[, 3:4, drop = FALSE]) <
        rowSums(median_group[, 1:2, drop = FALSE])
    } else {
      p <- matrix(model$cp[c, right + 1L] - model$cp[c, left], nr, k)
      n <- matrix(model$cn[c, right + 1L] - model$cn[c, left], nr, k)
      pt <- rowSums(p)
      nt <- rowSums(n)
      # Match pROC's trapezoid arithmetic, not just its algebraic value:
      # last-bit differences change OptimizeRank's ordering of nearly tied AUCs.
      sensitivity <- p / pt
      specificity <- (nt - n) / nt
      ovr <- ((1 + sensitivity) / 2) * specificity +
        (sensitivity / 2) * (1 - specificity)
      # A binary indicator has median 0, 0.5 or 1 depending on its count.
      binary_median <- function(x, total) (x >= total / 2) / 2 + (x > total / 2) / 2
      flip <- binary_median(n, nt) > binary_median(p, pt)
      reverse_sensitivity <- (pt - p) / pt
      reverse_specificity <- n / nt
      reverse_auc <- ((1 + reverse_sensitivity) / 2) * reverse_specificity +
        (reverse_sensitivity / 2) * (1 - reverse_specificity)
      ovr[flip] <- reverse_auc[flip]
      first <- matrix(model$first[c, , ][as.vector(interval)], nr, k)
      # Absent groups never occur in unique(predictor) in the reference.
      ovr[p + n == 0] <- Inf
      ordering <- order(rep(seq_len(nr), k), as.vector(ovr), as.vector(first))
      p <- matrix(p[ordering], nr, k, byrow = TRUE)
      n <- matrix(n[ordering], nr, k, byrow = TRUE)
      cumulative_n <- n
      if (k > 1L) for (j in 2:k) cumulative_n[, j] <- cumulative_n[, j - 1L] + n[, j]
      auc <- rowSums(p * (cumulative_n - n / 2)) / (pt * nt)
      reverse <- norm_median_bin(n) > norm_median_bin(p)
    }
    result[, c] <- ifelse(reverse, 1 - auc, auc)
  }
  result
}

norm_info <- function(cuts, ranked) {
  k <- length(cuts) + 1L
  groups <- k - 1L - findInterval(seq_len(nrow(ranked)) - 1L, cuts)
  paste(vapply(0:(k - 1L), function(g) {
    paste0(g, "=", paste(ranked$Subtype[groups == g], collapse = ","))
  }, character(1)), collapse = "; ")
}

subtype_norm_dp <- function(df, ranked, k, method, started, verbose, objective = "median") {
  # Test
  if (FALSE) {
    # The synthetic generator and complete runnable examples live in
    # test/normCCS-5-more/; no patient data are needed.
    model <- norm_prepare_counts(df, ranked)
    norm_solve_dp(model$loss[1L, , ], k)
  }
  objective <- norm_validate_objective(objective)
  if (verbose) luckyBase::LuckyVerbose(paste0("subtypeNorm: DP starts and ", objective, " refinement..."), levels = 2)
  model <- norm_prepare_counts(df, ranked)
  nc <- length(model$cohorts)
  weights <- rbind(rep(1 / nc, nc), c(0, rep(1 / (nc - 1L), nc - 1L)), diag(nc))
  # Also seed reverse directions because legacy pROC auto-direction can flip.
  weights <- rbind(weights, -weights)
  loss_matrix <- matrix(model$loss, nrow = nc)
  starts <- lapply(seq_len(nrow(weights)), function(i) {
    cost <- matrix(as.vector(crossprod(weights[i, ], loss_matrix)), model$m, model$m)
    norm_solve_dp(cost, k)$cuts
  })
  starts <- unique(c(starts, list(as.integer(floor(seq_len(k - 1L) * model$m / k)))))
  score <- function(cuts) norm_objective_score(norm_auc_batch(cuts, model, method), model$cohorts, objective)
  candidate_order <- function(values, cuts) {
    if (objective == "median") return(order(-values, seq_along(values)))
    auc <- norm_auc_batch(cuts, model, method)
    means <- norm_objective_score(auc, model$cohorts, "mean")
    infos <- apply(cuts, 1L, norm_info, ranked = ranked)
    order(-values, -means, infos)
  }
  proposal_winner <- function(values, proposals) {
    if (objective == "median") return(which.max(values))
    ties <- which(values == max(values))
    if (length(ties) == 1L) return(ties)
    ties[candidate_order(values[ties], proposals[ties, , drop = FALSE])[1L]]
  }
  start_scores <- score(do.call(rbind, starts))
  selected <- candidate_order(start_scores, do.call(rbind, starts))[seq_len(min(6L, length(starts)))]
  candidates <- starts[selected]
  evaluations <- length(starts)
  refinements <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    cuts <- starts[[selected[i]]]
    best <- start_scores[selected[i]]
    converged <- FALSE
    pair_moves <- 0L
    for (pass in seq_len(10L)) {
      improved <- FALSE
      for (j in seq_along(cuts)) {
        lower <- if (j == 1L) 1L else cuts[j - 1L] + 1L
        upper <- if (j == length(cuts)) model$m - 1L else cuts[j + 1L] - 1L
        proposals <- matrix(rep(cuts, each = upper - lower + 1L), ncol = length(cuts))
        proposals[, j] <- lower:upper
        values <- score(proposals)
        evaluations <- evaluations + length(values)
        winner <- proposal_winner(values, proposals)
        if (values[winner] > best + 1e-12) {
          cuts <- as.integer(proposals[winner, ])
          best <- values[winner]
          improved <- TRUE
        }
      }
      # A single-boundary optimum can still improve by moving two neighboring
      # boundaries together. Search that neighborhood only after singles stall.
      if (!improved && length(cuts) > 1L) {
        for (j in seq_len(length(cuts) - 1L)) {
          lower <- if (j == 1L) 1L else cuts[j - 1L] + 1L
          upper <- if (j + 1L == length(cuts)) model$m - 1L else cuts[j + 2L] - 1L
          if (upper - lower < 2L) next
          first <- rep(lower:(upper - 1L), times = (upper - lower):1L)
          second <- unlist(lapply(lower:(upper - 1L), function(a) (a + 1L):upper), use.names = FALSE)
          pair_best <- best
          pair_cuts <- cuts
          for (offset in seq.int(1L, length(first), by = 512L)) {
            ids <- offset:min(offset + 511L, length(first))
            proposals <- matrix(rep(cuts, each = length(ids)), ncol = length(cuts))
            proposals[, j] <- first[ids]
            proposals[, j + 1L] <- second[ids]
            values <- score(proposals)
            evaluations <- evaluations + length(values)
            winner <- proposal_winner(values, proposals)
            if (values[winner] > pair_best + 1e-12) {
              pair_cuts <- as.integer(proposals[winner, ])
              pair_best <- values[winner]
            } else if (objective != "median" && values[winner] == pair_best &&
                       pair_best > best + 1e-12) {
              tied <- rbind(pair_cuts, proposals[winner, ])
              preferred <- proposal_winner(c(pair_best, values[winner]), tied)
              pair_cuts <- as.integer(tied[preferred, ])
            }
          }
          if (pair_best > best + 1e-12) {
            cuts <- pair_cuts
            best <- pair_best
            improved <- TRUE
            pair_moves <- pair_moves + 1L
          }
        }
      }
      if (!improved) {
        converged <- TRUE
        break
      }
    }
    candidates[[length(candidates) + 1L]] <- cuts
    refinements[[i]] <- list(initial_auc = start_scores[selected[i]], final_auc = best,
                             passes = pass, converged = converged, pair_moves = pair_moves)
  }
  candidates <- unique(candidates)
  infos <- vapply(candidates, norm_info, character(1), ranked = ranked)
  raw <- lapply(seq_along(candidates), function(i) {
    cuts <- candidates[[i]]
    group <- k - 1L - findInterval(match(as.character(df$Subtype), as.character(ranked$Subtype)) - 1L, cuts)
    do.call(rbind, lapply(seq_along(model$rows), function(c) {
      ids <- model$rows[[c]]
      data.frame(Cohort = model$cohorts[c], Info = infos[i], nSample = length(ids),
                 get_perform_markers(df$response[ids], group[ids], k, method),
                 stringsAsFactors = FALSE)
    }))
  })
  raw <- do.call(rbind, raw)
  final_scores <- vapply(split(raw$ROCAUC, raw$Info), stats::median, numeric(1), na.rm = TRUE)
  fast_auc <- norm_auc_batch(do.call(rbind, candidates), model, method)
  reference_auc <- matrix(raw$ROCAUC, nrow = length(candidates), byrow = TRUE)
  if (any(!is.finite(fast_auc)) || any(!is.finite(reference_auc)) ||
      any(abs(reference_auc - fast_auc) > 1e-10)) {
    stop("Fast normCCS scoring disagrees with the reference ROC evaluator; results were not returned.")
  }
  if (any(abs(norm_objective_score(fast_auc, model$cohorts, objective) -
              norm_objective_score(reference_auc, model$cohorts, objective)) > 1e-10)) {
    stop("Fast normCCS objective disagrees with reference ROC objective; results were not returned.")
  }
  summary <- data.frame(Info = names(final_scores), medianROCAUC = unname(final_scores), stringsAsFactors = FALSE)
  summary <- summary[order(-summary$medianROCAUC, summary$Info), , drop = FALSE]
  rownames(summary) <- NULL
  if (objective != "median") summary <- norm_objective_summary(raw, objective)
  winner <- match(summary$Info[1L], infos)
  result <- list(Raw = raw, Summary = summary, SearchDiagnostics = list(
    method = "dp", objective = "legacy_median_auc_including_All",
    roc_method = method, refinement_start_limit = 6L, refinement_pass_limit = 10L,
    refinement_neighborhood = "single_and_adjacent_pair_boundaries",
    refinement_pair_batch_size = 512L,
    surrogate = "fixed_direction_weighted_mean_raw_auc",
    exact_global_optimum = FALSE, dp_subproblems_exact = TRUE,
    M = model$m, K = k, n_dp_starts = nrow(weights), n_unique_starts = length(starts),
    n_unique_candidates = length(candidates), n_scored_proposals = evaluations,
    n_reference_verifications = length(reference_auc),
    selected_cutpoints = candidates[[winner]], ranked_subtypes = as.character(ranked$Subtype),
    best_initial_median_auc = max(start_scores), refinements = refinements,
    cohort_only_median_auc = stats::median(raw$ROCAUC[raw$Info == summary$Info[1L] & raw$Cohort != "All"]),
    workers = 1L, interval_storage_bytes = as.numeric(utils::object.size(model)),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started)))
  if (objective != "median") {
    diagnostics <- result$SearchDiagnostics
    diagnostics$objective <- objective
    diagnostics$objective_scope <- "cohorts_only"
    diagnostics$n_real_cohorts <- nc - 1L
    diagnostics$best_initial_median_auc <- NULL
    diagnostics$best_initial_objective_score <- max(start_scores)
    diagnostics$best_final_objective_score <- summary$objectiveScore[1L]
    diagnostics$refinements <- lapply(refinements, function(x) {
      names(x)[1:2] <- c("initial_objective_score", "final_objective_score")
      x
    })
    settings <- norm_objective_settings(objective)
    diagnostics <- c(diagnostics, settings[setdiff(names(settings), "norm_objective")])
    result$SearchDiagnostics <- diagnostics
  }
  result
}


# Internal state ordering and inherited normCCS search. No files or backends
# are changed; the continuous engine is intentionally single-core.
norm_validate_mode <- function(mode) {
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% c("discrete", "continuous")) {
    stop("norm_mode must be 'discrete' or 'continuous'.")
  }
  mode
}

norm_validate_top_n <- function(top_n) {
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n) || top_n > .Machine$integer.max) {
    stop("norm_top_n must be a positive integer no larger than .Machine$integer.max.")
  }
  invisible(top_n)
}

norm_validate_rank_marker <- function(marker) {
  choices <- c("NRR", "RR", "NRR+RR", "RR+NRR", "MeanNRR", "LowerTailNRR",
               "ShrunkLogOR", "RankFusion")
  if (!is.character(marker) || length(marker) != 1L || is.na(marker) ||
      !marker %in% choices) {
    stop("norm_rank_marker must be one of: ", paste(choices, collapse = ", "), ".")
  }
  invisible(marker)
}

# Return the ordered state table. Supplying the legacy summary preserves its
# original factor/tie order exactly; added methods use deterministic tie breaks.
norm_rank_subtypes <- function(df, norm_rank_marker = "NRR", summary = NULL) {
  # Test
  if (FALSE) {
    df <- expand.grid(Cohort = c("A", "B"), Subtype = 1:4, response = 0:1)
    df$SampleIDs <- seq_len(nrow(df))
    df$tumor_type <- "Lung"
    norm_rank_subtypes(df, "ShrunkLogOR")
  }
  norm_validate_rank_marker(norm_rank_marker)
  legacy <- norm_rank_marker %in% c("NRR", "RR", "NRR+RR", "RR+NRR")
  if (is.null(summary) || !legacy) {
    counts <- plyr::ddply(df, c("Cohort", "Subtype"), function(x) {
      data.frame(size = nrow(x), nResponse = sum(as.character(x$response) == "1"))
    })
    counts$response_rate <- counts$nResponse / counts$size
    totals <- plyr::ddply(counts, "Cohort", function(x) {
      data.frame(total = sum(x$size), positives = sum(x$nResponse))
    })
    index <- match(counts$Cohort, totals$Cohort)
    baseline <- totals$positives[index] / totals$total[index]
    counts$normalized_response_rate <- (counts$response_rate - baseline) / baseline
    if (is.null(summary)) {
      summary <- plyr::ddply(counts, "Subtype", function(x) {
        data.frame(medianNRR = stats::median(x$normalized_response_rate),
                   medianRR = stats::median(x$response_rate))
      })
    }
  }
  if (legacy) {
    if (norm_rank_marker == "NRR") return(dplyr::arrange(summary, dplyr::desc(medianNRR)))
    if (norm_rank_marker == "RR") return(dplyr::arrange(summary, dplyr::desc(medianRR)))
    if (norm_rank_marker == "NRR+RR") {
      return(dplyr::arrange(summary, dplyr::desc(medianNRR), dplyr::desc(medianRR)))
    }
    return(dplyr::arrange(summary, dplyr::desc(medianRR), dplyr::desc(medianNRR)))
  }
  if (norm_rank_marker == "RankFusion") {
    score <- -(rank(-summary$medianNRR, ties.method = "average") +
                 rank(-summary$medianRR, ties.method = "average")) / 2
  } else {
    counts$centered_log_odds <- log((counts$nResponse + 0.5) /
      (counts$size - counts$nResponse + 0.5)) -
      log((totals$positives[index] + 0.5) /
        (totals$total[index] - totals$positives[index] + 0.5))
    scores <- plyr::ddply(counts, "Subtype", function(x) {
      value <- switch(norm_rank_marker,
        MeanNRR = mean(x$normalized_response_rate),
        LowerTailNRR = {
          ordered <- sort(x$normalized_response_rate)
          mass <- length(ordered) * 0.25
          whole <- floor(mass)
          fraction <- mass - whole
          (sum(ordered[seq_len(whole)]) +
             if (fraction > 0) fraction * ordered[whole + 1L] else 0) / mass
        },
        ShrunkLogOR = mean(x$centered_log_odds))
      data.frame(rank_score = value)
    })
    score <- scores$rank_score[match(summary$Subtype, scores$Subtype)]
  }
  summary$rankScore <- score
  summary <- summary[order(-score, -summary$medianNRR, -summary$medianRR,
                           as.character(summary$Subtype)), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

# Bellman recurrence across mandatory parent blocks: allocate exactly one extra
# grade. Each block may stay intact or split once. The fixed-direction additive
# AUC surrogate is exact here; the public nonadditive objective is scored later.
norm_inherited_dp <- function(cost, parent_cuts) {
  edges <- c(0L, parent_cuts, nrow(cost))
  value <- c(0, Inf)
  selected <- NA_integer_
  for (j in seq_len(length(edges) - 1L)) {
    left <- edges[j] + 1L
    right <- edges[j + 1L]
    intact <- cost[left, right]
    next_value <- c(value[1L] + intact, value[2L] + intact)
    if (left < right) {
      splits <- seq.int(left, right - 1L)
      loss <- cost[cbind(left, splits)] + cost[cbind(splits + 1L, right)]
      best <- which.min(loss)
      alternative <- value[1L] + loss[best]
      if (alternative < next_value[2L]) {
        next_value[2L] <- alternative
        selected <- splits[best]
      }
    }
    value <- next_value
  }
  if (is.na(selected)) stop("No legal inherited split remains.")
  list(cuts = sort(c(parent_cuts, selected)), loss = value[2L])
}

# Original ROC evaluator, including legacy binary performance columns.
norm_reference_candidates <- function(df, ranked, candidates, model, method) {
  rows <- lapply(candidates, function(cuts) {
    k <- length(cuts) + 1L
    group <- k - 1L - findInterval(
      match(as.character(df$Subtype), as.character(ranked$Subtype)) - 1L, cuts)
    info <- norm_info(cuts, ranked)
    do.call(rbind, lapply(seq_along(model$rows), function(c) {
      ids <- model$rows[[c]]
      metrics <- if (k == 1L) {
        data.frame(ROCAUC = as.numeric(pROC::roc(df$response[ids], group[ids],
                                               quiet = TRUE)$auc))
      } else {
        get_perform_markers(df$response[ids], group[ids], k, method)
      }
      data.frame(Cohort = model$cohorts[c], Info = info, nSample = length(ids),
                 metrics, stringsAsFactors = FALSE)
    }))
  })
  plyr::rbind.fill(rows)
}

subtypeNorm_continuous <- function(
    df, n_norm_subtype = 2, norm_roc_method = "Raw", norm_rank_marker = "RR",
    numCores = NULL, verbose = TRUE,
    norm_search_method = c("bruteforce", "dp"), norm_objective = "median",
    norm_top_n = 10L
) {
  # Test
  if (FALSE) {
    df <- expand.grid(Cohort = c("A", "B"), Subtype = 1:5, response = 0:1)
    df$SampleIDs <- paste0("synthetic_", seq_len(nrow(df)))
    df$tumor_type <- "Lung"
    subtypeNorm_continuous(df, 4, norm_top_n = 3, verbose = FALSE)
  }
  objective <- norm_validate_objective(norm_objective)
  norm_validate_top_n(norm_top_n)
  norm_validate_rank_marker(norm_rank_marker)
  norm_search_method <- match.arg(norm_search_method)
  if (!is.numeric(n_norm_subtype) || length(n_norm_subtype) != 1L ||
      !is.finite(n_norm_subtype) || n_norm_subtype < 1 ||
      n_norm_subtype != floor(n_norm_subtype)) {
    stop("n_norm_subtype must be a single integer >= 1 in continuous mode.")
  }
  norm_validate_input(df, max(2, n_norm_subtype), norm_roc_method, numCores)
  started <- proc.time()[["elapsed"]]
  ranked <- norm_rank_subtypes(df, norm_rank_marker)
  if (n_norm_subtype > nrow(ranked)) {
    stop("n_norm_subtype must not exceed the number of observed subtypes.")
  }
  # Bound both count-tensor and result growth before allocating either.
  nc <- length(unique(df$Cohort)) + 1L
  interval_bytes <- if (n_norm_subtype > 1L) 64 * nc * nrow(ranked)^2 else 0
  estimated_bytes <- interval_bytes +
    as.double(norm_top_n) * n_norm_subtype * nc * 192
  if (!is.finite(estimated_bytes) || estimated_bytes > 2 * 1024^3) {
    stop("Continuous search exceeds the interval/result-memory budget; reduce norm_top_n or the state count.")
  }
  costs <- list()
  if (n_norm_subtype == 1L) {
    cohorts <- sort(unique(as.character(df$Cohort)))
    rows <- c(list(seq_len(nrow(df))), lapply(cohorts, function(c) which(df$Cohort == c)))
    model <- list(m = nrow(ranked), cohorts = c("All", cohorts), rows = rows)
  } else {
    model <- norm_prepare_counts(df, ranked)
    weights <- rbind(rep(1 / nc, nc), c(0, rep(1 / (nc - 1L), nc - 1L)), diag(nc))
    weights <- rbind(weights, -weights)
    loss <- matrix(model$loss, nrow = nc)
    costs <- lapply(seq_len(nrow(weights)), function(i) {
      matrix(as.vector(crossprod(weights[i, ], loss)), model$m, model$m)
    })
  }
  frontier <- list(integer())
  raw_layers <- layers <- vector("list", n_norm_subtype)
  evaluations <- dp_problems <- reference_checks <- 0L
  for (k in seq_len(n_norm_subtype)) {
    if (verbose) luckyBase::LuckyVerbose(paste0("subtypeNorm_continuous: K = ", k), levels = 2)
    parents <- frontier
    dp_starts <- list()
    if (k == 1L) {
      candidates <- list(integer())
    } else {
      dp_starts <- unlist(lapply(parents, function(cuts) {
        lapply(costs, function(cost) norm_inherited_dp(cost, cuts)$cuts)
      }), recursive = FALSE)
      dp_problems <- dp_problems + length(dp_starts)
      # A parent has M-K+1 legal splits. Complete this small neighborhood so
      # median/OptimizeRank candidates cannot be missed by the DP surrogate.
      extensions <- unlist(lapply(parents, function(cuts) {
        lapply(setdiff(seq_len(model$m - 1L), cuts), function(split) sort(c(cuts, split)))
      }), recursive = FALSE)
      candidates <- unique(c(dp_starts, extensions))
    }
    cut_matrix <- matrix(unlist(candidates, use.names = FALSE),
                         nrow = length(candidates), ncol = k - 1L, byrow = TRUE)
    auc <- if (k == 1L) matrix(0.5, 1L, nc) else norm_auc_batch(cut_matrix, model, norm_roc_method)
    scores <- norm_objective_score(auc, model$cohorts, objective)
    means <- norm_objective_score(auc, model$cohorts, "mean")
    infos <- vapply(candidates, norm_info, character(1), ranked = ranked)
    ordering <- if (objective == "median") order(-scores, infos) else order(-scores, -means, infos)
    retained <- head(ordering, norm_top_n)
    frontier <- candidates[retained]
    reference <- norm_reference_candidates(df, ranked, frontier, model, norm_roc_method)
    actual <- matrix(reference$ROCAUC, nrow = length(retained), ncol = nc, byrow = TRUE)
    if (any(!is.finite(actual)) || any(!is.finite(auc)) ||
        any(abs(actual - auc[retained, , drop = FALSE]) > 1e-10)) {
      stop("Fast continuous normCCS scoring disagrees with the reference ROC evaluator; results were not returned.")
    }
    evaluations <- evaluations + length(candidates)
    reference_checks <- reference_checks + length(actual)
    raw_layers[[k]] <- reference
    layers[[k]] <- list(K = k, method = if (k == 1L) "bruteforce" else "inherited_dp",
      parent_cutpoints = if (k == 1L) list() else parents,
      parent_infos = if (k == 1L) character() else vapply(parents, norm_info, character(1), ranked = ranked),
      dp_start_cutpoints = unique(dp_starts), n_evaluated = length(candidates),
      n_retained = length(frontier), retained_cutpoints = frontier,
      retained_infos = infos[retained], objective_scores = scores[retained])
  }
  raw <- plyr::rbind.fill(raw_layers)
  if (objective == "median") {
    values <- vapply(split(raw$ROCAUC, raw$Info), stats::median, numeric(1), na.rm = TRUE)
    summary <- data.frame(Info = names(values), medianROCAUC = unname(values), stringsAsFactors = FALSE)
    summary <- summary[order(-summary$medianROCAUC, summary$Info), , drop = FALSE]
    rownames(summary) <- NULL
  } else {
    summary <- norm_objective_summary(raw, objective)
  }
  list(Raw = raw, Summary = summary, SearchDiagnostics = list(
    method = "dp", mode = "continuous", rank_marker = norm_rank_marker,
    objective = objective, roc_method = norm_roc_method,
    M = model$m, K = n_norm_subtype, top_n = norm_top_n, layers = layers,
    ranked_subtypes = as.character(ranked$Subtype),
    inheritance = "preserve_all_parent_cutpoints_and_add_one",
    surrogate = "fixed_direction_weighted_mean_raw_auc",
    exact_global_optimum = FALSE, dp_subproblems_exact = TRUE,
    complete_single_split_neighborhood = TRUE,
    n_dp_subproblems = dp_problems, n_scored_proposals = evaluations,
    n_unique_candidates = nrow(summary), n_reference_verifications = reference_checks,
    workers = 1L, interval_storage_bytes = as.numeric(utils::object.size(model)),
    elapsed_seconds = unname(proc.time()[["elapsed"]] - started)))
}
