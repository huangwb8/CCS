


#' @title Downstream analysis: Molecular subtypes & Clinical features
#' @description Downstream analysis: Molecular subtypes & Clinical features
#' @inheritParams subtypeEffect
#' @inheritParams CCSPublicParams
#' @param cohort_level The level of cohort used in scater+bar plot
#' @param color_tissue a named vector like \code{c(Brain = "#984EA3",Oral = "#BF5B17")}
#' @param color_cohort Default is \code{NULL}, which means that all cohort names were colored as black. You can also use a named vector like \code{c('Kim2018_PD-1_gastric' = "#E31A1C")}.
#' @param color_subtype Default is \code{NULL}, which means that all raw subtype names were colored as \code{scales::hue_pal()(n)}. You can also use a named vector like \code{c('0' = "#E31A1C")}.
#' @param n_norm_subtype The number of normalized subtypes you want to explore
#' @param norm_roc_method The method to call ROC AUC. One of \code{'OptimizeRank'} and \code{'Raw'}.
#' @param norm_rank_marker The marker for normalization rank. Interpretation:
#' \itemize{
#'   \item \code{NRR} Ranked only by normalized response rate.
#'   \item \code{RR} Ranked only by raw response rate.
#'   \item \code{NRR+RR} Ranked by normalized response rate and then by raw response rate.
#'   \item \code{RR+NRR} Ranked by raw response rate and then by normalized response rate.
#' }
#' @param plot_layout_widths \code{layout_widths} for \code{\link[patchwork]{plot_layout}}
#' @param plot_ROCCutoff The cut-off value for visualization of ROC AUC
#' @param plot_silence_cohort Whether to silence cohort names in the plot
#' @param plot_silence_subtype Whether to silence subtype names in the plot
#' @importFrom tidyr `%>%`
#' @importFrom dplyr summarize arrange desc filter
#' @importFrom plyr ddply
#' @return A list
#' \itemize{
#'   \item \code{Repeat} Data or parameters for repeatability
#'   \item \code{Plot} 1. Scatter plot + bar plot for RR & NRR; 2. Forest plot for RR.
#'   \item \code{Data} 1. ROC data; 2. Normalized subtype results 3. Clinic utility from \code{\link{subtypeEffect}}.
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
    norm_rank_marker = c('NRR', 'RR', 'NRR+RR', 'RR+NRR')[1],
    plot_ROCCutoff = 0.8,
    plot_layout_widths = c(9, 1),
    plot_silence_cohort = F,
    plot_silence_subtype = F,
    numCores = NULL,
    verbose = TRUE
){

  # Test
  if(F){

    library(luckyBase)
    Plus.library(c('CCS', 'GSClassifier', 'plotly','cowplot','tidyr','ggplot2','ggpubr','purrr','furrr','stringi','digest', 'pROC','ComplexHeatmap','scales','plyr','dplyr','forestplot','ggrepel','writexl','readxl','patchwork','gtable','grid',"reshape2","circlize","parallel","foreach","doParallel","pROC"))
    source('./R/ccs_base.R',encoding = 'utf-8')


    data = read_xlsx('E:/iProjects/RCheck/GSClassifier/routine01/test/DataFrame_CCS+ClinicFeature_PanIMTv20240726+CDSDBv20240726.xlsx')
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
      cohort_level = intersect(cohort_level, as.character(unique(df$Cohort)))
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
    data_norm <- subtypeNorm(
      df,
      n_norm_subtype = n_norm_subtype,
      norm_roc_method = norm_roc_method,
      norm_rank_marker = norm_rank_marker,
      numCores = numCores,
      verbose = verbose
    )
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
      n_norm_subtype = n_norm_subtype
    ),
    Plot = list(
      ScatterBarPlot = plot_r,
      ForestPlot = plot_f
    ),
    Data = list(
      ROC = data_roc,
      Normalization = data_norm,
      ClinicUtility = data_utility
    )
  )
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
        if (length(unique(sub_df$response)) < 2) {
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
    'Merge' = p1234
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
    numCores = NULL, verbose = TRUE
){

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
    if(norm_rank_marker == 'RR+NRR'){
      df_cutoff <- arrange(df_cutoff, desc(medianRR), desc(medianNRR))
    } else if(norm_rank_marker == 'NRR+RR'){
      df_cutoff <- arrange(df_cutoff, desc(medianNRR), desc(medianRR))
    } else if(norm_rank_marker == 'RR'){
      df_cutoff <- arrange(df_cutoff, desc(medianRR))
    } else if(norm_rank_marker == 'NRR'){
      df_cutoff <- arrange(df_cutoff, desc(medianNRR))
    } else {
      stop('Not available rank marker. Please use one of "RR" and "NRR"!')
    }

    combinations <- combn(1:(nrow(df_cutoff)-1), n_norm_subtype - 1)
  }

  # Parallel processing: Find best combination
  if(T){

    # Parallel: Start
    if(is.null(numCores)){
      numCores <- detectCores() - 1  # Use all but one core
    }
    registerDoParallel(cores = numCores)

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


