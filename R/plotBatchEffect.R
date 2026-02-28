setGeneric("plotBatchEffect", function(object, ...) {
  standardGeneric("plotBatchEffect")
})

#' @rdname CCS-method.plotBatchEffect
#' @title CCS method: plotBatchEffect
#' @description \code{plotBatchEffect} method for \code{CCS} class
#' @inheritParams CCSPublicParams
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @import ggplot2
#' @seealso \code{\link{ccs}}.
#' @exportMethod plotBatchEffect
setMethod(
  "plotBatchEffect",
  signature(object = "CCS"),
  function(object,
           data) {
    # Test
    if (F) {
      library(luckyBase)
      Plus.library(c("tidyr", "dplyr", "plyr", "reshape2", "purrr", "GSClassifier", "ggplot2", "cowplot"))
      # 跨平台路径适配
      root_path <- ifelse(Sys.info()["sysname"] == "Darwin", "/Volumes/2T01/winE", "E:")
      object <- readRDS(file.path(root_path, "iProjects/RCheck/GSClassifier/test02/ccs/PADv20240810/prediction/PanIMTv20240726_legacy_02/resCCS_PADv20240810.rds"))
      mt <- object@Data$Probability$d1
      data <- readRDS(file.path(root_path, "iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810.rds"))[c(1, 2)]
      # data[['STAD']][['GSE21983']] <- NULL
    }

    # Matrix 01: Raw gene expression
    mt_1 <- getGeneMatrix(object, data)


    # Matrix 02: d1 matrix
    mt_2 <- object@Data$Probability$d1

    # Alignment
    coSample <- intersect(rownames(mt_1), rownames(mt_2))
    mt_l <- list()
    mt_l[[1]] <- mt_1[coSample, ]
    mt_l[[2]] <- mt_2[coSample, ]

    # Plot data
    mt_l_plot <- llply(mt_l, function(mt) plotOneBatchEffect(object, data, mt))

    # Result
    # p <- plot_grid(
    #   mt_l_plot[[1]],
    #   mt_l_plot[[2]] + theme(axis.title.y = element_blank()),
    #   align = "h", ncol = 2, rel_widths = c(0.5, 0.5)
    # )

    # Plot
    p <- plot_grid(
      mt_l_plot[[1]][["Tissue"]] + labs(title = "Tissue-level") + theme(axis.title = element_blank()),
      mt_l_plot[[1]][["Cohort"]] + labs(title = "Cohort-level") + theme(axis.title = element_blank()),
      mt_l_plot[[2]][["Tissue"]] + theme(axis.title = element_blank()),
      mt_l_plot[[2]][["Cohort"]] + theme(axis.title = element_blank()),
      align = "h", ncol = 2, rel_widths = c(0.5, 0.5), rel_heights = c(0.5, 0.5)
    )
    # 添加坐标轴标题（Nature 级别样式）
    # 使用 cowplot::ggdraw + draw_label 实现更精确的控制
    base_font_size <- 11

    x_axis_title <- ggdraw() +
      draw_label(
        "Sample",
        fontface = "bold",
        size = base_font_size,
        hjust = 0.5,
        vjust = 0.5
      )
    y_axis_title <- ggdraw() +
      draw_label(
        "Expression",
        fontface = "bold",
        size = base_font_size,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5
      )

    # 组合图形：使用更协调的比例
    p2 <- plot_grid(
      y_axis_title,
      plot_grid(p, x_axis_title, ncol = 1, rel_heights = c(0.92, 0.08)),
      ncol = 2, rel_widths = c(0.03, 0.97)
    )
    # tiff(file.path(root_path, '/iProjects/RCheck/GSClassifier/routine01/test/PADv20240810_plotBatchEffect_test01.tif'), width = 40, height = 22, res = 300, units = "cm"); print(p2); dev.off()
    return(p2)
  }
)


#### Assistant functions ####

#' @description Get gene expression matrix from a \code{CCSDataList}
#' @importFrom tidyr `%>%`
#' @importFrom plyr ldply
#' @importFrom dplyr filter
#' @importFrom purrr flatten
#' @importFrom GSClassifier geneMatch
getGeneMatrix <- function(object, data) {
  data_f <- flatten(data)

  data_f2 <- ldply(data_f, function(X) {
    if (F) {
      X <- data_f[[1]]
    }

    X_subset <- geneMatch(
      X[["expr"]],
      geneAnnotation = object@Repeat$geneAnnotation,
      geneid = object@Repeat$geneid,
      matchmode = c("fix", "free")[1]
    )[["Subset"]]

    if (is.one.true(colMeans(X_subset) > 100)) {
      LuckyVerbose("Outlier samples: ", paste(colnames(X_subset), collapse = "; "))
    }

    return(cbind(SampleIDs = colnames(X_subset), as.data.frame(t(X_subset))))
  }) %>% filter(!duplicated(SampleIDs))

  data_f3 <- as.matrix(data_f2[, -c(1, 2)])
  rownames(data_f3) <- as.character(data_f2[, 2])

  return(data_f3)
}


#' @description Plot batch effect
#' @param mt a matrix with feature cols and sample rows like \code{d1} matrix in \code{ccs} object.
#' @inheritParams ccs
#' @inheritParams cluster
#' @importFrom reshape2 melt
#' @importFrom purrr flatten
#' @importFrom plyr llply
#' @importFrom dplyr arrange filter left_join
#' @import ggplot2
plotOneBatchEffect <- function(object, data, mt) {
  # Reshape data
  mt <- mt[!duplicated(rownames(mt)), ]
  mt_df <- cbind(ID = rownames(mt), as.data.frame(mt))
  mt_df2 <- mt_df %>% melt(id.vars = "ID", measure.vars = colnames(mt), variable.name = "variable", value.name = "value")

  # Annotation
  if (T) {
    data <- data %>%
      flatten() %>%
      llply(., function(x) colnames(x[["expr"]]))
    annot <- data.frame()
    for (i in 1:length(data)) { # i=1
      annot <- rbind(annot, cbind(ID = data[[i]], cohort = names(data)[i]))
    }

    # Alignment
    coSample <- Reduce(intersect, list(as.character(annot$ID), rownames(mt), names(object@Data$CancerType)))
    annot <- annot[match(coSample, annot$ID), ]
    annot$tissue <- object@Data$CancerType[coSample]
    annot <- annot %>% arrange(tissue, cohort)
    mt <- mt[coSample, ]
    mt_df2 <- mt_df2 %>% filter(ID %in% coSample)
  }

  # Plot data
  mt_df3 <- left_join(mt_df2, annot, by = "ID")
  mt_df3$ID <- factor(mt_df3$ID, levels = rev(annot$ID))
  nSample <- length(unique(mt_df3$ID))
  β <- 5 / nSample

  # Here I don't want to use facet() series (～￣▽￣)～

  # Plot - Tissue-level
  if (T) {
    p_batch1 <- ggplot(mt_df3, aes(x = ID, y = value, color = tissue, fill = "white")) +
      geom_boxplot(
        # outlier.shape = NA,  # Hide outliers for clarity
        outlier.size = 0.001 * β,
        width = 0.006 * β, linewidth = 0.001 * β
      ) +
      stat_summary(fun = median, geom = "point", shape = 23, size = 40 * β, color = "black", fill = "black") +
      labs(x = "Sample", y = "Tissue-level") +
      guides(color = "none", fill = "none") +
      theme_bw() +
      theme(
        # Nature 级别样式
        text = element_text(family = "sans", colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        # 边框与网格
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        # 隐藏 x 轴刻度（样本过多）
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

    # tiff(paste0('./test/',project,'_plotBatch_mt-Matrix_ByTissue.tif'), width = 20, height = 10, res = 600, units = "cm"); print(p_batch1); dev.off()
  }

  # Plot - Cohort-level
  if (T) {
    p_batch2 <- ggplot(mt_df3, aes(x = ID, y = value, color = cohort, fill = "white")) +
      geom_boxplot(
        # outlier.shape = NA,  # Hide outliers for clarity
        outlier.size = 0.001 * β,
        width = 0.006 * β, linewidth = 0.001 * β
      ) +
      stat_summary(fun = median, geom = "point", shape = 23, size = 40 * β, color = "black", fill = "black") +
      labs(x = "Sample", y = "Cohort-level") +
      guides(color = "none", fill = "none") +
      theme_bw() +
      theme(
        # Nature 级别样式
        text = element_text(family = "sans", colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        # 边框与网格
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        # 隐藏 x 轴刻度（样本过多）
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  # Merge
  # p_batch <- plot_grid(
  #   p_batch1 + theme(axis.title.x = element_blank()),
  #   p_batch2,
  #   align = "v", ncol = 1, rel_heights = c(0.5, 0.5)
  # )
  p_batch <- list()
  p_batch[["Tissue"]] <- p_batch1
  p_batch[["Cohort"]] <- p_batch2
  # tiff(paste0('./test/',project,'_plotBatch_mt-Matrix_Tissue+Cohort.tif'), width = 20, height = 20, res = 300, units = "cm"); print(p_batch); dev.off()
  return(p_batch)
}
