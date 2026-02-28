

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
  signature(object='CCS'),
  function(
    object,
    data
  ){

    # Test
    if(F){
      library(luckyBase)
      Plus.library(c('tidyr','dplyr','plyr','reshape2','purrr','GSClassifier','ggplot2','cowplot'))
      object <- readRDS('E:/iProjects/RCheck/GSClassifier/routine01/ccs/PADv20240810_variant_01/resCCS_PADv20240810.rds')
      mt <- object@Data$Probability$d1
      data <- readRDS('E:/iProjects/CCS_Data/report/DataListForCCS_GEO+cBioPortal+UCXCXenav20240809_PADv20240810.rds')[c(1,2)]
      # data[['STAD']][['GSE21983']] <- NULL
    }

    # Matrix 01: Raw gene expression
    mt_1 <- getGeneMatrix(object, data)


    # Matrix 02: d1 matrix
    mt_2 <- object@Data$Probability$d1

    # Alignment
    coSample <- intersect(rownames(mt_1), rownames(mt_2))
    mt_l <- list()
    mt_l[[1]] <- mt_1[coSample,]
    mt_l[[2]] <- mt_2[coSample,]

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
      mt_l_plot[[1]][['Tissue']] + labs(title = 'Tissue-level') + theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5)),
      mt_l_plot[[1]][['Cohort']] + labs(title = 'Cohort-level') + theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.5)),
      mt_l_plot[[2]][['Tissue']] + theme(axis.title = element_blank()),
      mt_l_plot[[2]][['Cohort']] + theme(axis.title = element_blank()),
      align = "h", ncol = 2, rel_widths = c(0.5, 0.5),rel_heights = c(0.5,0.5)
    )
    x_axis_title <- ggdraw() + draw_label("Sample", x = 0.5, y = 0.5) # fontface = "bold",
    p2 <-  plot_grid(
      p, x_axis_title,
      ncol = 1,
      rel_heights = c(0.95, 0.05)
    )
    # tiff(paste0('E:/iProjects/RCheck/GSClassifier/routine01/test/PADv20240810_plotBatchEffect_test01.tif'), width = 40, height = 22, res = 300, units = "cm"); print(p2); dev.off()
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
getGeneMatrix <- function(object, data){

  data_f <- flatten(data)

  data_f2 <- ldply(data_f, function(X){

    if(F){
      X <- data_f[[1]]
    }

    X_subset <- geneMatch(
      X[["expr"]],
      geneAnnotation = object@Repeat$geneAnnotation,
      geneid = object@Repeat$geneid,
      matchmode = c("fix", "free")[1]
    )[['Subset']]

    if(is.one.true(colMeans(X_subset)>100)){
      LuckyVerbose('Outlier samples: ', paste(colnames(X_subset), collapse = '; '))
    }

    return(cbind(SampleIDs = colnames(X_subset),as.data.frame(t(X_subset))))

  }) %>% filter(!duplicated(SampleIDs))

  data_f3 <- as.matrix(data_f2[,-c(1,2)]); rownames(data_f3) <- as.character(data_f2[,2])

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
plotOneBatchEffect <- function(object,data,mt){

  # Reshape data
  mt <- mt[!duplicated(rownames(mt)),]
  mt_df <- cbind(ID = rownames(mt), as.data.frame(mt))
  mt_df2 <- mt_df %>% melt(id.vars = 'ID', measure.vars = colnames(mt), variable.name = 'variable', value.name = "value")

  # Annotation
  if(T){
    data <- data %>% flatten() %>% llply(.,function(x)colnames(x[['expr']]))
    annot <- data.frame()
    for(i in 1:length(data)){ # i=1
      annot <- rbind(annot,cbind(ID=data[[i]],cohort = names(data)[i]))
    }

    # Alignment
    coSample <- Reduce(intersect, list(as.character(annot$ID), rownames(mt), names(object@Data$CancerType)))
    annot <- annot[match(coSample, annot$ID),]
    annot$tissue <- object@Data$CancerType[coSample]
    annot <- annot %>% arrange(tissue, cohort)
    mt <- mt[coSample,]
    mt_df2 <- mt_df2 %>% filter(ID %in% coSample)
  }

  # Plot data
  mt_df3 <- left_join(mt_df2, annot, by='ID')
  mt_df3$ID <- factor(mt_df3$ID, levels = rev(annot$ID))
  nSample <- length(unique(mt_df3$ID)); β = 5/nSample

  # Here I don't want to use facet() series (～￣▽￣)～

  # Plot - Tissue-level
  if(T){

    p_batch1 <- ggplot(mt_df3,aes(x=ID, y=value, color=tissue, fill='white'))+
      geom_boxplot(
        # outlier.shape = NA,  # Hide outliers for clarity
        outlier.size = 0.001*β,
        width = 0.006*β, linewidth = 0.001*β) +
      stat_summary(fun=median, geom="point", shape=23, size=40*β,color='black',fill='black') +
      labs(x = 'Sample',y = 'Tissue-level') +
      guides(color = "none", fill = "none") +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", linewidth=1.5),
        axis.line = element_line(colour = "black", linewidth=0, linetype = 1),
        axis.title = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

    # tiff(paste0('./test/',project,'_plotBatch_mt-Matrix_ByTissue.tif'), width = 20, height = 10, res = 600, units = "cm"); print(p_batch1); dev.off()
  }

  # Plot - Cohort-level
  if(T){
    p_batch2 <- ggplot(mt_df3,aes(x=ID, y=value, color=cohort, fill='white'))+
      geom_boxplot(
        # outlier.shape = NA,  # Hide outliers for clarity
        outlier.size = 0.001*β,
        width = 0.006*β, linewidth = 0.001*β) +
      stat_summary(fun=median, geom="point", shape=23, size=40*β,color='black',fill='black') +
      labs(x = 'Sample',y = 'Cohort-level') +
      guides(color = "none", fill = "none") +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", linewidth=1.5),
        axis.line = element_line(colour = "black", linewidth=0, linetype = 1),
        axis.title = element_text(size = 10),
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
  p_batch[['Tissue']] <- p_batch1
  p_batch[['Cohort']] <- p_batch2
  # tiff(paste0('./test/',project,'_plotBatch_mt-Matrix_Tissue+Cohort.tif'), width = 20, height = 20, res = 300, units = "cm"); print(p_batch); dev.off()
  return(p_batch)

}


