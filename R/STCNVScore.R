# source('R/boundary_define_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
# OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
# assay = "Spatial"
# Sample = "CRC1"

#' Title Score ST data with infercnv result
#'
#' add cnv score to ST Seurat object
#'
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save figures and processed data
#' @param Sample Name of your sample
#' @param assay Name of assay uesed in infercnv
#'
#' @return A Seurat object of ST data with CNV sub-clusters and CNV scores
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
#' OutDir <- "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#' Sample <- "CRC1"
#' assay = "Spatial"
#' TumorST <- STCNVScore(TumorST = TumorST, assay = assay, OutDir = OutDir, Sample = Sample)
#'
STCNVScore <- function(TumorST = TumorST,
                       assay = c("Saptial","Morph"),
                       OutDir = NULL,
                       Sample = Sample) {

  if (is.null(OutDir) == TRUE) {
    OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }

  cnv_outdir = paste(OutDir, "InferCNV/output_", assay, sep = "")

  ## CNV Label
  cell_groupings <- read.table(paste(cnv_outdir,"/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",sep = ""),header = T)
  cell_groupings$class <- gsub("all_observations.all_observations.","",cell_groupings$cell_group_name)

  # set cnv score of normal cluster == 0
  NormalCluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(
    split(
      TumorST@meta.data[, c("seurat_clusters", "NormalScore")],
      TumorST@meta.data[, c("seurat_clusters", "NormalScore")]$seurat_clusters
    ),
    function(test) mean(test$NormalScore)
  )), decreasing = T)[1]]

  cell_groupings <- cell_groupings[!cell_groupings$cell %in% rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == NormalCluster,]),]

  class_df <- data.frame(cell_groups = rownames(as.data.frame.array(table(cell_groupings$class))),
                         cnv_label = c(1:8))

  cell_groupings$cnv_label <- class_df$cnv_label[match(cell_groupings$class,class_df$cell_groups)]
  rownames(cell_groupings) <- NULL

  infercnv.label <- cell_groupings[,c("cell","cnv_label")] %>% tibble::column_to_rownames(.,var = 'cell')

  infercnv.label <- rbind(infercnv.label, data.frame(
    row.names = rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == NormalCluster,]),
    cnv_label = rep("Normal", length(rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == NormalCluster,])))
  ))

  TumorST@meta.data$CNVLabel <- infercnv.label$cnv_label[match(rownames(TumorST@meta.data), rownames(infercnv.label))]

  .cluster_cols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
  )

  ## Plot CNV Label
  pdf(paste(OutDir, Sample, "_cnv_label.pdf", sep = ""), width = 7, height = 7)
  p <- SpatialDimPlot(TumorST, group.by = "CNVLabel", cols = .cluster_cols, pt.size.factor = 1, alpha = 0.6) +
    scale_fill_manual(values = .cluster_cols)
  print(p)
  dev.off()

  pdf(paste(OutDir, Sample, "_reduction_cnvlabel.pdf", sep = ""), width = 7, height = 7)
  p <- DimPlot(TumorST, group.by = "CNVLabel", cols = .cluster_cols) +
    scale_fill_manual(values = .cluster_cols)
  print(p)
  dev.off()

  ## CNV Score
  cell_gene_state <- read.table(paste(cnv_outdir,'/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_genes.dat',sep = ""),header = T)
  cell_gene_state$class <- gsub("all_observations.all_observations.","",cell_gene_state$cell_group_name)
  cell_gene_state$cnv_label <- class_df$cnv_label[match(cell_gene_state$class,class_df$cell_groups)]

  cnv_label_score <- data.frame(cnv_label = unique(cell_gene_state$cnv_label),
                                cnv_score = unlist(lapply(unique(cell_gene_state$cnv_label),function(label){

                                  cell_gene_state_sub <- cell_gene_state[cell_gene_state$cnv_label == label,]
                                  sub_label_score <- sum(abs(cell_gene_state_sub$state-3))
                                  return(sub_label_score)

                                }))
  )

  cell_scores_CNV <- rbind(cnv_label_score,data.frame(cnv_label="Normal",cnv_score = 0))
  TumorST@meta.data$cnv_score <- cell_scores_CNV$cnv_score[match(TumorST@meta.data$CNVLabel,cell_scores_CNV$cnv_label)]
  cell_scores_CNVA <- TumorST@meta.data[,c("CNVLabel","cnv_score")]

  ##plot cnv score
  pdf(paste(OutDir,Sample, "_cnv_observation_vlnplot.pdf", sep = ""), width = 6, height = 4)
  p <- ggplot(cell_scores_CNVA, aes(x = CNVLabel, y = cnv_score, fill = CNVLabel)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(stat = "boxplot", alpha = 1, width = .5, outlier.size = 0.5) +
    labs(y = "CNV_scores") +
    ggpubr::stat_compare_means() +
    scale_fill_manual(values = .cluster_cols) +
    theme(
      axis.text = element_text(colour = "black"),
      panel.background = element_blank(), panel.grid = element_blank(), legend.title = element_text(face = "bold"),
      axis.line = element_line(colour = "black"), legend.text = element_text(size = 10), title = element_text(face = "bold")
    ) +
    labs(title = "CNV Scores") +
    NoLegend()
  print(p)
  dev.off()

  return(TumorST)
}
