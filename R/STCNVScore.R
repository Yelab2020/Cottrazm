#source('R/transition_define_settings.R')
#infercnv.dend <- read.tree(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",package = "Cottrazm"))
#cnv_table <- read.table(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt",package = "Cottrazm"))
#TumorST <- readr::read_rds("YourPath/TumorTransition/1.TransitionDefine/CRC1/TumorSTClustered.rds.gz")
#OutDir = "YourPath/TumorTransition/Fig/1.TransitionDefine/CRC1/"
#Sample = "CRC1"

#' Title Score ST data with infercnv result
#'
#' add cnv score to ST seurat object
#'
#' @param infercnv.dend InferCNV random tree subclustered result of ST data
#' @param cnv_table CNV scores of observations in ST data
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save figures and processed data
#' @param Sample Name of your sample
#'
#' @return A Seurat object of ST data with CNV subclusters and CNV scroes
#' @export
#'
#' @examples
#' infercnv.dend <- read.tree(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",package = "Cottrazm"))
#' cnv_table <- read.table(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt",package = "Cottrazm"))
#' TumorST <- readr::read_rds("YourPath/TumorTransition/1.TransitionDefine/CRC1/TumorSTClustered.rds.gz")
#' OutDir = "YourPath/TumorTransition/Fig/1.TransitionDefine/CRC1/"
#' Sample = "CRC1"
#' TumorST <- STCNVScore(infercnv.dend = infercnv.dend,cnv_table = cnv_table,TumorST = TumorST,OutDir = OutDir,Sample = Sample)

STCNVScore <- function(infercnv.dend = infercnv.dend,
                        cnv_table = cnv_table,
                        TumorST = TumorST,
                        OutDir = NULL,
                        Sample = Sample){

  if (is.null(OutDir) == TRUE){
    OutDir = paste(getwd(),"/",Sample,"/",sep = "")
    dir.create(OutDir)
  }

  ##CNV Cluster
  infercnv.label <- dendextend::cutree(infercnv.dend,k=8)
  infercnv.label <- as.data.frame(infercnv.label)
  infercnv.label <- rbind(infercnv.label,data.frame(row.names = rownames(TumorST@meta.data)[!rownames(TumorST@meta.data) %in% infercnv.label$row.names],
                                                    infercnv.label=rep("Normal",length(rownames(TumorST@meta.data)[!rownames(TumorST@meta.data) %in% infercnv.label$row.names]))))

  TumorST@meta.data$CNVLabel <- infercnv.label$infercnv.label[match(rownames(TumorST@meta.data),rownames(infercnv.label))]

  pdf(paste(OutDir,Sample,"_cnv_label.pdf",sep = ""),width = 7,height = 7)
  p <- SpatialDimPlot(TumorST,group.by = "CNVLabel",cols = .cluster_cols,pt.size.factor = 1,alpha = 0.6)
  print(p)
  dev.off()

  pdf(paste(OutDir,Sample,"_reduction_cnvlabel.pdf",sep = ""),width = 7,height = 7)
  p <- DimPlot(TumorST,group.by = "CNVLabel",cols = .cluster_cols)
  print(p)
  dev.off()

  #CNV Score
  colnames(cnv_table) <- gsub("\\.1","-1",colnames(cnv_table))
  # Replicate the table
  cnv_score_table <- as.matrix(cnv_table)
  cnv_score_table_pts <- abs(cnv_score_table - 3)

  # Scores are stored in “cnv_score_table_pts”. Use colSums to add up scores for each cell and store as vector
  cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
  colnames(cell_scores_CNV) <- "cnv_score"
  cell_scores_CNV$CNVLabel <- TumorST@meta.data$CNVLabel[match(rownames(cell_scores_CNV),rownames(TumorST@meta.data))]

  TumorST@meta.data$CNVScores <- cell_scores_CNV$cnv_score[match(rownames(TumorST@meta.data),rownames(cell_scores_CNV))]

  pdf(paste(OutDir,'cnv_observation_vlnplot.pdf',sep = ""),width = 6,height = 4)
  p <- ggplot(cell_scores_CNV,aes(x=CNVLabel,y=cnv_score,fill=CNVLabel))+
    geom_violin(alpha=0.5)+
    geom_boxplot(stat = "boxplot",alpha=1,width=.5,outlier.size = 0.5)+
    labs(y='CNV_scores')+
    ggpubr::stat_compare_means()+
    scale_fill_manual(values = .cluster_cols)+
    theme(axis.text = element_text(colour = 'black'),
          panel.background = element_blank(),panel.grid = element_blank(),legend.title = element_text(face = 'bold'),
          axis.line = element_line(colour = 'black'),legend.text = element_text(size = 10),title = element_text(face = 'bold'))+labs(title = 'CNV Scores')+NoLegend()
  print(p)
  dev.off()

  return(TumorST)
}
