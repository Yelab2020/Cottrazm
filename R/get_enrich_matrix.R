# source('R/spatial_decon_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# sig_exp <- readr::read_rds("YourPath/sig_exp.rds.gz")
# clustermarkers_list <- readr::read_rds("YourPath/clustermarkers_list.rds.gz")
# #get st and sc filtered expr data
# TumorST <- NormalizeData(TumorST,assay = "Spatial")
# TumorST@meta.data$Decon_topics <- paste(TumorST@meta.data$Location,TumorST@meta.data$seurat_clusters,sep = "_")
# expr_values = as.matrix(TumorST@assays$Spatial@data) #ST expr log
# nolog_expr = 2^(expr_values)-1 #ST expr nolog
# meta_data <- TumorST@meta.data[,c("nCount_Spatial","Decon_topics","Location")]
#
# #Signature score
# for(cluster in names(clustermarkers_list)){
#   cluster_markers = clustermarkers_list[[cluster]][1:25]
#   cluster_score <- apply(TumorST@assays$Spatial@data[rownames(TumorST@assays$Spatial@data) %in% cluster_markers,],2,mean)
#   meta_data <- cbind(meta_data,cluster_score)
# }
# colnames(meta_data) <- c("nCount_Spatial,"Decon_topics","Location",names(clustermarkers_list))
#
# #filter st and sc feature
# intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
# filter_sig = sig_exp[intersect_gene,]

#' Title get enrich matrix
#'
#' Construct binary matrix for PAGE enrichment
#'
#' @param filter_sig  The sig_exp filtered by features
#' @param clustermarkers_list A list of markers of each cell type of scRNAseq dataset
#'
#' @return A binary matrix features as row, cell types as column
#' @export
#'
#' @examples
#' enrich_matrix <- get_enrich_matrix(filter_sig = filter_sig,clustermarkers_list = clustermarkers_list)
#'
get_enrich_matrix <- function(filter_sig = filter_sig, clustermarkers_list = clustermarkers_list) {

  # construct enrich matrix initial
  enrich_matrix <- matrix(0, nrow = dim(filter_sig)[1], ncol = dim(filter_sig)[2])
  rownames(enrich_matrix) <- rownames(filter_sig)
  colnames(enrich_matrix) <- colnames(filter_sig)

  for (i in 1:ncol(enrich_matrix)) {
    cluster <- colnames(enrich_matrix)[i]
    feature <- intersect(clustermarkers_list[[cluster]], rownames(enrich_matrix))
    enrich_matrix[feature, cluster] <- 1
  }

  return(enrich_matrix)
}
