# source('R/spatial_decon_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# sig_exp <- readr::read_rds("YourPath/sig_exp.rds.gz")
# clustermarkers_list <- readr::read_rds("YourPath/clustermarkers_list.rds.gz")
#
# #get st and sc filtered expr data
# TumorST <- NormalizeData(TumorST,assay = "Spatial")
# TumorST@meta.data$Decon_topics <- paste(TumorST@meta.data$Location,TumorST@meta.data$seurat_clusters,sep = "_")
# expr_values = as.matrix(TumorST@assays$Spatial@data) #ST expr log
# nolog_expr = 2^(expr_values)-1 #ST expr nolog
#
# #filter st and sc feature
# intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
# filter_sig = sig_exp[intersect_gene,]
# filter_log_expr = expr_values[intersect_gene,]
# #get enrich matrix
# enrich_matrix <- get_enrich_matrix(filter_sig = filter_sig,clustermarkers_list = clustermarkers_list)

#' Title enrich analysis
#'
#' PAGE enrich analysis
#'
#' @param filter_log_expr A matrix of filtered log transformed ST data expression
#' @param enrich_matrix A binary matrix features as row, cell types as column
#'
#' @return A matrix of enrich result, cell types as row, spatial spots as column
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
#' sig_exp <- readr::read_rds("YourPath/sig_exp.rds.gz")
#' clustermarkers_list <- readr::read_rds("YourPath/clustermarkers_list.rds.gz")
#' # get st and sc filtered expr data
#' TumorST <- NormalizeData(TumorST,assay = "Spatial")
#' TumorST@meta.data$Decon_topics <- paste(TumorST@meta.data$Location,TumorST@meta.data$seurat_clusters,sep = "_")
#' expr_values = as.matrix(TumorST@assays$Spatial@data) #ST expr log
#' nolog_expr = 2^(expr_values)-1 #ST expr nolog
#' # filter st and sc feature
#' intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
#' filter_sig = sig_exp[intersect_gene,]
#' filter_log_expr = expr_values[intersect_gene,]
#' # get enrich matrix
#' enrich_matrix <- get_enrich_matrix(filter_sig = filter_sig,clustermarkers_list = clustermarkers_list)
#' enrich_result <- enrich_analysis(filter_log_expr = filter_log_expr,enrich_matrix = enrich_matrix)
#'
enrich_analysis <- function(filter_log_expr,
                            enrich_matrix) {


  # calculate mean gene expression
  mean_gene_expr <- log2(rowMeans(2^filter_log_expr - 1, dims = 1) + 1)

  geneFold <- filter_log_expr - mean_gene_expr
  # calculate sample/spot mean and sd
  cellColMean <- apply(geneFold, 2, mean)
  cellColSd <- apply(geneFold, 2, stats::sd)
  # get enrichment scores
  enrichment <- matrix(data = NA, nrow = dim(enrich_matrix)[2], ncol = length(cellColMean))
  for (i in (1:dim(enrich_matrix)[2])) {
    signames <- rownames(enrich_matrix)[which(enrich_matrix[, i] == 1)]
    sigColMean <- apply(geneFold[signames, ], 2, mean)
    m <- length(signames)
    vectorX <- NULL
    for (j in (1:length(cellColMean))) {
      Sm <- sigColMean[j]
      u <- cellColMean[j]
      sigma <- cellColSd[j]
      zscore <- (Sm - u) * m^(1 / 2) / sigma
      vectorX <- append(vectorX, zscore)
    }
    enrichment[i, ] <- vectorX
  }
  rownames(enrichment) <- colnames(enrich_matrix)
  colnames(enrichment) <- names(cellColMean)
  return(enrichment)
}
