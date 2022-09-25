#' Title Get the mean expression of significant features across scRNAseq cell types
#'
#' Get mean expression of scRNAseq significant features across each cell type
#'
#' @param se.obj A scRNAseq Seurat object of corresponding tissue
#' @param DefineTypes The column in scRNAseq Surat object metadata: cell types defined in scRNAseq data (usually including: T, B, Myeloid cells, Fibroblast cells, Endothelial cells, Malignant cells, Adjacent normal cells and .etc)
#' @param sig_scran Markers of each scRNAseq cell type
#'
#' @return A matrix features as row, scRNAseq cell types as column
#' @export
#'
#' @examples
#' sig_exp <- get_sig_exp(se.obj = sc_obj,DefineTypes = "MajorTypes",sig_scran = unique(unlist(clustermarkers_list)))

#'
get_sig_exp <- function(se.obj = WholeTissueSC,
                        DefineTypes = "Majortypes",
                        sig_scran = sig_scran # unique genes
) {
  norm_exp <- 2^(se.obj@assays$RNA@data) - 1
  id <- se.obj@meta.data[, DefineTypes]
  ExprSubset <- norm_exp[sig_scran, ]
  Sig_exp <- NULL
  for (i in unique(id)) {
    Sig_exp <- cbind(Sig_exp, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
  }
  colnames(Sig_exp) <- unique(id) # mean expression of each cell type
  return(Sig_exp)
}
