# source('R/spatial_decon_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# DeconData <- openxlsx::read.xlsx(system.file("extdata/DeconData.xlsx",package = "Cottrazm"))
# plot_col <- colnames(DeconData)[2:ncol(DeconData)]

#' Title Bar plot of deconvolution result
#'
#' Bar plot of deconvolution result in location groups
#'
#' @param DeconData A data frame of deconvolution result of each spatial spot
#' @param TumorST A Seurat object with defined malignant (Mal) spots, boundary (Bdy) spots, and non-malignnat (nMal) spots
#' @param plot_col Names of cell types to plot
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
#' DeconData <- oepnxlsx::read.xlsx(system.file("extdata/DeconData.xlsx", package = "Cottrazm"))
#' plot_col <- colnames(DeconData)[2:ncol(DeconData)]
#' DeconBarplot(DeconData = DeconData, TumorST = TumorST, plot_col = plot_col)
#'
DeconBarplot <- function(DeconData = DeconData,
                         TumorST = TumorST,
                         plot_col = plot_col) {
  metadata <- DeconData
  metadata[, "Location"] <- TumorST@meta.data$Location[match(metadata$cell_ID, rownames(TumorST@meta.data))]

  Location <- do.call(rbind, lapply(unique(metadata$Location), function(x) {
    metadata_split <- metadata[metadata$Location == x, ] %>% as.data.frame()
    sub <- colSums(metadata_split[, plot_col], na.rm = T) %>%
      data.frame() %>%
      tibble::rownames_to_column() %>%
      set_colnames(., c("Types", "Sum")) %>%
      dplyr::mutate(Per = 100 * Sum / sum(Sum))
  }))
  Location$Group <- rep(unique(metadata$Location), each = length(plot_col))
  Location$Types <- factor(Location$Types, levels = plot_col)
  Barplot <- ggplot(Location, aes(x = Group, y = Per, fill = Types)) +
    geom_bar(stat = "identity") +
    theme_bw()
  return(Barplot)
}
