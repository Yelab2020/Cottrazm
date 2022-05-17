# source('R/spatial_dwls_settings.R')
# TumorST <- readr::read_rds("YourPath/Tumorboundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# DeconData <- read.table(system.file("extdata/DeconData.csv",package = "Cottrazm"),sep = ",",header = T,row.names = "X")
# plot_col <- colnames(DeconData)[2:ncol(DeconData)]

#' Title Barplot of deconvolution result
#'
#' Barplot of deconvolution result in boundary define groups
#'
#' @param DeconData A data frame of deconvolution result of each spatial spot
#' @param TumorST A Seurat object with defined tumor spots, boundary/transition zone spots, and out spots
#' @param plot_col Names of cell types to plot
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/Tumorboundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
#' DeconData <- read.table(system.file("extdata/DeconData.csv",package = "Cottrazm"),sep = ",",header = T,row.names = "X")
#' plot_col <- colnames(DeconData)[2:ncol(DeconData)]
#' DeconBarplot(DeconData = DeconData, TumorST = TumorST, plot_col = plot_col)
#'
DeconBarplot <- function(DeconData = DeconData,
                         TumorST = TumorST,
                         plot_col = plot_col){

  metadata = DeconData
  metadata[,'BoundaryDefine'] = TumorST@meta.data$Location[match(metadata$cell_ID,rownames(TumorST@meta.data))]

  Boundary <- do.call(rbind,lapply(unique(metadata$BoundaryDefine),function(x){
    metadata_split <- metadata[metadata$BoundaryDefine == x,] %>% as.data.frame()
    sub <- colSums(metadata_split[,plot_col],na.rm = T) %>% data.frame()  %>% tibble::rownames_to_column() %>% set_colnames(.,c("Types","Sum")) %>%
      dplyr::mutate(Per=100*Sum/sum(Sum))
  }))
  Boundary$Group <- rep(unique(metadata$BoundaryDefine),each=length(plot_col))
  Boundary$Types <- factor(Boundary$Types,levels = plot_col)
  Barplot <- ggplot(Boundary,aes(x=Group,y=Per,fill=Types))+
    geom_bar(stat = "identity")+
    theme_bw()
  return(Barplot)
}
