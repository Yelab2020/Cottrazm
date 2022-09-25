# source('R/boundary_define_settings.R')
# InDir = paste(system.file("extdata/outs",package = "Cottrazm"),"/",sep = "")
# Sample = "CRC1"
# OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"

#' Title ST data preprocess
#'
#' read ST raw data and quality control
#'
#' @param InDir Path to file contain Spaceranger out result (at least including 'filtered_feature_bc_matrix','filtered_feature_bc_matrix.h5','spatial')
#' @param Sample Name of your sample
#' @param OutDir Path to file save figures and processed data
#'
#' @return A Seurat object
#' @export
#'
#' @examples
#' InDir <- paste(system.file("extdata/outs", package = "Cottrazm"), "/", sep = "")
#' Sample <- "CRC1"
#' OutDir <- "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#' TumorST <- STPreProcess(InDir = InDir, OutDir = OutDir, Sample = Sample)
#'
STPreProcess <- function(InDir = InDir, Sample = Sample, OutDir = NULL) {
  if (is.null(OutDir) == TRUE) {
    OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }

  # read files------
  aa_try <- try(
    Xdata <- Seurat::Read10X(data.dir = paste(InDir, "filtered_feature_bc_matrix", sep = "")),
    silent = T
  )
  if (is(aa_try, "try-error")) {
    Xdata <- Seurat::Read10X_h5(filename = paste(InDir, "filtered_feature_bc_matrix.h5", sep = ""))
  } else {
    Xdata <- Xdata
  }

  XF <- CreateSeuratObject(counts = Xdata, project = Sample, min.spots = 0, assay = "Spatial")
  # read image filesa
  Ximage <- Read10X_Image(image.dir = paste(InDir, "spatial", sep = ""))
  Seurat::DefaultAssay(Ximage) <- "Spatial"
  # link matrix and image file
  Ximage <- Ximage[colnames(XF)]
  XF[["image"]] <- Ximage
  TumorST <- XF

  # QC----
  dir.create(paste(OutDir, "QC", sep = ""))
  TumorST[["Mito.percent"]] <- PercentageFeatureSet(TumorST, pattern = "^MT-")

  pdf(paste(OutDir, "QC/Vlnplot.pdf", sep = ""), width = 6, height = 4)
  p <- VlnPlot(TumorST, features = c("nFeature_Spatial", "nCount_Spatial", "Mito.percent"), pt.size = 0, combine = F)
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoLegend() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
  }
  p <- cowplot::plot_grid(plotlist = p, ncol = 3)
  print(p)
  dev.off()

  pdf(paste(OutDir, "QC/featurplot.pdf", sep = ""), width = 7, height = 7)
  p <- SpatialFeaturePlot(TumorST, features = c("nFeature_Spatial", "nCount_Spatial", "Mito.percent"), combine = F)
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
  }
  print(cowplot::plot_grid(plotlist = p, ncol = 3))
  dev.off()

  QCData <- TumorST@meta.data[, c("nCount_Spatial", "nFeature_Spatial", "Mito.percent")]
  openxlsx::write.xlsx(QCData, paste(OutDir, "QC/QCData.xlsx", sep = ""), overwrite = T)

  return(TumorST)
}
