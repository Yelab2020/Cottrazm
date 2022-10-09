# source('R/boundary_define_settings.R')
# source('R/boundary_define_utility.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
# Sample = "CRC1"
# OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
# TumorST <- STCNVScore(TumorST = TumorST, assay = "Spatial", OutDir = OutDir, Sample = Sample)
# MalLabel = c(1,2)
# TumorSTn <- BoundaryDefine(TumorST = TumorST, MalLabel = MalLabel, OutDir = OutDir, Sample = Sample)

#' Title Boundary plot
#'
#' Visualize boundary define results of ST data
#'
#' @param TumorSTn A subset Seurat object of ST data with defined malignant (Mal) spots, boundary (Bdy) spots, and non-malignant (nMal) spots
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save figures and processed data
#' @param Sample Name of your sample
#'
#' @return A Seurat object with defined malignant (Mal) spots, boundary (Bdy) spots, and non-malignant (nMal) spots
#' @export
#'
#' @examples
#'
#' Sample <- "CRC1"
#' OutDir <- "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
#' TumorST <- STCNVScore(TumorST = TumorST, assay = "Spatial", OutDir = OutDir, Sample = Sample)
#' MalLabel <- c(1, 2)
#' TumorSTn <- BoundaryDefine(TumorST = TumorST, MalLabel = MalLabel, OutDir = OutDir, Sample = Sample)
#' TumorST <- BoundaryPlot(TumorSTn = TumorSTn, TumorST = TumorST, OutDir = OutDir, Sample = Sample)
#'
BoundaryPlot <- function(TumorSTn = TumorSTn,
                         TumorST = TumorST,
                         OutDir = OutDir,
                         Sample = Sample) {
  if (is.null(OutDir) == TRUE) {
    OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }

  # get position
  position <- TumorST@images$image@coordinates

  # get spots neighbors
  dists <- compute_interspot_distances(position = position, scale.factor = 1.05)
  df_j <- find_neighbors(position = position, radius = dists$radius, method = "manhattan")

  # get barcode
  Mal_barcode <- rownames(TumorSTn@meta.data)[grep("Mal", TumorSTn@meta.data$LabelNew)]
  Bdy_barcode <- rownames(TumorSTn@meta.data)[grep("Bdy", TumorSTn@meta.data$LabelNew)]

  Normal_Bdy_barcode <- unique(unlist(nbrs(df_j = df_j, MalCellIDAdd = Mal_barcode, CellIDRaw = c(Bdy_barcode, Mal_barcode))))

  nMal_barcode <- rownames(TumorST@meta.data) %>% .[!. %in% c(Mal_barcode, Bdy_barcode, Normal_Bdy_barcode)]

  Barcode_Ann <- data.frame(
    barcode = c(Mal_barcode, Bdy_barcode, Normal_Bdy_barcode, nMal_barcode),
    Location = c(
      rep("Mal", length(Mal_barcode)),
      rep("Bdy", length(c(Bdy_barcode, Normal_Bdy_barcode))),
      rep("nMal", length(nMal_barcode))
    )
  )

  TumorST@meta.data$Location <- Barcode_Ann$Location[match(rownames(TumorST@meta.data), Barcode_Ann$barcode)]
  TumorST_DefineColors <- data.frame(
    Location = c("Mal", "Bdy", "nMal"),
    colors = c("#CB181D", "#1f78b4", "#fdb462")
  )
  TumorST@meta.data$Location <- factor(TumorST@meta.data$Location, levels = TumorST_DefineColors$Location)
  TumorST_DefineColors$colors <- as.character(TumorST_DefineColors$colors)

  pdf(paste(OutDir, Sample, "_BoundaryDefine.pdf", sep = ""), width = 7, height = 7)
  p <- SpatialDimPlot(TumorST, group.by = "Location", cols = TumorST_DefineColors$colors) +
    scale_fill_manual(values = TumorST_DefineColors$colors )
  print(p)
  dev.off()

  readr::write_rds(TumorST, paste(OutDir, Sample, "_BoundaryDefine.rds.gz", sep = ""), compress = "gz")

  return(TumorST)
}
