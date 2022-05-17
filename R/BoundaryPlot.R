#source('R/boundary_define_settings.R')
#source('R/boundary_define_utility.R')
#infercnv.dend <- read.tree(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",package = "Cottrazm"))
#cnv_table <- read.table(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt",package = "Cottrazm"))
#TumorST <- readr::read_rds("YourPath/Tumorboundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
#Sample = "CRC1"
#OutDir = "YourPath/Tumorboundary/Fig/1.BoundaryDefine/CRC1/"
#TumorST <- STCNVScore(infercnv.dend = infercnv.dend,cnv_table = cnv_table,TumorST = TumorST,OutDir = OutDir)
#TumorLabel = c(1,2)
#TumorSTn <- BoundaryDefine(TumorST = TumorST, TumorLabel = TumorLabel, OutDir = OutDir, Sample = Sample)

#' Title Boundary plot
#'
#' Visualize boundary define results of ST data
#'
#' @param TumorSTn A subset Seurat object of ST data with defined tumor, boundary/transition zone, and normal spots
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save figures and processed data
#' @param Sample Name of your sample
#'
#' @return A Seurat object with defined tumor spots, boundary/transition zone spots, and out spots
#' @export
#'
#' @examples
#'
#' Sample = "CRC1"
#' OutDir = "YourPath/Tumorboundary/Fig/1.BoundaryDefine/CRC1/"
#' infercnv.dend <- read.tree(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",package = "Cottrazm"))
#' cnv_table <- read.table(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt",package = "Cottrazm"))
#' TumorST <- readr::read_rds("YourPath/Tumorboundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
#' TumorST <- STCNVScore(infercnv.dend = infercnv.dend,cnv_table = cnv_table,TumorST = TumorST,OutDir = OutDir)
#' TumorLabel = c(1,2)
#' TumorSTn <- BoundaryDefine(TumorST = TumorST, TumorLabel = TumorLabel, OutDir = OutDir, Sample = Sample)
#' TumorST <- BoundaryPlot(TumorSTn = TumorSTn, TumorST = TumorST, OutDir = OutDir,Sample = Sample)
#'
BoundaryPlot <- function(TumorSTn = TumorSTn,
                         TumorST = TumorST,
                         OutDir = OutDir,
                         Sample = Sample){


  if (is.null(OutDir) == TRUE)
    OutDir = paste(getwd(),"/",Sample,"/",sep = "")
  dir.create(OutDir)

  #get position
  position <- TumorST@images$image@coordinates
  position$spot.ids <- seq_len(nrow(position))

  #get spots neighbors
  dists <- compute_interspot_distances(position = position,scale.factor = 1.05)
  df_j <- find_neighbors(position = position,radius = dists$radius,method = "manhattan")
  names(df_j) <- position$spot.ids

  #get barcode
  Tumor_barcode = rownames(TumorSTn@meta.data)[grep("Tumor",TumorSTn@meta.data$LabelNew)]
  Trans_barcode = rownames(TumorSTn@meta.data)[grep("Trans",TumorSTn@meta.data$LabelNew)]

  Tumor_ID <- position[Tumor_barcode,]$spot.ids
  Trans_ID <- position[Trans_barcode,]$spot.ids

  Normal_Trans_ID = unique(unlist(nbrs(df_j = df_j,TumorSpotIDAdd = Tumor_ID,SpotIDRaw = c(Trans_ID,Tumor_ID))))
  Normal_Trans_barcode <- rownames(position[position$spot.ids %in% Normal_Trans_ID,])

  Out_barcode = rownames(TumorST@meta.data) %>% .[! .%in% c(Tumor_barcode,Trans_barcode,Normal_Trans_barcode)]

  Barcode_Ann <- data.frame(barcode=c(Tumor_barcode,Trans_barcode,Normal_Trans_barcode,Out_barcode),
                            Location= c(rep("Tumor",length(Tumor_barcode)),
                                           rep("Trans",length(c(Trans_barcode,Normal_Trans_barcode))),
                                           rep("OutSpots",length(Out_barcode))))
  #Barcode_Ann$Location <- gsub("1|2|3|4|5","",Barcode_Ann$Location)
  TumorST@meta.data$Location <- Barcode_Ann$Location[match(rownames(TumorST@meta.data),Barcode_Ann$barcode)]
  TumorST_DefineColors <- data.frame(Location=c("Tumor","Trans","OutSpots"),
                                     colors=c("#CB181D","#1f78b4","#fdb462"))
  TumorST@meta.data$Location <- factor(TumorST@meta.data$Location,levels=TumorST_DefineColors$Location)
  TumorST_DefineColors$colors <- as.character(TumorST_DefineColors$colors)

  pdf(paste(OutDir,Sample,"BoundaryDefine.pdf",sep = ""),width = 7,height = 7)
  p <- SpatialDimPlot(TumorST,group.by = "Location",cols = TumorST_DefineColors$colors)
  print(p)
  dev.off()

  readr::write_rds(TumorST,paste(OutDir,Sample,"BoundaryDefine.rds.gz",sep = ""),compress = "gz")

  return(TumorST)
}
