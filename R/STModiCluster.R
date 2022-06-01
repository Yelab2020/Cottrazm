#source('R/boundary_define_settings.R')
#InDir = paste(system.file("extdata/outs",package = "Cottrazm"),"/",sep = "")
#Sample = "CRC1"
#OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#TumorST <- STPreProcess(InDir = InDir,OutDir = OutDir,Sample = Sample)

#' Title Morphological adjusted cluster determination
#'
#' Grep HE staining information as Morphological matrix to adjusted gene expression matrix and use adjusted assay to cluster ST data
#'
#' @param TumorST A Seurat object
#' @param res Value of the resolution parameter
#' @param InDir Path to file contain Spaceranger out result
#' @param Sample Name of your sample
#' @param OutDir Path to file save figures and processed data
#'
#' @return A Seurat object with morphological adjusted expression matrix and determined clusters
#' @export
#'
#' @examples
#' InDir = paste(system.file("extdata/outs",package = "Cottrazm"),"/",sep = "")
#' Sample = "CRC1"
#' OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#' TumorST <- STPreProcess(InDir = InDir,OutDir = OutDir,Sample = Sample)
#' TumorST <- STModiCluster(InDir = InDir,Sample = Sample,OutDir = OutDir,TumorST = TumorST, res = 1.5)

STModiCluster <- function(InDir = InDir,
                          Sample = Sample,
                          OutDir = NULL,
                          TumorST = TumorST,
                          res = 1.5){

  if (is.null(OutDir) == TRUE){
    OutDir = paste(getwd(),"/",Sample,"/",sep = "")
    dir.create(OutDir)
  }

  #Adjusted_expr_mtx
  use_condaenv("TumorBoundary", required = TRUE)
  source_python(system.file("python/Rusedtile.py",package = "Cottrazm"))
  Adjusted_expr_mtx <- ME_normalize(inDir=InDir,outDir=OutDir,sample=Sample)

  #Create Morph seu.obj
  aa_try <- try(
    rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts),silent = T
  )
  if(is(aa_try,'try-error')){
    library(Matrix)
    Adjusted_expr_mtx <-  Matrix::readMM(paste(OutDir,Sample,"_raw_SME_normalizeA.mtx",sep = ""))
  }else{
    rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts)
    colnames(Adjusted_expr_mtx) <- rownames(TumorST@assays$Spatial@counts)
  }

  Adjusted_expr_mtxF <- t(as.matrix(Adjusted_expr_mtx))
  MorphMatirxSeurat <- CreateSeuratObject(counts = as.matrix(Adjusted_expr_mtxF))

  #Add morph feature as assay to TumorST
  MorphMatirxSeurat <- subset(MorphMatirxSeurat,cells = rownames(TumorST@meta.data))
  TumorST@assays$Morph <- MorphMatirxSeurat@assays$RNA

  #Use Morph as assay Cluster
  TumorST <- NormalizeData(TumorST,assay = "Morph")
  TumorST <- FindVariableFeatures(object =TumorST, mean.function = ExpMean, dispersion.function = LogVMR,assay = "Morph")
  TumorST <- ScaleData(object = TumorST,assay = "Morph")#,vars.to.regress = c('Mito.percent','Ribo.percent'))
  TumorST <- RunPCA(object = TumorST,  npcs = 50, verbose = FALSE,assay = "Morph")
  TumorST <- FindNeighbors(TumorST, reduction = "pca", dims = 1:50,assay = "Morph")
  TumorST <- RunUMAP(object = TumorST, dims = 1:50,assay = "Morph")
  TumorST <- FindClusters(TumorST, resolution = res, algorithm = 1,graph.name = "Morph_snn")

  TumorST@meta.data$seurat_clusters <- TumorST@meta.data[, paste("Morph_snn_res.",res,sep = "")]
  #plot cluster result
  pdf(paste(OutDir,Sample,"_Spatial_SeuratCluster.pdf",sep = ""),width = 7,height = 7)
  p <- SpatialDimPlot(TumorST,group.by = "seurat_clusters",cols = .cluster_cols,pt.size.factor = 1,alpha = 0.8)+labs(title = paste('Resolution = ',res,sep = ""))
  print(p)
  dev.off()

  pdf(paste(OutDir,Sample,"_UMAP_SeuratCluster.pdf",sep = ""),width = 7,height = 7)
  p <- DimPlot(TumorST,group.by = "seurat_clusters",cols = .cluster_cols)+labs(title = paste('Resolution = ',res,sep = ""))
  print(p)
  dev.off()

  #Add ImmuneScore
  Normalfeatures <- c("NKG7","KLRC1","CCR7","FOXP3","CTLA4","CD8B","CD8A","CXCR6","CD3D","CD3E","PTPRC","MZB1","CD79A","SELL","CD19","AICDA")
  TumorST@meta.data$NormalScore <- apply(TumorST@assays$Spatial@data[rownames(TumorST@assays$Spatial@data) %in% Normalfeatures,],2,mean)

  pdf(paste(OutDir,Sample,"_NormalScore.pdf",sep = ""),width = 6,height = 4)
  p <- VlnPlot(TumorST,features = "NormalScore",pt.size = 0,group.by = "seurat_clusters",cols = .cluster_cols)+
    geom_boxplot()+
    geom_hline(yintercept = max(unlist(lapply(split(TumorST@meta.data[,c("seurat_clusters","NormalScore")],TumorST@meta.data[,c("seurat_clusters","NormalScore")]$seurat_clusters),
                                              function(test)median(test$NormalScore)))),linetype="dashed")+
    ggpubr::stat_compare_means()+NoLegend()
  print(p)
  dev.off()

  NormalCluster = levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[,c("seurat_clusters","NormalScore")],
                                                                            TumorST@meta.data[,c("seurat_clusters","NormalScore")]$seurat_clusters),
                                                                      function(test)mean(test$NormalScore))),decreasing = T)[1]]
  print(paste("NormalCluster = ",NormalCluster,sep = ""))

  #save CNV annotation file
  cellAnnotation <- data.frame(CellID=rownames(TumorST@meta.data),DefineTypes=TumorST@meta.data[,'seurat_clusters'])
  dir.create(paste(OutDir,"InferCNV",sep = ""))
  write.table(cellAnnotation,paste(OutDir,'InferCNV/CellAnnotation.txt',sep=""),sep = '\t',row.names = F,col.names = F,quote = F)

  return(TumorST)
}
