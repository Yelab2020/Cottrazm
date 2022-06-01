#source('R/boundary_define_settings.R')
#InDir = paste(system.file("extdata/outs",package = "Cottrazm"),"/",sep = "")
#Sample = "CRC1"
#OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#TumorST <- STPreProcess(InDir = InDir,OutDir = OutDir,Sample = Sample)
#TumorST <- STModiCluster(InDir = InDir,Sample = Sample,OutDir = OutDir,TumorST = TumorST, res = 1.5)

#' Title ST data InferCNV
#'
#' Run infercnv of clustered Seurat object of tumor ST data
#'
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save infercnv results
#' @param assay The name of assay used in InferCNV (Morph: Morphological adjusted gene expression, Spatial: gene expression)
#'
#' @return A large CNV object
#' @export
#'
#' @examples
#' InDir = paste(system.file("extdata/outs",package = "Cottrazm"),"/",sep = "")
#' Sample = "CRC1"
#' OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#' TumorST <- STPreProcess(InDir = InDir,OutDir = OutDir,Sample = Sample)
#' TumorST <- STModiCluster(InDir = InDir,Sample = Sample,OutDir = OutDir,TumorST = TumorST, res = 1.5)
#' STInferCNV <- STCNV(TumorST = TumorST,OutDir = OutDir,assay = "Spatial")

STCNV <- function(TumorST = TumorST,
                  assay = c("Morph","Spatial"),
                  OutDir = NULL){

  if (is.null(OutDir) == TRUE){
    OutDir = paste(getwd(),"/",Sample,"/",sep = "")
    dir.create(OutDir)
  }

  #Run CNV random tree model
  #generate matrix
  matrix <- Seurat::GetAssayData(TumorST,slot = 'counts',assay = assay) %>% as.matrix()

  #Creat infercnv objcter
  annotation_file = paste(OutDir,"InferCNV/CellAnnotation.txt",sep = "")

  NormalCluster = levels(TumorST$seurat_clusters)[order(unlist(lapply(split(TumorST@meta.data[,c("seurat_clusters","NormalScore")],
                                                                            TumorST@meta.data[,c("seurat_clusters","NormalScore")]$seurat_clusters),
                                                                      function(test)mean(test$NormalScore))),decreasing = T)[1]]

  ref_cluster = NormalCluster
  gene_order_file = system.file("extdata/gencode_v38_gene_pos.txt",package = "Cottrazm")

  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = matrix,
                                       annotations_file = annotation_file,
                                       delim = '\t',
                                       gene_order_file = gene_order_file,
                                       ref_group_names = ref_cluster) #set reference group, if NULL, take mean of all spots as reference

  # perform infercnv operations to reveal cnv signal
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir = paste(OutDir,"InferCNV/output_",assay,sep = ""),  # output files
                                cluster_by_groups=F,   # cluster or not
                                analysis_mode = "subclusters",
                                denoise=T,
                                HMM=T,
                                tumor_subcluster_partition_method = "random_trees",
                                HMM_type = "i6",
                                BayesMaxPNormal = 0,
                                num_threads = 30,
                                plot_probabilities = F,
                                save_rds = F,
                                save_final_rds = F,
                                no_plot = F,
                                output_format = NA,
                                useRaster = T,
                                up_to_step = 17
                                )
  return(infercnv_obj)
}
