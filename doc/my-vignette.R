## ----echo = TRUE--------------------------------------------------------------
# library(Seurat)
# library(magrittr)
# library(dplyr)
# library(Matrix)
# library(ggplot2)
# library(stringr)
# library(RColorBrewer)
# library(patchwork)
# library(ggtree)
# library(BiocGenerics)
# library(readr)
# library(rtracklayer)
# library(infercnv)
# library(phylogram)
# library(utils)
# library(dendextend)
# library(assertthat)
# library(reticulate)
# library(openxlsx)
# library(scatterpie)
# library(cowplot)
# library(stats)
# library(quadprog)
# library(data.table)
# library(Rfast)
# library(ggrepel)
# library(tibble)
# library(clusterProfiler)
# library(utils)
# library(org.Hs.eg.db)

## ----echo = TRUE--------------------------------------------------------------
#devtools::install_github('Yelab2020/Cottrazm')
#install.packages("Cottrazm.tar.gz", repos = NULL, type = "source")
#library(Cottrazm)

## -----------------------------------------------------------------------------
# print('STPreProcess')
# InDir = "Spaceranger/outs/"
# Sample = "YourSampleName"
# OutDir = "YourOutDir/"
# TumorST <-
#   STPreProcess(InDir = InDir,
#                OutDir = OutDir,
#                Sample = Sample)

## -----------------------------------------------------------------------------
# print('STModiCluster')
# res = 1.5
# TumorST <-
#   STModiCluster(
#     InDir = InDir,
#     Sample = Sample,
#     OutDir = OutDir,
#     TumorST = TumorST,
#     res = res
#   )

## -----------------------------------------------------------------------------
# print('STCNV')
# STInferCNV <-
#   STCNV(TumorST = TumorST,
#         OutDir  = OutDir,
#         assay = "Spatial")

## -----------------------------------------------------------------------------
# print('STCNVScore')
# 
# infercnv.dend <-
#   read.tree(
#     file =
#       "YourPath/infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt"
#   )
# 
# cnv_table <-
#   read.table(
#     file =
#       'YourPath/infercnv.17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt'
#   )
# 
# TumorST <-
#   STCNVScore(
#     infercnv.dend = infercnv.dend,
#     cnv_table = cnv_table,
#     TumorST = TumorST,
#     Sample = Sample,
#     OutDir = OutDir
#   )

## -----------------------------------------------------------------------------
# TumorSTn <-
#   TransitionDefine(
#     TumorST = TumorST,
#     TumorLabel = c(1,2),
#     OutDir = OutDir,
#     Sample = Sample
#   )

## -----------------------------------------------------------------------------
# TumorST <-
#   TransitionPlot(
#     TumorSTn = TumorSTn,
#     TumorST = TumorST,
#     OutDir = OutDir,
#     Sample = Sample
#   )

## -----------------------------------------------------------------------------
#clustermarkers_list
# sc_obj <- readr::read_rds("YourPath/scRNAseqRef.rds.gz")
# clustermarkers <-
#   Seurat::FindAllMarkers(object = sc_obj,
#                          logfc.threshold = 0.25,
#                          only.pos = TRUE)
# 
# clustermarkers_list <- split(clustermarkers, clustermarkers$cluster)
# clustermarkers_list <-
#   lapply(names(clustermarkers_list), function(cluster) {
#     sub_markers <- clustermarkers_list[[cluster]]$gene
#   })
# names(clustermarkers_list) <-
#   names(split(clustermarkers, clustermarkers$cluster))
# 
# #sig_exp
# sig_exp <-
#   get_sig_exp(
#     se.obj = sc_obj,
#     DefineTypes = "MajorTypes",
#     sig_scran = unique(unlist(clustermarkers_list))
#   )

## -----------------------------------------------------------------------------
# #Processed ST data
# TumorST <- NormalizeData(TumorST, assay = "Spatial")
# TumorST@meta.data$Decon_topics <-
#   paste(TumorST@meta.data$Location,
#         TumorST@meta.data$seurat_clusters,
#         sep = "_")
# expr_values = as.matrix(TumorST@assays$Spatial@data) #ST expr log
# nolog_expr = 2 ^ (expr_values) - 1 #ST expr nolog
# meta_data <-
#   TumorST@meta.data[, c("nCount_Spatial", "Decon_topics", "Location")]
# 
# #Signature score
# for (cluster in names(clustermarkers_list)) {
#   cluster_markers = clustermarkers_list[[cluster]][1:25]
#   cluster_score <-
#     apply(TumorST@assays$Spatial@data[rownames(TumorST@assays$Spatial@data) %in% cluster_markers, ], 2, mean)
#   meta_data <- cbind(meta_data, cluster_score)
# }
# colnames(meta_data) <-
#   c("nCount_Spaital",
#     "Decon_topics",
#     "Location",
#     names(clustermarkers_list))
# 
# #filter st and sc feature
# intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
# filter_sig = sig_exp[intersect_gene, ]
# filter_expr = nolog_expr[intersect_gene, ]
# filter_log_expr = expr_values[intersect_gene, ]
# 
# #enrichment analysis
# enrich_matrix <-
#   get_enrich_matrix(filter_sig = filter_sig,
#                     clustermarkers_list = clustermarkers_list)
# enrich_result <-
#   enrich_analysis(filter_log_expr = filter_log_expr,
#                   enrich_matrix = enrich_matrix)

## -----------------------------------------------------------------------------
# DeconData <- SpatialDecon(enrich_matrix = enrich_matrix,
#                          enrich_result = enrich_result,
#                          filter_expr = filter_expr,
#                          filter_sig = filter_sig,
#                          clustermarkers_list = clustermarkers_list,
#                          meta_data = meta_data,
#                          malignant_cluster = "Malignant epithelial cells",
#                          tissue_cluster = "Epithelial cells",
#                          stromal_cluster = "Fibroblast cells")
# 
# print(head(DeconData))


## -----------------------------------------------------------------------------
# plot_col <- colnames(DeconData)[2:ncol(DeconData)]
# img_path = "inst/extdata/outs/spatial/tissue_lowres_image.png"
# DeconPieplot(
#   DeconData = DeconData,
#   TumorST = TumorST,
#   plot_col = plot_col,
#   img_path = img_path,
#   pie_scale = 0.4,
#   scatterpie_alpha = 0.8,
#   border_color = "grey"
# )
# 
# DeconBarplot(DeconData = DeconData,
#              TumorST = TumorST,
#              plot_col = plot_col)

## -----------------------------------------------------------------------------
# TumorSTRecon <-
#   SpatialRecon(
#     TumorST = TumorST,
#     sig_exp = sig_exp,
#     clustermarkers_list = clustermarkers_list,
#     DeconData = DeconData,
#     Location = c("Trans")
#   )

## -----------------------------------------------------------------------------
# DiffGenes <- FindDiffGenes(TumorST = TumorST, assay = "Spatial")

## -----------------------------------------------------------------------------
# DiffVolcanoplot(
#   DiffGenes = DiffGenes,
#   Location = "Trans",
#   cut_off_pvalue = 2,
#   cut_off_logFC = 0.25,
#   n = 15
# )

## -----------------------------------------------------------------------------
# transition_enrich <-
#   FeatureEnrichment(
#     DiffGenes = DiffGenes,
#     cut_off_logFC = 0.25,
#     cut_off_pvalue = 0.05,
#     Location = "Trans"
#   )

## -----------------------------------------------------------------------------
sessionInfo()

