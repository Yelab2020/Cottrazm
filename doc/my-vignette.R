## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE,warning=FALSE---------------------------------------------
library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(patchwork)
library(ggtree)
library(BiocGenerics)
library(readr)
library(rtracklayer)
library(infercnv)
library(phylogram)
library(utils)
library(dendextend)
library(assertthat)
library(reticulate)
library(openxlsx)
library(scatterpie)
library(cowplot)
library(stats)
library(quadprog)
library(data.table)
library(Rfast)
library(ggrepel)
library(tibble)
library(clusterProfiler)
library(utils)
library(org.Hs.eg.db)

## ----eval=FALSE---------------------------------------------------------------
#  #from local path
#  devtools::install_github('Yelab2020/Cottrazm')
#  #from github
#  install.packages("Cottrazm.tar.gz", repos = NULL, type = "source")
#  library(Cottrazm)

## ----eval = FALSE-------------------------------------------------------------
#  print('STPreProcess')
#  InDir = paste(system.file("extdata/outs",package = "Cottrazm"),"/",sep = "")
#  Sample = "CRC1"
#  OutDir = paste(getwd(),"/",Sample,"/",sep = "")
#  TumorST <-
#    STPreProcess(
#      InDir = InDir,
#      OutDir = OutDir,
#      Sample = Sample
#      )

## ----eval = FALSE-------------------------------------------------------------
#  print('STModiCluster')
#  res = 1.5
#  TumorST <-
#    STModiCluster(
#      InDir = InDir,
#      Sample = Sample,
#      OutDir = OutDir,
#      TumorST = TumorST,
#      res = res
#    )

## ----eval = FALSE-------------------------------------------------------------
#  print('STCNV' )
#  STInferCNV <-
#    STCNV(
#      TumorST = TumorST,
#      OutDir  = OutDir,
#      assay = "Spatial"
#      )

## ----eval = FALSE-------------------------------------------------------------
#  print('STCNVScore')
#  
#  TumorST <-
#    STCNVScore(
#      TumorST = TumorST,
#      assay = "Spatial",
#      Sample = Sample,
#      OutDir = OutDir
#      )

## ----eval = FALSE-------------------------------------------------------------
#  TumorSTn <-
#    BoundaryDefine(
#      TumorST = TumorST,
#      MalLabel = c(7,8),
#      OutDir = OutDir,
#      Sample = Sample
#    )

## ----eval = FALSE-------------------------------------------------------------
#  TumorST <-
#    BoundaryPlot(
#      TumorSTn = TumorSTn,
#      TumorST = TumorST,
#      OutDir = OutDir,
#      Sample = Sample
#    )

## ----echo=FALSE---------------------------------------------------------------
TumorST <- readr::read_rds('/work/xzz123/Project/2021/TumorBoundary/package/CRC1/CRC1_BoundaryDefine.rds.gz')
Seurat::SpatialDimPlot(TumorST,group.by = 'Location',cols = c("#CB181D", "#1f78b4", "#fdb462"))

## ----eval = FALSE-------------------------------------------------------------
#  #clustermarkers_list
#  sc_obj <- readr::read_rds("YourPath/scRNAseqRef.rds.gz")
#  clustermarkers <-
#    Seurat::FindAllMarkers(object = sc_obj,
#                           logfc.threshold = 0.25,
#                           only.pos = TRUE)
#  
#  clustermarkers_list <- split(clustermarkers, clustermarkers$cluster)
#  clustermarkers_list <-
#    lapply(names(clustermarkers_list), function(cluster) {
#      sub_markers <- clustermarkers_list[[cluster]]$gene
#    })
#  names(clustermarkers_list) <-
#    names(split(clustermarkers, clustermarkers$cluster))
#  
#  #sig_exp
#  sig_exp <-
#    get_sig_exp(
#      se.obj = sc_obj,
#      DefineTypes = "MajorTypes",
#      sig_scran = unique(unlist(clustermarkers_list))
#    )
#  
#  #In this vignette we read processed clustermarkers_list and sig_exp
#  clustermarkers_list <- readr::read_rds(system.file("inst/extdata/clustermarkers_list.rds.gz",package = "Cottrazm"))
#  sig_exp <- readr::read_rds(system.file("inst/extdata/sig_exp.rds.gz",package = "Cottrazm"))

## ----eval = FALSE-------------------------------------------------------------
#  #Processed ST data
#  TumorST <- NormalizeData(TumorST, assay = "Spatial")
#  TumorST@meta.data$Decon_topics <-
#    paste(TumorST@meta.data$Location,
#          TumorST@meta.data$seurat_clusters,
#          sep = "_")
#  expr_values = as.matrix(TumorST@assays$Spatial@data) #ST expr log
#  nolog_expr = 2 ^ (expr_values) - 1 #ST expr nolog
#  meta_data <-
#    TumorST@meta.data[, c("nCount_Spatial", "Decon_topics", "Location")]
#  
#  #Signature score
#  for (cluster in names(clustermarkers_list)) {
#    cluster_markers = clustermarkers_list[[cluster]][1:25]
#    cluster_score <-
#      apply(TumorST@assays$Spatial@data[rownames(TumorST@assays$Spatial@data) %in% cluster_markers, ], 2, mean)
#    meta_data <- cbind(meta_data, cluster_score)
#  }
#  colnames(meta_data) <-
#    c("nCount_Spatial",
#      "Decon_topics",
#      "Location",
#      names(clustermarkers_list))
#  
#  #filter st and sc feature
#  intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
#  filter_sig = sig_exp[intersect_gene, ]
#  filter_expr = nolog_expr[intersect_gene, ]
#  filter_log_expr = expr_values[intersect_gene, ]
#  
#  #enrichment analysis
#  enrich_matrix <-
#    get_enrich_matrix(filter_sig = filter_sig,
#                      clustermarkers_list = clustermarkers_list)
#  enrich_result <-
#    enrich_analysis(filter_log_expr = filter_log_expr,
#                    enrich_matrix = enrich_matrix)

## ----eval = FALSE-------------------------------------------------------------
#  DeconData <- SpatialDecon(
#    enrich_matrix = enrich_matrix,
#    enrich_result = enrich_result,
#    filter_expr = filter_expr,
#    filter_sig = filter_sig,
#    clustermarkers_list = clustermarkers_list,
#    meta_data = meta_data,
#    malignant_cluster = "Malignant epithelial cells",
#    tissue_cluster = "Epithelial cells",
#    stromal_cluster = "Fibroblast cells"
#    )

## ----eval=FALSE---------------------------------------------------------------
#  DeconData <- openxlsx::read.xlsx(system.file("inst/extdata/DeconData.xlsx",package = "Cottrazm"))

## ----echo=FALSE---------------------------------------------------------------
DeconData <- openxlsx::read.xlsx('/work/xzz123/Project/2021/TumorBoundary/package/Cottrazm/inst/extdata/DeconData.xlsx')
print(head(DeconData))

## ----eval=FALSE---------------------------------------------------------------
#  plot_col <- colnames(DeconData)[2:ncol(DeconData)]
#  img_path = system.file("inst/extdata/outs/spatial/tissue_lowres_image.png",package = "Cottrazm")
#  
#  DeconPieplot(
#    DeconData = DeconData,
#    TumorST = TumorST,
#    plot_col = plot_col,
#    img_path = img_path,
#    pie_scale = 0.4,
#    scatterpie_alpha = 0.8,
#    border_color = "grey"
#    )
#  
#  DeconBarplot(
#    DeconData = DeconData,
#    TumorST = TumorST,
#    plot_col = plot_col
#    )

## ----echo=FALSE---------------------------------------------------------------
plot_col <- colnames(DeconData)[2:ncol(DeconData)]
img_path = "/work/xzz123/Project/2021/TumorBoundary/package/Cottrazm/inst/extdata/outs/spatial/tissue_lowres_image.png"

  DeconData = DeconData;
  TumorST = TumorST;
  plot_col = plot_col;
  img_path = img_path;
  pie_scale = 0.4;
  scatterpie_alpha = 0.8;
  border_color = "grey"

DeconData <- DeconData[DeconData$cell_ID %in% rownames(TumorST@meta.data), ]

  ## Preprocess data
  slice <- names(TumorST@images)[1]

  spatial_coord <- data.frame(TumorST@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("cell_ID") %>%
    dplyr::mutate(
      imagerow_scaled =
        imagerow * TumorST@images[[slice]]@scale.factors$lowres,
      imagecol_scaled =
        imagecol * TumorST@images[[slice]]@scale.factors$lowres
    ) %>%
    dplyr::inner_join(DeconData, by = "cell_ID")

  ### Load histological image into R
  img_path <- img_path # lowers image png(input dir)
  img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))

  if (img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  } else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }

  # Convert image to grob object
  img_grob <- grid::rasterGrob(img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc")
  )


  ## Plot spatial scatterpie plot
  scatterpie_plt <- ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = 0,
      ymax = -nrow(img)
    ) +
    scatterpie::geom_scatterpie(
      data = spatial_coord,
      ggplot2::aes(
        x = imagecol_scaled,
        y = imagerow_scaled
      ),
      cols = plot_col,
      color = border_color,
      alpha = scatterpie_alpha,
      pie_scale = pie_scale,
      lwd = 0.1
    ) +
    ggplot2::scale_y_reverse() +
    ggplot2::ylim(nrow(img), 0) +
    ggplot2::xlim(0, ncol(img)) +
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_void() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = NULL,
      ylim = NULL,
      expand = TRUE,
      clip = "on"
    )
  
 print(scatterpie_plt)
 
 ##Decon Bar
  DeconData = DeconData;
  TumorST = TumorST;
  plot_col = plot_col

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
  print(Barplot)


## ----eval=FALSE---------------------------------------------------------------
#  TumorSTRecon <-
#    SpatialRecon(
#      TumorST = TumorST,
#      sig_exp = sig_exp,
#      clustermarkers_list = clustermarkers_list,
#      DeconData = DeconData,
#      Location = "Bdy"
#    )

## ----eval=FALSE---------------------------------------------------------------
#  DiffGenes <- FindDiffGenes(TumorST = TumorST, assay = "Spatial")

## ----eval=FALSE---------------------------------------------------------------
#  DiffVolcanoplot(
#    DiffGenes = DiffGenes,
#    Location = "Bdy",
#    cut_off_pvalue = 2,
#    cut_off_logFC = 0.25,
#    n = 15
#  )

## ----echo=FALSE---------------------------------------------------------------
DiffGenes <- readr::read_rds('/work/xzz123/Project/2021/TumorBoundary/Fig/3.TransComposition/CRC1/1.SpatialDiffGenes.rds')
names(DiffGenes) <- c("nMal","Bdy","Mal")
  DiffGenes = DiffGenes;
  Location = "Bdy";
  cut_off_pvalue = 2;
  cut_off_logFC = 0.25;
  n = 15

# get Location DEG
  LocationDiff <- DiffGenes[[Location]]
  LocationDiff <- LocationDiff[LocationDiff$Symbol %in% grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", LocationDiff$Symbol, invert = T, value = T), ]

  dataset <- LocationDiff %>%
    tibble::rownames_to_column() %>%
    set_colnames(., c("gene", colnames(LocationDiff)))
  dataset$color <- ifelse(-log10(LocationDiff$FDR) < cut_off_pvalue | abs(LocationDiff$Diff) <= cut_off_logFC, "grey",
    ifelse(LocationDiff$Diff > 0, "Bdy", "Other")
  )
  datasetN <- dataset %>%
    dplyr::group_by(color) %>%
    dplyr::top_n(n, abs(Diff))
  datasetN <- datasetN[datasetN$color != "grey", ]

  # plot
  p <- ggplot(dataset, aes(x = Diff, y = (-log10(pvalue)), color = color)) +
    geom_point(aes(fill = color), size = 1) +
    scale_color_manual(values = c("red", "grey", "blue")) +
    ggrepel::geom_text_repel(
      data = datasetN, aes(label = datasetN$Symbol), max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
      size = 3, point.padding = unit(0.35, "lines"), segment.color = "black", show.legend = F, color = "black", fontface = "bold"
    ) +
    geom_vline(xintercept = c(-cut_off_logFC, cut_off_logFC), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = cut_off_pvalue, lty = 4, col = "black", lwd = 0.8) +
    labs(x = "log2(FoldChange)", y = "-log10(pvalueue)") +
    theme(panel.background = element_rect(color = "black", fill = NA), legend.position = "bottom", legend.title = element_blank())
  print(p)

## ----eval=FALSE---------------------------------------------------------------
#  boundary_enrich <-
#    FeatureEnrichment(
#      DiffGenes = DiffGenes,
#      cut_off_logFC = 0.25,
#      cut_off_pvalue = 0.05,
#      Location = "Bdy"
#    )

## -----------------------------------------------------------------------------
sessionInfo()

