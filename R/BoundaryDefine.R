# source('R/boundary_define_settings.R')
# source('R/boundary_define_utility.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
# Sample = "CRC1"
# OutDir = "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
# TumorST <- STCNVScore(TumorST = TumorST,assay = "Spatial",OutDir = OutDir,Sample = Sample)
# MalLabel = c(1,2)

#' Title Boundary define
#'
#' Find boundary of tumor in ST data
#'
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save figures and processed data
#' @param Sample Name of your sample
#' @param MalLabel CNV labels with high cnv scores which could be defined as malignant label. If NULL, Cottrazm will take the two CNVLabels with the highest average CNVScores.
#'
#' @return A subset Seurat object of ST data with defined malignant spots (Mal), boundary spots (Bdy), and non-malignant spots (nMal)
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
#' Sample <- "CRC1"
#' OutDir <- "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
#' TumorST <- STCNVScore(TumorST = TumorST, assay = "Spatial", OutDir = OutDir, Sample = Sample)
#' MalLabel <- c(1, 2)
#' TumorSTn <- BoundaryDefine(TumorST = TumorST, MalLabel = MalLabel, OutDir = OutDir, Sample = Sample)
#'
BoundaryDefine <- function(TumorST = TumorST,
                           MalLabel = NULL,
                           OutDir = NULL,
                           Sample = Sample) {
  if (is.null(OutDir) == TRUE) {
    OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }

  # get UMAPembeddings
  UMAPembeddings <- as.data.frame(TumorST@reductions$umap@cell.embeddings %>% set_colnames(., c("x", "y")))
  # get position
  slice <- names(TumorST@images)[1]
  position <- data.frame(TumorST@images[[slice]]@coordinates)
  position$spot.ids <- seq_len(nrow(position))

  # get spots neighbors
  dists <- compute_interspot_distances(position = position, scale.factor = 1.05)
  df_j <- find_neighbors(position = position, radius = dists$radius, method = "manhattan")

  # set primiary MalCellID by CNV
  if (is.null(MalLabel) == TRUE) {
    MalLabel <-
      order(unlist(lapply(
        split(
          #TumorST@meta.data[, c("CNVLabel", "CNVScores")],
          #TumorST@meta.data[, c("CNVLabel", "CNVScores")]$CNVLabel
          TumorST@meta.data[, c("CNVLabel", "cnv_score")],
          TumorST@meta.data[, c("CNVLabel", "cnv_score")]$CNVLabel
        ),
        function(test) {
          #median(test$CNVScores)
          median(test$cnv_scores)
        }
      )), decreasing = T)[1:2]
  }

  MalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$CNVLabel %in% MalLabel, ])

  NormalCluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(
    split(
      TumorST@meta.data[, c("seurat_clusters", "NormalScore")],
      TumorST@meta.data[, c("seurat_clusters", "NormalScore")]$seurat_clusters
    ),
    function(test) mean(test$NormalScore)
  )), decreasing = T)[1]]
  NormalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == NormalCluster, ])

  # get better MalCellID by umap position
  CNV_seurat_df <- as.data.frame.array(table(TumorST@meta.data$CNVLabel, TumorST@meta.data$seurat_clusters))[MalLabel, ]
  ClusterID <- c()
  for (cluster in levels(TumorST@meta.data$seurat_clusters)) {
    sub_cluster_colsum <- sum(CNV_seurat_df[, cluster])
    if (sub_cluster_colsum > table(TumorST@meta.data$seurat_clusters)[cluster] * 0.5) {
      ClusterID <- c(ClusterID, cluster)
    }
  }

  # generate MalCluster tibble
  CiMal <- tibble::tibble(cluster = ClusterID)
  CiMal <- CiMal %>%
    dplyr::mutate(sub_MalCellID = purrr::map(.x = cluster, .f = function(.x) {
      sub_MalCellID <- intersect(MalCellID, rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == .x, ]))
      return(sub_MalCellID)
    })) %>%
    dplyr::mutate(sub_CiMal = purrr::map(.x = sub_MalCellID, .f = function(.x) {
      sub_CiMal <- apply(UMAPembeddings[.x, ], 2, mean)
      return(sub_CiMal)
    }))

  CiNormal <- apply(UMAPembeddings[NormalCellID, ], 2, mean)

  # filter MalCell
  MalCellIDsi <- purrr::map2(.x = CiMal$sub_MalCellID, .y = CiMal$sub_CiMal, function(.x, .y) {
    sub_MalCellsi <- lapply(.x, function(id) {
      pos <- UMAPembeddings[id, ]
      rt <- sqrt(sum((pos - .y)^2))
      rn <- sqrt(sum((pos - CiNormal)^2))
      if (rt < 1 / 3 * rn) {
        return(id)
      }
    }) %>% unlist()
    return(sub_MalCellsi)
  }) %>% unlist()

  # deal with lonely spots
  MalCellIDL <- lapply(MalCellIDsi, function(name) {
    ifelse(length(df_j[[name]][df_j[[name]] %in% MalCellIDsi]) == 0, name, NA)
  }) %>%
    unlist() %>%
    na.omit()

  BdyCellID <- NULL

  nbrs_of_MalL <- nbrs(df_j = df_j, MalCellIDAdd = MalCellIDL, CellIDRaw = c(MalCellIDsi, NormalCellID, BdyCellID))
  ClusterL <- do.call(rbind, lapply(names(nbrs_of_MalL), function(celll) {
    sub <- TumorST@meta.data[celll, ]$seurat_clusters
    nbrs_of_celll <- nbrs_of_MalL[[celll]]
    cluster <- do.call(rbind, lapply(nbrs_of_celll, function(idl) {
      pos <- UMAPembeddings[idl, ]
      rt <- sqrt(sum((pos - unlist(CiMal[CiMal$cluster == sub, ]$sub_CiMal))^2))
      rn <- sqrt(sum((pos - CiNormal)^2))
      cluster <- data.frame(CellID = idl, Location = ifelse(rt < 1 / 3 * rn, "Mal", "Bdy"))
    }))
    return(cluster)
  }))

  ClusterL <- as.data.frame.array(table(ClusterL$CellID, ClusterL$Location))
  aa_try <- try(
    aa <- as.character(lapply(rownames(ClusterL), function(spotID) {
      ifelse(ClusterL[spotID, ]$Mal > 0, "Mal", "Bdy")
    })),
    silent = T
  )
  if (!is(aa_try, "try-error")) {
    ClusterL$Location <- aa
  } else {
    ClusterL$Location <- rep("Bdy", nrow(ClusterL))
  }
  ClusterL <- data.frame(CellID = rownames(ClusterL), Location = as.character(ClusterL$Location))

  # prepare initial IDs for repeat
  MalCellID <- c(MalCellIDsi, ClusterL[ClusterL$Location == "Mal", ]$CellID)
  BdyCellID <- as.character(ClusterL[ClusterL$Location == "Bdy", ]$CellID)
  MalCellIDN <- MalCellID
  n <- 1
  TumorSTn <- TumorST
  Clustern <- rbind(
    data.frame(CellID = NormalCellID, Location = rep("Normal", length(NormalCellID))),
    data.frame(CellID = MalCellID, Location = rep("Mal", length(MalCellID))),
    data.frame(CellID = BdyCellID, Location = rep("Bdy", length(BdyCellID)))
  )
  repeat {
    if (length(MalCellIDN) < 3) {
      print(paste("Stop at n = ", (n - 1), sep = ""))
      return(TumorSTn)
    } else {

      nbrs_of_Mal <- nbrs(df_j = df_j, MalCellIDAdd = MalCellIDN, CellIDRaw = c(MalCellID, NormalCellID, BdyCellID))

      if (length(unique(unlist(nbrs_of_Mal))) < 3) {
        print(paste("Stop at n = ", (n - 1), sep = ""))
        return(TumorSTn)
      } else {
        TumorSTn <- subset(TumorST, cells = c(unique(unlist(nbrs_of_Mal)), MalCellID, NormalCellID, BdyCellID))
        TumorSTn@meta.data$Label <- Clustern$Location[match(rownames(TumorSTn@meta.data), Clustern$CellID)]

        if (n == 1){
          TumorSTn@meta.data$Label <- factor(TumorSTn@meta.data$Label,levels = c("Normal","Bdy","Mal"))
        }else{
          TumorSTn@meta.data$Label <- factor(TumorSTn@meta.data$Label, levels = c("Normal", "Bdy", "Mal", paste("Mal", c(1:n-1), sep = "")))
        }

        #nbr location define
        ClusterAdd <- ClusterUpdate(
          x = n,
          position = position,
          df_j = df_j,
          UMAPembeddings = UMAPembeddings,
          MalCellIDN = MalCellIDN,
          BdyCellID = BdyCellID,
          NormalCellID = NormalCellID,
          MalCellID = MalCellID
        )
        Clustern <- rbind(Clustern, ClusterAdd)
        TumorSTn@meta.data$LabelNew <- Clustern$Location[match(rownames(TumorSTn@meta.data), as.character(Clustern$CellID))]
        TumorSTn@meta.data$LabelNew <- factor(TumorSTn@meta.data$LabelNew, levels = c("Normal", "Bdy", "Mal", paste("Mal", c(1:n), sep = "")))


        pdf(paste(OutDir, Sample, "_with_nbrs", n, ".pdf", sep = ""), width = 7, height = 7)
        p1 <- SpatialDimPlot(TumorSTn, cols = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:n]), group.by = "Label")+
          scale_fill_manual(values = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:n]))
        print(p1)
        dev.off()

        pdf(paste(OutDir, Sample, "_reduction_with_nbrs", n, ".pdf", sep = ""), width = 7, height = 7)
        p2 <- DimPlot(TumorSTn, cols = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:n]), group.by = "Label") +
          scale_fill_manual(values = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:n]))+
          theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) + labs(title = NULL)
        print(p2)
        dev.off()

        pdf(paste(OutDir, Sample, "_out_", n, ".pdf", sep = ""), width = 7, height = 7)
        p3 <- SpatialDimPlot(TumorSTn, group.by = "LabelNew", cols = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:(n + 1)]))+
          scale_fill_manual(values = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:(n + 1)]))
        print(p3)
        dev.off()

        pdf(paste(OutDir, Sample, "_reduction_out_", n, ".pdf", sep = ""), width = 7, height = 7)
        p4 <- DimPlot(TumorSTn, group.by = "LabelNew", cols = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:(n + 1)])) +
          scale_fill_manual(values = c("#33a02c", "#1f78b4", rev(c("#fef0d9", "#fdd49e", "#fdbb84", "#fc8d59", "#ef6548", "#d7301f", "#990000"))[1:(n + 1)]))+
          theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) + labs(title = NULL)
        print(p4)
        dev.off()

        MalCellID <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% c("Mal", paste("Mal", c(1:n), sep = "")), ])
        MalCellIDN <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% paste("Mal", n, sep = ""), ])
        BdyCellID <- rownames(TumorSTn@meta.data[TumorSTn$LabelNew == "Bdy", ])

        n <- n + 1
      }
    }
    if (n > 6) {
      print(paste("Stop at n = ", (n - 1), sep = ""))
      return(TumorSTn)
    }
  }
  return(TumorSTn)
}
