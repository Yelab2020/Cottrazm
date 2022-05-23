#source('R/transition_define_settings.R')
#source('R/transition_define_utility.R')
#infercnv.dend <- read.tree(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",package = "Cottrazm"))
#cnv_table <- read.table(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt",package = "Cottrazm"v))
#TumorST <- readr::read_rds("YourPath/TumorTransition/1.TransitionDefine/CRC1/TumorSTClustered.rds.gz")
#Sample = "CRC1"
#OutDir = "YourPath/TumorTransition/Fig/1.TransitionDefine/CRC1/"
#TumorST <- STCNVScore(infercnv.dend = infercnv.dend,cnv_table = cnv_table,TumorST = TumorST,OutDir = OutDir)
#TumorLabel = c(1,2)

#' Title transition define
#'
#' Find boudnary/transition zone of tumor in ST data
#'
#' @param TumorST A Seurat object with morphological adjusted expression matrix and determined clusters
#' @param OutDir Path to file save figures and processed data
#' @param Sample Name of your sample
#' @param TumorLabel CNV labels with high cnv scores which could be defined as tumor. If NULL, Cottrazm will take the two CNVLabels with the highest average CNVScores.
#'
#' @return A subset Seurat object of ST data with defined tumor, transition zone, and normal spots
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorTransition/1.TransitionDefine/CRC1/TumorSTClustered.rds.gz")
#' Sample = "CRC1"
#' OutDir = "YourPath/TumorTransition/Fig/1.TransitionDefine/CRC1/"
#' infercnv.dend <- read.tree(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations_dendrogram.txt",package = "Cottrazm"))
#' cnv_table <- read.table(file = system.file("extdata/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.observations.txt",package = "Cottrazm"))
#' TumorST <- STCNVScore(infercnv.dend = infercnv.dend,cnv_table = cnv_table,TumorST = TumorST,OutDir = OutDir)
#' TumorLabel = c(1,2)
#' TumorSTn <- TransitionDefine(TumorST = TumorST, TumorLabel = TumorLabel, OutDir = OutDir, Sample = Sample)
#'
TransitionDefine <- function(TumorST = TumorST,
                           TumorLabel = NULL,
                           OutDir = NULL,
                           Sample = Sample){

  if (is.null(OutDir) == TRUE){
    OutDir = paste(getwd(),"/",Sample,"/",sep = "")
    dir.create(OutDir)
  }

  #get UMAPembeddings
  UMAPembeddings <- as.data.frame(TumorST@reductions$umap@cell.embeddings %>% set_colnames(.,c("x","y")))
  #get position
  slice <- names(TumorST@images)[1]
  position <- data.frame(TumorST@images[[slice]]@coordinates)
  position$spot.ids <- seq_len(nrow(position))

  #get spots neighbors
  dists <- compute_interspot_distances(position = position,scale.factor = 1.05)
  df_j <- find_neighbors(position = position,radius = dists$radius,method = "manhattan")
  names(df_j) <- position$spot.ids

  #set primiary tumorcellID by CNV
  TumorST@meta.data$CNVLabel <- factor(TumorST@meta.data$CNVLabel,levels = c(1:8,"Normal"))

  if (is.null(TumorLabel) == TRUE)
    TumorLabel <-
    order(unlist(lapply(split(TumorST@meta.data[, c("CNVLabel", "CNVScores")],
                              TumorST@meta.data[, c("CNVLabel", "CNVScores")]$CNVLabel),
                        function(test)
                          median(test$CNVScores))), decreasing = T)[1:2]

  TumorCellID <- rownames(TumorST@meta.data[TumorST@meta.data$CNVLabel %in% TumorLabel,])
  TumorSpotID <- position[TumorCellID,]$spot.ids

  NormalCellID <- rownames(TumorST@meta.data[TumorST@meta.data$CNVLabel == "Normal",])
  NormalSpotID <- position[rownames(position) %in% NormalCellID,]$spot.ids

  #get better tumorcellID by umap position
  #intergrate CNV tumor cluster and suerat TumorLabel
  CNV_seurat_df <- as.data.frame.array(table(TumorST@meta.data$CNVLabel,TumorST@meta.data$seurat_clusters))[TumorLabel,]
  ClusterID <- c()
  for(cluster in levels(TumorST@meta.data$seurat_clusters)){
    sub_cluster_colsum <- sum(CNV_seurat_df[,cluster])
    if(sub_cluster_colsum > table(TumorST@meta.data$seurat_clusters)[cluster]*0.5){
      ClusterID <- c(ClusterID,cluster)
    }
  }

  #generate CiTumor tibble
  CiTumor <- tibble::tibble(cluster= ClusterID)
  CiTumor <- CiTumor %>% dplyr::mutate(sub_TumorCellID = purrr::map(.x = cluster,.f = function(.x){
    sub_TumorCellID <- intersect(TumorCellID,rownames(TumorST@meta.data[TumorST@meta.data$seurat_clusters == .x,]))
    return(sub_TumorCellID)
  })) %>% dplyr::mutate(sub_Citumor = purrr::map(.x = sub_TumorCellID,.f = function(.x){
    sub_CiTumor <- apply(UMAPembeddings[.x,],2,mean)
    return(sub_CiTumor)
  }))

  CiNormal <- apply(UMAPembeddings[NormalCellID,],2,mean)

  TumorCellIDsi <- purrr::map2(.x = CiTumor$sub_TumorCellID,.y = CiTumor$sub_Citumor,function(.x,.y){
    sub_TumorCellsi <- lapply(.x, function(id){
      pos <- UMAPembeddings[id,]
      rt <- sqrt(sum((pos-.y)^2))
      rn <- sqrt(sum((pos-CiNormal)^2))
      if(rt < 1/3 * rn){return(id)}
    }) %>% unlist()
    return(sub_TumorCellsi)
  }) %>% unlist()
  #TumorCellIDsi <- TumorCellID
  TumorSpotIDsi <- position[TumorCellIDsi,]$spot.ids

  #deal with lonely spots
  TumorSpotIDL <- lapply(TumorSpotIDsi, function(name){ifelse(length(df_j[[name]][df_j[[name]] %in% TumorSpotIDsi]) == 0, name, NA)}) %>% unlist() %>% na.omit()
  TumorCellIDL <- rownames(position[position$spot.ids %in% TumorSpotIDL,])

  TransCellID <- NULL
  TransSpotID <- position[TransCellID,]$spot.ids

  nbrs_of_tumorL <- nbrs(df_j = df_j,TumorSpotIDAdd = TumorSpotIDL,SpotIDRaw = c(TumorSpotIDsi,NormalSpotID,TransSpotID))
  ClusterL <- do.call(rbind,lapply(names(nbrs_of_tumorL),function(spotl){
    sub <- TumorST@meta.data[rownames(position[position$spot.ids == spotl,]),]$seurat_clusters
    nbrs_of_spotl <- rownames(position[position$spot.ids %in% nbrs_of_tumorL[[spotl]],])
    cluster <- do.call(rbind,lapply(nbrs_of_spotl, function(idl){
      pos <- UMAPembeddings[idl,]
      rt <- sqrt(sum((pos-unlist(CiTumor[CiTumor$cluster == sub,]$sub_Citumor))^2))
      rn <- sqrt(sum((pos-CiNormal)^2))
      cluster <- data.frame(CellID = idl, Location= ifelse(rt < 1/3*rn, "Tumor","Trans"))
    }))
    return(cluster)
  }))

  ClusterL <- as.data.frame.array(table(ClusterL$CellID,ClusterL$Location))
  aa_try <- try(
    aa <- as.character(lapply(rownames(ClusterL),function(spotID){ifelse(ClusterL[spotID,]$Tumor > 0,"Tumor","Trans")})),silent = T
  )
  if(!is(aa_try,'try-error')){
    ClusterL$Location <- aa
  }else{
    ClusterL$Location <- rep("Trans",nrow(ClusterL))
  }
  ClusterL <- data.frame(CellID=rownames(ClusterL),Location=as.character(ClusterL$Location))

  #prepare initial IDs for repeat
  TumorCellID <- c(TumorCellIDsi,ClusterL[ClusterL$Location == "Tumor",]$CellID)
  TransCellID <- as.character(ClusterL[ClusterL$Location == "Trans",]$CellID)
  TumorCellIDN <- TumorCellID
  n=1
  TumorSTn <- TumorST
  Clustern <- rbind(data.frame(CellID=NormalCellID,Location=rep("Normal",length(NormalCellID))),
                    data.frame(CellID=TumorCellID,Location=rep("Tumor",length(TumorCellID))),
                    data.frame(CellID=TransCellID,Location=rep("Trans",length(TransCellID))))
  repeat {

    if (length(TumorCellIDN) < 3){
      print(paste("Stop at n = ",(n-1),sep = ""))
      return(TumorSTn)

    }else{

      TransSpotID <- position[TransCellID,]$spot.ids
      TumorSpotIDN <- position[TumorCellIDN,]$spot.ids
      TumorSpotID <- position[TumorCellID,]$spot.ids
      nbrs_of_tumor <- nbrs(df_j = df_j,TumorSpotIDAdd = TumorSpotIDN,SpotIDRaw = c(TumorSpotID,NormalSpotID,TransSpotID))

      if(length(unique(unlist(nbrs_of_tumor))) < 3){
        print(paste("Stop at n = ",(n-1),sep = ""))
        return(TumorSTn)

      }else{

        TumorSTn <- subset(TumorST,cells= c(rownames(position[position$spot.ids %in% unique(unlist(nbrs_of_tumor)),]),TumorCellID,NormalCellID,TransCellID))
        TumorSTn@meta.data$Label <- Clustern$Location[match(rownames(TumorSTn@meta.data),Clustern$CellID)]

        ClusterAdd <- ClusterUpdate(x = n,
                                    position = position,
                                    df_j = df_j,
                                    UMAPembeddings = UMAPembeddings,
                                    TumorCellIDN = TumorCellIDN,
                                    TransCellID = TransCellID,
                                    NormalCellID = NormalCellID,
                                    TumorCellID = TumorCellID)
        Clustern <- rbind(Clustern,ClusterAdd)
        TumorSTn@meta.data$LabelNew <- Clustern$Location[match(rownames(TumorSTn@meta.data),as.character(Clustern$CellID))]
        TumorSTn@meta.data$LabelNew <- factor(TumorSTn@meta.data$LabelNew,levels = c("Normal","Trans","Tumor",paste("Tumor",c(1:n),sep = "")))


        pdf(paste(OutDir,Sample,"_with_nbrs",n,".pdf",sep = ""),width = 7,height = 7)
        p1 <- SpatialDimPlot(TumorSTn,cols = c("#33a02c","#1f78b4",rev(c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'))[1:n]),group.by = "Label")
        print(p1)
        dev.off()

        pdf(paste(OutDir,Sample,"_reduction_with_nbrs",n,".pdf",sep = ""),width = 7,height = 7)
        p2 <- DimPlot(TumorSTn,cols= c("#33a02c","#1f78b4",rev(c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'))[1:n]),group.by = "Label")+
          theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),axis.text = element_blank())+labs(title = NULL)
        print(p2)
        dev.off()

        pdf(paste(OutDir,Sample,"_out_",n,".pdf",sep = ""),width = 7,height = 7)
        p3 <- SpatialDimPlot(TumorSTn,group.by = "LabelNew",cols = c("#33a02c","#1f78b4",rev(c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'))[1:(n+1)]))
        print(p3)
        dev.off()

        pdf(paste(OutDir,Sample,"_reduction_out_",n,".pdf",sep = ""),width = 7,height = 7)
        p4 <- DimPlot(TumorSTn,group.by = "LabelNew",cols = c("#33a02c","#1f78b4",rev(c('#fef0d9','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#990000'))[1:(n+1)]))+
          theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),axis.text = element_blank())+labs(title = NULL)
        print(p4)
        dev.off()
        readr::write_rds(TumorSTn,paste(OutDir,Sample,"_out_",n,".rds.gz",sep = ""),compress = "gz")

        TumorCellID <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% c("Tumor",paste("Tumor",c(1:n),sep = "")),])
        TumorCellIDN <- rownames(TumorSTn@meta.data[TumorSTn@meta.data$LabelNew %in% paste("Tumor",n,sep = ""),])
        TransCellID <- rownames(TumorSTn@meta.data[TumorSTn$LabelNew == "Trans",])

        n <- n+1

      }
    }
    if(n > 5){
      print(paste("Stop at n = ",n,sep = ""))
      return(TumorSTn)

    }

  }
  return(TumorSTn)
}
