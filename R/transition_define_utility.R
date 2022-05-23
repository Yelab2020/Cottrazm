
#TransitionDefine utility


#' Title get neighbors of spot
#'
#' get neighbor spots' id of TumorSpotIDAdd
#'
#' @param df_j <list> contain spaital spot and its' neighbors
#' @param TumorSpotIDAdd <character> spatial spot id to find neighbors
#' @param SpotIDRaw <character> spatial spot id already defined
#'
#' @return <list> contain spatial spot and its' neighbors
#' @export
#'
nbrs <- function(df_j = df_j,TumorSpotIDAdd = TumorSpotIDAdd,SpotIDRaw = SpotIDRaw){
  nbrs_of_Tumor <- lapply(TumorSpotIDAdd,function(id){
    nbs = df_j[[id]]
    nbs = nbs[!nbs %in% SpotIDRaw]
  })
  names(nbrs_of_Tumor) <- TumorSpotIDAdd
  return(nbrs_of_Tumor)
}


#' Title compute_interspots_distance
#'
#' @param position <data.frame> spatial spots information(spot id, row, cols ,imagerow, imagecol)
#' @param scale.factor <numeric> scale factor
#'
#' @return
#' @export
#'
compute_interspot_distances <- function(position, scale.factor = 1.05) {
  cols <- c("row", "col", "imagerow", "imagecol")
  assertthat::assert_that(all(cols %in% colnames(position)))
  dists <- list()
  dists$xdist <- coef(lm(position$imagecol ~ position$col))[2]
  dists$ydist <- coef(lm(position$imagerow ~ position$row))[2]
  dists$radius <- (abs(dists$xdist) + abs(dists$ydist)) * scale.factor


  dists
}

#' Title find neighbors
#'
#' @param position <data.frame> spatial spots information(spot id, row, cols ,imagerow, imagecol)
#' @param radius <numeric> in dists
#' @param method <charactor> method to calculate distance
#'
#' @return <list> of all spatial spots' neighbors
#' @export
#'
find_neighbors <- function(position, radius,method = c("manhattan", "euclidean")) {
  method <- match.arg(method)
  pdist <- as.matrix(stats::dist(as.matrix(position[,c("imagecol","imagerow")]), method=method))
  neighbors <- (pdist <= radius & pdist > 0)  # <= less than or equal to;
  df_j2 <- sapply(seq_len(nrow(position)),function(x) as.vector(which(neighbors[x, ])))
  df_j2
}


#' Title updat cluster dataframe update
#'
#' @param TumorCellIDN <character> new defined tumor spots ids
#' @param position <data.frame> spatial spots information(spot id, row, cols ,imagerow, imagecol)
#' @param df_j <list> contain spaital spot and its' neighbors
#' @param UMAPembeddings <data.frame> spatial spots umap reduction coordinates
#' @param NormalCellID <character> normal spots ids
#' @param TransCellID <character> trans spots ids
#' @param TumorCellID <characoter> all tumor spots ids
#' @param x <interger> repeat times
#'
#' @return <data.frame> neighbors of new defined tumor spots ids and their define types
#' @export
#'
ClusterUpdate <- function(TumorCellIDN = TumorCellIDN,
                          position = position,
                          df_j = df_j,
                          UMAPembeddings = UMAPembeddings,
                          NormalCellID = NormalCellID,
                          TransCellID = TransCellID,
                          TumorCellID = TumorCellID,
                          x=n){

  P <- lapply(TumorCellIDN,function(i){
    cellspot <- position[i,]$spot.ids #tumor center
    ncellID <- rownames(position[position$spot.ids %in% df_j[[cellspot]],]) #neighbors of tumor center
    cpos <- UMAPembeddings[i,] #umap embedding of tumor center
    npos <- UMAPembeddings[ncellID,] #umap embedding of tumor neighbor

    ntumorID <- ncellID[ncellID %in% TumorCellID] #tumor cellid in 7spots systerm
    Citumor <- data.frame(t(apply(rbind(npos[ntumorID,],cpos),2,mean))) #Ci of tumor in 7spots system
    rtumor <- lapply(c(ntumorID,i), function(id){ #radius of nbrs of tumor to Ci of tumor
      sqrt(sum((rbind(npos,cpos)[id,]-Citumor)^2))
    }) %>% unlist()

    ntransID <- ncellID[ncellID %in% TransCellID]
    nbrsID <- ncellID[!ncellID %in% TumorCellID & !ncellID %in% TransCellID & !ncellID %in% NormalCellID]

    if (length(ntransID) <= 1){
      p <- lapply(nbrsID,function(id){
        sqrt(sum((npos[id,]-Citumor)^2))
      })
      d <- data.frame(cellID=nbrsID,p1= unlist(p),p2 = unlist(p))
      d$cluster <- lapply(nbrsID,function(id){
        ifelse(d[d$cellID == id,]$p1 <= 0.8*max(rtumor), "Tumor","Trans")
      }) %>% unlist()

    }else{
      Citrans <- data.frame(t(apply(npos[ntransID,],2,mean)))
      rtrans <- lapply(ntransID,function(id){
        sqrt(sum((npos[id,]-Citrans)^2))
      }) %>% unlist()

      p1 <- lapply(nbrsID,function(id){
        sqrt(sum((npos[id,]-Citumor)^2))
      })
      p2 <- lapply(nbrsID, function(id){
        sqrt(sum((npos[id,]-Citrans)^2))
      })
      d <- data.frame(cellID=nbrsID,p1= unlist(p1),p2=unlist(p2))
      d$cluster <- lapply(nbrsID,function(id){
        ifelse(d[d$cellID == id,]$p1 <= 0.8*max(rtumor) & d[d$cellID == id,]$p2 > max(rtrans), "Tumor","Trans")
      }) %>% unlist()
    }

    return(d)
  })
  names(P) <- TumorCellIDN

  PDF <- do.call(rbind,lapply(P, data.frame))
  rownames(PDF) <- NULL
  Cluster <- as.data.frame.array(table(PDF$cellID,PDF$cluster))
  aa_try <- try(
    aa <- as.character(lapply(rownames(Cluster),function(spotID){ifelse(Cluster[spotID,]$Tumor > 0,paste("Tumor",x,sep = ""),"Trans")})),silent = T
  )
  if(!is(aa_try,'try-error')){
    Cluster$Location <- aa
  }else{
    Cluster$Location <- rep("Trans",nrow(Cluster))
  }
  ClusterNew <- data.frame(CellID=rownames(Cluster),Location=as.character(Cluster$Location))
  return(ClusterNew)
}

