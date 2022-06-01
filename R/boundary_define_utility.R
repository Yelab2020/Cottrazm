
#BoundaryDefine utility


#' Title get neighbors of spot
#'
#' get neighbor spots' id of MalSpotIDAdd
#'
#' @param df_j <list> contain spatial spot and its' neighbors
#' @param MalSpotIDAdd <character> spatial malignant spots' id to find neighbors
#' @param SpotIDRaw <character> spatial spots' id already defined
#'
#' @return <list> contain spatial spot and its' neighbors
#' @export
#'
nbrs <- function(df_j = df_j,MalSpotIDAdd = MalSpotIDAdd,SpotIDRaw = SpotIDRaw){
  nbrs_of_Mal <- lapply(MalSpotIDAdd,function(id){
    nbs = df_j[[id]]
    nbs = nbs[!nbs %in% SpotIDRaw]
  })
  names(nbrs_of_Mal) <- MalSpotIDAdd
  return(nbrs_of_Mal)
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
#' @param method <character> method to calculate distance
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


#' Title update cluster data frame
#'
#' @param MalCellIDN <character> new defined malignant (Mal) spots' ids
#' @param position <data.frame> spatial spots information(spot id, row, cols ,imagerow, imagecol)
#' @param df_j <list> contain spaital spot and its' neighbors
#' @param UMAPembeddings <data.frame> spatial spots' umap reduction coordinates
#' @param NormalCellID <character> normal spots' ids
#' @param BdyCellID <character> boundary (Bdy) spots' ids
#' @param MalCellID <character> all malignant spots' ids
#' @param x <integer> repeat times
#'
#' @return <data.frame> neighbors of new defined malignant spots' ids and their define types
#' @export
#'
ClusterUpdate <- function(MalCellIDN = MalCellIDN,
                          position = position,
                          df_j = df_j,
                          UMAPembeddings = UMAPembeddings,
                          NormalCellID = NormalCellID,
                          BdyCellID = BdyCellID,
                          MalCellID = MalCellID,
                          x=n){

  P <- lapply(MalCellIDN,function(i){
    cellspot <- position[i,]$spot.ids #center malignant spot of hexagon system
    ncellID <- rownames(position[position$spot.ids %in% df_j[[cellspot]],]) #neighbors of center malignant spot
    cpos <- UMAPembeddings[i,] #umap coordinate of center malignant spot
    npos <- UMAPembeddings[ncellID,] #umap coordinate of neighbors of center malignant spot

    nMalID <- ncellID[ncellID %in% MalCellID] #malignant spot id in hexagon system
    CiMal <- data.frame(t(apply(rbind(npos[nMalID,],cpos),2,mean))) #Ci of malignant spots in hexagon system
    rMal <- lapply(c(nMalID,i), function(id){ #radius of malignant nbrs of CiMal
      sqrt(sum((rbind(npos,cpos)[id,]-CiMal)^2))
    }) %>% unlist()

    nBdyID <- ncellID[ncellID %in% BdyCellID]
    nbrsID <- ncellID[!ncellID %in% MalCellID & !ncellID %in% BdyCellID & !ncellID %in% NormalCellID]

    if (length(nBdyID) <= 1){
      p <- lapply(nbrsID,function(id){
        sqrt(sum((npos[id,]-CiMal)^2))
      })
      d <- data.frame(cellID=nbrsID,p1= unlist(p),p2 = unlist(p))
      d$cluster <- lapply(nbrsID,function(id){
        ifelse(d[d$cellID == id,]$p1 <= 0.8*max(rMal), "Mal","Bdy")
      }) %>% unlist()

    }else{
      CiBdy <- data.frame(t(apply(npos[nBdyID,],2,mean)))
      rBdy <- lapply(nBdyID,function(id){
        sqrt(sum((npos[id,]-CiBdy)^2))
      }) %>% unlist()

      p1 <- lapply(nbrsID,function(id){
        sqrt(sum((npos[id,]-CiMal)^2))
      })
      p2 <- lapply(nbrsID, function(id){
        sqrt(sum((npos[id,]-CiBdy)^2))
      })
      d <- data.frame(cellID=nbrsID,p1= unlist(p1),p2=unlist(p2))
      d$cluster <- lapply(nbrsID,function(id){
        ifelse(d[d$cellID == id,]$p1 <= 0.8*max(rMal) & d[d$cellID == id,]$p2 > max(rBdy), "Mal","Bdy")
      }) %>% unlist()
    }

    return(d)
  })
  names(P) <- MalCellIDN

  PDF <- do.call(rbind,lapply(P, data.frame))
  rownames(PDF) <- NULL
  Cluster <- as.data.frame.array(table(PDF$cellID,PDF$cluster))
  aa_try <- try(
    aa <- as.character(lapply(rownames(Cluster),function(spotID){ifelse(Cluster[spotID,]$Mal > 0,paste("Mal",x,sep = ""),"Bdy")})),silent = T
  )
  if(!is(aa_try,'try-error')){
    Cluster$Location <- aa
  }else{
    Cluster$Location <- rep("Bdy",nrow(Cluster))
  }
  ClusterNew <- data.frame(CellID=rownames(Cluster),Location=as.character(Cluster$Location))
  return(ClusterNew)
}

