
# BoundaryDefine utility


#' Title get neighbors of spot
#'
#' get neighbor spots' id of MalSpotIDAdd
#'
#' @param df_j a list contain spatial spot and its' neighbors
#' @param MalCellIDAdd  barcodes of spatial malignant spots to find neighbors
#' @param CellIDRaw barcodes of spatial spots already defined location
#'
#' @return a list contain spatial spot and its' neighbors
#' @export
#'
nbrs <- function(df_j = df_j, MalCellIDAdd = MalCellIDAdd, CellIDRaw = CellIDRaw) {
  nbrs_of_Mal <- lapply(MalCellIDAdd, function(id) {
    nbs <- df_j[[id]]
    nbs <- nbs[!nbs %in% CellIDRaw]
  })
  names(nbrs_of_Mal) <- MalCellIDAdd
  return(nbrs_of_Mal)
}


#' Title compute_interspots_distance
#'
#' @param position a data frame contain spatial spots' information(spot id, row, cols ,imagerow, imagecol)
#' @param scale.factor scale factor
#'
#' @return a list contain predicted radius
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
#' @param radius predicted radius of spatial spots
#' @param method method to calculate distance
#'
#' @return a list contain of all spatial spots' neighbors
#' @export
#'
find_neighbors <- function(position, radius, method = c("manhattan", "euclidean")) {
  method <- match.arg(method)
  pdist <- as.matrix(stats::dist(as.matrix(position[, c("imagecol", "imagerow")]), method = method))
  neighbors <- (pdist <= radius & pdist > 0) # <= less than or equal to;
  df_j2 <- sapply(seq_len(nrow(position)), function(x) as.vector(names(neighbors[x,])[which(neighbors[x,])]))
  names(df_j2) <- rownames(position)
  return(df_j2)
}


#' Title update cluster data frame
#'
#' @param MalCellIDN barcodes of newlydefined malignant (Mal) spots
#' @param position <data.frame> spatial spots information(spot id, row, cols ,imagerow, imagecol)
#' @param df_j a list contain spatial spots and thier neighbors
#' @param UMAPembeddings <data.frame> spatial spots' umap reduction coordinates
#' @param NormalCellID barcodes of normal spots
#' @param BdyCellID barcodes of boundary (Bdy) spots
#' @param MalCellID barcodes of all malignant spots
#' @param x repeat times
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
                          x = n) {
  P <- lapply(MalCellIDN, function(i) {
    #cellspot <- position[i, ]$spot.ids # center malignant spot of hexagon system
    ncellID <- df_j[[i]] # neighbors of center malignant spot
    cpos <- UMAPembeddings[i, ] # umap coordinate of center malignant spot
    npos <- UMAPembeddings[ncellID, ] # umap coordinate of neighbors of center malignant spot

    nMalID <- ncellID[ncellID %in% MalCellID] # malignant spot id in hexagon system
    CiMal <- data.frame(t(apply(rbind(npos[nMalID, ], cpos), 2, mean))) # Ci of malignant spots in hexagon system
    rMal <- lapply(c(nMalID, i), function(id) { # radius of malignant nbrs of CiMal
      sqrt(sum((rbind(npos, cpos)[id, ] - CiMal)^2))
    }) %>% unlist()

    nBdyID <- ncellID[ncellID %in% BdyCellID]
    nbrsID <- ncellID[!ncellID %in% MalCellID & !ncellID %in% BdyCellID & !ncellID %in% NormalCellID]

    if (length(nBdyID) <= 1) {
      p <- lapply(nbrsID, function(id) {
        sqrt(sum((npos[id, ] - CiMal)^2))
      })
      d <- data.frame(cellID = nbrsID, p1 = unlist(p), p2 = unlist(p))
      d$cluster <- lapply(nbrsID, function(id) {
        ifelse(d[d$cellID == id, ]$p1 <= 0.8 * max(rMal), "Mal", "Bdy")
      }) %>% unlist()
    } else {
      CiBdy <- data.frame(t(apply(npos[nBdyID, ], 2, mean)))
      rBdy <- lapply(nBdyID, function(id) {
        sqrt(sum((npos[id, ] - CiBdy)^2))
      }) %>% unlist()

      p1 <- lapply(nbrsID, function(id) {
        sqrt(sum((npos[id, ] - CiMal)^2))
      })
      p2 <- lapply(nbrsID, function(id) {
        sqrt(sum((npos[id, ] - CiBdy)^2))
      })
      d <- data.frame(cellID = nbrsID, p1 = unlist(p1), p2 = unlist(p2))
      d$cluster <- lapply(nbrsID, function(id) {
        ifelse(d[d$cellID == id, ]$p1 <= 0.8 * max(rMal) & d[d$cellID == id, ]$p2 > max(rBdy), "Mal", "Bdy")
      }) %>% unlist()
    }

    return(d)
  })
  names(P) <- MalCellIDN

  PDF <- do.call(rbind, lapply(P, data.frame))
  rownames(PDF) <- NULL
  Cluster <- as.data.frame.array(table(PDF$cellID, PDF$cluster))
  aa_try <- try(
    aa <- as.character(lapply(rownames(Cluster), function(cellID) {
      ifelse(Cluster[cellID, ]$Mal > 0, paste("Mal", x, sep = ""), "Bdy")
    })),
    silent = T
  )
  if (!is(aa_try, "try-error")) {
    Cluster$Location <- aa
  } else {
    Cluster$Location <- rep("Bdy", nrow(Cluster))
  }
  ClusterNew <- data.frame(CellID = rownames(Cluster), Location = as.character(Cluster$Location))
  return(ClusterNew)
}
