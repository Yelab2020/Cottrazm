

# title:get_feature_weight
# function: from filter_sig get feature weight of each cell type
get_feature_weight <- function(filter_sig = filter_sig) {
  W_st <- do.call(rbind, lapply(rownames(filter_sig), function(feature) {
    a <- filter_sig[feature, ]
    scale_a <- a / sum(a)
    return(scale_a)
  }))
  rownames(W_st) <- rownames(filter_sig)

  return(W_st)
}


# title: get_obs
# function: get each Boundary spot deconvolution data
get_obs <- function(DeconData = DeconData, obs_ID = SubID) {
  obs <- tibble::tibble(cell_ID = DeconData[DeconData$cell_ID %in% obs_ID, ]$cell_ID) %>%
    dplyr::mutate(Decon = purrr::map(.x = cell_ID, .f = function(.x) {
      sub_Decon <- DeconData[DeconData$cell_ID == .x, ][, 2:ncol(DeconData)]
      sub_DeconF <- as.data.frame(sub_Decon[, colnames(sub_Decon)[which(sub_Decon > 0)]])
      if (ncol(sub_DeconF) == 1) {
        colnames(sub_DeconF) <- colnames(sub_Decon)[which(sub_Decon > 0)]
      }
      return(sub_DeconF)
    }))

  return(obs)
}



# title:get_sub_mtx
# function: from Bdy_expr to sub_Bdy_expr
get_recon_mtx <- function(TumorST = TumorST,
                          sig_exp = sig_exp,
                          clustermarkers_list = clustermarkers_list,
                          DeconData = DeconData,
                          Location = Location) {
  star_time <- Sys.time()

  # get st expr
  SubID <- rownames(TumorST@meta.data[TumorST@meta.data$Location %in% Location, ])
  expr_values <- as.matrix(TumorST@assays$Spatial@counts)
  sub_expr <- expr_values[, SubID]

  # get intersect genes
  ClusterMarkers <- unlist(clustermarkers_list)[grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", unlist(clustermarkers_list), invert = T)]
  intersect_genes <- intersect(rownames(expr_values), unique(ClusterMarkers))

  # filter st-expr and sig
  filter_sig <- sig_exp[intersect(rownames(sig_exp), intersect_genes), ]
  expr_values_filter <- expr_values[intersect(rownames(sig_exp), intersect_genes), ]
  sub_expr_filter <- sub_expr[intersect(rownames(sig_exp), intersect_genes), ]

  # get w_st
  W_st <- get_feature_weight(filter_sig = filter_sig)

  # get obs data
  obs <- get_obs(DeconData = DeconData, obs_ID = SubID)

  # subset mtx
  mtx <- c()
  pb <- txtProgressBar(style = 3)
  for (i in 1:length(obs$cell_ID)) {
    spot <- obs$cell_ID[i]
    # get fraction of sub_matrix
    decon <- as.data.frame(t(unlist(obs[obs$cell_ID == spot, ]$Decon)))
    sub_matrix <- matrix(0, nrow = nrow(sub_expr_filter), ncol = ncol(decon))
    colnames(sub_matrix) <- paste(colnames(decon), spot, sep = "_")
    rownames(sub_matrix) <- rownames(sub_expr_filter)

    # calculate feature
    for (feature in rownames(sub_matrix)) {
      sub_wt <- W_st[feature, colnames(decon)]
      deno <- sum(sub_wt * decon)

      if (deno == 0) {
        feature_modi <- rep(sub_expr_filter[feature, spot], ncol(sub_matrix))
      } else {
        feature_modi <- sub_expr_filter[feature, spot] / deno * sub_wt
      }
      sub_matrix[feature, ] <- feature_modi
    }

    mtx <- cbind(mtx, sub_matrix)
    # print(i)
    setTxtProgressBar(pb, i / length(obs$cell_ID))
  }

  end_time <- Sys.time()
  close(pb)
  run_time <- end_time - star_time

  print(sprintf("Time to get reconstructed matrix was %shours", run_time))

  return(mtx)
}
