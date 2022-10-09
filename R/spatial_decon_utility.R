
# title: optimize_deconvolute_dwls
optimize_deconvolute_dwls <- function(exp,
                                      Signature) {
  ###### overlap signature with spatial genes
  Genes <- intersect(rownames(Signature), rownames(as.matrix(exp)))
  S <- Signature[Genes, ]
  S <- Matrix::as.matrix(S)
  Bulk <- Matrix::as.matrix(exp)
  subBulk <- as.matrix(Bulk[Genes, ])
  allCounts_DWLS <- NULL
  all_exp <- rowMeans(as.matrix(exp))
  solution_all_exp <- solve_OLS_internal(S, all_exp[Genes])
  constant_J <- find_dampening_constant(S, all_exp[Genes], solution_all_exp)
  # print(constant_J)
  for (j in 1:(dim(subBulk)[2])) {
    B <- subBulk[, j]
    solDWLS <- optimize_solveDampenedWLS(S, B, constant_J)
    allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
  }
  colnames(allCounts_DWLS) <- colnames(exp)
  return(allCounts_DWLS)
}

# title: solve_OLS_internal
solve_OLS_internal <- function(S,
                               B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")

  pd_D_mat <- nearPD(D / sc)
  solution <- quadprog::solve.QP(as.matrix(pd_D_mat$mat), d / sc, A, bzero)$solution
  #solution <- quadprog::solve.QP(D, d, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# title: find_dampening_constant
find_dampening_constant <- function(S,
                                    B,
                                    goldStandard) {
  solutionsSd <- NULL
  # goldStandard is used to define the weights
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsScaledMinusInf <- wsScaled
  # ignore infinite weights
  if (max(wsScaled) == "Inf") {
    wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
  }
  # try multiple values of the dampening constant (multiplier)
  # for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
    multiplier <- 1 * 2^(j - 1)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    solutions <- NULL
    seeds <- c(1:100)
    for (i in 1:100) {
      set.seed(seeds[i]) # make nondeterministic
      subset <- sample(length(ws), size = length(ws) * 0.5) # randomly select half of gene set
      # solve dampened weighted least squares for subset
      fit <- stats::lm(B[subset] ~ -1 + S[subset, ], weights = wsDampened[subset])
      sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
      solutions <- cbind(solutions, sol)
    }
    solutionsSd <- cbind(solutionsSd, apply(solutions, 1, stats::sd))
  }
  # choose dampening constant that results in least cross-validation variance
  j <- which.min(colMeans(solutionsSd^2))
  return(j)
}


# title: optimize_solveDampenedWLS
optimize_solveDampenedWLS <- function(S,
                                      B,
                                      constant_J) {
  # first solve OLS, use this solution to find a starting point for the weights
  solution <- solve_OLS_internal(S, B)
  # now use dampened WLS, iterate weights until convergence
  iterations <- 0
  changes <- c()
  # find dampening constant for weights using cross-validation
  j <- constant_J
  change <- 1
  while (change > .01 & iterations < 1000) {
    newsolution <- solve_dampened_WLSj(S, B, solution, j)
    # decrease step size for convergence
    solutionAverage <- rowMeans(cbind(newsolution, matrix(solution, nrow = length(solution), ncol = 4)))
    change <- norm(Matrix::as.matrix(solutionAverage - solution))
    solution <- solutionAverage
    iterations <- iterations + 1
    changes <- c(changes, change)
  }
  # print(round(solution/sum(solution),5))
  return(solution / sum(solution))
}

# title: solve_dampened_WLSj
solve_dampened_WLSj <- function(S,
                                B,
                                goldStandard,
                                j) {
  multiplier <- 1 * 2^(j - 1)
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsDampened <- wsScaled
  wsDampened[which(wsScaled > multiplier)] <- multiplier
  W <- diag(wsDampened)
  D <- t(S) %*% W %*% S
  d <- t(S) %*% W %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")

  pd_D_mat <- nearPD(D / sc)
  solution <- quadprog::solve.QP(as.matrix(pd_D_mat$mat), d / sc, A, bzero)$solution
  #solution <- quadprog::solve.QP(D / sc, d / sc, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# title: spot_propotion_initial
spot_proportion_initial <- function(enrich_matrix = enrich_matrix,
                                    enrich_result = enrich_result,
                                    filter_expr = filter_expr,
                                    filter_sig = filter_sig,
                                    clustermarkers_list = clustermarkers_list,
                                    meta_data = meta_data,
                                    malignant_cluster = "Malignant epithelial cells",
                                    tissue_cluster = "Epithelial cells",
                                    stromal_cluster = "Fibroblast cells") {

  #### select  cell types####
  # select topic from module score
  topic <- meta_data$Decon_topics
  topic_sort <- sort(unique(topic))

  # signature score cutoff and signature score selected cell types
  ct_sig0 <- lapply(names(clustermarkers_list), function(cluster) {
    cutoff_sig <- quantile(meta_data[, cluster], 0.75)
    test <- as.data.frame(do.call(rbind, lapply(topic_sort, function(topic) {
      topic_qun <- median(meta_data[meta_data$Decon_topics == topic, ][, cluster])
    })))

    ct <- as.character(topic_sort[which(test - cutoff_sig > 0)])
    return(ct)
  })

  names(ct_sig0) <- names(clustermarkers_list)

  # calculate enrichment scroe cutoff
  cutoff_enrich <- lapply(rownames(enrich_result), function(cluster) {
    cluster_enrich <- enrich_result[cluster, ]
    cluster_enrichF <- quantile(cluster_enrich[cluster_enrich > 0], 0.75)
    return(cluster_enrichF)
  }) %>% unlist()
  names(cutoff_enrich) <- rownames(enrich_result)

  # create dwls_init matrix
  dwls_results <- matrix(0, nrow = dim(enrich_matrix)[2], ncol = dim(filter_expr)[2])
  rownames(dwls_results) <- colnames(enrich_matrix)
  colnames(dwls_results) <- colnames(filter_expr)

  # start select cell type of each decon-topic and initial deconvolution
  for (i in 1:length(topic_sort)) {

    # signature score ct
    ct_sig <- c()
    for (cluster in names(clustermarkers_list)) {
      if (topic_sort[i] %in% ct_sig0[[cluster]]) {
        ct_sig <- c(ct_sig, cluster)
      }
    }

    # split enrich result of topic i
    cluster_i_enrich <- as.matrix(enrich_result[, which(topic == topic_sort[i])]) # enrichment result of topic i

    # select top 2 cluster
    row_i_max <- Rfast::rowMaxs(cluster_i_enrich, value = TRUE)
    names(row_i_max) <- rownames(cluster_i_enrich)
    ct_sort <- sort(row_i_max[which(row_i_max - cutoff_enrich > 0)], decreasing = T)[2]
    ct_enrich <- names(row_i_max[which(row_i_max - cutoff_enrich > 0)])[which(row_i_max[which(row_i_max - cutoff_enrich > 0)] >= ct_sort)]


    if (malignant_cluster %in% ct_sig & tissue_cluster %in% ct_sig) {
      ct_sig <- ct_sig[ct_sig != malignant_cluster]
    } # signature score can't tell malignant form epithelial

    # combine ct_sig & ct_enrich
    if (length(ct_sig) == 0) {
      ct <- ct_enrich
    } else {
      if (length(ct_sig) == 1) {
        if (ct_sig %in% ct_enrich) {
          ct <- ct_enrich
        } else {
          ct <- c(ct_sig, ct_enrich[1])
        }
      } else {
        ct <- unique(c(ct_sig, ct_enrich))
      }
    }

    if (length(ct) == 0) {
      ct <- names(sort(row_i_max - cutoff_enrich, decreasing = T))[1:2]
    }

    # in case tissue/malignant cluster cover change/expression of other cluster (except for low feature clusters)
    if (malignant_cluster %in% ct | tissue_cluster %in% ct) {
      ct_add <- names(sort(row_i_max - cutoff_enrich, decreasing = T)[1:3])
      ct <- unique(c(ct, ct_add))
    }

    if (unlist(strsplit(topic_sort[i], "_"))[1] == "Mal") {
      ct <- unique(c(ct, malignant_cluster))
    }
    if (median(meta_data[meta_data$Decon_topics == topic_sort[i], ]$nCount_Spatial) < 5000) {
      ct <- unique(c(ct, stromal_cluster))
    }

    # na.omit and print-check
    ct <- as.character(na.omit(ct))
    print(as.character(topic_sort[i]))
    print(ct)

    ## select cluster unique gene
    ct_gene <- c()
    for (j in 1:length(ct)) {
      sig_gene_j <- rownames(enrich_matrix)[which(enrich_matrix[, ct[j]] == 1)]
      ct_gene <- c(ct_gene, sig_gene_j)
    }

    # intersect sc cluster unique gene with st filter expr
    uniq_ct_gene <- intersect(unique(ct_gene), rownames(filter_expr))
    select_enrichig_exp <- as.matrix(filter_sig[uniq_ct_gene, ct]) # select sc cluster mean ref
    cluster_i_cell <- which(topic == topic_sort[i])
    cluster_cell_exp <- as.matrix(filter_expr[uniq_ct_gene, cluster_i_cell]) # spots in topic-i; genes from sc contribute to topic-i

    cluster_i_dwls <- optimize_deconvolute_dwls(cluster_cell_exp, select_enrichig_exp)
    dwls_results[ct, cluster_i_cell] <- cluster_i_dwls
  }

  ##### remove negative values
  for (i in dim(dwls_results)[1]) {
    negtive_index <- which(dwls_results[i, ] < 0)
    dwls_results[i, negtive_index] <- 0
  }

  return(dwls_results)
}


# title: spot_deconvolution
spot_deconvolution <- function(expr,
                               meta_data,
                               ct_exp,
                               enrich_matrix,
                               binary_matrix) {

  # topic information
  topic <- meta_data$Decon_topics
  topic_sort <- sort(unique(topic))
  #### initialize dwls matrix
  dwls_results <- matrix(0, nrow = dim(ct_exp)[2], ncol = dim(expr)[2])
  rownames(dwls_results) <- colnames(ct_exp)
  colnames(dwls_results) <- colnames(expr)
  # print(binary_matrix)
  for (i in 1:length(topic_sort)) {
    cluster_i_matrix <- as.matrix(binary_matrix[, which(topic == topic_sort[i])])
    row_i_max <- Rfast::rowMaxs(cluster_i_matrix, value = TRUE)
    ct_i <- rownames(cluster_i_matrix)[which(row_i_max == 1)]
    ######## calculate proportion based on binarized deconvolution results at first step
    if (length(ct_i) == 1) {
      dwls_results[ct_i[1], which(topic == topic_sort[i])] <- 1
    } else {
      ct_gene <- c()
      for (j in 1:length(ct_i)) {
        sig_gene_j <- rownames(enrich_matrix)[which(enrich_matrix[, ct_i[j]] == 1)]
        ct_gene <- c(ct_gene, sig_gene_j)
      }
      uniq_ct_gene <- intersect(rownames(expr), unique(ct_gene))
      select_enrichig_exp <- ct_exp[uniq_ct_gene, ct_i]
      cluster_i_cell <- which(topic == topic_sort[i])
      cluster_cell_exp <- as.matrix(expr[uniq_ct_gene, cluster_i_cell])
      colnames(cluster_cell_exp) <- colnames(expr)[cluster_i_cell]
      ###### calculate
      ###### overlap signature with spatial genes
      all_exp <- rowMeans(cluster_cell_exp)
      solution_all_exp <- solve_OLS_internal(select_enrichig_exp, all_exp)
      constant_J <- find_dampening_constant(select_enrichig_exp, all_exp, solution_all_exp)
      ###### deconvolution for each spot
      for (k in 1:(dim(cluster_cell_exp)[2])) {
        B <- Matrix::as.matrix(cluster_cell_exp[, k])
        ct_enrichpot_k <- rownames(cluster_i_matrix)[which(cluster_i_matrix[, k] == 1)]
        if (length(ct_enrichpot_k) == 1) {
          dwls_results[ct_enrichpot_k[1], colnames(cluster_cell_exp)[k]] <- 1
        } else {
          ct_k_gene <- c()
          for (m in 1:length(ct_enrichpot_k)) {
            sig_gene_k <- rownames(enrich_matrix)[which(enrich_matrix[, ct_enrichpot_k[m]] == 1)]
            ct_k_gene <- c(ct_k_gene, sig_gene_k)
          }
          uniq_ct_k_gene <- intersect(rownames(ct_exp), unique(ct_k_gene))
          S_k <- Matrix::as.matrix(ct_exp[uniq_ct_k_gene, ct_enrichpot_k])
          solDWLS <- optimize_solveDampenedWLS(S_k, B[uniq_ct_k_gene, ], constant_J)
          dwls_results[names(solDWLS), colnames(cluster_cell_exp)[k]] <- solDWLS
        }
      }
    }
  }
  ##### remove negative values
  for (i in dim(dwls_results)[1]) {
    negtive_index <- which(dwls_results[i, ] < 0)
    dwls_results[i, negtive_index] <- 0
  }
  return(dwls_results)
}
