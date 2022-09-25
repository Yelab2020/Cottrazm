# source('R/spatial_recon_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# DiffGenes <- FindDiffGenes(TumorST = TumorST,assay = "Spatial")

#' Title Differential expressed genes volcano plot
#'
#' Draw a volcano plot of features in boundary defined ST data
#'
#' @param DiffGenes A list of data frame contain differential expressed genes, p value, log2FC and adjusted p value in Mal, Bdy and nMal spots
#' @param Location Name in location defined ST data (Mal, Bdy, and nMal)
#' @param cut_off_pvalue Cut off of p value, features p value < 10^(-cut_off_pvalue) will be selected as significant
#' @param cut_off_logFC Cut off of log2FC, features log2FC > cut_off_logFC will be selected as significant
#' @param n Number of features labeled in volcano plot, default is 10
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
#' DiffGenes <- FindDiffGenes(TumorST = TumorST, assay = "Spatial")
#' DiffVolcanoplot(DiffGenes = DiffGenes, Location = "Bdy", cut_off_pvalue = 2, cut_off_logFC = 0.25, n = 10)
#'
DiffVolcanoplot <- function(DiffGenes = DiffGenes,
                            Location = c("Mal", "Bdy", "nMal"),
                            cut_off_pvalue = 2,
                            cut_off_logFC = 0.25,
                            n = NULL) {
  if (is.null(n) == TRUE) {
    n <- 10
  }

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
}
