# source('R/spatial_recon_settings.R')
# TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# assay = "Spaital"

#' Title Find differential expressed genes
#'
#' Find differential expressed genes in malignant (Mal) spots, boundary (Bdy) spots, and  non-malignant (nMal) spots
#'
#' @param TumorST A Seurat object with defined  malignant (Mal) spots, boundary (Bdy) spots, and  non-malignant (nMal) spots
#' @param assay The name of assay used to find DEG (Morph: Morphological adjusted gene expression, Spatial: gene expression)
#'
#' @return A list of data frame contain differential expressed genes, p-value, log2FC and adjusted p value in Mal, Bdy, and nMal spots
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
#' assay <- "Spatial"
#' DiffGenes <- FindDiffGenes(TumorST = TumorST, assay = assay)
#'
FindDiffGenes <- function(TumorST = TumorST,
                          assay = c("Morph", "Spatial")) {
  if (assay == "Spatial") {
    TumorST <- NormalizeData(TumorST, assay = assay)
  }

  my.t.test <- function(...) {
    obj <- try(t.test(...), silent = TRUE)
    if (is(obj, "try-error")) {
      return(NA)
    } else {
      return(obj$p.value)
    }
  }

  # prepare gene x cell dataframe
  TumorSTm <- GetAssayData(TumorST, assay = assay, slot = "data")
  # TumorSTm <- apply(TumorSTm,2,function(x){signif(x,digits = 2)})
  TumorSTm <- data.frame(TumorSTm)
  TumorSTm <- TumorSTm %>% tibble::rownames_to_column(., var = "gene")
  colnames(TumorSTm) <- c("gene", as.character(TumorST$Location))

  # prepare majortypes-colnames df
  MajorTypes <- data.frame(
    MajorTypes = data.frame(do.call(rbind, strsplit(colnames(TumorSTm[, 2:ncol(TumorSTm)]), "\\.")))$X1,
    SampleID = colnames(TumorSTm[2:ncol(TumorSTm)])
  )

  # find Trans diff genes vs other
  DiffGenes <- lapply(split(MajorTypes[, "SampleID"], MajorTypes$MajorTypes), function(ID) {
    OtherSampleID <- MajorTypes$SampleID[!(MajorTypes$SampleID) %in% ID]
    DiffGenes <- apply(TumorSTm[, colnames(TumorSTm) %in% MajorTypes$SampleID], 1, function(x) {
      CanExp <- x[ID] %>% as.numeric()
      CanExp <- CanExp[!is.na(CanExp)] %>% as.numeric()
      OtherExp <- x[OtherSampleID] %>% as.numeric()
      OtherExp <- OtherExp[!is.na(OtherExp)] %>% as.numeric()
      c(Diff = mean(CanExp) - mean(OtherExp), pvalue = my.t.test(CanExp, OtherExp))
    }) %>%
      t() %>%
      data.frame() # Welch t-test
    DiffGenes$Symbol <- TumorSTm$gene
    DiffGenes <- DiffGenes[!is.na(DiffGenes$pvalue), ]
    DiffGenes$FDR <- p.adjust(DiffGenes$pvalue, method = "fdr")
    return(DiffGenes)
  })
  return(DiffGenes)
}
