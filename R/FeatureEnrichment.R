# source('R/recon_spatial_TME_settings.R')
# TumorST <- readr::read_rds("YourPath/Tumorboundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
# DiffGenes <- FindDiffGenes(TumorST = TumorST,assay = "Spatial")

#' Title Enrichment analysis of features of boundary location
#'
#' Use differentially expressed genes of boundary location to find the KEGG and GO pathway these locations DEG enriched.
#'
#' @param DiffGenes A list of dataframe contain differential expressed genes and its p value, log2FC and adjusted p value in tumor, boundary/transition zone and outspots
#' @param cut_off_logFC Cut off of log2FC, features log2FC > cut_off_logFC will be selected as significant
#' @param cut_off_pvalue Cut off of p value, features p value < cut_off_pvalue will be selected as significant
#' @param Location Names of boundary Location to analyse feature enrichment.
#'
#' @return A tibble data frame contain KEGG and GO enrichment result
#' @export
#'
#' @examples
#' TumorST <- readr::read_rds("YourPath/Tumorboundary/1.BoundaryDefine/CRC1/TumorSTBoundaryDefine.rds.gz")
#' DiffGenes <- FindDiffGenes(TumorST = TumorST,assay = "Spatial")
#' TransFeatureEnrich <- FeatureEnrichment(DiffGenes = DiffGenes, cut_off_logFC = 0.25, cut_off_pvalue = 0.05, Location = "Trans")
#'
FeatureEnrichment <- function(DiffGenes = DiffGenes,
                               cut_off_logFC = 0.25,
                               cut_off_pvalue = 0.05,
                               Location = c("Trans","Tumor","OutSpots")){

  #get DiffGeneSig
  DiffGenesSig <- lapply(DiffGenes,function(sub){
    sub <- sub[sub$Diff >= cut_off_logFC  & sub$FDR <= cut_off_pvalue ,]
    sub <- sub[sub$Symbol %in% grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL",sub$Symbol,invert = T,value = T),]
    sub <- sub[!is.na(sub$Diff),]
  })

  Enrichment <- tibble::tibble(Location = Location)
  Enrichment <- Enrichment %>% dplyr::mutate(LocationDiffFeatures =
                                               purrr::map(.x = Location,.f = function(.x){DiffGenesSig[[.x]]$Symbol}))

  #Goenrichment
  Enrichment <- Enrichment %>% dplyr::mutate(GO = purrr::map(.x = LocationDiffFeatures,.f = function(.x){
    GO <- clusterProfiler::enrichGO(gene= .x,
                                    keyType = "SYMBOL",
                                    OrgDb         = "org.Hs.eg.db",
                                    ont           = "BP",
                                    pAdjustMethod = "fdr",
                                    pvalueCutoff  = 0.2)
    return(GO)
  }))

  #KEGGenrichment
  Enrichment <- Enrichment %>% dplyr::mutate(KEGG = purrr::map(.x = LocationDiffFeatures,.f = function(.x){
    geneL <-  clusterProfiler::bitr(.x,fromType = "SYMBOL",
                                    toType = c("ENTREZID", "ENSEMBL"), OrgDb = "org.Hs.eg.db")

    KEGG <- clusterProfiler::enrichKEGG(gene = geneL$ENTREZID,
                                        organism = "hsa",
                                        keyType = "kegg",
                                        pAdjustMethod = "BH",
                                        minGSSize  = 3,
                                        pvalueCutoff = 0.2,
                                        qvalueCutoff  = 0.2)
    return(KEGG)
  }))

  return(Enrichment)
}
