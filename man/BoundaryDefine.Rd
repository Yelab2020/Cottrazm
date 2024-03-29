% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BoundaryDefine.R
\name{BoundaryDefine}
\alias{BoundaryDefine}
\title{Title Boundary define}
\usage{
BoundaryDefine(
  TumorST = TumorST,
  MalLabel = NULL,
  OutDir = NULL,
  Sample = Sample
)
}
\arguments{
\item{TumorST}{A Seurat object with morphological adjusted expression matrix and determined clusters}

\item{MalLabel}{CNV labels with high cnv scores which could be defined as malignant label. If NULL, Cottrazm will take the two CNVLabels with the highest average CNVScores.}

\item{OutDir}{Path to file save figures and processed data}

\item{Sample}{Name of your sample}
}
\value{
A subset Seurat object of ST data with defined malignant spots (Mal), boundary spots (Bdy), and non-malignant spots (nMal)
}
\description{
Find boundary of tumor in ST data
}
\examples{
TumorST <- readr::read_rds("YourPath/TumorBoundary/1.BoundaryDefine/CRC1/TumorSTClustered.rds.gz")
Sample <- "CRC1"
OutDir <- "YourPath/TumorBoundary/Fig/1.BoundaryDefine/CRC1/"
TumorST <- STCNVScore(TumorST = TumorST, assay = "Spatial", OutDir = OutDir, Sample = Sample)
MalLabel <- c(1, 2)
TumorSTn <- BoundaryDefine(TumorST = TumorST, MalLabel = MalLabel, OutDir = OutDir, Sample = Sample)

}
