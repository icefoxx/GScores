% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calScores.R
\name{run_JASMINE}
\alias{run_JASMINE}
\title{Run JASMINE scoring on a Seurat object}
\usage{
run_JASMINE(
  object,
  genesets,
  assay = "RNA",
  slot = "counts",
  method = "likelihood",
  methods = "JASMINE"
)
}
\arguments{
\item{object}{A Seurat object containing the single-cell data.}

\item{genesets}{A list of gene sets to be used for JASMINE scoring}

\item{method}{'oddsratio' or 'likelihood'.}

\item{ncores}{Integer specifying the number of cores to be used for parallel computation. Default is 1.}
}
\value{
The input Seurat object with the JASMINE scores added as a new assay called "JASMINE".
}
\description{
This function performs JASMINE scoring on a given Seurat object using the specified gene sets.
It calculates the JASMINE scores and stores them in a new assay called "JASMINE".
}
\author{
Qiong Zhang
}
