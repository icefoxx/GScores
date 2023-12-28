#' Extract expression matrix from Seurat object or matrix
#'
#' This function extracts the expression matrix from either a Seurat object or a matrix object. If a Seurat object is provided, the function can optionally specify the assay and slot to be used for extracting the expression matrix. If a matrix is provided, the function will create a Seurat object with the provided matrix and default metadata fields. The function can also specify the minimum number of cells and features required to keep genes in the expression matrix.
#'
#' @name assay2mtx
#' @author Qiong Zhang
#' @param object A Seurat object or matrix containing expression data
#' @param assay A character string specifying the assay name to extract expression data from the Seurat object
#' @param slot A character string specifying the slot name to extract expression data from the Seurat object
#' @param min.cells An integer specifying the minimum number of cells required to keep genes in the expression matrix
#' @param min.features An integer specifying the minimum number of features required to keep genes in the expression matrix
#' @return A matrix containing the expression data, with gene names as row names and cell barcodes as column names
#' @examples
#' expmtx <- assay2mtx(object = seurat_object, assay = "RNA", slot = "counts", min.cells = 10, min.features = 5)
#' expmtx <- assay2mtx(object = expression_matrix, min.cells = 10, min.features = 5)
#'
assay2mtx <- function(object, assay = NULL, slot = "counts", min.cells = 3, min.features = 1, update = F) {
  if (any(methods::is(object) %in% "Seurat")) {
    assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
    object <- if (isTRUE(update)) Seurat::UpdateSeuratObject(object) else object
  } else {
    object <- Seurat::CreateSeuratObject(counts = object, project = "GSES", assay = assay, min.cells = min.cells, min.features = min.features)
    assay <- "RNA"
  }
  Seurat::GetAssayData(object, assay = assay, slot = slot)
}

#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @param y vector with 2 values with up and down boundaries
#' @name scaleXY
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scaleXY <- function(x, y) {
  (x - y[1]) / (y[2] - y[1])
}


#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @param y vector with 2 values with max and min value
#' @name scale01
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scale01 <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' Transform a list comprised of character vectors with equal lengths to matrix
#'
#' @param vectorList Vector list
#' @name sameSizeVectorList2Matrix
#' @return Matrix contains character
#' @author Qiong Zhang
#'
sameSizeVectorList2Matrix <- function(vectorList){
  require('Matrix', quietly = T)
  sm_i<-NULL
  sm_j<-NULL
  sm_x<-NULL
  for (k in 1:length(vectorList)) {
    sm_i <- c(sm_i,rep(k,length(vectorList[[k]]@i)))
    sm_j <- c(sm_j,vectorList[[k]]@i)
    sm_x <- c(sm_x,vectorList[[k]]@x)
  }
  return (Matrix::sparseMatrix(i=sm_i, j=sm_j, x=sm_x, dims=c(length(vectorList), vectorList[[1]]@length)))
}


#' Calculate the up- and down- range of ranking for gene sets with gene background
#'
#' @param geneSets Genesets
#' @param geneAll Background genes
#' @name calRange
#' @return Up- and down- range of ranking for gene sets
#' @author Qiong Zhang
#'
calRange <- function(geneSets, geneAll) {
  .geneSets.name <- names(geneSets)
  .ranges <- do.call(cbind, lapply(geneSets, function(x){
    .gene.num.gs <- length(x)
    .gene.num.al <- length(geneAll)
    .z <- ceiling((.gene.num.gs + 1)/2)
    return(c(.z,  .gene.num.al - .z + 1 ))
  }))
  .ranges
}

