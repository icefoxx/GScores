#' Calculate Signature score for genes from single-cell data
#'
#' @name GScores
#' @param object Seurat object.
#' @param genesets A list of genesets, for example:
#'     \code{list(CD8CTL = c("CD3D", "CD3E", "CD4-", "CD8A", "CD8B", "KLRF1", "GZMA", "CCL5"),
#'     Bcell = c("CD19", "CD79A", "CD3D-", "CD3E-"))}
#'     Users can also specify genes with negative effects by suffixing the '-', see examples of CD8CTL
#' @param assay The sce object assay where the data is to be found
#' @param slot Which slot used in the analysis
#' @param methods methods used to calculate gene set signatures
#' @param scale scale gene set signatures (from 0 to 1)
#' @examples
#' library(GScores)
#' genesets <- list(CD8CTL = c("CD3D", "CD3E", "CD4-", "CD8A", "CD8B", "KLRF1", "GZMA", "CCL5"),
#'     Bcell = c("CD19", "CD79A", "CD3D-", "CD3E-"))
#' methods <-  c('Gsignatures', 'AUCell', 'UCell', 'AddModuleScore','ssgsea', 'VAM', 'plage', 'pagoda2', 'VISION')
#' scores <- GScores(pbmc3k.final, genesets, assay = 'RNA', slot = 'counts', methods = methods, scale = T)
#' head(scores)
#' @export
#'
GScores <- function(object, genesets, assay = "RNA", slot = "counts", methods = c("GSignatures"), scale = F, return.type = 'obj') {

  .genesets <- getExistGenes(object, genesets)
  .paralst <- list(methods = methods, slot = slot, assya = assay)
  .para.run <- checkPara(.paralst)
  .run.fun <- function(funct, object, genesets, assay, slot, methods) {
    out <- tryCatch(
      {
        print_message(stringr::str_glue('Run {methods} for gene set scoring'))
        object <- funct(object, genesets = genesets, assay = assay, slot = slot, methods = methods)
        return(object)
      },
      error=function(cond) {
        print_message(paste0("Scoring failure for the method: ", methods))
        print_message("Here's the original error message: ")
        print_message(cond)
      },
      finally={
        print_message(paste0("Processed gene set scoring by: ", methods))
      }
    )
  }
#  methods = .out.m, callfun = .out.r, slots = .out.s
  for( i in seq_along(.para.run$methods)){
    .cal.method <- .para.run$methods[[i]]
    .cal.fun <- .para.run$callfun[[i]]
    .cal.slot <- .para.run$slots[[i]]
    .cal.run <- get(paste0('run_', .cal.fun))
    object <- .run.fun(funct = .cal.run, object = object, genesets = .genesets, assay = assay, slot = .cal.slot, methods = .cal.method)
  }
  if (isTRUE(scale)){
    for(i in seq_along(.para.run$methods)){
      .cal.method <- .para.run$methods[[i]]
      .scale.score <- apply(object[[.cal.method]]@data, 1, scale01)
      object <- SeuratObject::SetAssayData(object, slot = "scale.data",
                                           new.data = t(.scale.score), assay = .cal.method)
    }
  }
  if (return.type == 'scores') {
    .slot.target <- 'data'
    if (isTRUE(scale)) .slot.target <- 'scale.data'
    .s.o <- fetchFeasMeta(object, names(.genesets), slot = .slot.target, assay = methods)
    return(.s.o)
  } else {
    return(object)
  }
}
