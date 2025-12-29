
#' Print message
#'
#' @param ... Character string to print
#' @param time Logical to print time
#' @param verbose Logical to print message
#' @name print_message
#' @author Qiong Zhang
print_message <- function(..., time = T, verbose = T){
  if (verbose){
    if (time){
      message(Sys.time() , ": ",  ...)
    } else {
      message(...)
    }
  }
}

#' Extract positive and negative genes and seperate them to 2 lists from genesets
#'
#' This function seperates genes to 2 parts, positive and negative parts, based on their suffix label
#'
#' @name getPosNegSets
#' @author Qiong Zhang
#' @param genesets A list of gene sets to be used for gene signature scoring.
#' @return gene list contains 2 elements: positive and negative genesets.
#' @export
#'
#'@examples
#' \dontrun{
#' # Load example data and gene sets
#' genesets <- list(CD4 = c("CD3D", "CD3E", "CD4", "CD8A-", "CD8B-"), B = c("CD19+", "CD79A", "CD3D-"))
#'
#' genesets.pn <- getPosNegSets(genesets = genesets)
#' }
#'
splitGeneSet2PosNeg<<- function(genesets, geneall, ignore = T, discard = F){
  .geneSet.diff <- setdiff(gsub("-$", "", unlist(genesets)), geneall)
  if (length(.geneSet.diff) > 0 & isFALSE(ignore)){
    .geneSet.diff.string <- paste(.geneSet.diff, collapse = ', ')
    .geneSet.diff.string.err <- paste("Genes missing: ", .geneSet.diff.string)
    stop( .geneSet.diff.string.err )
  }
  .genes.l <- lapply(seq_along(genesets), function(x, geneall, discard){
    .gene.r <- genesets[[x]]
    .gene.g <- gsub("[-+]$", "", .gene.r)
    .gene.a <- .gene.r[which(.gene.g %in% geneall)]
    .gene.rev.index <- grepl('-$', .gene.a, perl=TRUE)
    .gene.index <- !.gene.rev.index
    .gene.p <- .gene.a[which(.gene.index)]
    .gene.p <- gsub("[-+]$", "", .gene.p)
    .gene.n <- .gene.a[which(.gene.rev.index)]
    return(list(pos = .gene.p, neg = .gene.n))
  }, geneall = geneall, discard = discard)
  names(.genes.l) <- names(genesets)
  return(.genes.l)
}

#' Fetch features and metadata
#'
#' This function fetches features and metadata based on the key words.
#'
#' @name fetchFeasMeta
#' @author Qiong Zhang
#' @param seu Seurat object
#' @param features which feature(s)
#' @param keywords A list corresponds to the feature parameter.
#' @param assay target assay in seu object
#' @param slot target slot in assay
#' @return A data table fits the value
#' @export
#'
#'@examples
#' \dontrun{
#' # Load example data and gene sets
#' filter <- list(seurat_annotations = c('Memory CD4 T', 'Naive CD4 T'), CD4 = c('>1', '<1.6'))
#'
#' fetchFeasMeta(seu = pbmc3k.final, features = c("seurat_annotations", 'CD4'), keywords = filter)
#' }
#'
fetchFeasMeta <- function(seu, features, slot = 'data', assay = "RNA", keywords = NULL, cell.ident = T){
  .asn <- Seurat::Assays(seu)
  .assay.n <- length(assay)
  if(!all(assay %in% .asn)) stop(stringr::str_glue('{assay} does not exist in {seu}'))
  if (isTRUE(cell.ident)) {
    .data.f <- data.frame(seu.ident = Seurat::Idents(seu))
  }else{
    .data.f <- data.frame()
  }
  .seu.meta <- seu@meta.data
  .seu.meta.name <- colnames(.seu.meta)
  if(any(features %in% .seu.meta.name)) {
    .f.i.m <- match(features, .seu.meta.name)
    .f.i.d <- .seu.meta[, .f.i.m]
    colnames(.f.i.d) <- intersect(features, .seu.meta.name)
    .data.f <- cbind(.data.f, .f.i.d)
  }
  .feas <- setdiff(features, .seu.meta.name)
  if (.assay.n > 1){
    .data.feas <- do.call(cbind, lapply(seq_along(assay), function(x){
      Seurat::DefaultAssay(seu) <- assay[x]
      .d <- Seurat::FetchData(object = seu, vars = .feas, slot = slot)
      .d.n <- paste0(assay[x], '.', colnames(.d))
      colnames(.d) <- .d.n
      return(.d)
      }))
  } else {
    Seurat::DefaultAssay(seu) <- assay
    .data.feas <- Seurat::FetchData(object = seu, vars = .feas, slot = slot)
  }
  .data.f <- data.table::data.table(cbind(.data.f, .data.feas), keep.rownames = 'cell.barcodes.selected')
  if (is.null(keywords) | (.assay.n > 1)) return(.data.f)
  .keywords.col <- names(keywords)
  if (all(.keywords.col %in% features)){
    .filters <- sapply(seq_along(.keywords.col), function(x){
      .n <- .keywords.col[x]
      .v <- keywords[[x]]
      .c <- ' %in% '
      .filter <- .filter <- paste0(.n, .c, keywords[x])
      if(stringr::str_starts(.v[1], '[><=]')) {
        .c <- ' '
        .filter <- paste0(.n, .c, keywords[[x]], collapse = ' & ')
      }
      return(.filter)
    }
    )
    .filters <- paste('(',.filters, ')', collapse = ' & ')
    print(.filters)
    .data.f <- .data.f[eval(parse(text = .filters))]
    return(.data.f)
  } else{
    .f.lost <- setdiff(features, .keywords.col)
    stop(stringr::str_glue('Features missing: {.f.lost}'))
  }
}

#' Fetch msigdb gene sets
#'
#' This function fetches getsets from msigdb.
#'
#' @name fetchSetsMsigdb
#' @author Qiong Zhang
#' @param seu Seurat object
#' @param features which feature(s)
#' @param keywords A list corresponds to the feature parameter.
#' @param assay target assay in seu object
#' @param slot target slot in assay
#' @return A data table fits the value
#' @export
#'
#'@examples
#' \dontrun{
#' # Load example data and gene sets
#' filter <- list(seurat_annotations = c('Memory CD4 T', 'Naive CD4 T'), CD4 = c('>1', '<1.6'))
#'
#' fetchFeasMeta(seu = pbmc3k.final, features = c("seurat_annotations", 'CD4'), keywords = filter)
#' }
#'
fetchSetsMsigdb <- function(species = 'Homo sapiens', category = 'H', subcategory = NULL, formatId = 'symbol', localdb = F){
  .formatId <- switch(formatId,
                      symbol = "gene_symbol",
                      geneid = "entrez_gene",
                      ensembl = "ensembl_gene",
                      "gene_symbol")
  if(isTRUE(localdb)) {
    .msigdbr.path.local <- system.file("extdata", 'msigdbr.rds', package = "GSignatures")
    .msigdbr.data <- readr::read_rds(.msigdbr.path.local)
    .h.dt <- .msigdbr.data[species]
    .h.dt <- .h.dt[gs_cat %in% category]
    if (!is.null(subcategory)) .h.dt <- .h.dt[grepl(subcategory, gs_subcat)]
  } else {
  .h.dt <- data.table::data.table(msigdbr::msigdbr(species = species, category = category, subcategory = subcategory))
  }
  .h.dt.s <- .h.dt[,.(gs_name, get(.formatId))]
  .h.dt.s <- unique(.h.dt.s)
  .h.dt.s.lst <- lapply(split(.h.dt.s, .h.dt.s$gs_name), function(x) x$V2)
  return(.h.dt.s.lst)
}


#' Fetch signal genes existed in expression matrix
#'
#'
#' @name getExistGenes
#' @author Qiong Zhang
#' @param expM Seurat object or expression matrix
#' @param geneSets which feature(s)
#' @return GeneSets contain genes in the expression matrix
#' @export
#'
getExistGenes <- function(expm, geneSets, trim = T, ratio_warn = T, topoint = T){

  .genes.all <- rownames(expm)
  .genes.inq <- unlist(geneSets)
  if (any(stringr::str_detect(.genes.inq, '-$'))) .genes.inq <- gsub("-$", "", .genes.inq)
  if (isTRUE(trim)) .genes.inq <- gsub('\\.\\d+', '', .genes.inq)
  .genes.inq.f <- .genes.inq %in% .genes.all
  if (!all(.genes.inq.f)) {
    if (isTRUE(ratio_warn)) {
      .ratio.g <- sum(.genes.inq.f) / length(.genes.inq.f) * 100
      if (.ratio.g < 50) warning(stringr::str_glue('{.ratio.g}% genes in the geneSets do NOT exist in the expression matrix'))
    }
  }
  .geneSets.loci.end <- cumsum(sapply(geneSets, length))
  .geneSets.loci.start <- c(1, .geneSets.loci.end[-length(.geneSets.loci.end)]+1)
  .geneSets.f <- lapply(seq_along(.geneSets.loci.end), function(x) {
    .gene.region <- c(.geneSets.loci.start[x] : .geneSets.loci.end[x])
    .genes <- .genes.inq[.gene.region][.genes.inq.f[.gene.region]]
    })
  if (isTRUE(topoint)) {
    names(geneSets) <- stringr::str_replace_all(names(geneSets), pattern = '_', '.')
    message("Feature names with underscores ('_') have been replaced with dot ('.')")
  }
  names(.geneSets.f) <- names(geneSets)
  .geneSets.f.f <- Filter(length, .geneSets.f)
  if (isTRUE(ratio_warn)) {
    .geneSets.missing <- setdiff(names(geneSets), names(.geneSets.f.f))
    if (length(.geneSets.missing) > 0 ) warning(stringr::str_glue('Genes of {.geneSets.missing} do NOT exist in the expression matrix'))
  }
  return(.geneSets.f.f)
}

#' Validate whether users had installed selected methods
#'
#' @name vldMethods
#' @author Qiong Zhang
#' @param expM Seurat object or expression matrix
#' @param geneSets which feature(s)
#' @return GeneSets contain genes in the expression matrix
#' @export
#'
vldMethods <- function(method.use, install = T){
  print_message(stringr::str_glue('Checck availability for methods'))
  .methods.avail <- c("GSignatures","VAM", "VISION", "pagoda2",
                   "AUCell", "UCell", "singscore", "ssgsea",
                   "JASMINE", "scSE", "wmean", "wsum",
                   "mdt", "viper", "gsva", "zscore",
                   "plage", "AddModuleScore", "scGSEA")
  .methods.pacs <- c("Gsignatures", "VAM", "VISION", "pagoda2+scde",
                     'AUCell', 'UCell', 'singscore', 'GSVA',
                     'build.in', 'build.in', 'decoupleR', 'decoupleR',
                     'decoupleR', 'decoupleR', 'GSVA', 'GSVA',
                     'GSVA', 'Seurat', 'gficf+pointr+RcppML')
  .methods.source <- list(github = c(VISION = "YosefLab/VISION",
                                     gficf = "gambalab/gficf",
                                     Gsignatures = "icefoxx/Gsignatures/"
                                     ))
  .method.diff.index <- match(method.use, .methods.avail)
  .method.diff <- method.use[is.na(.method.diff.index)]
  if (length(.method.diff) > 0) stop(stringr::str_glue('{.method.diff} do NOT exist in our package'))
  .pac.need.index <- .method.diff.index[!is.na(.method.diff.index)]
  .pac.need <- .methods.pacs[.pac.need.index]
  .package.lack <- .pac.need[ ! .pac.need == 'build.in']
  .ps <- unlist(stringr::str_split(.package.lack, '\\+'))
  .package.lack <- setdiff(.ps, rownames(utils::installed.packages()))
  if (length(.package.lack) > 0) {
    if (isTRUE(install)){
      for(i in seq_along(.package.lack)){
        .p <- .package.lack[i]
        print_message(stringr::str_glue('Missing packages {.p}. Automatical installation starts'))
        pacman::p_install(.p, character.only = T)
        if (!requireNamespace(.p, quietly = TRUE)) pacman::p_install_gh(.methods.source[['github']][.p])
        if (!requireNamespace(.p, quietly = TRUE)) stop(stringr::str_glue('Installation {.p} failure,  should be manually installed by users'))
        }
    } else {
      .ps <- paste0(.package.lack, collapse = ',')
      stop(stringr::str_glue('Missing package {.ps}, should be manually installed by users'))
    }
  }
  print_message(stringr::str_glue('Methods checking completed'))
}

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

#' Construct a weight vector or data.frame for gene sets
#'
#' @name geneset2w
#' @author Qiong Zhang
#' @param genesets Seurat object or expression matrix
#' @param outformat 'd' (data.frame) or l (list)
#' @return weight for genesets
#' @export
#'
geneset2w <- function(genesets, outformat = 'l'){
  .g.w <- lapply(seq_along(genesets), function(x){
    .g <- genesets[[x]]
    .g.s <- rep(1, length(genesets[[x]]))
    .g.i <- grepl('-$', .g)
    .g.s[which(.g.i)] <- -1
    .g.n <- stringr::str_remove(.g, "\\+$|-$")
    names(.g.s) <- .g.n
    return(.g.s)
  })
  names(.g.w) <- names(genesets)
  .out <- .g.w
  if (identical(outformat, 'd')) {
    .out <- do.call(rbind, lapply(seq_along(.g.w), function(x){
      data.frame(source = names(.g.w[x]),
                 target = names(.g.w[[x]]),
                 weight = .g.w[[x]]
                 )
    }))
    rownames(.out) <- NULL
  }
  return(.out)
}

#' Check parameters for functions
#'
#' @name checkPara
#' @author Qiong Zhang
#' @param paralst parameters list
#' @return Validate parameters in a ENV
#' @export
#'
checkPara <- function(paralst, auto = T){
  print_message(stringr::str_glue('Checking parameters'))
  .slots <- c('data', 'counts')
  vldMethods(method.use = paralst$methods)
  .methods <- c(GSignatures = 'GSignatures',
               AUCell = 'AUCell',
               UCell = 'UCell',
               singscore = 'singscore',
               GSXX = 'ssgsea',
               VAM = 'VAM',
               scSE = 'scSE',
               decoupleR = 'wmean',
               decoupleR = 'wsum',
               decoupleR = 'mdt',
               decoupleR = 'viper',
               JASMINE = 'JASMINE',
               GSXX = 'gsva',
               GSXX = 'zscore',
               GSXX = 'plage',
               AddModuleScore = 'AddModuleScore',
               scGSEA = 'scGSEA',
               pagoda2 = 'pagoda2',
               VISION = 'VISION')
  .adpt.slot <- rep(.slots, c(16,3))
  .out.i <- match(paralst$methods, .methods)
  .out.m <- .methods[.out.i]
  .out.s <- .adpt.slot[.out.i]
  .out.r <- names(.out.m)
  .diff.slot <- paralst$slots != .out.s
  if (any(.diff.slot)) {
    .diff.slot.methods <- paralst$methods[which(.diff.slot)]
    stop(stringr::str_glue('Change slot for methods {.diff.slot.methods}'))
  }
  print_message(stringr::str_glue('Parameters checking completed'))
  return(list(methods = .out.m, callfun = .out.r, slots = .out.s))
}

#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @param y vector with 2 values with up and down boundaries
#' @name scaleXY
#' @return rescale vector with value between y[2] and y[1]
#' @author Qiong Zhang
#'
scaleXY <- function(x, y) {
  (x - y[1]) / (y[2] - y[1])
}



