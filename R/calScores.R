#' Run UCell analysis on a Seurat object
#'
#' This function performs UCell analysis on a given Seurat object using the specified gene sets.
#' It calculates the UCell scores and stores them in a new assay called "UCell".
#'
#' @name run_UCell
#' @author Qiong Zhang
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for UCell scoring.
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @param MaxRank Integer specifying the maximum rank for UCell scoring. Default is 1500.
#'
#' @return The input Seurat object with the UCell scores added as a new assay called "UCell".
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data and gene sets
#' library(Seurat)
#' pbmc_small <- Seurat::pbmc_small()
#' example_genesets <- list(set1 = c("CD3D", "CD3E"), set2 = c("CD19", "CD79A"))
#'
#' # Run UCell analysis
#' pbmc_small_ucell <- runUCell(pbmc_small, genesets = example_genesets)
#' }
#'
#'
run_UCell <- function(object, genesets, assay = 'RNA', slot = 'counts', methods = 'UCell', MaxRank = 1500){
  print_message("UCell scoring starts")
  .expM <- assay2mtx(object, assay = assay, slot = 'counts')
  .us <- UCell::ScoreSignatures_UCell(matrix = .expM,
                                      features = genesets,
                                      maxRank = MaxRank,
                                      w_neg = 1)
  .us <- t(.us)
  rownames(.us) <- names(genesets)
  object[["UCell"]] <- SeuratObject::CreateAssayObject(counts = .us)
  print_message("Finish calculate UCell scores")
  return(object)
}


#' Run AUCell analysis on a Seurat object or expression matrix
#'
#' This function performs AUCell analysis on a given Seurat object or expression matrix using the specified gene sets.
#' It calculates the AUCell scores and stores them in a new assay called "AUCell".
#'
#' @name run_AUCell
#' @author Qiong Zhang
#' @param object A Seurat object containing the single-cell data or an expression matrix.
#' @param genesets A list of gene sets to be used for AUCell scoring.
#' @param aucratio A numeric value specifying the ratio for determining the maximum rank in AUCell scoring.
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#'
#' @return The input Seurat object or expression matrix with the AUCell scores added as a new assay called "AUCell".
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data and gene sets
#' library(AUCell)
#' pbmc_small <- Seurat::pbmc_small()
#' example_genesets <- list(set1 = c("CD3D", "CD3E"), set2 = c("CD19", "CD79A"))
#'
#' # Run AUCell analysis
#' pbmc_small_aucell <- run_AUCell(pbmc_small, genesets = example_genesets, aucratio = 0.05)
#' }
#'
run_AUCell <- function(object, genesets, assay = 'RNA', slot = 'counts', methods = 'AUCell', aucratio = 0.05){
  print_message("Start AUCell scoring")
  .expMat <- assay2mtx(object, assay = 'RNA', slot = 'counts')
#  .expMat <- ifelse(is.data.frame(object), object, tryCatch(as.data.frame(object), error = identity))
  .calAUCell <- function(cellindex, expM, genesets, aucratio, ncore){
    .gene.n <- nrow(expM)
    .exp <- expM[, cellindex]
    .aucell.rank <- AUCell::AUCell_buildRankings(.exp,
                                                 plotStats = F,
                                                 verbose = F,
                                                 splitByBlocks = TRUE)
    .genesets4aucell <- genesets |> purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
    .aucMaxRank <- ifelse (is.null(aucratio), ceiling(0.05 * .gene.n), ceiling(aucratio * .gene.n))
    .aucell.scores <- AUCell::AUCell_calcAUC(.genesets4aucell,
                                             .aucell.rank,
                                             aucMaxRank = .aucMaxRank,
                                             verbose = F)
    .aucell.scores <- SummarizedExperiment::assay(.aucell.scores)
  }
  if (methods::is(.expMat, "error")) {
    .expMat.lst <- intervalCell(.expMat)
    .aucell.scores.lst <- do.call(rbind, lapply(.expMat.lst, .calAUCell, .expMat, genesets, aucratio, ncore))
  } else {
    .aucell.scores.lst <- .calAUCell(1:ncol(.expMat), .expMat, genesets, aucratio, ncore)
  }
  object[["AUCell"]] <- SeuratObject::CreateAssayObject(counts = .aucell.scores.lst)
  print_message("AUCell scoring completed")
  return(object)
}

#' Run ssGSEA or c("gsva", "zscore", "plage") analysis on a Seurat object
#'
#' This function performs ssGSEA (single sample Gene Set Enrichment Analysis) on a given Seurat object using the specified gene sets.
#' It calculates the ssGSEA scores and stores them in a new assay called "ssgsea".
#' @author Qiong Zhang
#' @name run_GSXX
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for ssGSEA scoring.
#' @param methods could be one of c("gsva", "ssgsea", "zscore", "plage")
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#'
#' @return The input Seurat object with the ssGSEA scores added as a new assay called "ssgsea".
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data and gene sets
#' library(Seurat)
#' pbmc_small <- Seurat::pbmc_small()
#' example_genesets <- list(set1 = c("CD3D", "CD3E"), set2 = c("CD19", "CD79A"))
#'
#' # Run ssGSEA analysis
#' pbmc_small_ssgsea <- runssGSEA(pbmc_small, genesets = example_genesets, aucratio = 0.05)
#' }
#'
run_GSXX <- function(object, genesets, assay = 'RNA', slot = 'data', methods = 'ssgsea', kcdf = 'Gaussian', ssgsea.norm = F){
  print_message(stringr::str_glue("start {methods} scoring"))
  .expMat <- assay2mtx(object, assay = assay, slot = slot)
  .genesets4ss <- genesets |> purrr::discard(.p = function(x){all(stringr::str_detect(x, pattern = "\\+$|-$"))})
  if (methods %in% c('plage')) {
    .n <- ncol(.expMat)
    .rmean <- Matrix::rowSums(.expMat) / .n
    .rvar  <- Matrix::rowSums(.expMat^2) / n - .rmean ^ 2
    .expMat   <- .expMat[.rvar > 0, , drop = FALSE]
    .genesets2 <- lapply(.genesets4ss, function(gs) intersect(gs, rownames(.expMat)))
    .genesets2 <- .genesets2[lengths(.genesets2) >= 2]
    .genes_in_sets <- unique(unlist(.genesets2, use.names = FALSE))
    .expMat <- as.matrix(.expMat[.genes_in_sets, , drop = FALSE])
    storage.mode(.expMat) <- "double"
    .genesets4ss <- .genesets2
  }
  .ss.scores <- GSVA::gsva(.expMat,
                           .genesets4ss,
                           method = methods,
                           kcdf = kcdf,
                           ssgsea.norm = ssgsea.norm,
                           verbose = F)
  object[[methods]] <- SeuratObject::CreateAssayObject(counts = as.matrix(.ss.scores))
  print_message(stringr::str_glue("{methods} scoring completed"))
  return(object)
}

#' Run singscore analysis on a Seurat object
#'
#' This function performs singscore analysis on a given Seurat object using the specified gene sets.
#' It calculates the singscore scores and stores them in a new assay called "singscore".
#'
#' @name run_singscore
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for singscore analysis.
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the singscore scores added as a new assay called "singscore".
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data and gene sets
#' library(Seurat)
#' pbmc_small <- Seurat::pbmc_small()
#' example_genesets <- list(set1 = c("CD3D", "CD3E"), set2 = c("CD19", "CD79A"))
#'
#' # Run singscore analysis
#' pbmc_small_singscore <- runsingscore(pbmc_small, genesets = example_genesets, aucratio = 0.05)
#' }
#'
run_singscore <- function(object, genesets, assay = 'RNA', slot = 'counts', methods = 'singscore' ){
  print_message("Start singscore scoring")
  .calSingscore <- function(cellindex, expMat, genesets){
    .ss.scores <- list()
    .ss.rank <- tryCatch(singscore::rankGenes(as.data.frame(expMat[, cellindex])), error = identity)
    for (i in seq_along(genesets)){
      .genes <- genesets[[i]]
      .geneset.pos <- NULL
      .geneset.neg <- NULL
      if (any(stringr::str_detect(.genes, pattern = "\\+$|-$"))) {
        .geneset.pos <- stringr::str_match(.genes, pattern = "(.+)\\+$")[,2] |> purrr::discard(is.na)
        .geneset.neg <- stringr::str_match(.genes, pattern = "(.+)\\-$")[,2] |> purrr::discard(is.na)
      } else {
        .geneset.pos <- .genes |> purrr::discard(is.na)
      }
      if (length(.geneset.neg) > 0){
        if (length(.geneset.pos) > 0){
          .ss.scores[[i]] <- singscore::simpleScore(.ss.rank, upSet = .geneset.pos, downSet = .geneset.neg, centerScore = F)
        } else {
          .ss.scores[[i]] <- singscore::simpleScore(.ss.rank, upSet = .geneset.neg, centerScore = F)
        }
      } else {
        .ss.scores[[i]] <- singscore::simpleScore(.ss.rank, upSet = .geneset.pos, centerScore = F)
      }
      .ss.scores[[i]] <- .ss.scores[[i]]$TotalScore
    }
    names(.ss.scores) <- names(genesets)
    .ss.scores <- do.call(rbind, .ss.scores)
    colnames(.ss.scores) <- colnames(expMat)
    return(.ss.scores)
  }
  .ss.scores.all <- NULL
  .expMat <- assay2mtx(object, assay = assay, slot = slot)
  if (methods::is(.expMat, "error")) {
    .expMat.lst <- intervalCell(.expMat)
    .ss.scores.all <- do.call(rbind, lapply(.expMat.lst, .calSingscore, .expMat, genesets))
  } else {
    .ss.scores.all <- .calSingscore(1:ncol(.expMat), .expMat, genesets)
  }
  object[["singscore"]] <- SeuratObject::CreateAssayObject(counts = .ss.scores.all)
  print_message("singscore scoring completed")
  return(object)
}

#' Run JASMINE scoring on a Seurat object
#'
#' This function performs JASMINE scoring on a given Seurat object using the specified gene sets.
#' It calculates the JASMINE scores and stores them in a new assay called "JASMINE".
#'
#' @name run_JASMINE
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for JASMINE scoring
#' @param method 'oddsratio' or 'likelihood'.
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the JASMINE scores added as a new assay called "JASMINE".
#' @export
#'
run_JASMINE <- function(object, genesets, assay = 'RNA', slot = 'counts', method = 'likelihood', methods = 'JASMINE'){

  print_message(stringr::str_glue("Start JASMINE scoring"))
  #all codes (annotation removal) were from JASMINE github except for the Multicore calculation

  RankCalculation <- function(x,genes){
    subdata = x[x!=0]
    DataRanksUpdated=rank(subdata)
    DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]
    CumSum = ifelse(length(DataRanksSigGenes),mean(DataRanksSigGenes,na.rm = TRUE),0 )
    FinalRawRank = CumSum/length(subdata)
    return(FinalRawRank)
  }

  ORCalculation <- function(data, genes){
    GE = data[which(rownames(data) %in% genes),]
    NGE = data[-which(rownames(data) %in% genes),]
    SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))
    NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))
    SigGenesNE = nrow(GE) - SigGenesExp
    SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)
    NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)
    NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)
    NSigGenesNE = NSigGenesNE - SigGenesNE
    OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)
    return(OR)
  }

  LikelihoodCalculation <- function(data,genes){
    GE = data[which(rownames(data) %in% genes),]
    NGE = data[-which(rownames(data) %in% genes),]
    SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))
    NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))
    SigGenesNE = nrow(GE) - SigGenesExp
    SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)
    NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)
    NSigGenesNE = nrow(data) - (NSigGenesExp + SigGenesExp)
    NSigGenesNE = NSigGenesNE - SigGenesNE
    LR1 = SigGenesExp*(NSigGenesExp + NSigGenesNE)
    LR2 = NSigGenesExp * (SigGenesExp + SigGenesNE)
    LR = LR1/LR2
    return(LR)
  }

  NormalizationJAS <- function(JAS_Scores){
    JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores)- min(JAS_Scores))
    return(JAS_Scores)
  }

  .JASMINE <- function(data,genes,method){
    idx = match(genes,rownames(data))
    idx = idx[!is.na(idx)]
    if(length(idx)> 1){
      RM = apply(data,2,function(x) RankCalculation(x,genes))
      RM = NormalizationJAS(RM)
      if(method == "oddsratio"){
        OR = ORCalculation(data,genes)
        OR = NormalizationJAS(OR)
        JAS_Scores = (RM + OR)/2
      }else if(method == "likelihood"){
        LR = LikelihoodCalculation(data,genes)
        LR = NormalizationJAS(LR)
        JAS_Scores = (RM + LR)/2
      }
      FinalScores = data.frame(JAS_Scores)
      return(FinalScores)
    }
  }
  .expM <- assay2mtx(object, assay = 'RNA', slot = 'counts')
  .score <- do.call(cbind, lapply(seq_along(genesets), function(x) {
    .JASMINE(.expM, genes = genesets[[x]], method = method)
  }))
  colnames(.score) <- names(genesets)
  print_message(stringr::str_glue("JASMINE scoring completed"))
  View(.score)
  object[['JASMINE']] <- SeuratObject::CreateAssayObject(counts = t(.score))
  return(object)
}

#' Run VAM scoring on a Seurat object
#'
#' This function performs VAM scoring on a given Seurat object using the specified gene sets.
#' It calculates the VAM scores and stores them in a new assay called "VAM".
#'
#' @name run_VAM
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for VAM scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the VAM scores added as a new assay called "VAM".
#' @export
#'
run_VAM <- function(object, genesets, assay = 'RNA', slot = 'data', methods = 'VAM'){
  Seurat::DefaultAssay(object) <- assay
  print_message('Start VAM scoring')
  .genesets.vam <- VAM::createGeneSetCollection(gene.ids = rownames(object),
                                                      gene.set.collection = genesets)
  object <- VAM::vamForSeurat(seurat.data = object,
                               gene.set.collection = .genesets.vam,
                               center = F,
                               gamma = T,
                               sample.cov = F,
                               return.dist = F)
  object[["VAM"]] <- object[['VAMcdf']]
  object[['VAMcdf']] <- NULL
  print_message('VAM scoring completed')
  return(object)
}

#' Run scSE scoring on a Seurat object
#'
#' This function performs scSE scoring on a given Seurat object using the specified gene sets.
#' It calculates the scSE scores and stores them in a new assay called "scSE".
#'
#' @name run_scSE
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for scSE scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the scSE scores added as a new assay called "scSE".
#' @export
#'
run_scSE <- function(object, genesets, assay = 'RNA', slot = "counts", methods = 'scSE'){
  print_message('Start scSE scoring')
  .expMat <- assay2mtx(object, assay = assay, slot = slot)
  .expMat.umi.sum <- Matrix::colSums(.expMat)
  .genesets.umi.sum <- do.call(rbind, lapply(seq_along(genesets), function(x) Matrix::colSums(.expMat[genesets[[x]], ])))
  .genesets.umi.score <- .genesets.umi.sum/.expMat.umi.sum * 100
  colnames(.genesets.umi.score) <- colnames(.expMat)
  rownames(.genesets.umi.score) <- names(genesets)
  print_message('scSE scoring completed')
  object[["scSE"]] <- SeuratObject::CreateAssayObject(counts = .genesets.umi.score)
  return(object)
}

#' Run VISION scoring on a Seurat object
#'
#' This function performs VISION scoring on a given Seurat object using the specified gene sets.
#' It calculates the VISION scores and stores them in a new assay called "VISION".
#'
#' @name run_VISION
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for VISION scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the VISION scores added as a new assay called "VISION".
#' @export
#'
run_VISION <- function(object, genesets, assay = 'RNA', slot = "counts", methods = 'VISION'){
  print_message('Start VISION scoring')
  .genesets.weight <- geneset2w(genesets, 'l')
  .genesets.vision <- lapply(seq_along(.genesets.weight), function(x) VISION::createGeneSignature(name = names(.genesets.weight[x]), sigData = .genesets.weight[[x]]))
  .vision.obj <- VISION::Vision(object,
                               signatures = .genesets.vision,
                               projection_methods = NULL,
                               assay = assay)
  .vision.obj <- VISION::analyze(.vision.obj)
  .vision.score <- VISION::getSignatureScores(.vision.obj)
  object[['VISION']] <- SeuratObject::CreateAssayObject(counts = t(.vision.score))
  print_message('VISION scoring completed')
  return(object)
}

#' Run Seurat scoring on a Seurat object
#'
#' This function performs Seurat scoring on a given Seurat object using the specified gene sets.
#' It calculates the Seurat scores and stores them in a new assay called "Seurat".
#'
#' @name run_AddModuleScore
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for Seurat scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the Seurat scores added as a new assay called "Seurat".
#' @export
#'
run_AddModuleScore <- function(object, genesets, assay = 'RNA', slot = "counts", methods = 'AddModuleScore', seed = 10086, ncores = 1){
  print_message('Start Seurat scoring')
  .g.n <- names(genesets)
  .tmp.name <- stringr::str_c(paste0(sample(letters, 10), collapse = ''), '_')
  object <- Seurat::AddModuleScore(object,
                                   features = genesets,
                                   name = .tmp.name,
                                   seed = seed,
                                   assay = assay)
  .object.meta.new <- object@meta.data
  .object.meta.new.name <- colnames(.object.meta.new)
  .genesets.score.name <- paste0(.tmp.name, 1:length(.g.n))
  .genesets.score.index <- match(.genesets.score.name, .object.meta.new.name)
  .genesets.score.matrix <- .object.meta.new[, .genesets.score.index, drop = F]
  colnames(.genesets.score.matrix) <- .g.n
  object[["AddModuleScore"]] <- SeuratObject::CreateAssayObject(counts = t(.genesets.score.matrix))
  object@meta.data <- object@meta.data[, -.genesets.score.index]
  print_message('AddModuleScore scoring completed')
  return(object)
}

#' Run pagoda2 scoring on a Seurat object
#'
#' This function performs pagoda2 scoring on a given Seurat object using the specified gene sets.
#' It calculates the pagoda2 scores and stores them in a new assay called "pagoda2".
#'
#' @name run_pagoda2
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for pagoda2 scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the pagoda2 scores added as a new assay called "pagoda2".
#' @export
#'
run_pagoda2 <- function(object, genesets, assay = 'RNA', slot = "counts", methods = 'pagoda2', seed = 10086, ncores = 1){
  print_message('Start pagoda2 scoring')
  .genesets.pagoda2.env <- list2env(genesets)
  .expM <- assay2mtx(object, assay = assay, slot = "counts")
  .expM.n <- ncol(.expM)
  # The following codes are from pagoda2's documents which are too complicated to understands. I just modified the n.odgens and max.pathway.size and min.pathway.size and n.randomizations to let the function fit several specific conditions.
  .data.pagoda2 <- pagoda2::basicP2proc(.expM,
                                       n.cores = ncores,
                                       n.odgenes=3e3,
                                       get.largevis=FALSE,
                                       make.geneknn=FALSE,
                                       get.tsne = F)
  .newtestPathwayOverdispersion <- function(setenv,
                                       type='counts',
                                       max.pathway.size=6e3,
                                       min.pathway.size=1,
                                       n.randomizations=3,
                                       verbose=FALSE,
                                       n.cores=self$n.cores,
                                       score.alpha=0.05,
                                       plot=FALSE,
                                       cells=NULL,
                                       adjusted.pvalues=TRUE,
                                       z.score = stats::qnorm(0.05/2, lower.tail = FALSE),
                                       use.oe.scale = FALSE,
                                       return.table=FALSE,
                                       name='pathwayPCA',
                                       correlation.distance.threshold=0.2,
                                       loading.distance.threshold=0.01,
                                       top.aspects=Inf,
                                       recalculate.pca=FALSE,
                                       save.pca=TRUE) {

    if (!requireNamespace("scde", quietly=TRUE)){
      stop("You need to install package 'scde' to be able to use testPathwayOverdispersion().")
    }

    nPcs <- 1
    if (type=='counts') {
      x <- self$counts
      # apply scaling if using raw counts
      x@x <- x@x*rep(self$misc[['varinfo']][colnames(x),'gsf'],diff(x@p))
    } else {
      if (!type %in% names(self$reductions)) { stop("Reduction ",type,' not found')}
      x <- self$reductions[[type]]
    }
    if (!is.null(cells)) {
      x <- x[cells,]
    }

    proper.gene.names <- colnames(x)

    if (is.null(self$misc[['pwpca']]) || recalculate.pca) {
      if (verbose) {
        message("determining valid pathways")
      }

      # def function
      sn <- function(x){
        names(x) <- x
        return(x)
      }


      # determine valid pathways
      gsl <- ls(envir = setenv)
      gsl.ng <- unlist(parallel::mclapply(sn(gsl), function(go) sum(unique(get(go, envir = setenv)) %in% proper.gene.names),mc.cores=n.cores,mc.preschedule=TRUE))
      gsl <- gsl[gsl.ng >= min.pathway.size & gsl.ng<= max.pathway.size]
      names(gsl) <- gsl

      if (verbose) {
        message("processing ", length(gsl), " valid pathways")
      }

      cm <- Matrix::colMeans(x)
      # def function
      papply <-utils::getFromNamespace("papply", "pagoda2")
      pwpca <- papply(gsl, function(sn) {
        lab <- proper.gene.names %in% get(sn, envir = setenv)
        if (sum(lab)<1) {
          return(NULL)
        }
        pcs <- irlba::irlba(x[,lab], nv=nPcs, nu=0, center=cm[lab])
        pcs$d <- pcs$d/sqrt(nrow(x))
        pcs$rotation <- pcs$v
        pcs$v <- NULL

        # get standard deviations for the random samples
        ngenes <- sum(lab)
        z <- do.call(rbind,lapply(seq_len(n.randomizations), function(i) {
          si <- sample(ncol(x), ngenes)
          pcs <- irlba::irlba(x[,si], nv=nPcs, nu=0, center=cm[si])$d
        }))
        z <- z/sqrt(nrow(x))

        # local normalization of each component relative to sampled PC1 sd
        avar <- pmax(0, (pcs$d^2-mean(z[, 1]^2))/stats::sd(z[, 1]^2))

        if (T) {
          # flip orientations to roughly correspond with the means
          pcs$scores <- as.matrix(t(x[,lab] %*% pcs$rotation) - as.numeric((cm[lab] %*% pcs$rotation)))
          cs <- unlist(lapply(seq_len(nrow(pcs$scores)), function(i) sign(stats::cor(pcs$scores[i,], colMeans(t(x[, lab, drop = FALSE])*abs(pcs$rotation[, i]))))))
          pcs$scores <- pcs$scores*cs
          pcs$rotation <- pcs$rotation*cs
          rownames(pcs$rotation) <- colnames(x)[lab]
        } # don't bother otherwise - it's not significant
        return(list(xp=pcs,z=z,n=ngenes))
      }, n.cores = n.cores,mc.preschedule=TRUE)
      if (save.pca) {
        self$misc[['pwpca']] <- pwpca
      }
    } else {
      if (verbose) {
        message("reusing previous overdispersion calculations")
        pwpca <- self$misc[['pwpca']]
      }
    }

    if (verbose) {
      message("scoring pathway od signifcance")
    }

    # score overdispersion
    true.n.cells <- nrow(x)

    pagoda.effective.cells <- function(pwpca, start = NULL) {
      n.genes <- unlist(lapply(pwpca, function(x) rep(x$n, nrow(x$z))))
      var <- unlist(lapply(pwpca, function(x) x$z[, 1]))
      if (is.null(start)) { start <- true.n.cells*2 } # start with a high value
      of <- function(p, v, sp) {
        sn <- p[1]
        vfit <- (sn+sp)^2/(sn*sn+1/2) -1.2065335745820*(sn+sp)*((1/sn + 1/sp)^(1/3))/(sn*sn+1/2)
        residuals <- (v-vfit)^2
        return(sum(residuals))
      }
      x <- stats::nlminb(objective = of, start = c(start), v = var, sp = sqrt(n.genes-1/2), lower = c(1), upper = c(true.n.cells))
      return((x$par)^2+1/2)
    }
    n.cells <- pagoda.effective.cells(pwpca)

    vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {
      vars <- as.numeric((pwpca[[i]]$xp$d))
      cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$rotation)))
    })))

    # fix p-to-q mistake in qWishartSpike
    qWishartSpikeFixed <- function (q, spike, ndf = NA, pdim = NA, var = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)  {
      params <- RMTstat::WishartSpikePar(spike, ndf, pdim, var, beta)
      stats::qnorm(q, mean = params$centering, sd = params$scaling, lower.tail, log.p)
    }

    # add right tail approximation to ptw, which gives up quite early
    pWishartMaxFixed <- function (q, ndf, pdim, var = 1, beta = 1, lower.tail = TRUE) {
      params <- RMTstat::WishartMaxPar(ndf, pdim, var, beta)
      q.tw <- (q - params$centering)/(params$scaling)
      p <- RMTstat::ptw(q.tw, beta, lower.tail, log.p = TRUE)
      p[p == -Inf] <- stats::pgamma((2/3)*q.tw[p == -Inf]^(3/2), 2/3, lower.tail = FALSE, log.p = TRUE) + lgamma(2/3) + log((2/3)^(1/3))
      p
    }

    vshift <- 0
    ev <- 0

    vdf$var <- vdf$var-(vshift-ev)*vdf$n
    basevar <- 1
    vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
    #vdf$z <- qnorm(pWishartMax(vdf$var, n.cells, vdf$n, log.p = TRUE, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
    vdf$z <- stats::qnorm(pWishartMaxFixed(vdf$var, n.cells, vdf$n, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
    # def function
    bh.adjust <-utils::getFromNamespace("bh.adjust", "pagoda2")
    vdf$cz <- stats::qnorm(bh.adjust(stats::pnorm(as.numeric(vdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
    vdf$ub <- RMTstat::qWishartMax(score.alpha/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
    vdf$ub.stringent <- RMTstat::qWishartMax(score.alpha/nrow(vdf)/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)

    if (plot) {
      test_pathway_par <- graphics::par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
      on.exit(graphics::par(test_pathway_par))
      un <- sort(unique(vdf$n))
      on <- order(vdf$n, decreasing = FALSE)
      pccol <- grDevices::colorRampPalette(c("black", "grey70"), space = "Lab")(max(vdf$npc))
      plot(vdf$n, vdf$var/vdf$n, xlab = "gene set size", ylab = "PC1 var/n", ylim = c(0, max(vdf$var/vdf$n)), col = grDevices::adjustcolor(pccol[vdf$npc],alpha=0.1),pch=19)
      graphics::lines(vdf$n[on], (vdf$exp/vdf$n)[on], col = 2, lty = 1)
      graphics::lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on], col = 2, lty = 2)
    }

    rs <- (vshift-ev)*vdf$n
    vdf$oe <- (vdf$var+rs)/(vdf$exp+rs)
    vdf$oec <- (vdf$var+rs)/(vdf$ub+rs)

    df <- data.frame(name = names(pwpca)[vdf$i], npc = vdf$npc, n = vdf$n, score = vdf$oe, z = vdf$z, adj.z = vdf$cz, stringsAsFactors = FALSE)
    if (adjusted.pvalues) {
      vdf$valid <- vdf$cz  >=  z.score
    } else {
      vdf$valid <- vdf$z  >=  z.score
    }

    if (!any(vdf$valid)) {
      stop("No significantly overdispersed pathways found at z.score threshold of ",z.score)
    }

    # apply additional filtering based on >0.5 sd above the local random estimate
    vdf$valid <- vdf$valid & unlist(lapply(pwpca,function(x) !is.null(x$xp$scores)))
    vdf$name <- names(pwpca)[vdf$i]

    if (return.table) {
      df <- df[vdf$valid, ]
      df <- df[order(df$score, decreasing = TRUE), ]
      return(df)
    }
    if (verbose) {
      message("compiling pathway reduction")
    }
    # calculate pathway reduction matrix

    # return scaled patterns
    xmv <- do.call(rbind, lapply(pwpca[vdf$valid], function(x) {
      xm <- x$xp$scores
    }))

    if (use.oe.scale) {
      xmv <- (xmv -rowMeans(xmv))* (as.numeric(vdf$oe[vdf$valid])/sqrt(apply(xmv, 1, stats::var)))
      vdf$sd <- as.numeric(vdf$oe)
    } else {
      # chi-squared
      xmv <- (xmv-rowMeans(xmv)) * sqrt((stats::qchisq(stats::pnorm(vdf$z[vdf$valid], lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells)/apply(xmv, 1, stats::var))
      vdf$sd <- sqrt((stats::qchisq(stats::pnorm(vdf$z, lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells))

    }
    rownames(xmv) <- paste("#PC", vdf$npc[vdf$valid], "# ", names(pwpca)[vdf$i[vdf$valid]], sep = "")
    rownames(vdf) <- paste("#PC", vdf$npc, "# ", vdf$name, sep = "")
    self$misc[['pathwayODInfo']] <- vdf

    # collapse gene loading
    if (verbose) {
      message("clustering aspects based on gene loading ... ",appendLF=FALSE)
    }
    tam2 <- pagoda2::pagoda.reduce.loading.redundancy(list(xv=xmv,xvw=matrix(1,ncol=ncol(xmv),nrow=nrow(xmv))),pwpca,NULL,plot=FALSE,distance.threshold=loading.distance.threshold,n.cores=n.cores)
    if (verbose) {
      message(nrow(tam2$xv)," aspects remaining")
    }
    if (verbose) {
      message("clustering aspects based on pattern similarity ... ",appendLF=FALSE)
    }
    tam3 <- pagoda2::pagoda.reduce.redundancy(tam2, distance.threshold=correlation.distance.threshold,top=top.aspects)
    if (verbose) {
      message(nrow(tam3$xv)," aspects remaining\n")
    }
    tam2$xvw <- tam3$xvw <- NULL # to save space
    tam3$env <- setenv

    # clean up aspect names, as GO ids are meaningless
    names(tam3$cnam) <- rownames(tam3$xv) <- paste0('aspect',1:nrow(tam3$xv))

    self$misc[['pathwayOD']] <- tam3
    self$reductions[[name]] <- tam3$xv
    invisible(tam3)
  }

  .data.pagoda2$.newtestPathwayOverdispersion <- .newtestPathwayOverdispersion
  environment(.data.pagoda2$.newtestPathwayOverdispersion) <- .data.pagoda2$.__enclos_env__
  .data.pagoda2$.newtestPathwayOverdispersion(.genesets.pagoda2.env,
                                           verbose=T,
                                           recalculate.pca=F,
                                           top.aspects=15, score.alpha =0)
  ## complicated codes finished
  .ss <- do.call(rbind, lapply(seq_along(.data.pagoda2$misc$pwpca), function(x, nc) {
    .s <- .data.pagoda2$misc$pwpca[[x]]$xp$score
    `if`(is.null(.s), rep(NULL, nc), .s)
    }, nc = .expM.n))
  colnames(.ss) <- colnames(.expM)
  rownames(.ss) <- names(genesets)
  object[["pagoda2"]] <- SeuratObject::CreateAssayObject(counts = .ss)
  print_message('pagoda2 scoring completed')
  return(object)
}


#' Run scGSEA scoring on a Seurat object
#'
#' This function performs scGSEA scoring on a given Seurat object using the specified gene sets.
#' It calculates the scGSEA scores and stores them in a new assay called "scGSEA".
#'
#' @name run_scGSEA
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for scGSEA scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the scGSEA scores added as a new assay called "scGSEA".
#' @export
#'
run_scGSEA <- function(object, genesets, assay = 'RNA', slot = "counts", methods = 'scGSEA', ncores = 1){
  print_message('Start scGSEA scoring')
  .expM <- assay2mtx(object, assay = assay, slot = slot)
  .expM.data <- gficf::gficf(M = .expM)
  .expM.data <- gficf::runPCA(data = .expM.data, dim = 10, use.odgenes = T)

  # slightly modified from "runScGSEA" function https://github.com/gambalab/gficf/blob/master/R/pathwayAnalisys.R
  .runScGSEA.m <- function (data, geneID, species, category, subcategory = NULL,
                           pathway.list = NULL, nsim = 10000, nt = 0, minSize = 15,
                           maxSize = Inf, verbose = TRUE, seed = 10086, nmf.k = 100,
                           fdr.th = 0.05, gp = 0, rescale = "none", normalization = "gficf")
  {
    if (nt == 0) nt = parallel::detectCores()
    geneID = base::match.arg(arg = geneID, choices = c("ensamble",
                                                       "symbol"), several.ok = F)
    rescale = base::match.arg(arg = rescale, choices = c("none",
                                                         "byGS", "byCell"), several.ok = F)
    normalization = base::match.arg(arg = normalization, choices = c("gficf",
                                                                     "cpm"), several.ok = F)
    options(RcppML.threads = nt)
    set.seed(seed)
    # def function
    tsmessage <-utils::getFromNamespace("tsmessage", "gficf")

    if (is.null(data$scgsea)) {
      data$scgsea = list()
      if (normalization == "gficf") {
        if (!is.null(data$pca) && data$pca$type == "NMF") {
          if (data$dimPCA < nmf.k || data$pca$use.odgenes) {
            tsmessage("... Performing NMF", verbose = verbose)
            tmp = RcppML::nmf(A = data$gficf, k = nmf.k)
            data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,
                                                sparse = T,
                                                dimnames = list(rownames(data$gficf),NULL))
            data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,
                                                  sparse = T,
                                                  dimnames = list(NULL,colnames(data$gficf))))
            rm(tmp)
            gc()
          }else {
            tsmessage(paste0("Found NMF reduction with k greaten or equal to ",
                             nmf.k), verbose = T)
            pointr::ptr("tmp", "data$pca$genes")
            data$scgsea$nmf.w = tmp
            pointr::ptr("tmp2", "data$pca$cells")
            tmp2 <- NULL
            data$scgsea$nmf.h = tmp2
            rm(tmp, tmp2)
            gc()
          }
        }else {
          tsmessage("... Performing NMF", verbose = verbose)
          tmp = RcppML::nmf(data$gficf, k = nmf.k)
          data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,
                                              sparse = T,
                                              dimnames = list(rownames(data$gficf),NULL))
          data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,
                                                sparse = T,
                                                dimnames = list(NULL,colnames(data$gficf))))
          rm(tmp)
          gc()
        }
      }else {
        tsmessage("... Performing NMF", verbose = verbose)
        tmp = RcppML::nmf(log1p(data$rawCounts), k = nmf.k)
        data$scgsea$nmf.w <- Matrix::Matrix(data = tmp$w,
                                            sparse = T)
        data$scgsea$nmf.h <- t(Matrix::Matrix(data = tmp$h,
                                              sparse = T))
        rm(tmp)
        gc()
      }
    }else {
      tsmessage("Found a previous scGSEA, thus the already computed NMF will be used",
                verbose = T)
      tsmessage("If you want to recompute NMF, please call resetScGSEA first",
                verbose = T)
      data$scgsea$es <- NULL
      data$scgsea$nes <- NULL
      data$scgsea$pval <- NULL
      data$scgsea$fdr <- NULL
      data$scgsea$pathways <- NULL
      data$scgsea$x <- NULL
      data$scgsea$stat <- NULL
    }
    tsmessage("Loading pathways...", verbose = verbose)
    if (!is.list(pathway.list)) {
      gs = msigdbr::msigdbr(species = species, category = category,
                            subcategory = subcategory)
      if (geneID == "symbol") {
        data$scgsea$pathways = split(x = gs$gene_symbol,
                                     f = gs$gs_name)
      }else {
        data$scgsea$pathways = split(x = gs$ensembl_gene,
                                     f = gs$gs_name)
      }
    }else {
      data$scgsea$pathways = pathway.list
      rm(pathway.list)
    }
    data$scgsea$es = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                    ncol = ncol(data$scgsea$nmf.w))
    data$scgsea$nes = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                     ncol = ncol(data$scgsea$nmf.w))
    data$scgsea$pval = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                      ncol = ncol(data$scgsea$nmf.w))
    data$scgsea$fdr = Matrix::Matrix(data = 0, nrow = length(data$scgsea$pathways),
                                     ncol = ncol(data$scgsea$nmf.w))
    rownames(data$scgsea$es) = rownames(data$scgsea$nes) = rownames(data$scgsea$pval) = rownames(data$scgsea$fdr) = names(data$scgsea$pathways)
    tsmessage("Performing GSEA...", verbose = verbose)
    oldw <- getOption("warn")
    options(warn = -1)
    pb = utils::txtProgressBar(min = 0, max = ncol(data$scgsea$nmf.w), initial = 0, style = 3)
    nt_fgsea <- ceiling(length(data$scgsea$pathways)/100)
    nt_fgsea <- ifelse(nt_fgsea > nt, nt, nt_fgsea)
    bpparameters <- BiocParallel::SnowParam(nt_fgsea)
    for (i in 1:ncol(data$scgsea$nmf.w)) {
      df = as.data.frame(fgsea::fgseaMultilevel(pathways = data$scgsea$pathways,
                                                stats = data$scgsea$nmf.w[, i], nPermSimple = nsim,
                                                gseaParam = gp, BPPARAM = bpparameters, minSize = minSize,
                                                maxSize = maxSize))[, 1:7]
      data$scgsea$es[df$pathway, i] = df$ES
      data$scgsea$nes[df$pathway, i] = df$NES
      data$scgsea$pval[df$pathway, i] = df$pval
      data$scgsea$fdr[df$pathway, i] = df$padj
      utils::setTxtProgressBar(pb, i)
    }
    base::close(pb)
    on.exit(options(warn = oldw))
    ix = is.na(data$scgsea$nes)
    if (sum(ix) > 0) {
      data$scgsea$nes[ix] = 0
      data$scgsea$pval[ix] = 1
      data$scgsea$fdr[ix] = 1
    }
    data$scgsea$x = data$scgsea$nes
    data$scgsea$x[data$scgsea$x < 0 | data$scgsea$fdr >= fdr.th] = 0
    data$scgsea$x = Matrix::Matrix(data = data$scgsea$nmf.h %*% t(data$scgsea$x), sparse = T)
    data$scgsea$stat = df[, c("pathway", "size")]

    # def function
    armaColSum <-utils::getFromNamespace("armaColSum", "gficf")

    data$scgsea$x = data$scgsea$x[, armaColSum(data$scgsea$x) > 0]
    if (rescale != "none") {
      if (rescale == "byGS") {
        data$scgsea$x = t(data$scgsea$x)
        data$scgsea$x = t((data$scgsea$x - rowMeans(data$scgsea$x))/apply(data$scgsea$x, 1, stats::sd))
      }
      if (rescale == "byCell") {
        data$scgsea$x = (data$scgsea$x - rowMeans(data$scgsea$x))/apply(data$scgsea$x, 1, stats::sd)
      }
    }
    return(data)
  }

  if (ncol(.expM.data) > 10000) {
    nmf.k <- 100
  }else{
    nmf.k <- 50
  }

  # perform
  .expM.data <- .runScGSEA.m(data = .expM.data,
                      geneID = "symbol",
                      minSize = minGSSize,
                      pathway.list = genesets,
                      seed = seeds,
                      nmf.k = nmf.k,
                      fdr.th = 1,
                      rescale = "none",
                      verbose = T,
                      nt = ncores)

  object[["scGSEA"]] <- SeuratObject::CreateAssayObject(counts = t(data$scgsea$x))
  print_message('scGSEA scoring completed')
}


#' Run several methods in decoupleR package to perform gene set scoring on a Seurat object
#'
#' This function performs several scoring on a given Seurat object using the specified gene sets.
#' It calculates the several scores and stores them in a new assay called "several".
#'
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for several scoring
#' @param ncores Integer specifying the number of cores to be used for parallel computation. Default is 1.
#' @author Qiong Zhang
#' @return The input Seurat object with the several scores added as a new assay named by related method.
#' @export
#'
run_decoupleR <- function(object, genesets, assay = 'RNA', slot = 'counts', methods = 'viper', scale = T){
  print_message(stringr::str_glue('Start {methods} scoring'))
  .minGS <- min(sapply(genesets, length))
  .call.method <- get(paste0('run_', methods), asNamespace('decoupleR'))
  .genesets <-  genesets
  .expM <- assay2mtx(object, assay = assay, slot = slot)
  if (methods == 'viper') {
    .expM <- as.matrix(.expM)
    if (length(.genesets) < 2) {
    .genesets <- rep(.genesets, 3)
    names(.genesets) <- paste0(names(genesets), c('', '_dup1', '_dup2'))
    .minGS <- 1
    }
  }
  .genesets.weight <- geneset2w(.genesets, 'd')
  .acts <- .call.method(mat = .expM,
                        .genesets.weight,
                        .mor='weight',
                        minsize = 0)
  .key.statistic <- c(paste0('corr_', methods), methods)
  .acts <- data.table::data.table(.acts)
  .acts <- .acts[statistic %in% .key.statistic & source %in% names(genesets), 2:4]
  .scores <- data.table::dcast(.acts, formula = source~condition, value.var='score', fill=0)
  .scores.names <- .scores$source
  .scores <- as.data.frame(.scores[, -1])
  rownames(.scores) <- .scores.names
  object[[methods]] <- SeuratObject::CreateAssayObject(counts = .scores)
  print_message(stringr::str_glue('{methods} scoring completed'))
  return(object)
}


#' Run GSignatures scoring on a Seurat object
#'
#' This function performs GSignatures scoring on a given Seurat object using the specified gene sets.
#' It calculates the GSignatures scores and stores them in a new assay called "GSignatures".
#'
#' @param object A Seurat object containing the single-cell data.
#' @param genesets A list of gene sets to be used for GSignatures scoring
#' @author Qiong Zhang
#' @return The input Seurat object with the GSignatures scores added as a new assay called "GSignatures".
#' @export
#'
run_GSignatures <- function(object, genesets, assay = "RNA", slot = "counts"){
  print_message('Start GSignatures scoring')
  object <- Gsignatures::GSignatures(object, genesets, assay = assay)
  print_message('GSignatures scoring completed')
  return(object)
}

