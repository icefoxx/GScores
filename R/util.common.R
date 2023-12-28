
#' Simulate the matrix based on the sparsity
#'
#' @param Mat Expression or binary matrix
#' @param prop The ratio of simulation cells
#' @param seed Seed used
#' @param addname Whether to add the prefix of 'pseudo' to simulation cells
#' @param merge Whether to merge with original matrix
#' @param random Use global sparsity or random select the cell sparsity
#' @name simSparseMat
#' @return Simulation matrix with similar sparsity as original one
#' @author Qiong Zhang
#' @export
#'
simSparseMat <- function(Mat, prop = 0.1, seed = 10086, addname = T, merge = T, random = T) {
  require(Matrix, quietly = T)
  set.seed(seed)
  .mat.spa <- diff(Mat@p)
  .ncol <- length(.mat.spa)
  .nrow <- nrow(Mat)
  .sample.n <- .ncol * prop
  .sample.n <- ifelse(.sample.n < 300, 300, .sample.n)
  .sparsity <- .mat.spa/.nrow
  if (isFALSE(random)) .sparsity <- rep(sum(.mat.spa)/(.ncol * .nrow), .sample.n)
  .col.sample <- sample(.ncol, .sample.n)
  .s <- lapply(1:.sample.n, function(x){
    .col.sparsity <- .sparsity[x]
    .vs <- as(sample(c(1,0), .nrow, prob = c(.col.sparsity, 1-.col.sparsity), replace = T), 'sparseVector')
  })
  .ms <- list2SparseMat(.s)
  if(isTRUE(addname)){
    rownames(.ms) <- rownames(Mat)
    colnames(.ms) <- paste0("pseudo", .col.sample)
  }
  if(isTRUE(merge)) .ms <- cbind(Mat, .ms)
  return(.ms)
}


#' Transfer the sparseVector list to Matrix format
#'
#' @param list sparseVector list
#' @name list2SparseMat
#' @return Sparse matrix
#' @author Qiong Zhang
#'
list2SparseMat <- function (list) {
  require(Matrix, quietly = T)
  .vl <- lapply(list, as, 'sparseVector')
  .vln <- unique(sapply(.vl, length))
  stopifnot( length(.vln) == 1 )
  return( Matrix::sparseMatrix(
    x = unlist(lapply(.vl, slot, "x")),
    i = unlist(lapply(.vl, slot, "i")),
    p = c(0,cumsum(sapply(.vl, function(x){length(x@x)}))),
    dims=c(.vln,length(.vl))
  ))
}


#' Random drop n percent of expressed gene for cells and return the matrix
#'
#' @param m Expression matrix of gene(row) X cell (column)
#' @param n The percentage of genes should be dropped
#'
#' @return Expression matrix with dropped genes
#' @author Qiong Zhang
#' @export
#'
get_random_loss_sparse_matrix <- function(m, n) {
  require(Matrix)
  .s <- sample(1:1e9, 1)
  set.seed(.s)
  .g <- nrow(m)
  .m <- lapply(1:ncol(m), function(c, g){
    .exped.g <- which(m[, c] > 0)
    .g.n <- length(.exped.g)
    .s.n <- as.integer((1-n/100) * .g.n)
    .sampled_g <- sample(.exped.g, .s.n)
    Matrix::sparseVector( x=rep(1,.s.n), i=.sampled_g, length=g )
  }, .g)
  .ms <- list2SparseMat(.m)
  rownames(.ms) <- rownames(m)
  colnames(.ms) <- paste("pseudo",n, .s , colnames(m), sep = '.')
  return(.ms*m)
}
