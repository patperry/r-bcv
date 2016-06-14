Rotate <- function(Data){
  data = Data
  #generate a random rotation
  dim <- ncol(data)
  z <- matrix(rnorm(dim * dim), dim, dim)
  qr <- qr(z)
  q <- qr.Q(qr)
  sign <- sample(c(-1, 1), dim, replace=TRUE)
  rot <- q %*% diag(sign, dim)
  #rotate the columns of the data matrix
  x_rot <- data %*% rot
  return(x_rot)
}

Uncorrelate2 <- function(Data,K){
  Fit <- stats::kmeans(Data, K, nstart = 100)
  Pool <- Data - Fit$center[Fit$cluster,]

  SVD <- svd(Pool)
  E = 10^(-5)
  result <- Data%*%SVD$v%*%(diag(1/(SVD$d+E)))
  result <- Rotate(result)

  return(result)
}

#' Choosing the number of clusters of a matrix using correlation-corrected Gabriel CV method.
#'
#'
#' Perform correlation-corrected version of Gabriel cross-validation for determining the
#' number of clusters of a matrix
#'
#' @param x The matrix to find the number of clusters.
#' @param krow	The number of row folds.
#' @param kcol	The number of column folds.
#' @param maxcenters The upper bound of searching range, i.e k searched in 1:maxcenters
#' @param classify.method The classifier used inside the algorithm,
#' with possible choice "nearest", "lda-equal" and "lda-proportions".
#' @details This function is a correlation-corrected verison of \code{\link{cv.kmeans.gabriel}}.
#' The \code{\link{cv.kmeans.gabriel}} works well when the correlation between the dimensions of
#' matrix \eqn{x} is low. It tends to overestimate the number of clusters \eqn{k} when the
#' correlation between dimension is high. This function overcomes such difficulty by
#' rotating the original matrix so that the correlation between dimensions is diminished.
#' It then calls the \code{\link{cv.kmeans.gabriel}} on the rotated matrix to find the number
#' \eqn{k}. Simulation shows this function works well when all the clusters have similar
#' covariance structure.
#'
#' @return The number of clusters chosen by the algorithm.
#'
#'
#' @author Wei Fu, Patrick O. Perry
#'
#' @examples
#' ## generate a 100x2 matrix with single cluster centered at (2,1)
#' ## the correlation rho between dimensions is 0.6
#' library(MASS)
#' sigma = matrix(c(1,0.6,0.6,1),nrow = 2, ncol=2)
#' x <- mvrnorm(n = 100, mu = c(2,1), sigma)
#'
#' ## the k returned by cv.kmeans.gabriel function
#' cv.kmeans.gabriel(x, 5, 2, maxcenters=10, classify.method="nearest")$$centers
#'
#' ## the k returned by correlation-corrected gabriel.cor.correct function
#' gabriel.cor.correct(x, 5, 2, maxcenters=10, classify.method="nearest")
#'
#' @export
gabriel.cor.correct <- function(x, krow = 5, kcol = 2, maxcenters = 10, classify.method = "nearest"){

  K.original <- bcv::cv.kmeans.gabriel(x, krow, kcol, maxcenters, classify.method)$centers

  DATA <- Uncorrelate2(x,K.original)

  K2 <- bcv::cv.kmeans.gabriel(DATA, krow, kcol, maxcenters, classify.method)$centers

  return(K2)
}
