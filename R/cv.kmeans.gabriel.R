cluster_kmeans <- function(x, centers, seed = 2651513, iter.max = 100,
                           nstart = 100, ...)
{
  if (!is.na(seed)) {
    if ((exists0 <- exists(".Random.seed", envir=globalenv()))) {
      seed0 <- .Random.seed
    }
    set.seed(seed)
  }

  cl <- stats::kmeans(x, centers, iter.max=iter.max, nstart=nstart, ...)

  if (!is.na(seed)) {
    if (exists0) {
      assign(".Random.seed", seed0, envir=globalenv())
    } else if (exists(".Random.seed", envir=globalenv())) {
      rm(".Random.seed", envir=globalenv())
    }
  }

  cl
}

classify_lda <- function(x, grouping, prior=c("proportions", "equal"))
{
  g <- as.factor(grouping)
  prior <- match.arg(prior)

  n <- length(grouping)
  ng <- nlevels(g)

  if (prior == "equal") {
    prior <- rep(1/ng, ng)
  } else {
    counts <- as.vector(table(g))
    prior <- counts / n
  }

  MASS::lda(x, g, prior)
}

classify_nearest <- function(x, grouping)
{
  g <- as.factor(grouping)
  nc <- tabulate(g)

  gmat <- model.matrix( ~ g - 1)
  sums <- t(gmat) %*% x
  rownames(sums) <- levels(g)
  means <- sums / nc

  structure(list(means=means, levels=levels(g)),
            class="classify_nearest")
}


predict.classify_nearest <- function(object, x, ...)
{
  means <- object$means
  levels <- object$levels

  d <- scale(x %*% t(means), center=0.5 * rowSums(means * means),
             scale=FALSE)
  class <- factor(levels[apply(d, 1, nnet::which.is.max)],
                  levels=levels)
  list(class=class)
}

#' Choosing the number of clusters of a matrix using Gabriel CV method.
#'
#'
#' Perform Gabriel-style cross-validation for determining the
#' number of clusters of a matrix.
#'
#' @param x The matrix to find the number of clusters.
#' @param krow	The number of row folds.
#' @param kcol	The number of column folds.
#' @param maxcenters The upper bound of searching range, i.e k searched in 1:maxcenters
#' @param classify.method The classifier used inside the algorithm,
#' with possible choice "nearest", "lda-equal" and "lda-proportions".
#'
#' @details Through gabriel cross-validation, the data matrix has been decomposed
#' into four blocks:
#' \deqn{\mathbf{X} = \begin{bmatrix} X_{train} & Y_{train} \\ X_{test}  & Y_{test} \end{bmatrix},}
#' where \eqn{X} and \eqn{Y} are formed by CV on the column. Note here \eqn{Y} is not the
#' real responses since the problem is a unsupervised learning problem. Instead, \eqn{Y}
#' is the leaf-out part of each CV (on the column ) replicate while \eqn{X} is the leaf-in part.
#'
#' For each CV replicate, the algorithm applies kmeans with parameter k = 1, ..., maxcenters on
#' \eqn{Y_{train}}, which yields k clusters' assigment; builds a classifier using the training
#' data (i.e. \eqn{X_{train}} and \eqn{Y_{train}}) where \eqn{Y_{train}} has been replaced by
#' the cluster assigments; predicts the cluster assigments of \eqn{X_{test}} using the
#' classifier just trained; and calculates the prediction error by measuring the difference
#' between \eqn{Y_{test}} and the predicted cluster assigments (each cluster
#' assigment is replaced by the correspoding cluster center). The overall prediction
#' error is average across all CV replicates and the algorithm picks the k
#' that minimizes such error.
#'
#' @return
#' \item{msep}{The mean square error of prediction (MSEP); this is a
#'  matrix whose columns contain the mean square errors in the predictions of the holdout
#'  sets for k = 1, ..., maxcenters across the different replicates.}
#'\item{centers}{The number of clusters chosen by the algorithm.}
#'
#' @author Wei Fu, Patrick O. Perry
#'
#' @examples
#' ## generate a 40x2 matrix with two clusters centered at (4,4) and (0,0)
#' cluster1 <- matrix(c(4,4),nrow = 20, ncol = 2)+ matrix(rnorm(40),nrow = 20, ncol = 2)
#' cluster2 <- matrix(c(0,0),nrow = 20, ncol = 2)+ matrix(rnorm(40),nrow = 20, ncol = 2)
#' x = rbind(cluster1,cluster2 )
#'
#' ## cluster the matrix with gabriel CV method with 5 by 2 cross-validation
#' cv.kmeans.gabriel(x, 5, 2, maxcenters=10, classify.method="nearest")
#'
#' @export
cv.kmeans.gabriel <- function(x, krow = 5, kcol = 2, maxcenters = 10,
                              classify.method = c("nearest", "lda-equal", "lda-proportions"))
{
  classify.method <- match.arg(classify.method)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (n < 2)
    stop("x should have at least two rows")
  if (p < 2)
    stop("x should have at least two columns")
  if ((krow > n) || (krow <= 1))
    stop("krow outside allowable range")
  if ((kcol > p) || (kcol <= 1))
    stop("kcol outside allowable range")
  if (maxcenters <= 0)
    stop("maxcenters should be positive")

  krow.o <- krow
  krow <- bcv:::round.fold(n, krow)
  kcol.o <- kcol
  kcol <- bcv:::round.fold(p, kcol)
  if (krow != krow.o)
    warning("krow has been set to ", krow)
  if (kcol != kcol.o)
    warning("kcol has been set to ", kcol)

  s.r <- bcv:::choose.sets(n, krow)
  s.c <- bcv:::choose.sets(p, kcol)
  n0 <- n - max(table(s.r))
  p0 <- p - max(table(s.c))
  maxcenters.o <- maxcenters
  maxcenters <- min(n0, round(maxcenters.o))
  if (!missing(maxcenters) && maxcenters != maxcenters.o)
    warning("maxcenters has been set to ", maxcenters)

  if (classify.method == "nearest") {
    classify <- classify_nearest
  } else if (classify.method == "lda-equal") {
    classify <- function(x, grouping)
      classify_lda(x, grouping, prior="equal")
  } else if (classify.method == "lda-proportions") {
    classify <- function(x, grouping)
      classify_lda(x, grouping, prior="proportions")
  } else if (classify.method == "svm") {
    stop("not yet implemented")
  }

  msep <- matrix(NA, krow * kcol, maxcenters)

  for (k in seq_len(krow)) {
    for (l in seq_len(kcol)) {
      test <- s.r == k
      train <- !test
      response <- s.c == l
      predictor <- !response

      for (centers in seq_len(maxcenters)) {
        if (centers == 1) {
          fit <- colMeans(x[train,response,drop=FALSE])
          err <- scale(x[test,response,drop=FALSE], center=fit, scale=FALSE)
        } else {
          cl <- cluster_kmeans(x[train,response,drop=FALSE], centers=centers)
          cluster <- factor(cl$cluster, levels=seq_len(centers))
          fit <- classify(x[train,predictor,drop=FALSE], cluster)
          pred <- predict(fit, x[test,predictor,drop=FALSE])$class
          err <- x[test,response,drop=FALSE] - cl$centers[pred,,drop=FALSE]
        }

        msep[k + (l - 1) * krow, centers] <- mean(err^2)
      }
    }
  }

  msep.mean <- colMeans(msep)
  centers <- which.min(msep.mean)
  list(msep = msep, centers=centers)
}
