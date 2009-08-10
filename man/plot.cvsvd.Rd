\name{plot.cvsvd}
\alias{plot.cvsvd}
\title{
  Plot the Result of an SVD Cross-Validation
}
\description{
  Plot the result of \code{\link{cv.svd.gabriel}} or 
  \code{\link{cv.svd.wold}}, optionally with error bars.
}
\usage{
  \method{plot}{cvsvd}(x, errorbars = TRUE, add = FALSE,
       xlab = "Rank", ylab = "Mean Sq. Prediction Error", 
       col = "blue", col.errorbars = "gray50", ...)
}
\arguments{
  \item{x}{the result of a \code{\link{cv.svd.gabriel}} or 
           \code{link{cv.svd.wold}} computation.}
  \item{errorbars}{indicates whether or not to add error bars.}
  \item{col}{the color to use for showing prediction error.}
  \item{col.errorbars}{the color to use for the error bars.}
  \item{add}{indicates whether or not to add to the current plot.}
  \item{xlab}{the label for the x axis.}
  \item{ylab}{the label for the y axis.}
  \item{\dots}{additional arguments for \code{plot}.}
}
\details{
  Plot the result of \code{\link{cv.svd.gabriel}} or 
  \code{\link{cv.svd.wold}}.  This plots a the estimated prediction
  error as a function of rank, optionally with error bars.
  
  If \code{add} is \code{TRUE}, the current plot is not cleared.
}
\author{Patrick O. Perry}

\seealso{
  \code{\link{cv.svd.gabriel}}, 
  \code{\link{cv.svd.wold}},
  \code{\link{print.cvsvd}}
  \code{\link{summary.cvsvd}}
}
\examples{
  # generate a rank-2 matrix plus noise
  n <- 50; p <- 20; k <- 2
  u <- matrix( rnorm( n*k ), n, k )
  v <- matrix( rnorm( p*k ), p, k )
  e <- matrix( rnorm( n*p ), n, p )
  x <- u \%*\% t(v) + e
  
  # perform 5-fold Wold-style cross-validtion
  cvw <- cv.svd.wold( x, 5, maxrank=10 )
  
  # perform (2,2)-fold Gabriel-style cross-validation
  cvg <- cv.svd.gabriel( x, 2, 2, maxrank=10 )
  
  # plot the results
  par( mfrow=c(2,1) )
  plot( cvw, main="Wold-style CV")
  plot( cvg, main="Gabriel-style CV")
}
