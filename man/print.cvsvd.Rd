\name{print.cvsvd}
\alias{print.cvsvd}
\title{
  Print the Result of an SVD Cross-Validation
}
\description{
  Print the result of \code{\link{cv.svd.gabriel}} or 
  \code{\link{cv.svd.wold}}.
}
\usage{
  \method{print}{cvsvd}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{the result of a \code{\link{cv.svd.gabriel}} or 
           \code{\link{cv.svd.wold}} computation.}
  \item{digits}{the digits of precision to show in the output.}
  \item{\dots}{additional arguments to \code{print}.}
}
\details{
  Print a table of the estimated prediction errors and the standard errors
  of the estimate.  Put an asterisk (\code{*}) next to the minimum and a
  plus (\code{+}) next to the "one standard error rule" choice.
}
\author{Patrick O. Perry}
\seealso{
  \code{\link{cv.svd.gabriel}}, 
  \code{\link{cv.svd.wold}},
  \code{\link{plot.cvsvd}}
  \code{\link{summary.cvsvd}}
}
