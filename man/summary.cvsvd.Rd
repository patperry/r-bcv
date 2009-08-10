\name{summary.cvsvd}
\alias{summary.cvsvd}
\title{
  Summarize the Result of an SVD Cross-Validation
}
\description{
  Summarize the result of \code{\link{cv.svd.gabriel}} or 
  \code{\link{cv.svd.wold}}.
}
\usage{
  \method{summary}{cvsvd}(object, ...)
}
\arguments{
  \item{object}{the result of a \code{\link{cv.svd.gabriel}} or 
        \code{\link{cv.svd.wold}} computation.}
  \item{\dots}{additional arguments to \code{summary}.}
}
\details{
  Print a table of the estimated prediction errors and the standard errors
  of the estimate.  Put an asterisk (\code{*}) next to the minimum and a
  plus (\code{+}) next to the "one standard error rule" choice.
}
\value{
  \item{nfolds }{the number of cross-validation folds}
  \item{maxrank }{the maximum rank for which prediction error is estimated.} 
  \item{msep.mean }{the average mean square error of prediction (MSEP)
    accross all folds for ranks 0, 1, ...,  \code{maxrank}.}
  \item{msep.se }{the standard errors of the \code{msep.mean} estimates.}
  \item{rank.best }{the rank with the minimum \code{msep.mean} value.}
  \item{rank.1se }{the smallest rank within one standard error of the
                   minimum \code{msep.mean} value.}
}
\author{Patrick O. Perry}
\seealso{
  \code{\link{cv.svd.gabriel}}, 
  \code{\link{cv.svd.wold}},
  \code{\link{plot.cvsvd}}
  \code{\link{print.cvsvd}}
}
