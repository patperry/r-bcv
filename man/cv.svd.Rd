\name{cv.svd}
\alias{cv.svd}
\alias{cv.svd.gabriel}
\alias{cv.svd.wold}
\title{
  Cross-Validation for choosing the rank of an SVD approximation.
}
\description{
  Perform Wold- or Gabriel-style cross-validation for determining the
  appropriate rank SVD approximation of a matrix.
}
\usage{
  cv.svd.gabriel(x, krow = 2, kcol = 2, 
                 maxrank = floor(min(n - n/krow, p - p/kcol)))
                 
  cv.svd.wold(x, k = 5, maxrank = 20, tol = 1e-4, maxiter = 20)
}
\arguments{
  \item{x}{the matrix to cross-validate.}
  \item{k}{the number of folds (for Wold-style CV).}
  \item{krow}{the number of row folds (for Gabriel-style CV).}
  \item{kcol}{the number of column folds (for Gabriel-style CV).}
  \item{maxrank}{the maximum rank to cross-validate up to.}
  \item{tol}{the convergence tolerance for \code{\link{impute.svd}}.}
  \item{maxiter}{the maximum number of iterations for 
                 \code{\link{impute.svd}}.}
}
\details{
  These functions are for cross-validating the SVD of a matrix.  They
  assume a model $X = U D V' + E$ with the terms being signal
  and noise, and try to find the best rank to truncate the SVD of \code{x}
  at for minimizing prediction error.  Here, prediction error is measured
  as sum of squares of residuals between the truncated SVD and the
  signal part.
  
  For both types of cross-validation, in each replicate we leave out part
  of the matrix, fit an SVD approximation to the left-in part, and measure
  prediction error on the left-out part.  
  
  In Wold-style cross-validation, the holdout set is "speckled", a random set 
  of elements in the matrix.  The missing elements are predicted using 
  \code{\link{impute.svd}}.
  
  In Gabriel-style cross-validation, the holdout set is "blocked".  We 
  permute the rows and columns of the matrix, and leave out the lower-right
  block.  We use a modified Schur-complement to predict the held-out block.
  In Gabriel-style, there are \code{krow*kcol} total folds.
}
\value{
  \item{call }{the function call}
  \item{msep }{the mean square error of prediction (MSEP); this is a matrix 
    whose columns contain the mean square errors in the predictions of the
    holdout sets for ranks 0, 1, ...,  \code{maxrank} across the  different 
    replicates.}
  \item{maxrank }{the maximum rank for which prediction error is estimated; 
    this is equal to \code{nrow(msep)+1}.}
  
  \item{krow }{the number of row folds (for Gabriel-style only).}
  \item{kcol }{the number of column folds (for Gabriel-style only).}
  \item{rowsets }{the partition of rows into \code{krow} holdout sets
                  (for Gabriel-style only).}
  \item{colsets }{the partition of the columns into \code{kcol} holdout sets
                  (for Gabriel-style only).}
  
  \item{k }{the number of folds (for Wold-style only).}
  \item{sets }{the partition of indices into \code{k} holdout sets
               (for Wold-style only).}
}
\references{
  Gabriel, K.R. (2002). Le biplot--outil d'explaration de 
    \enc{données}{donnees} multidimensionelles.
    \emph{J. Roy. Stat. Soc. Series B} \bold{40} 186--196.
  
  Owen, A.B. and Perry, P.O. (2009). Bi-cross-validation of the SVD and the 
     non-negative matrix factorization.
     \emph{Annals of Applied Statistics} \bold{3}(2) 564--594.
    
  Wold, S. (1978).  Cross-validitory estimation of the number of components
    in factor and principal components models.  
    \emph{Technometrics} \bold{20}(4) 397--405.
}
\author{Patrick O. Perry}
\note{
  Gabriel's version of cross-validation was for leaving out a single element 
  of the matrix, which corresponds to n-by-p-fold.  Owen and Perry generalized 
  Gabriel's idea to larger holdouts, showing that 2-by-2-fold cross-validation 
  often works better.
  
  Wold's original version of cross-validation did not use the EM algorithm
  to estimate the SVD.  He recommend using the NIPALS algorithm instead, which
  has since faded into obscurity.
  
  Wold-style cross-validation takes a lot more computation than Gabriel-style.  
  The \code{maxrank}, \code{tol}, and \code{maxiter} have been chosen to give 
  up some accuracy in the name of expediency.  They may need to be adjusted to 
  get the best results.
}

\seealso{
    \code{\link{impute.svd}}, 
    \code{\link{plot.cvsvd}},
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
  (cvw <- cv.svd.wold( x, 5, maxrank=10 ))
  
  # perform (2,2)-fold Gabriel-style cross-validation
  (cvg <- cv.svd.gabriel( x, 2, 2, maxrank=10 ))
}
