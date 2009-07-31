\name{impute.svd}
\alias{impute.svd}
\title{
  Missing value imputation via the SVDImpute algorithm
}
\description{
  Given a matrix with missing values, impute the missing entries using
  a low-rank SVD approximation estimated by the EM algorithm.
}
\usage{
  impute.svd(x, k = min(n, p), tol = max(n, p) * 1e-10, maxiter = 100)
}
\arguments{
  \item{x}{a matrix to impute the missing entries of.}
  \item{k}{the rank of the SVD approximation.}
  \item{tol}{the convergence tolerance for the EM algorithm.}
  \item{maxiter}{the maximum number of EM steps to take.}
}
\details{
  Impute the missing values of \code{x} as follows:  First, initialize
  all \code{NA} values to the column means, or \code{0} if all entries
  in the column are missing.  Then, until convergence, compute the
  first \code{k} terms of the SVD of the completed matrix.  Replace the
  previously missing values with their approximations from the SVD, and 
  compute the RSS between the \emph{non-missing} values and the SVD.
  
  Declare convergence if 
    \code{
        abs(rss0 - rss1) / (.Machine$double.eps + rss1) < tol
    },
  where \code{rss0} and \code{rss1} are the RSS values computed from
  successive iterations.  Stop early after \code{maxiter} iterations
  and issue a warning.
}
\value{
  \item{x}{the completed version of the matrix.}
  \item{rss}{the sum of squares between the SVD approximation and the 
    non-missing values in \code{x}.}
  \item{iter}{the number of EM iterations before algorithm stopped.}
}
\references{
  Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., 
    Tibshirani, R., Botstein, D. and Altman, R.B. (2001).  
    Missing value estimation methods for DNA microarrays.
    \emph{Bioinformatics} \bold{17}(6), 520--525.
}
\author{Patrick O. Perry}

\seealso{
    \code{\link{cv.svd.wold}}
}

\examples{
  # Generate a matrix with missing entries    
  n <- 20
  p <- 10
  u <- rnorm( n )
  v <- rnorm( p )
  xfull <- u \%*\% rbind( v ) + rnorm( n*p )
  miss  <- sample( seq_len( n*p ), n )
  x       <- xfull
  x[miss] <- NA
      
  # impute the missing entries with a rank-1 SVD approximation
  xhat <- impute.svd( x, 1 )$x   
  
  # compute the prediction error for the missing entries
  sum( ( xfull-xhat )^2 )
}
