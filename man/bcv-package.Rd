\name{bcv-package}
\alias{bcv-package}
\alias{bcv}
\docType{package}
\title{
    Cross-Validation for the SVD (Bi-Cross-Validation)
}
\description{
    This package implements methods for choosing the rank of an SVD 
    approximation via cross validation.  It provides both Gabriel-style 
    "block" holdouts and Wold-style "speckled" holdouts.  Also included is an 
    implementation of the SVDImpute algorithm.
}
\details{
    \tabular{ll}{
    Package: \tab bcv\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2009-08-15\cr
    License: \tab BSD3\cr
    }
    
    Basic usage is to call either \code{\link{cv.svd.gabriel}} or
    \code{\link{cv.svd.wold}}.
}

\author{Patrick O. Perry <patperry@gmail.com>}

\keyword{package}
\seealso{
    \code{\link{impute.svd}},
    \code{\link{cv.svd.gabriel}},
    \code{\link{cv.svd.wold}},
    \code{\link{plot.cvsvd}},
    \code{\link{print.cvsvd}}
    \code{\link{summary.cvsvd}}
}
