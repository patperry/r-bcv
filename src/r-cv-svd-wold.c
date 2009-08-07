 
#include <assert.h>
 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "bcv-svd-wold.h"
#include "r-cv-svd-wold.h"


SEXP 
R_cv_svd_wold (SEXP xx, SEXP kk, SEXP maxrankmaxrank, SEXP toltol, 
               SEXP maxitermaxiter, SEXP setssets)
{
    bcv_error_t err = 0;
    bcv_index_t m, n, i, k, maxiter, maxrank;
    bcv_svd_wold_t *wold = NULL;
    double tol, *msep;
    SEXP msepmsep, dimdim;

    m       = INTEGER (getAttrib (xx, R_DimSymbol))[0];
    n       = INTEGER (getAttrib (xx, R_DimSymbol))[1];
    k       = asInteger (kk);
    maxrank = asInteger (maxrankmaxrank);
    tol     = asReal (toltol);
    maxiter = asInteger (maxitermaxiter);

    PROTECT (msepmsep = allocVector (REALSXP, (maxrank + 1) * k));
    PROTECT (dimdim   = allocVector (INTSXP, 2));
    INTEGER (dimdim) [0] = maxrank + 1;
    INTEGER (dimdim) [1] = k;
    setAttrib (msepmsep, R_DimSymbol, dimdim);
    msep = NUMERIC_POINTER (msepmsep);

    bcv_matrix_t x       = { m, n, NUMERIC_POINTER (xx), BCV_MAX (m,1) };
    bcv_partition_t part = { m*n, k, INTEGER_POINTER (setssets) };

    wold = bcv_svd_wold_alloc (m*n, m, n);
    if (!wold)
        error ("could not allocate enough memory for Wold "
               " cross-validation of a %d-by-%d matrix", m, n);

    bcv_svd_wold_init (wold, &x, &part);
    
    for (i = 0; i < k; i++)
    {
        R_CheckUserInterrupt ();
        err = bcv_svd_wold_get_msep (wold, i, tol, maxiter, msep, maxrank);
            
        if (err)
            error ("the SVD algorithm did not converge for the (%d)"
                   " holdout", i);
        
        msep += maxrank + 1;
    }
    
    bcv_svd_wold_free (wold);
    
    UNPROTECT (2);
    return msepmsep;
}
