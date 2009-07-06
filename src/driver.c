/* driver.c
 */
 
#include <assert.h>
#include <string.h>
 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "bcv-partition.h"
#include "bcv-svd-gabriel.h"
#include "driver.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))


SEXP 
driver_svd (SEXP xx, SEXP KK, SEXP LL, SEXP max_rank, SEXP s_r, SEXP s_c)
{
    bcv_index_t M, N, K, L, i, j, kmax;
    bcv_error_t bcv_error;
    SEXP rss_R, dim;
    double *rss;
    
    if (!isMatrix (xx) || !isNumeric (xx))
        error ("x should be a matrix");

    M    = INTEGER (getAttrib (xx, R_DimSymbol))[0];
    N    = INTEGER (getAttrib (xx, R_DimSymbol))[1];
    K    = asInteger (KK);
    L    = asInteger (LL);
    kmax = asInteger (max_rank);

    if (kmax < 0)
        error ("max_rank should be non-negative");

    PROTECT (rss_R = allocVector (REALSXP, (kmax + 1) * K * L));
    rss = NUMERIC_POINTER (rss_R);

    bcv_matrix_t x           = { M, N, NUMERIC_POINTER (xx), MAX (M,1) };
    bcv_partition_t row_part = { M, K, INTEGER_POINTER (s_r) };
    bcv_partition_t col_part = { N, L, INTEGER_POINTER (s_c) };
    bcv_gabriel_holdin_t max_holdin = { M, N };
    
    bcv_svd_gabriel_t *bcv = bcv_svd_gabriel_alloc (max_holdin, M, N);
    
    if (!bcv)
        error ("Could not allocate enough memory for Gabriel "
               " cross-validation of a %d-by-%d matrix", M, N);

    bcv_svd_gabriel_init (bcv, &x, &row_part, &col_part);

    for (j = 0; j < L; j++)
    {
        for (i = 0; i < K; i++)
        {
            R_CheckUserInterrupt ();
            bcv_error = bcv_svd_gabriel_get_rss (bcv, i, j, rss, kmax);
            
            if (bcv_error)
                error ("The SVD algorithm did not converge for the (%d,%d)"
                       " holdin", i, j);
            
            rss += kmax + 1;
        }
    }

    bcv_svd_gabriel_free (bcv);

    PROTECT (dim = allocVector (INTSXP, 2));
    INTEGER (dim) [0] = (kmax + 1);
    INTEGER (dim) [1] = K * L;
    setAttrib (rss_R, R_DimSymbol, dim);
    
    UNPROTECT (2);
    return rss_R;
}
