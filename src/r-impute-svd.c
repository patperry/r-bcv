#include <assert.h>
 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "bcv-svd-impute.h"
#include "r-impute-svd.h"

static bcv_index_t
bcv_count_missing (const bcv_matrix_t *a);

static void
bcv_find_missing (const bcv_matrix_t *a, bcv_index_t *indices);

static bcv_index_t
bcv_svd_impute (bcv_svd_impute_t *impute,
                bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                const bcv_index_t *indices, bcv_index_t num_indices,
                bcv_index_t k, double tol, bcv_index_t maxiter);


SEXP 
R_impute_svd (SEXP xx, SEXP kk, SEXP toltol, SEXP maxitermaxiter)
{
    bcv_index_t m, n, k, maxiter, num_indices, iter;
    bcv_index_t *indices = NULL;
    bcv_svd_impute_t *impute = NULL;
    double tol, rss;
    size_t impute_size;
    SEXP xhatxhat, dimdim, rssrss, iteriter, res;

    m       = INTEGER (getAttrib (xx, R_DimSymbol))[0];
    n       = INTEGER (getAttrib (xx, R_DimSymbol))[1];
    k       = asInteger (kk);
    tol     = asReal (toltol);
    maxiter = asInteger (maxitermaxiter);

    PROTECT (xhatxhat = allocVector (REALSXP, m * n));

    bcv_matrix_t x    = { m, n, NUMERIC_POINTER (xx),       BCV_MAX (m,1) };
    bcv_matrix_t xhat = { m, n, NUMERIC_POINTER (xhatxhat), BCV_MAX (m,1) };

    num_indices = bcv_count_missing (&x);
    if (num_indices > 0)
    {
        indices = (void *) R_alloc (num_indices, sizeof (bcv_index_t));
        bcv_find_missing (&x, indices);
    }
        
    impute_size = bcv_svd_impute_size (m, n);

    if (!impute_size)
        error ("could not allocate enough memory to impute missing values"
               " for a %d-by-%d matrix", m, n);

    impute = (void *) R_alloc (impute_size, 1);
    iter   = bcv_svd_impute (impute, &xhat, &x, indices, num_indices, k, tol,
                             maxiter);
    rss    = bcv_svd_impute_get_rss (impute);

    PROTECT (dimdim = allocVector (INTSXP, 2));
    INTEGER (dimdim) [0] = m;
    INTEGER (dimdim) [1] = n;
    setAttrib (xhatxhat, R_DimSymbol, dimdim);

    PROTECT (rssrss = allocVector (REALSXP, 1));
    REAL (rssrss) [0] = rss;
    
    PROTECT (iteriter = allocVector (INTSXP, 1));
    INTEGER (iteriter) [0] = iter;

    PROTECT (res = allocVector (VECSXP, 3));
    SET_VECTOR_ELT (res, 0, xhatxhat);
    SET_VECTOR_ELT (res, 1, rssrss);
    SET_VECTOR_ELT (res, 2, iteriter);
    
    UNPROTECT (5);
    return res;
}


bcv_index_t
bcv_count_missing (const bcv_matrix_t *a)
{
    bcv_index_t i, m, n, lda, result = 0;
    double *data;
    
    _bcv_assert_valid_matrix (a);
    
    m    = a->m;
    n    = a->n;
    data = a->data;
    lda  = a->lda;
    
    /* not implemented for non-contiguous storage */
    assert (lda == BCV_MAX (m,1)); 
    
    for (i = 0; i < m*n; i++)
    {
        if (ISNA (data[i]))
            result++;
    }
    
    return result;
}


void
bcv_find_missing (const bcv_matrix_t *a, bcv_index_t *indices)
{
    bcv_index_t i, m, n, lda;
    double *data;

    _bcv_assert_valid_matrix (a);
    assert (indices);

    m    = a->m;
    n    = a->n;
    data = a->data;
    lda  = a->lda;
    
    /* not implemented for non-contiguous storage */
    assert (lda == BCV_MAX (m,1)); 

    for (i = 0; i < m*n; i++)
    {
        if (ISNA (data[i]))
            *indices++ = i;
    }
}


bcv_index_t
bcv_svd_impute (bcv_svd_impute_t *impute,
                bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                const bcv_index_t *indices, bcv_index_t num_indices,
                bcv_index_t k, double tol, bcv_index_t maxiter)
{
    bcv_error_t err = 0;
    bcv_index_t iter = 0;
    double rss0, rss1 = BCV_DBL_POSINF, delta;
    
    bcv_svd_impute_init (impute, xhat, x, indices, num_indices);

    /* quick return for empty matrices */
    if (x->m == 0 || x->n == 0)
        return 0;
        
    do
    {
        R_CheckUserInterrupt ();
        
        rss0 = rss1;
        iter++;
        
        err = bcv_svd_impute_step (impute, xhat, x, indices, num_indices, k);
        rss1  = bcv_svd_impute_get_rss (impute);
        delta = fabs (rss1 - rss0) / (BCV_DBL_EPSILON + rss1);
    }
    while (!err && delta > tol && iter < maxiter);
    
    if (err)
        error ("Error computing the SVD of the imputed matrix.");
    
    return iter;
}
