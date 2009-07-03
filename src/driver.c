/* driver.c
 */
 
#include <string.h>
 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "bcv-partition.h"
#include "bcv-svd.h"
#include "driver.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef struct _perm_t
{
    bcv_partition_t *row_part;
    bcv_partition_t *col_part;
    int m, n;
    int *ir, *jc;
} perm_t;


void 
perm_init (perm_t *perm, int M, int N, int K, int L, int *s_r, int *s_c)
{
    perm->row_part = bcv_partition_alloc (M, K);
    perm->col_part = bcv_partition_alloc (N, L);
    if (!perm->row_part || !perm->col_part) 
        error ("Could not allocate enough memory to BCV"
               " an %d-by-%d matrix.", M, N);

    bcv_partition_init (perm->row_part, M, K, s_r);
    bcv_partition_init (perm->col_part, N, L, s_c);

    perm->ir  = (int *) R_alloc (M, sizeof (int));
    perm->jc  = (int *) R_alloc (N, sizeof (int));
}

void
perm_deinit (perm_t *perm)
{
    bcv_partition_free (perm->row_part);
    bcv_partition_free (perm->col_part);
}

void
perm_select (perm_t *perm, int k, int l)
{
    perm->m = bcv_partition_get_perm (perm->row_part, k, perm->ir);
    perm->n = bcv_partition_get_perm (perm->col_part, l, perm->jc);
}


SEXP 
driver_svd (SEXP xx, SEXP KK, SEXP LL, SEXP max_rank, SEXP s_r, SEXP s_c)
{
    double *x_data;
    bcv_index_t M, N, K, L, k, i, j, kmax;
    bcv_svd_t *bcv;
    bcv_holdin_t holdin;
    SEXP rss_R, dim;
    double *rss;
    perm_t perm;
    
    if (!isMatrix (xx) || !isNumeric (xx))
        error ("x should be a matrix");


    x_data = NUMERIC_POINTER (xx);
    M      = INTEGER (getAttrib (xx, R_DimSymbol))[0];
    N      = INTEGER (getAttrib (xx, R_DimSymbol))[1];
    K      = asInteger (KK);
    L      = asInteger (LL);
    kmax   = asInteger (max_rank);

    bcv_matrix_t x = { M, N, x_data, MAX (M,1) };
    perm_init (&perm, M, N, K, L, INTEGER_POINTER (s_r), INTEGER_POINTER (s_c));
    
    PROTECT (rss_R = allocVector (REALSXP, (kmax + 1) * K * L));
    rss = NUMERIC_POINTER (rss_R);
    bcv_holdin_t max_holdin = { M, N };
    bcv = bcv_svd_alloc (max_holdin, M, N);
    
    if (!bcv)
        error ("Could not allocate bcv_svd_t for size (%d,%d)", M, N);

    for (j = 0; j < L; j++)
    {
        for (i = 0; i < K; i++)
        {
            R_CheckUserInterrupt ();
            
            perm_select (&perm, i, j);
            if (perm.n <= 0 || perm.n >= N)
            {
                REprintf ("L: %d  N: %d  n: %d\n", L, N, perm.n);
            }
            holdin.m = perm.m;
            holdin.n = perm.n;
            bcv_svd_initp (bcv, holdin, &x, perm.ir, perm.jc);
            /* TODO: check for error return */

            *rss++ = bcv_svd_get_resid_rss (bcv);
            for (k = 0; k < kmax; k++)
            {
                bcv_svd_update_resid (bcv, 1.0, k);
                *rss++ = bcv_svd_get_resid_rss (bcv);            
            }
        }
    }

    bcv_svd_free (bcv);
    perm_deinit (&perm);

    PROTECT (dim = allocVector (INTSXP, 2));
    INTEGER (dim) [0] = (kmax + 1);
    INTEGER (dim) [1] = K * L;
    setAttrib (rss_R, R_DimSymbol, dim);
    
    UNPROTECT (2);
    return rss_R;
}
