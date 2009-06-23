/* driver.c
 */
 
#include <string.h>
 
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include "bcv-svd.h"
#include "driver.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef struct _perm_t
{
    int M, N, K, L, m, n;
    int *s_r, *s_c;
    int *ir, *jc;
} perm_t;

void
perm_debug (perm_t *perm)
{
    int i, j, M = perm->M, N = perm->N;
    
    REprintf ("perm_t {\n"
              "    M   : %d\n"
              "    N   : %d\n"
              "    K   : %d\n"
              "    L   : %d\n"
              "    m   : %d\n"
              "    n   : %d\n",
              perm->M, perm->N, perm->K, perm->L, perm->m, perm->n);
    
    REprintf ("    s_r : { ");
    for (i = 0; i < M; i++) REprintf ("%d ", perm->s_r[i]);
    REprintf ("}\n");

    REprintf ("    s_c : { ");
    for (j = 0; j < N; j++) REprintf ("%d ", perm->s_c[j]);
    REprintf ("}\n");

    REprintf ("    ir  : { ");
    for (i = 0; i < M; i++) REprintf ("%d ", perm->ir[i]);
    REprintf ("}\n");

    REprintf ("    jc  : { ");
    for (j = 0; j < N; j++) REprintf ("%d ", perm->jc[j]);
    REprintf ("}\n");
    
    REprintf ("}\n");
}

void 
perm_init (perm_t *perm, int M, int N, int K, int L, int *s_r, int *s_c)
{
    perm->M = M;
    perm->N = N;
    perm->K = K;
    perm->L = L;
    
    perm->s_r = s_r;
    perm->s_c = s_c;
    perm->ir  = (int *) R_alloc (M, sizeof (int));
    perm->jc  = (int *) R_alloc (N, sizeof (int));
}

void
perm_select (perm_t *perm, int k, int l)
{
    int *ir, *jc, *s_r, *s_c;
    int i, j, M, N;
    int m = 0, n = 0, mc, nc;
    
    M   = perm->M;
    N   = perm->N;
    s_r = perm->s_r;
    s_c = perm->s_c;
    ir  = perm->ir;
    jc  = perm->jc;
    
    mc = M;
    for (i = 0; i < M; i++)
    {
        if (s_r[i] != k)
            ir[i] = m++;
        else
            ir[i] = --mc;
    }
    for (i = 0; i < M; i++)
        if (ir[i] >= m)
            ir[i] = (M - 1) - (ir[i] - m);
    
    nc = N;
    for (j = 0; j < N; j++)
    {
        if (s_c[j] != l)
            jc[j] = n++;
        else
            jc[j] = --nc;
    }
    for (j = 0; j < N; j++)
        if (jc[j] >= n)
            jc[j] = (N - 1) - (jc[j] - n);
    
    perm->m = m;
    perm->n = n;
}


SEXP 
driver_svd (SEXP xx, SEXP KK, SEXP LL, SEXP max_rank, SEXP s_r, SEXP s_c)
{
    double *x;
    int M, N, K, L, m, n, r, mn, k, i, j, kmax;
    bcv_svd_t *bcv;
    SEXP mse_R, dim;
    double *mse;
    perm_t perm;
    
    if (!isMatrix (xx) || !isNumeric (xx))
        error ("x should be a matrix");

    x    = NUMERIC_POINTER (xx);
    M    = INTEGER (getAttrib (xx, R_DimSymbol))[0];
    N    = INTEGER (getAttrib (xx, R_DimSymbol))[1];
    K    = asInteger (KK);
    L    = asInteger (LL);
    kmax = asInteger (max_rank);

    perm_init (&perm, M, N, K, L, INTEGER_POINTER (s_r), INTEGER_POINTER (s_c));
    
    PROTECT (mse_R = allocVector (REALSXP, (kmax + 1) * K * L));
    mse = NUMERIC_POINTER (mse_R);
    bcv = bcv_svd_alloc (M, N);

    for (j = 1; j <= L; j++)
    {
        for (i = 1; i <= K; i++)
        {
            R_CheckUserInterrupt ();
            
            perm_select (&perm, i, j);
            if (perm.n <= 0 || perm.n >= N)
            {
                REprintf ("L: %d  N: %d  n: %d\n", L, N, perm.n);
            }
            bcv_svd_initp (bcv, M, N, perm.m, perm.n, x, M, perm.ir, perm.jc);

            //bcv_svd_debug (bcv);
            //perm_debug (&perm);

            bcv_svd_decompose (bcv);
            //bcv_svd_debug (bcv);
        
            *mse++ = bcv_svd_get_resid_mse (bcv);
            for (k = 0; k < kmax; k++)
            {
                bcv_svd_update_resid (bcv, -1.0, k);
                *mse++ = bcv_svd_get_resid_mse (bcv);            
            }
        }
    }

    bcv_svd_free (bcv);

    PROTECT (dim = allocVector (INTSXP, 2));
    INTEGER (dim) [0] = (kmax + 1);
    INTEGER (dim) [1] = K * L;
    setAttrib (mse_R, R_DimSymbol, dim);
    
    UNPROTECT (2);
    return mse_R;
}
