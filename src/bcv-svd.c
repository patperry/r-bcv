
#include <assert.h>
#include <strings.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include "bcv-svd.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define BLOCKSIZE 64


struct _bcv_svd
{
    int M_max, N_max;
    int M, N, m, n, mn;
    double *x11, *x12, *x21, *x22, *d;
    int ld11, ld12, ld21, ld22;
};


static int decompose (int M, int N, int m, int n, double *x11, int ld11, 
                      double *x12, int ld12, double *x21, int ld21, double *d);

static void update (int M, int N, int m, int n, int k, double alpha,
                    const double *x21, int ld21, const double *vt, int ldvt,
                    const double *x12, int ld12, double *x22, int ld22);
                    
static void set_identity (int m, int n, double *a, int lda);

static void print_matrix (const char *name, int m, int n, const double *x, 
                          int ldx);



bcv_svd_t *
bcv_svd_alloc (int M_max, int N_max)
{
    bcv_svd_t *bcv;
    
    assert( M_max >= 1 );
    assert( N_max >= 1 );
    
    if (   (bcv      = malloc (sizeof (bcv_svd_t)))
        && (bcv->x11 = malloc (M_max * N_max * sizeof (double)))
        && (bcv->d   = malloc (MIN (M_max, N_max) *  sizeof (double)))
       )
    {
        bcv->M_max = M_max;
        bcv->N_max = N_max;
        bcv->ld11  = M_max;
        bcv->ld12  = M_max;
        bcv->ld21  = M_max;
        bcv->ld22  = M_max;
        bcv->m     = 0;
        bcv->n     = 0;
        bcv->mn    = 0;

        return bcv;
    } 
        
    error ("Could not allocate bcv_svd_t for size (%d,%d)", M_max, N_max);
}


void 
bcv_svd_init (bcv_svd_t *bcv, int M, int N, int m, int n, double *x, int ldx)
{
    bcv_svd_initp (bcv, M, N, m, n, x, ldx, NULL, NULL);
}

void 
bcv_svd_initp (bcv_svd_t *bcv, int M, int N, int m, int n, double *x, int ldx, 
               int *p, int *q)
{
    int one = 1;
    int i, j;
    
    assert (bcv);
    assert (x);
    assert (0 < M && M <= bcv->M_max);
    assert (0 < N && N <= bcv->N_max);
    assert (0 < m && m < M);
    assert (0 < n && n < N);
    assert (ldx >= M);
    
    bcv->M  = M;
    bcv->N  = N;
    bcv->m  = m;
    bcv->n  = n;
    bcv->mn = MIN (m, n);

    bcv->x21 = bcv->x11 + m;
    bcv->x12 = bcv->x11 + n * ldx;
    bcv->x22 = bcv->x11 + m + n * ldx;

    if (p == NULL && q == NULL)
    {
        F77_CALL (dlacpy) ("R", &M, &N, x, &ldx, bcv->x11, &(bcv->ld11));
    }
    else if (p == NULL) /* permute the columns only */
    {
        for (j = 0; j < N; j++)
        {
            F77_CALL (dcopy) (&M, x + j * ldx, &one, 
                              bcv->x11 + q[j] * (bcv->ld11), &one);
        }
    }
    else if (q == NULL) /* permute the rows only */
    {
        for (i = 0; i < M; i++)
        {
            F77_CALL (dcopy) (&N, x + i, &ldx, bcv->x11 + p[i], &(bcv->ld11));
        }
    }
    else /* permute both rows and columns */
    {
        double *src, *dst;

        for (j = 0; j < N; j++)
        {
            src = x + j * ldx;
            dst = bcv->x11 + q[j] * (bcv->ld11);
            
            for (i = 0; i < M; i++)
                dst[p[i]] = src[i];
        }
    }
}

void 
bcv_svd_free (bcv_svd_t *bcv)
{
    if (bcv)
    {
        free (bcv->d);
        free (bcv->x11);
        free (bcv);
    }
}


int 
bcv_svd_decompose (bcv_svd_t *bcv)
{
    int info;
    
    assert (bcv);
    
    info = decompose (bcv->M, bcv->N, bcv->m, bcv->n, bcv->x11, bcv->ld11, 
                      bcv->x12, bcv->ld12, bcv->x21, bcv->ld21, bcv->d);
                          
    return info;
}


void 
bcv_svd_get_resid (const bcv_svd_t *bcv, int *m2, int *n2, double **resid, 
                   int *ldr)
{
    assert (bcv);
    assert (m2);
    assert (n2);
    assert (resid);
    assert (ldr);
    
    *m2    = bcv->M - bcv->m;
    *n2    = bcv->N - bcv->n;
    *resid = bcv->x22;
    *ldr   = bcv->ld22;
}


int 
bcv_svd_get_max_rank (bcv_svd_t *bcv)
{
    assert (bcv);
    return bcv->mn;
}


void
bcv_svd_update_resid (bcv_svd_t *bcv, double scale, int k)
{
    double alpha;

    assert (bcv);
    assert (k < bcv->mn);
    
    alpha = scale / bcv->d[k];
    update (bcv->M, bcv->N, bcv->m, bcv->n, k, alpha, bcv->x21, bcv->ld21,
            bcv->x11, bcv->ld11, bcv->x12, bcv->ld12, bcv->x22, bcv->ld22);
}


double 
bcv_svd_get_resid_mse (const bcv_svd_t *bcv)
{
    int m2 = bcv->M - bcv->m;
    int n2 = bcv->N - bcv->n;
    double frob;
    double mse;
    
    frob = F77_CALL (dlange) ("F", &m2, &n2, bcv->x22, &(bcv->ld22), NULL);
    mse  = (frob * frob) / (m2 * n2);
    
    return mse;
}


void
bcv_svd_debug (const bcv_svd_t *bcv)
{
    assert (bcv);
    
    REprintf ("bcv_svd_t {\n"
              "    M_max : %d\n"
              "    N_max : %d\n"
              "    M     : %d\n"
              "    N     : %d\n"
              "    m     : %d\n"
              "    n     : %d\n"
              "    mn    : %d\n",
              bcv->M_max, bcv->N_max, bcv->M, bcv->N, bcv->m, bcv->n, bcv->mn);
    
    print_matrix ("x11", bcv->m, bcv->n, bcv->x11, bcv->ld11);
    print_matrix ("x12", bcv->m, bcv->N - bcv->n, bcv->x12, bcv->ld12);
    print_matrix ("x21", bcv->M - bcv->m, bcv->n, bcv->x21, bcv->ld21);
    print_matrix ("x22", bcv->M - bcv->m, bcv->N - bcv->n, bcv->x22, bcv->ld22);
    print_matrix ("  d", 1, bcv->mn, bcv->d, 1);

    REprintf ("}\n");
}


/*
 * Decompose x11 = Q U D V^T P^T
 * Set       x11 := V^T
 *           x12 := U^T Q^T x12
 *           x21 := x21 P
 *           d   := D
 */
static int
decompose (int M, int N, int m, int n, double *x11, int ld11, double *x12, 
           int ld12, double *x21, int ld21, double *d)
{
    int m2 = M - m;
    int n2 = N - n;
    int mn = MIN (m, n);
    double e[mn], tauq[mn], taup[mn];
    int lwork = (m + n) * BLOCKSIZE;
    double work[lwork];
    double _d;
    int info = 0;
    int zero = 0;
    int one = 1;
    double dzero = 0.0;
    int i;
    char *uplo;

    assert (m > 0);
    assert (n > 0);
    
    if (m == 1) /* x11 = d v^T */
    {
        d[0] = F77_CALL (dnrm2) (&n, x11, &ld11);

        _d = 1.0 / d[0];
        F77_CALL (dscal) (&n, &_d, x11, &ld11);
    }
    else if (n == 1) /* x11 = u d */
    {
        d[0] = F77_CALL (dnrm2) (&m, x11, &one);
        
        /* x12 := u^T x12 */
        _d = 1.0 / d[0];
        F77_CALL (dgemv) ("T", &m, &n2, &_d, x12, &ld12, x11, &one, &dzero,
                          work, &one);
        F77_CALL (dcopy) (&n2, work, &one, x12, &ld12);

        x11[0] = 1.0;
    }
    else 
    {
        /* decompose x11 := Q B P^T */ 
        F77_CALL (dgebrd) (&m, &n, x11, &ld11, d, e, tauq, taup, work, &lwork, 
                           &info);
    
        if (info == 0)
        {
            /* set x21 := x21 P */
            F77_CALL (dormbr) ("P", "R", "N", &m2, &n, &m, x11, 
                               &ld11, taup, x21, &ld21, work, &lwork, &info);
            assert (info == 0);
        
            /* set x12 := Q^T x12 */
            F77_CALL (dormbr) ("Q", "L", "T", &m, &n2, &n, x11, 
                               &ld11, tauq, x12, &ld12, work, &lwork, &info);
            assert (info == 0);

            /* decompose B = U S V^T
             *    set x11 := V^T 
             *        x12 := U^T x12
             */
            uplo = (m >= n) ? "U" : "L";
            set_identity (mn, n, x11, ld11);
            F77_CALL (dbdsqr) (uplo, &mn, &n, &zero, &n2, d, e, x11, &ld11, 
                               NULL, &one, x12, &ld12, work, &info);
        }
    }
    
    return info;
}

/*
 * We have x22_hat = x21 v d x12
 *                 = \sum_{k=1}^{K} d[k] * (x21 v(:,k)) (x12(k,:))^T
 *
 * This function updates x22 as:
 *         x22 := x22 + alpha * (x21 v(:,k)) (x12(k,:))^T
 */
static void
update (int M, int N, int m, int n, int k, double alpha,
        const double *x21, int ld21, const double *vt, int ldvt,
        const double *x12, int ld12, double *x22, int ld22)
{
    int m2 = M - m;
    int n2 = N - n;
    int mn = MIN (m, n);
    double u[m2];
    double done = 1.0;
    double dzero = 0.0;
    int one = 1;
    
    assert (k < mn);

    /* u := x21 * v(:,k)  ( = x21 : vt(k,:) ) */
    F77_CALL (dgemv) ("N", &m2, &n, &done, x21, &ld21, vt + k, &ldvt,
                      &dzero, u, &one);

    /* x22 := x22 + alpha * u x12(k,:)^T */
    F77_CALL (dger) (&m2, &n2, &alpha, u, &one, (double *) (x12 + k), &ld12, 
                     x22, &ld22);
}


static void
set_identity (int m, int n, double *a, int lda)
{
    int j, mn;
    
    if (m > 0 && n > 0)
    {
        if (lda == m)
        {
            bzero (a, m * n * sizeof (double));

            mn = MIN (m, n);
            for (j = 0; j < mn; j++, a += lda + 1)
                *a = 1.0;
        }
        else
        {
            for (j = 0; j < n; j++, a += lda)
            {
                bzero (a, m * sizeof (double));
                if (j < m) a[j] = 1.0;
            }
        }
    }
}

static void
print_matrix (const char *name, int m, int n, const double *x, int ldx)
{
    int i,j;
 
    REprintf ("    %s (0x%x) :\n", name, x );
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (x[i + j*ldx] >= 0) REprintf (" ");
            REprintf ("%4g    ", x[i + j*ldx]);
        }
        REprintf("\n          ");
    }
    REprintf("\n");
}
