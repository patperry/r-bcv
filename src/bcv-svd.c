
#include <assert.h>
#include <strings.h>
#include <R.h>
#include <R_ext/Lapack.h>
#include "bcv-svd.h"
#include "bcv-types.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define BLOCKSIZE 64


struct _bcv_svd
{
    int M_max, N_max;
    bcv_matrix_t *x11; bcv_matrix_t *x12;
    bcv_matrix_t *x21; bcv_matrix_t *x22;
    double *d;
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
    
    if (   (bcv            = malloc (sizeof (bcv_svd_t)))
        && (bcv->x11       = malloc (sizeof (bcv_matrix_t)))
        && (bcv->x21       = malloc (sizeof (bcv_matrix_t)))
        && (bcv->x12       = malloc (sizeof (bcv_matrix_t)))
        && (bcv->x22       = malloc (sizeof (bcv_matrix_t)))
        && (bcv->x11->data = malloc (M_max * N_max * sizeof (double)))
        && (bcv->d         = malloc (MIN (M_max, N_max) *  sizeof (double)))
       )
    {
        bcv->M_max = M_max;
        bcv->N_max = N_max;
        
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
    
    bcv->x11->m   = m;
    bcv->x11->n   = n;
    bcv->x11->lda = M;

    bcv->x21->m   = M - m;
    bcv->x21->n   = n;
    bcv->x21->data = bcv->x11->data + m;
    bcv->x21->lda = M;
    
    bcv->x12->m    = m;
    bcv->x12->n    = N - n;
    bcv->x12->data = bcv->x11->data + n * M;
    bcv->x12->lda  = M;

    bcv->x22->m    = M - m;
    bcv->x22->n    = N - n;
    bcv->x22->data = bcv->x11->data + m + n * M;
    bcv->x22->lda  = M;

    if (p == NULL && q == NULL)
    {
        F77_CALL (dlacpy) ("R", &M, &N, x, &ldx, bcv->x11->data, &(bcv->x11->lda));
    }
    else if (p == NULL) /* permute the columns only */
    {
        for (j = 0; j < N; j++)
        {
            F77_CALL (dcopy) (&M, x + j * ldx, &one, 
                              bcv->x11->data + q[j] * (bcv->x11->lda), &one);
        }
    }
    else if (q == NULL) /* permute the rows only */
    {
        for (i = 0; i < M; i++)
        {
            F77_CALL (dcopy) (&N, x + i, &ldx, bcv->x11->data + p[i], &(bcv->x11->lda));
        }
    }
    else /* permute both rows and columns */
    {
        double *src, *dst;

        for (j = 0; j < N; j++)
        {
            src = x + j * ldx;
            dst = bcv->x11->data + q[j] * (bcv->x11->lda);
            
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
    bcv_index_t M, N, m, n;
    
    assert (bcv);
    
    m = bcv->x11->m;
    n = bcv->x11->n;
    M = m + bcv->x22->m;
    N = n + bcv->x22->n;
    
    info = decompose (M, N, m, n, 
                      bcv->x11->data, bcv->x11->lda, 
                      bcv->x12->data, bcv->x12->lda, 
                      bcv->x21->data, bcv->x21->lda, 
                      bcv->d);
                          
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
    
    *m2    = bcv->x22->m;
    *n2    = bcv->x22->n;
    *resid = bcv->x22->data;
    *ldr   = bcv->x22->lda;
}


int 
bcv_svd_get_max_rank (bcv_svd_t *bcv)
{
    assert (bcv);
    return MIN (bcv->x11->m, bcv->x11->n);
}


void
bcv_svd_update_resid (bcv_svd_t *bcv, double scale, int k)
{
    double alpha;
    bcv_index_t M, N, m, n, mn;
    
    assert (bcv);
    
    m  = bcv->x11->m;
    n  = bcv->x11->n;
    M  = m + bcv->x22->m;
    N  = n + bcv->x22->n;
    mn = MIN (m,n);

    assert (k < mn);

    alpha = scale / bcv->d[k];
    update (M, N, m, n, k, alpha, 
            bcv->x21->data, bcv->x21->lda,
            bcv->x11->data, bcv->x11->lda, 
            bcv->x12->data, bcv->x12->lda, 
            bcv->x22->data, bcv->x22->lda);
}


double 
bcv_svd_get_resid_mse (const bcv_svd_t *bcv)
{
    int m2 = bcv->x22->m;
    int n2 = bcv->x22->n;
    double frob;
    double mse;
    
    frob = F77_CALL (dlange) ("F", &m2, &n2, 
                              bcv->x22->data, &(bcv->x22->lda), NULL);
    mse  = (frob * frob) / (m2 * n2);
    
    return mse;
}


void
bcv_svd_debug (const bcv_svd_t *bcv)
{
    bcv_index_t M, N, m, n, mn;
    
    assert (bcv);
    
    m = bcv->x11->m;
    n = bcv->x11->n;
    M = m + bcv->x22->m;
    N = n + bcv->x22->n;
    mn = MIN (m,n);

    
    REprintf ("bcv_svd_t {\n"
              "    M_max : %d\n"
              "    N_max : %d\n"
              "    M     : %d\n"
              "    N     : %d\n"
              "    m     : %d\n"
              "    n     : %d\n"
              "    mn    : %d\n",
              bcv->M_max, bcv->N_max, M, N, m, n, mn);
    
    print_matrix ("x11", bcv->x11->m, bcv->x11->n, bcv->x11->data, bcv->x11->lda);
    print_matrix ("x12", bcv->x12->m, bcv->x12->n, bcv->x12->data, bcv->x12->lda);
    print_matrix ("x21", bcv->x21->m, bcv->x21->n, bcv->x21->data, bcv->x21->lda);
    print_matrix ("x22", bcv->x22->m, bcv->x22->n, bcv->x22->data, bcv->x22->lda);
    print_matrix ("  d", 1, mn, bcv->d, 1);

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
