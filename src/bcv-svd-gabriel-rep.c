
#include <assert.h>
#include <stdint.h>
#include <strings.h>
#include "bcv-align-private.h"
#include "bcv-svd-gabriel-rep.h"
#include "bcv-matrix-private.h"


struct _bcv_svd_grep
{
    bcv_matrix_t *x11; bcv_matrix_t *x12;
    bcv_matrix_t *x21; bcv_matrix_t *x22;
    double *d;
    void *work;
};


static void
bcv_svd_grep_init_storage (bcv_svd_grep_t *bcv, bcv_gabriel_holdin_t holdin, 
                           bcv_index_t M, bcv_index_t N);

static bcv_error_t
bcv_svd_grep_decompose (bcv_svd_grep_t *bcv);

static bcv_index_t
bcv_svd_grep_decompose_work_len (bcv_gabriel_holdin_t holdin, bcv_index_t M, 
                                 bcv_index_t N);


size_t
bcv_svd_grep_size (bcv_gabriel_holdin_t holdin, bcv_index_t M, bcv_index_t N)
{
    size_t total, decompose_e_tauq_taup, update_u, work_len, result = 0;
    bcv_index_t m, n, mn;
    bcv_index_t decompose_lwork;
    
    _bcv_assert_valid_gabriel_holdin (&holdin, M, N);
    m  = holdin.m;
    n  = holdin.n;
    mn = BCV_MIN (m,n);
    
    /* space for the bcv_svd_grep_t and x11, x12, x21, x22 */
    total = (sizeof (bcv_svd_grep_t) 
             + (_bcv_alignof (bcv_matrix_t) - 1) 
             + 4 * sizeof (bcv_matrix_t)
             + (_bcv_alignof (double) - 1));
        
    /* space for the M*N data matrix */
    if (N == 0 
        || (M <= SIZE_MAX / N
            && M * N <= (SIZE_MAX - total) / sizeof (double)))
    {
        total +=  M * N * sizeof (double);
        
        /* space for d */
        if (mn <= (SIZE_MAX - total) / sizeof (double))
        {
            total += mn * sizeof (double);
            
            /* work space; the memory used by decompose() can be
             * re-used by update()  
             */
            decompose_e_tauq_taup = 3 * mn * sizeof (double);
            decompose_lwork = bcv_svd_grep_decompose_work_len (holdin, M, N);
            update_u  = M * sizeof (double);
            
            if (mn <= SIZE_MAX / 3
                && 3 * mn <= SIZE_MAX / sizeof (double)
                && (decompose_lwork > 0 || mn == 0)
                && decompose_lwork <= (SIZE_MAX - decompose_e_tauq_taup)
                                      / sizeof (double)
                && M <= SIZE_MAX / sizeof (double))
            {
                work_len = BCV_MAX (decompose_e_tauq_taup 
                                    + decompose_lwork * sizeof (double),
                                    update_u);

                if (work_len <= SIZE_MAX - total) 
                {
                    total += work_len;
                    result = total;
                }
            }
        }
    }
    
    return result;
}


size_t
bcv_svd_grep_align ()
{
    return _bcv_alignof (bcv_svd_grep_t);
}


bcv_svd_grep_t *
bcv_svd_grep_alloc (bcv_gabriel_holdin_t holdin, bcv_index_t M, bcv_index_t N)
{
    bcv_svd_grep_t *result = NULL;
    size_t size;
    
    assert (M >= 0);
    assert (N >= 0);
    _bcv_assert_valid_gabriel_holdin (&holdin, M, N);
    
    size = bcv_svd_grep_size (holdin, M, N);
    
    if (size > 0)
    {
        result = malloc (size);
    }

    return result;
}


void 
bcv_svd_grep_free (bcv_svd_grep_t *bcv)
{
    if (bcv)
        free (bcv);
}


void
bcv_svd_grep_init_storage (bcv_svd_grep_t *bcv, bcv_gabriel_holdin_t holdin, 
                           bcv_index_t M, bcv_index_t N)
{
    bcv_index_t m  = holdin.m;
    bcv_index_t n  = holdin.n;
    bcv_index_t mn = BCV_MIN (m,n);
    size_t matrix_align = _bcv_alignof (bcv_matrix_t);
    size_t double_align = _bcv_alignof (double);
    char *mem;
    
    assert (bcv);
    _bcv_assert_valid_gabriel_holdin (&holdin, M, N);
    
    mem            = (char *)bcv; mem += sizeof (bcv_svd_grep_t);
    
    mem = BCV_ALIGN_PTR (mem, matrix_align);
    bcv->x11       = (void *)mem; mem += sizeof (bcv_matrix_t);
    bcv->x21       = (void *)mem; mem += sizeof (bcv_matrix_t);
    bcv->x12       = (void *)mem; mem += sizeof (bcv_matrix_t);
    bcv->x22       = (void *)mem; mem += sizeof (bcv_matrix_t);
    
    mem = BCV_ALIGN_PTR (mem, double_align);
    bcv->x11->data = (void *)mem; mem += M * N * sizeof (double);
    bcv->d         = (void *)mem; mem += mn * sizeof (double);
    bcv->work      = (void *)mem;
    
    bcv->x11->m   = m;
    bcv->x11->n   = n;
    bcv->x11->lda = M;

    bcv->x21->m    = M - m;
    bcv->x21->n    = n;
    bcv->x21->data = bcv->x11->data + m;
    bcv->x21->lda  = M;
    
    bcv->x12->m    = m;
    bcv->x12->n    = N - n;
    bcv->x12->data = bcv->x11->data + n * M;
    bcv->x12->lda  = M;

    bcv->x22->m    = M - m;
    bcv->x22->n    = N - n;
    bcv->x22->data = bcv->x11->data + m + n * M;
    bcv->x22->lda  = M;
}


bcv_error_t
bcv_svd_grep_init (bcv_svd_grep_t *bcv, bcv_gabriel_holdin_t holdin,
                   const bcv_matrix_t *x)
{
    return bcv_svd_grep_init_with_perm (bcv, holdin, x, NULL, NULL);
}


bcv_error_t
bcv_svd_grep_init_with_perm (bcv_svd_grep_t *bcv, bcv_gabriel_holdin_t holdin,
                             const bcv_matrix_t *x, bcv_index_t *p, 
                             bcv_index_t *q)
{
    assert (bcv);
    _bcv_assert_valid_matrix (x);

    bcv_error_t result = 0;
    bcv_index_t M = x->m;
    bcv_index_t N = x->n;
    
    bcv_svd_grep_init_storage (bcv, holdin, M, N);

    bcv_matrix_t dst = { M, N, bcv->x11->data, M };
    _bcv_matrix_permute_copy (&dst, x, p, q);
    
    result = bcv_svd_grep_decompose (bcv);
    
    return result;
}


void 
bcv_svd_grep_get_resid (const bcv_svd_grep_t *bcv, bcv_matrix_t *resid)
{
    assert (bcv);
    assert (resid);
    
    resid->m    = bcv->x22->m;
    resid->n    = bcv->x22->n;
    resid->data = bcv->x22->data;
    resid->lda  = bcv->x22->lda;
}


bcv_index_t 
bcv_svd_grep_get_max_rank (const bcv_svd_grep_t *bcv)
{
    assert (bcv);
    return BCV_MIN (bcv->x11->m, bcv->x11->n);
}


void
bcv_svd_grep_get_holdout_sizes (const bcv_svd_grep_t *bcv, bcv_index_t *m2,
                                bcv_index_t *n2)
{
    assert (bcv);
    _bcv_assert_valid_matrix (bcv->x22);
    assert (m2);
    assert (n2);
    
    *m2 = bcv->x22->m;
    *n2 = bcv->x22->n;
}


double 
bcv_svd_grep_get_press (const bcv_svd_grep_t *bcv)
{
    double frob;
    double press;
    
    assert (bcv);
    
    frob  = _bcv_matrix_norm_frob (bcv->x22);
    press = (frob * frob);
    
    return press;
}


double
bcv_svd_grep_get_msep (const bcv_svd_grep_t *bcv)
{
    bcv_index_t m2, n2, holdout_size;
    double press, msep;
    
    m2 = bcv->x22->m;
    n2 = bcv->x22->n;
    holdout_size = m2 * n2;
    
    press = bcv_svd_grep_get_press (bcv);
    msep  = holdout_size == 0 ? 0.0 : press / holdout_size;
    
    return msep;
}


/*
 * Decompose x11 = Q Q1 D P1^T P^T
 * Set       x11  := P1^T
 *           x12  := Q1^T Q^T x12
 *           x21  := x21 P
 *           work := D
 */
static bcv_error_t
bcv_svd_grep_decompose (bcv_svd_grep_t *bcv)
{
    bcv_error_t result = 0; 
    bcv_index_t m, n, mn, m2, n2;
    bcv_matrix_uplo_t uplo;
    
    assert (bcv);
    assert (bcv->d);
    _bcv_assert_valid_matrix (bcv->x11);
    _bcv_assert_valid_matrix (bcv->x12);
    _bcv_assert_valid_matrix (bcv->x21);
    
    m  = bcv->x11->m;
    n  = bcv->x11->n;
    mn = BCV_MIN (m, n);
    m2 = bcv->x22->m;
    n2 = bcv->x22->n;

    if (mn > 0 && m2 > 0 && n2 > 0)
    {
        bcv_gabriel_holdin_t holdin = { m, n };
        bcv_index_t M = m + m2;
        bcv_index_t N = n + n2;
        bcv_index_t lwork = bcv_svd_grep_decompose_work_len (holdin, M, N);
    
        if (lwork > 0)
        {
            double *e    = bcv->work; 
            double *tauq = e + mn;
            double *taup = tauq + mn;
            double *work = taup + mn;

            /* decompose x11 := Q B P^T */ 
            _bcv_lapack_dgebrd (bcv->x11, bcv->d, e, tauq, taup, work, lwork);

            /* set x21 := x21 P */
            _bcv_lapack_dormbr (BCV_MATRIX_VECT_P, BCV_MATRIX_RIGHT, 
                                BCV_MATRIX_NOTRANS, bcv->x11, taup, bcv->x21, 
                                work, lwork);

            /* set x12 := Q^T x12 */
            _bcv_lapack_dormbr (BCV_MATRIX_VECT_Q, BCV_MATRIX_LEFT, 
                                BCV_MATRIX_TRANS, bcv->x11, tauq, bcv->x12, 
                                work, lwork);

            /* we can now drop the extra rows or columns; they never enter 
             * into the svd.
             */
            bcv->x11->m = mn;
            bcv->x11->n = mn;
            bcv->x12->m = mn;
            bcv->x21->n = mn;
    
            /* decompose B = Q1 S P1^T
             *    set x11 := P1^T 
             *        x12 := Q1^T x12
             */
            _bcv_matrix_set_identity (bcv->x11);
            uplo   = (m >= n) ? BCV_MATRIX_UPPER : BCV_MATRIX_LOWER;
            result = _bcv_lapack_dbdsqr (uplo, mn, bcv->d, e, 
                                         bcv->x11, NULL, bcv->x12, work);
        }
        else
        {
            result = 1; /* TODO: add code for no memory */
        }
    }
    
    return result;
}


static bcv_index_t
bcv_svd_grep_decompose_work_len (bcv_gabriel_holdin_t holdin, bcv_index_t M, 
                               bcv_index_t N)
{
    bcv_index_t result = 0;
    bcv_index_t m, n, mn;
    bcv_index_t dgebrd_lwork, dormbr_P_lwork, dormbr_Q_lwork, dbdsqr_lwork;
    _bcv_assert_valid_gabriel_holdin (&holdin, M, N);
    
    m  = holdin.m;
    n  = holdin.n;
    mn = BCV_MIN (m, n);

    /* We could be more precise below, replacing M with M - m and
     * N with N - n in the calls to _dormbr_work_len.  We prefer to
     * use the conservative values instead so that 
     * bcv_svd_grep_decompose_work_len() is monotonic in the holdin size. */
    dgebrd_lwork   = _bcv_lapack_dgebrd_work_len (m, n);
    dormbr_P_lwork = _bcv_lapack_dormbr_work_len (BCV_MATRIX_VECT_P, 
                                                  BCV_MATRIX_RIGHT, 
                                                  m, n, M, n);
    dormbr_Q_lwork = _bcv_lapack_dormbr_work_len (BCV_MATRIX_VECT_P, 
                                                  BCV_MATRIX_RIGHT, 
                                                  m, n, m, N);
    dbdsqr_lwork   = _bcv_lapack_dbdsqr_work_len (mn, BCV_FALSE);
    
    if (dgebrd_lwork > 0 
        && dormbr_P_lwork > 0 
        && dormbr_Q_lwork > 0 
        && dbdsqr_lwork > 0)
    {
        result = BCV_MAX (BCV_MAX (dgebrd_lwork, dormbr_P_lwork),
                          BCV_MAX (dormbr_Q_lwork, dbdsqr_lwork));
    }
    
    return result;
}


/*
 * This function updates x22 as:
 *         x22 := x22 - (scale/d[i]) * u[i] * v[i]^T
 */
void
bcv_svd_grep_update_resid (bcv_svd_grep_t *bcv, double scale, bcv_index_t i)
{
    assert (bcv);
    assert (0 <= i && i < bcv_svd_grep_get_max_rank (bcv));

    double alpha = -scale / bcv->d[i];
    
    /* x11 stores P1^T
     * Set p1 := ei^T P1^T = P1 ei 
     */
    bcv_vector_t p1 = { bcv->x11->n, bcv->x11->data + i, bcv->x11->lda };
    
    /* x21 stores x21 P
     * Set u := x21 P u = x21 P P1 ei
     *
     * In bcv_svd_grep_size, we have been careful to make sure that the extra work
     * space has room for at least M doubles.
     */
    bcv_index_t m2 = bcv->x21->m;
    bcv_vector_t u = { m2, bcv->work, 1 };
    _bcv_blas_dgemv (BCV_MATRIX_NOTRANS, 1.0, bcv->x21, &p1, 0.0, &u);

    /* x12 stores Q1^T Q^T x12
     * Set v := ei^T Q1^T Q^T x12
     */
    bcv_vector_t v = { bcv->x12->n, bcv->x12->data + i, bcv->x12->lda };
    
    /* Update x22 := x22 + alpha * u v^T */
    _bcv_blas_dger (alpha, &u, &v, bcv->x22);
}
