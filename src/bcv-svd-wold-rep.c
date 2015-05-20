
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include "bcv-align-private.h"
#include "bcv-svd-impute.h"
#include "bcv-svd-wold-rep.h"

struct _bcv_svd_wrep
{
    const bcv_matrix_t *x;
    bcv_wold_holdout_t holdout;

    bcv_svd_impute_t *impute;
    bcv_matrix_t *xhat;
};


bcv_svd_wrep_t *
bcv_svd_wrep_alloc (bcv_index_t M, bcv_index_t N)
{
    bcv_svd_wrep_t *result = NULL;
    size_t size = bcv_svd_wrep_size (M, N);

    if (size)
        result = malloc (size);
    
    return result;
}


size_t
bcv_svd_wrep_size (bcv_index_t M, bcv_index_t N)
{
    size_t result = 0, total, impute_size;
    
    total = (sizeof (bcv_svd_wrep_t)
             + (_bcv_alignof (bcv_matrix_t) - 1) 
             + sizeof (bcv_matrix_t)
             + (_bcv_alignof (double) - 1)
             + (bcv_svd_impute_align () - 1));
             
    /* space for xhat */
    if (N == 0
        || (M <= SIZE_MAX / N
            && M * N <= (SIZE_MAX - total) / sizeof (double)))
    {
        total += M * N * sizeof (double);
        
        /* space for the impute workspace */
        impute_size = bcv_svd_impute_size (M, N);
        if (impute_size > 0 &&
            impute_size <= SIZE_MAX - total)
        {
            total += impute_size;
            result = total;
        }
    }
    
    return result;
}


size_t
bcv_svd_wrep_align ()
{
    size_t result = _bcv_alignof (bcv_svd_wrep_t);
    
    return result;
}


void 
bcv_svd_wrep_free (bcv_svd_wrep_t *bcv)
{
    if (bcv)
        free (bcv);
}

void
bcv_svd_wrep_init (bcv_svd_wrep_t *bcv, bcv_wold_holdout_t holdout, 
                   const bcv_matrix_t *x)
{
    bcv_index_t m, n;
    char *mem;
    size_t matrix_align = _bcv_alignof (bcv_matrix_t);
    size_t impute_align = bcv_svd_impute_align ();
    size_t double_align = _bcv_alignof (double);
    
    assert (bcv);
    _bcv_assert_valid_matrix (x);
    
    m = x->m;
    n = x->n;
    _bcv_assert_valid_wold_holdout (&holdout, m, n);
    
    mem = (char *)bcv; mem += sizeof (bcv_svd_wrep_t);
    
    mem = BCV_ALIGN_PTR (mem, matrix_align);
    bcv->xhat = (void *)mem; mem += sizeof (bcv_matrix_t);

    mem = BCV_ALIGN_PTR (mem, double_align);
    bcv->xhat->data = (void *)mem; mem += m * n * sizeof (double);
    
    mem = BCV_ALIGN_PTR (mem, impute_align);
    bcv->impute = (void *)mem;
    
    bcv->xhat->m   = m;
    bcv->xhat->n   = n;
    bcv->xhat->lda = m;
    
    bcv->x       = x;
    bcv->holdout = holdout;
    
    bcv_svd_impute_init (bcv->impute, bcv->xhat, x, 
                         holdout.indices, holdout.num_indices);
}


double 
bcv_svd_wrep_get_press (const bcv_svd_wrep_t *bcv)
{
    bcv_index_t i, j, idx, ldx, ldxhat, m, n, num_indices, *indices;
    double *x, *xhat;
    double x_ij, xhat_ij, d_ij, abs_d_ij, press = 1.0, scale = 0.0;

    assert (bcv);
    _bcv_assert_valid_matrix (bcv->x);
    _bcv_assert_valid_matrix (bcv->xhat);
    
    indices     = bcv->holdout.indices;
    num_indices = bcv->holdout.num_indices;

    assert (indices || num_indices == 0);
    
    m      = bcv->x->m;
    n      = bcv->x->n;
    x      = bcv->x->data;
    ldx    = bcv->x->lda;
    xhat   = bcv->xhat->data;
    ldxhat = bcv->xhat->lda;
    
    assert (bcv->xhat->m == m);
    assert (bcv->xhat->n == n);
    
    if (m == 0 || n == 0)
        return 0.0;
    
    while (num_indices > 0)
    {
        idx = *indices++;
        num_indices--;
        
        if (ldx == m && ldxhat == m)
        {
            x_ij    = x[idx];
            xhat_ij = xhat[idx];
        }
        else
        {
            i       = idx % m;
            j       = idx / m;
            x_ij    = x[i + j * ldx];
            xhat_ij = xhat[i + j * ldxhat];
        }
        
        d_ij = x_ij - xhat_ij;
        if (d_ij != 0)
        {
            abs_d_ij = fabs (d_ij);
            if (scale < abs_d_ij)
            {
                double ratio = scale / abs_d_ij;
                press = 1.0 + press * (ratio * ratio);
                scale = abs_d_ij;
            }
            else
            {
                double ratio_inv = abs_d_ij / scale;
                press = press + ratio_inv * ratio_inv;
            }
        }
    }
    
    return press * (scale * scale);
}


double 
bcv_svd_wrep_get_msep (const bcv_svd_wrep_t *bcv)
{
    double press, msep;
    bcv_index_t holdout_size;
    
    assert (bcv);
    
    press        = bcv_svd_wrep_get_press (bcv);
    holdout_size = bcv->holdout.num_indices;
    msep         = holdout_size == 0 ? 0.0 : press / holdout_size;
    
    return msep;
}


bcv_index_t
bcv_svd_wrep_get_max_rank (const bcv_svd_wrep_t *bcv)
{
    bcv_index_t m, n, result;
    
    assert (bcv);
    _bcv_assert_valid_matrix (bcv->x);
    
    m      = bcv->x->m;
    n      = bcv->x->n;
    result = BCV_MIN (m, n);
    
    return result;
}


bcv_index_t
bcv_svd_wrep_get_holdout_size (const bcv_svd_wrep_t *bcv)
{
    bcv_index_t size;
    
    assert (bcv);
    
    size = bcv->holdout.num_indices;
    
    return size;
}


bcv_error_t 
bcv_svd_wrep_impute_step (bcv_svd_wrep_t *bcv, bcv_index_t k, 
                          double *rss)
{
    bcv_error_t err = 0;
    
    assert (bcv);
    assert (rss);
    
    err = bcv_svd_impute_step (bcv->impute, bcv->xhat, bcv->x,
                               bcv->holdout.indices, bcv->holdout.num_indices,
                               k);
    *rss = bcv_svd_impute_get_rss (bcv->impute);
    
    return err;
}
