
#include <assert.h>
#include <math.h>
#include <stdlib.h>
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
    
    result = calloc (1, sizeof (bcv_svd_wrep_t));
    
    return result;
}


size_t
bcv_svd_wrep_size (bcv_index_t M, bcv_index_t N)
{
    size_t result = 0;
    
    return result;
}


size_t
bcv_svd_wrep_align ()
{
    size_t result = __alignof__ (bcv_svd_wrep_t);
    
    return result;
}


void 
bcv_svd_wrep_free (bcv_svd_wrep_t *bcv)
{
    if (bcv)
    {
        if (bcv->impute)
            bcv_svd_impute_free (bcv->impute);
            
        if (bcv->xhat)
        {
            if (bcv->xhat->data)
                free (bcv->xhat->data);
                
            free (bcv->xhat);
        }
        
        free (bcv);
    }
}

void
bcv_svd_wrep_init (bcv_svd_wrep_t *bcv, bcv_wold_holdout_t holdout, 
                   const bcv_matrix_t *x)
{
    bcv_index_t m, n;
    
    assert (bcv);
    _bcv_assert_valid_matrix (x);
    
    m = x->m;
    n = x->n;

    _bcv_assert_valid_wold_holdout (&holdout, m, n);
    
    
    if ((bcv->xhat = calloc (1, sizeof (bcv_matrix_t)))
        && (bcv->xhat->data = malloc (m * n * sizeof (double)))
        && (bcv->impute = bcv_svd_impute_alloc (m, n)))
    {
        bcv->xhat->m   = m;
        bcv->xhat->n   = n;
        bcv->xhat->lda = m;
        
        bcv-> x      = x;
        bcv->holdout = holdout;
        
        bcv_svd_impute_init (bcv->impute, bcv->xhat, x, 
                             holdout.indices, holdout.num_indices);
    }
}


double 
bcv_svd_wrep_get_rss (const bcv_svd_wrep_t *bcv)
{
    bcv_index_t i, j, idx, ldx, ldxhat, m, n, num_indices, *indices;
    double *x, *xhat;
    double x_ij, xhat_ij, d_ij, abs_d_ij, rss = 1.0, scale = 0.0;

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
                rss   = 1.0 + rss * (ratio * ratio);
                scale = abs_d_ij;
            }
            else
            {
                double ratio_inv = abs_d_ij / scale;
                rss   = rss + ratio_inv * ratio_inv;
            }
        }
    }
    
    return rss * (scale * scale);
}


bcv_index_t
bcv_svd_wrep_get_max_rank (bcv_svd_wrep_t *bcv)
{
    bcv_index_t m, n, result;
    
    assert (bcv);
    _bcv_assert_valid_matrix (bcv->x);
    
    m      = bcv->x->m;
    n      = bcv->x->n;
    result = BCV_MIN (m, n);
    
    return result;
}


bcv_error_t 
bcv_svd_wrep_impute_step (bcv_svd_wrep_t *bcv, bcv_index_t k, 
                          double *train_rss)
{
    bcv_error_t err = 0;
    
    assert (bcv);
    assert (train_rss);
    
    err = bcv_svd_impute_step (bcv->impute, bcv->xhat, bcv->x,
                               bcv->holdout.indices, bcv->holdout.num_indices,
                               k);
    *train_rss = bcv_svd_impute_get_rss (bcv->impute);
    
    return err;
}
