
#include <math.h>
#include "bcv-svd-wold.h"
#include "bcv-svd-wold-rep.h"


struct _bcv_svd_wold 
{
    const bcv_matrix_t *x;
    const bcv_partition_t *part;
    bcv_svd_wrep_t *rep;
    bcv_index_t *holdout;
};


bcv_svd_wold_t *
bcv_svd_wold_alloc (bcv_index_t max_holdout, bcv_index_t M, bcv_index_t N)
{
    bcv_svd_wold_t *result = NULL;
    size_t size = bcv_svd_wold_size (M, N);
    
    if (size > 0)
    {
        if (!((result = calloc (1, size))
              && (result->holdout = 
                      malloc (max_holdout * sizeof (bcv_index_t)))))
        {
            bcv_svd_wold_free (result);
            result = NULL;
        }
    }
        
    return result;
}


size_t
bcv_svd_wold_size (bcv_index_t M, bcv_index_t N)
{
    size_t result = sizeof (bcv_svd_wold_t);
    
    return result;
}

size_t
bcv_svd_wold_align ()
{
    size_t result = __alignof__ (bcv_svd_wold_t);
    
    return result;
}


void
bcv_svd_wold_init (bcv_svd_wold_t *bcv, const bcv_matrix_t *x, 
                   const bcv_partition_t *part)
{
    bcv_index_t m, n;
    
    assert (bcv);
    _bcv_assert_valid_matrix (x);
    assert (part);
    
    m = x->m;
    n = x->n;
    
    bcv->rep  = bcv_svd_wrep_alloc (m, n);
    bcv->x    = x;
    bcv->part = part;
}


void
bcv_svd_wold_free (bcv_svd_wold_t *bcv)
{
    if (bcv) 
    {
        if (bcv->rep)
            bcv_svd_wrep_free (bcv->rep);
        
        free (bcv);
    }
}


bcv_error_t
bcv_svd_wold_get_rss (const bcv_svd_wold_t *bcv, bcv_index_t i,
                      double *rss, bcv_index_t max_rank,
                      double tol, bcv_index_t max_iter)
{
    bcv_error_t err = 0;
    bcv_index_t rank;
    bcv_wold_holdout_t holdout;
    
    assert (bcv);
    assert (0 <= i && i < bcv->part->k);
    
    holdout.indices     = bcv->holdout;
    holdout.num_indices = bcv_partition_get_set (bcv->part, i, 
                                                 holdout.indices);
    
    bcv_svd_wrep_init (bcv->rep, holdout, bcv->x);

    for (rank = 0; rank < max_rank; rank++)
    {
        bcv_index_t iter = 0;
        double rss0, rss1 = BCV_DBL_POSINF, delta = BCV_DBL_POSINF;

        for (iter = 0; !err && iter < max_iter && !(delta <= tol); iter++)
        {
            rss0 = rss1;
            iter++;
            
            err   = bcv_svd_wrep_impute_step (bcv->rep, rank, &rss1);
            delta = fabs (rss1 - rss0) / (BCV_DBL_EPSILON + rss1);
        }
        
        if (err)
            break;
            
        rss[rank] = bcv_svd_wrep_get_rss (bcv->rep);
    }
    
    return err;
}


bcv_index_t
bcv_svd_wold_get_max_rank (const bcv_svd_wold_t *bcv, bcv_index_t i)
{
    bcv_index_t m, n, mn;
    
    assert (bcv);
    _bcv_assert_valid_matrix (bcv->x);
    
    m  = bcv->x->m;
    n  = bcv->x->n;
    mn = BCV_MIN (m, n);
    
    return mn;
}
