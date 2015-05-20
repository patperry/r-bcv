
#include <math.h>
#include "bcv-align-private.h"
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
    size_t size = bcv_svd_wold_size (max_holdout, M, N);
    
    if (size > 0)
        result = malloc (size);
        
    return result;
}


size_t
bcv_svd_wold_size (bcv_index_t max_holdout, bcv_index_t M, bcv_index_t N)
{
    size_t result = 0, total;
    size_t wrep_size;
    
    total = (sizeof (bcv_svd_wold_t)
             + (bcv_svd_wrep_align() - 1)
             + (_bcv_alignof (bcv_index_t) - 1));
    
    /* space for the indices in the holdout set */
    if (max_holdout <= (SIZE_MAX - total) / sizeof (bcv_index_t))
    {
        total += max_holdout * sizeof (bcv_index_t);
        
        /* space for the replicate */
        wrep_size = bcv_svd_wrep_size (M, N);
        if (wrep_size > 0
            && wrep_size <= SIZE_MAX - total)
        {
            total += wrep_size;
            result = total;
        }
        
    }
    
    return result;
}

size_t
bcv_svd_wold_align ()
{
    size_t result = _bcv_alignof (bcv_svd_wold_t);
    
    return result;
}


void
bcv_svd_wold_init (bcv_svd_wold_t *bcv, const bcv_matrix_t *x, 
                   const bcv_partition_t *part)
{
    bcv_index_t m, n;
    size_t wrep_align  = bcv_svd_wrep_align ();
    size_t index_align = _bcv_alignof (bcv_index_t);
    char *mem;
    
    assert (bcv);
    _bcv_assert_valid_matrix (x);
    assert (part);
    
    m = x->m;
    n = x->n;
    
    assert (part->n == m*n);
    
    mem = (char *)bcv; mem += sizeof (bcv_svd_wold_t);
    
    mem = BCV_ALIGN_PTR (mem, wrep_align);
    bcv->rep = (void *)mem; mem += bcv_svd_wrep_size (m, n);
    
    mem = BCV_ALIGN_PTR (mem, index_align);
    bcv->holdout = (void *)mem;
    
    bcv->x    = x;
    bcv->part = part;
}


void
bcv_svd_wold_free (bcv_svd_wold_t *bcv)
{
    if (bcv) 
        free (bcv);
}


bcv_error_t
bcv_svd_wold_get_press (const bcv_svd_wold_t *bcv, bcv_index_t i,
                        double tol, bcv_index_t max_iter,
                        double *press, bcv_index_t max_rank)
{
    bcv_error_t err = 0;
    bcv_index_t rank;
    bcv_wold_holdout_t holdout;
    
    assert (bcv);
    assert (0 <= i && i < bcv->part->k);
    
    holdout.indices     = bcv->holdout;
    holdout.num_indices = bcv_partition_get_set (bcv->part, i, 
                                                 holdout.indices);

    for (rank = 0; rank <= max_rank; rank++)
    {
        bcv_index_t iter = 0;
        double rss0, rss1 = BCV_DBL_POSINF, delta = BCV_DBL_POSINF;

        bcv_svd_wrep_init (bcv->rep, holdout, bcv->x);

        for (iter = 0; !err && iter < max_iter && !(delta <= tol); iter++)
        {
            rss0  = rss1;
            err   = bcv_svd_wrep_impute_step (bcv->rep, rank, &rss1);
            delta = fabs (rss1 - rss0) / (BCV_DBL_EPSILON + rss1);
        }
        
        if (err)
            break;
            
        press[rank] = bcv_svd_wrep_get_press (bcv->rep);
    }
    
    return err;
}


bcv_error_t
bcv_svd_wold_get_msep (const bcv_svd_wold_t *bcv, bcv_index_t i,
                       double tol, bcv_index_t max_iter,
                       double *msep, bcv_index_t max_rank)
{
    bcv_error_t error;
    bcv_index_t rank, holdout_size;
    
    error = bcv_svd_wold_get_press (bcv, i, tol, max_iter, msep, max_rank);
    holdout_size = bcv_svd_wrep_get_holdout_size (bcv->rep);
    
    if (holdout_size > 0)
    {
        for (rank = 0; rank <= max_rank; rank++, msep++)
        {
            *msep = *msep / holdout_size;
        }
    }
    
    return error;    
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
