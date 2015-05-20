
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include "bcv-align-private.h"
#include "bcv-partition.h"
#include "bcv-svd-gabriel-rep.h"
#include "bcv-svd-gabriel.h"


struct _bcv_svd_gabriel
{
    const bcv_matrix_t *x;
    const bcv_partition_t *row_part;
    const bcv_partition_t *col_part;
    bcv_index_t *row_perm;
    bcv_index_t *col_perm;
    bcv_svd_grep_t *rep;
};


size_t
bcv_svd_gabriel_size (bcv_gabriel_holdin_t max_holdin, bcv_index_t M, 
                      bcv_index_t N)
{
    size_t result    = 0;
    size_t rep_size  = bcv_svd_grep_size (max_holdin, M, N);
    size_t rep_align = bcv_svd_grep_align ();    
    size_t total;
    
    total = (sizeof (bcv_svd_gabriel_t) 
             + (_bcv_alignof (bcv_index_t) - 1)
             + (rep_align - 1));
    
    if (M <= (SIZE_MAX - total) / sizeof (bcv_index_t))
    {
        total += M * sizeof (bcv_index_t);
        
        if (N <= (SIZE_MAX - total) / sizeof (bcv_index_t))
        {
            total += N * sizeof (bcv_index_t);
            
            if ((rep_size > 0 || M == 0 || N == 0)
                && rep_size <= SIZE_MAX - total)
            {
                total += rep_size;
                result = total;
            }
        }
    }

    return result;
}


size_t
bcv_svd_gabriel_align ()
{
    return _bcv_alignof (bcv_svd_gabriel_t);
}


bcv_svd_gabriel_t *
bcv_svd_gabriel_alloc (bcv_gabriel_holdin_t max_holdin, bcv_index_t M, 
                       bcv_index_t N)
{
    bcv_svd_gabriel_t *bcv = NULL;
    size_t size            = bcv_svd_gabriel_size (max_holdin, M, N);
    
    if (size > 0)
        bcv = malloc (size);

    return bcv;
}


void
bcv_svd_gabriel_free (bcv_svd_gabriel_t *bcv)
{
    if (bcv) 
        free (bcv);
}


void
bcv_svd_gabriel_init (bcv_svd_gabriel_t *bcv, const bcv_matrix_t *x, 
                      const bcv_partition_t *rows, 
                      const bcv_partition_t *cols)
{
    char *mem = (char *)bcv; mem += sizeof (bcv_svd_gabriel_t);
    size_t bcv_index_align = _bcv_alignof (bcv_index_t);
    size_t rep_align       = bcv_svd_grep_align ();
    
    assert (bcv);
    assert (x);
        
    bcv->x        = x;
    bcv->row_part = rows;
    bcv->col_part = cols;
    
    mem = BCV_ALIGN_PTR (mem, bcv_index_align);
    bcv->row_perm = (void *)mem; mem += (x->m) * sizeof (bcv_index_t);
    bcv->col_perm = (void *)mem; mem += (x->n) * sizeof (bcv_index_t);
    
    mem = BCV_ALIGN_PTR (mem, rep_align);
    bcv->rep = (void *)mem;
}


bcv_error_t
bcv_svd_gabriel_get_press (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                           bcv_index_t j, double *press, bcv_index_t max_rank)
{
    bcv_error_t error = 0;
    bcv_gabriel_holdin_t holdin;
    bcv_index_t rank = 0;
    
    assert (bcv);
    assert (press);
    
    holdin.m = bcv_partition_get_perm (bcv->row_part, i, bcv->row_perm);
    holdin.n = bcv_partition_get_perm (bcv->col_part, j, bcv->col_perm);
    
    error  = bcv_svd_grep_init_with_perm (bcv->rep, holdin, bcv->x,
                                          bcv->row_perm, bcv->col_perm);

    *press++ = bcv_svd_grep_get_press (bcv->rep);
    
    if (!error)
    {
        assert (0 <= max_rank 
                && max_rank <= bcv_svd_grep_get_max_rank (bcv->rep));
        
        for (rank = 0; rank < max_rank; rank++)
        {
            bcv_svd_grep_update_resid (bcv->rep, 1.0, rank);
            *press++ = bcv_svd_grep_get_press (bcv->rep);            
        }
    }
    
    return error;
}


bcv_error_t
bcv_svd_gabriel_get_msep (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                          bcv_index_t j, double *msep, bcv_index_t max_rank)
{
    bcv_error_t error;
    bcv_index_t m2, n2, size2, rank;
    
    error = bcv_svd_gabriel_get_press (bcv, i, j, msep, max_rank);
    bcv_svd_grep_get_holdout_sizes (bcv->rep, &m2, &n2);
    size2 = m2 * n2;
    
    if (size2 > 0)
    {
        for (rank = 0; rank <= max_rank; rank++, msep++)
        {
            *msep = *msep / size2;
        }
    }
    
    return error;
}


bcv_index_t
bcv_svd_gabriel_get_max_rank (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                              bcv_index_t j)
{
    assert (bcv);
    bcv_index_t m, n, mn;
    
    m  = bcv_partition_get_perm (bcv->row_part, i, bcv->row_perm);
    n  = bcv_partition_get_perm (bcv->col_part, j, bcv->col_perm);
    mn = BCV_MIN (m,n);
    
    return mn;
}
