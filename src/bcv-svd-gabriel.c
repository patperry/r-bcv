
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include "bcv-partition.h"
#include "bcv-svd-gabriel-rep.h"
#include "bcv-svd-gabriel.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))


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
             + (__alignof__ (bcv_index_t) - 1)
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
    return __alignof__ (bcv_svd_gabriel_t);
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
    void *mem = bcv; mem += sizeof (bcv_svd_gabriel_t);
    size_t bcv_index_align = __alignof__ (bcv_index_t);
    size_t rep_align       = bcv_svd_grep_align ();
    
    assert (bcv);
    assert (x);
        
    bcv->x        = x;
    bcv->row_part = rows;
    bcv->col_part = cols;
    
    mem += bcv_index_align - 1;
    mem = mem - ((uintptr_t) mem & (bcv_index_align - 1));
    
    bcv->row_perm = mem; mem += (x->m) * sizeof (bcv_index_t);
    bcv->col_perm = mem; mem += (x->n) * sizeof (bcv_index_t);
    
    mem += rep_align - 1;
    mem = mem - ((uintptr_t) mem & (rep_align - 1));
    
    bcv->rep = mem;
}


bcv_error_t
bcv_svd_gabriel_get_rss (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                         bcv_index_t j, double *rss, bcv_index_t max_rank)
{
    bcv_error_t error = 0;
    bcv_gabriel_holdin_t holdin;
    bcv_index_t rank = 0;
    
    assert (bcv);
    assert (rss);
    
    holdin.m = bcv_partition_get_perm (bcv->row_part, i, bcv->row_perm);
    holdin.n = bcv_partition_get_perm (bcv->col_part, j, bcv->col_perm);
    
    error  = bcv_svd_grep_init_with_perm (bcv->rep, holdin, bcv->x,
                                          bcv->row_perm, bcv->col_perm);

    *rss++ = bcv_svd_grep_get_rss (bcv->rep);
    
    if (!error)
    {
        bcv_index_t max_allowable_rank = bcv_svd_grep_get_max_rank (bcv->rep);
        assert (0 <= max_rank && max_rank <= max_allowable_rank);
        
        for (rank = 0; rank < max_rank; rank++)
        {
            bcv_svd_grep_update_resid (bcv->rep, 1.0, rank);
            *rss++ = bcv_svd_grep_get_rss (bcv->rep);            
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
    mn = MIN (m,n);
    
    return mn;
}
