
#include <assert.h>
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

bcv_svd_gabriel_t *
bcv_svd_gabriel_alloc (bcv_gabriel_holdin_t max_holdin, bcv_index_t M, 
                       bcv_index_t N)
{
    bcv_svd_gabriel_t *bcv = calloc (1, sizeof (bcv_svd_gabriel_t));
    
    if (bcv)
    {
        if (!(((bcv->row_perm = malloc (M * sizeof (bcv_index_t))) 
               || M == 0)
              && ((bcv->col_perm = malloc (N * sizeof (bcv_index_t)))
                  || N == 0)
              && ((bcv->rep      = bcv_svd_grep_alloc (max_holdin, M, N)))))
        {
            bcv_svd_gabriel_free (bcv);
            bcv = NULL;
        }
    }
        
    return bcv;
}

void
bcv_svd_gabriel_free (bcv_svd_gabriel_t *bcv)
{
    if (bcv)
    {
        if (bcv->row_perm) free (bcv->row_perm);
        if (bcv->col_perm) free (bcv->col_perm);
        if (bcv->rep)      bcv_svd_grep_free (bcv->rep);
        
        free (bcv);
    }
}


void
bcv_svd_gabriel_init (bcv_svd_gabriel_t *bcv, const bcv_matrix_t *x, 
                      const bcv_partition_t *rows, 
                      const bcv_partition_t *cols)
{
    assert (bcv);
    bcv->x = x;
    bcv->row_part = rows;
    bcv->col_part = cols;
}


bcv_error_t
bcv_svd_gabriel_get_rss (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                         bcv_index_t j, double *rss, bcv_index_t max_rank)
{
    bcv_error_t error = 0;
    bcv_gabriel_holdin_t holdin;
    bcv_index_t rank = 0;
    
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
