
#ifndef _BCV_SVD_WOLD_H
#define _BCV_SVD_WOLD_H

#include <assert.h>
#include <stdlib.h>
#include "bcv-partition.h"
#include "bcv-types.h"

typedef struct _bcv_svd_wold bcv_svd_wold_t;

typedef struct _bcv_wold_holdout {
    bcv_index_t *indices;
    bcv_index_t num_indices;
} bcv_wold_holdout_t;

#define _bcv_assert_valid_wold_holdout(x,M,N) \
    assert (x); \
    assert (0 <= (x)->num_indices && (x)->num_indices <= (M)*(N)); \
    assert ((x)->num_indices == 0 || (x)->indices)


bcv_svd_wold_t *
bcv_svd_wold_alloc (bcv_index_t M, bcv_index_t N);

size_t
bcv_svd_wold_size (bcv_index_t M, bcv_index_t N);

size_t
bcv_svd_wold_align ();

void
bcv_svd_wold_init (bcv_svd_wold_t *bcv, const bcv_matrix_t *x, 
                   const bcv_partition_t *indices);

void
bcv_svd_wold_free (bcv_svd_wold_t *bcv);

bcv_error_t
bcv_svd_wold_get_rss (const bcv_svd_wold_t *bcv, bcv_index_t i,
                      double *rss, bcv_index_t max_rank);

bcv_index_t
bcv_svd_wold_get_max_rank (const bcv_svd_wold_t *bcv, bcv_index_t i);

#endif /* _BCV_SVD_WOLD_H */
