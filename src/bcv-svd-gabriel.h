
#ifndef _BCV_SVD_GABRIEL_H
#define _BCV_SVD_GABRIEL_H

#include "bcv-partition.h"
#include "bcv-types.h"

typedef struct _bcv_svd_gabriel bcv_svd_gabriel_t;

/**
 * bcv_gabriel_holdin_t:
 * @m: the number of rows in the held-in set.
 * @n: the number of columns on the held-in set.
 *
 * A #bcv_gabriel_holdin_t specifies the dimensions of the held-in matrix.
 */
typedef struct _bcv_gabriel_holdin {
    bcv_index_t m;
    bcv_index_t n;
} bcv_gabriel_holdin_t;

#define _bcv_assert_valid_gabriel_holdin(x,M,N) \
    assert (x); \
    assert (0 <= (x)->m && (x)->m <= (M)); \
    assert (0 <= (x)->n && (x)->n <= (N));


bcv_svd_gabriel_t *
bcv_svd_gabriel_alloc (bcv_gabriel_holdin_t max_holdin, bcv_index_t M, 
                       bcv_index_t N);

void
bcv_svd_gabriel_init (bcv_svd_gabriel_t *bcv, const bcv_matrix_t *x, 
                      const bcv_partition_t *rows, 
                      const bcv_partition_t *cols);

void
bcv_svd_gabriel_free (bcv_svd_gabriel_t *bcv);

bcv_error_t
bcv_svd_gabriel_get_rss (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                         bcv_index_t j, double *rss, bcv_index_t max_rank);

bcv_index_t
bcv_svd_gabriel_get_max_rank (const bcv_svd_gabriel_t *bcv, bcv_index_t i,
                              bcv_index_t j);

#endif /* _BCV_SVD_GABRIEL_H */
