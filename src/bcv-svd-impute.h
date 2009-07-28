
#ifndef _BCV_SVD_IMPUTE_H
#define _BCV_SVD_IMPUTE_H

#include <stdlib.h>
#include "bcv-types.h"

typedef struct _bcv_svd_impute bcv_svd_impute_t;

bcv_svd_impute_t *
bcv_svd_impute_alloc (bcv_index_t m, bcv_index_t n);

void
bcv_svd_impute_free (bcv_svd_impute_t *impute);

size_t
bcv_svd_impute_size (bcv_index_t m, bcv_index_t n);

size_t
bcv_svd_impute_align ();

void
bcv_svd_impute_init (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices);

bcv_error_t
bcv_svd_impute_step (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices,
                     bcv_index_t k);

double
bcv_svd_impute_get_rss (const bcv_svd_impute_t *impute);

void
bcv_svd_impute_get_svd (const bcv_svd_impute_t *impute, bcv_matrix_t *udvt);

#endif /* _BCV_SVD_IMPUTE_H */
