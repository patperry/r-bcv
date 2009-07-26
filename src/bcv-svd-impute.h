
#ifndef _BCV_SVD_IMPUTE_H
#define _BCV_SVD_IMPUTE_H

#include <stdlib.h>
#include "bcv-types.h"

typedef struct _bcv_svd_impute bcv_svd_impute_t;

bcv_svd_impute_t *
bcv_svd_impute_alloc (bcv_index_t m, bcv_index_t n);

void
bcv_svd_impute_free (bcv_svd_impute_t *impute);

/*
size_t
bcv_svd_impute_size (bcv_index_t m, bcv_index_t n, bcv_index_t missing);

size_t
bcv_svd_impute_align ();
*/

bcv_error_t
bcv_svd_impute_init (bcv_svd_impute_t *impute, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices);



#endif /* _BCV_SVD_IMPUTE_H */
