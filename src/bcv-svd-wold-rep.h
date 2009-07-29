
#ifndef _BCV_SVD_WOLD_REP_H
#define _BCV_SVD_WOLD_REP_H

#include <stdlib.h>
#include "bcv-svd-wold.h"
#include "bcv-types.h"

typedef struct _bcv_svd_wrep bcv_svd_wrep_t;

bcv_svd_wrep_t *
bcv_svd_wrep_alloc (bcv_index_t M, bcv_index_t N);

size_t
bcv_svd_wrep_size (bcv_index_t M, bcv_index_t N);

size_t
bcv_svd_wrep_align ();

void 
bcv_svd_wrep_free (bcv_svd_wrep_t *bcv);

void
bcv_svd_wrep_init (bcv_svd_wrep_t *bcv, bcv_wold_holdout_t holdout, 
                   const bcv_matrix_t *x);

double 
bcv_svd_wrep_get_rss (const bcv_svd_wrep_t *bcv);

bcv_index_t
bcv_svd_wrep_get_max_rank (bcv_svd_wrep_t *bcv);

bcv_error_t 
bcv_svd_wrep_impute_step (bcv_svd_wrep_t *bcv, bcv_index_t k, 
                          double *train_rss);


#endif /* _BCV_SVD_WOLD_REP_H */
