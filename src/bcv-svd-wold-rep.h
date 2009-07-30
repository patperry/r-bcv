
#ifndef _BCV_SVD_WOLD_REP_H
#define _BCV_SVD_WOLD_REP_H

#include <stdlib.h>
#include "bcv-svd-wold.h"
#include "bcv-types.h"

typedef struct _bcv_svd_wrep bcv_svd_wrep_t;

/**
 * bcv_svd_wrep_alloc:
 * @M the number of rows
 * @N the number of columns
 *
 * Allocate enough memory for a #bcv_svd_wrep_t to performing a single
 * replicate of a Wold-style SVD cross-validation for matrix of the given
 * dimensions.  The memory should be freed by bcv_svd_wrep_free().
 */
bcv_svd_wrep_t *
bcv_svd_wrep_alloc (bcv_index_t M, bcv_index_t N);

/**
 * bcv_svd_wrep_size:
 * @M the number of rows
 * @N the number of columns
 *
 * Return the size (in bytes) for a wold-style CV replicate for a matrix of
 * the given dimensions.
 */
size_t
bcv_svd_wrep_size (bcv_index_t M, bcv_index_t N);

/**
 * bcv_svd_wrep_align:
 *
 * Return the alignment of a #bcv_svd_wrep_t.
 */
size_t
bcv_svd_wrep_align ();

/**
 * bcv_svd_wrep_free:
 * @bcv: a wold-style replicate workspace
 *
 * Free the memory allocated by bcv_svd_wrep_alloc().
 */
void 
bcv_svd_wrep_free (bcv_svd_wrep_t *bcv);

/**
 * bcv_svd_wrep_init:
 * @bcv: a wold-style replicate workspace
 * @holdout: the set of indices to hold out
 * @x: a matrix
 *
 * Initialize the workspace for performing a CV replicate of the SVD of
 * @x using thing given @holdout set.
 */
void
bcv_svd_wrep_init (bcv_svd_wrep_t *bcv, bcv_wold_holdout_t holdout, 
                   const bcv_matrix_t *x);

/**
 * bcv_svd_wrep_get_rss:
 * @bcv: an initialized wold-style replicate workspace
 *
 * Return the current RSS between the held-in set and its current estimate.
 */
double 
bcv_svd_wrep_get_rss (const bcv_svd_wrep_t *bcv);

/**
 * bcv_svd_wrep_get_max_rank:
 * @bcv: an initialized wold-style replicate workspace
 *
 * Return the maximum alowable SVD rank.  This is equal to the minimum of
 * the dimensions of the matrix that @bcv was initialized with.
 */
bcv_index_t
bcv_svd_wrep_get_max_rank (bcv_svd_wrep_t *bcv);

/**
 * bcv_svd_wrep_impute_step:
 * @bcv: an initialized wold-style replicate workspace
 * @k: the rank of the SVD approximation
 * @train_rss: (out) the rss between the training set and the estimate
 *
 * Perform one step of the rank @k SVD imputation algorithm, using the 
 * current estimate in @bcv for initial values.  Return the training
 * set RSS after the step in @train_rss.  Return 0 on success and non-zero
 * on failure to compute the SVD of the current estimate of @x.
 */
bcv_error_t 
bcv_svd_wrep_impute_step (bcv_svd_wrep_t *bcv, bcv_index_t k, 
                          double *train_rss);


#endif /* _BCV_SVD_WOLD_REP_H */
