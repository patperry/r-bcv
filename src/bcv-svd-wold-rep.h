
#ifndef _BCV_SVD_WOLD_REP_H
#define _BCV_SVD_WOLD_REP_H

#include <stdlib.h>
#include "bcv-svd-wold.h"
#include "bcv-types.h"

/**
 * bcv_svd_wrep_t:
 *
 * A #bcv_svd_wrep_t provides a workspace for performing a single replicate
 * of a Wold-style SVD cross-validation.
 */
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
 * bcv_svd_wrep_get_press:
 * @bcv: an initialized wold-style replicate workspace
 *
 * Return the current prediction error sum of squares (PRESS) between the
 * held-out set and its current estimate.
 */
double 
bcv_svd_wrep_get_press (const bcv_svd_wrep_t *bcv);

/**
 * bcv_svd_wrep_get_msep:
 * @bcv: an initialized wold-style replicate workspace
 *
 * Return the currentmean square error of prediction (MSEP) between the
 * held-out set and its current estimate.
 */
double 
bcv_svd_wrep_get_msep (const bcv_svd_wrep_t *bcv);


/**
 * bcv_svd_wrep_get_max_rank:
 * @bcv: an initialized wold-style replicate workspace
 *
 * Return the maximum alowable SVD rank.  This is equal to the minimum of
 * the dimensions of the matrix that @bcv was initialized with.
 */
bcv_index_t
bcv_svd_wrep_get_max_rank (const bcv_svd_wrep_t *bcv);

/**
 * bcv_svd_wrep_get_holdout_size:
 * @bcv: a Wold-style replicate workspace
 *
 * Return the number of elements in the holdout set.
 */
bcv_index_t
bcv_svd_wrep_get_holdout_size (const bcv_svd_wrep_t *bcv);

/**
 * bcv_svd_wrep_impute_step:
 * @bcv: an initialized wold-style replicate workspace
 * @k: the rank of the SVD approximation
 * @rss: (out) the RSS between the training set and the estimate
 *
 * Perform one step of the rank @k SVD imputation algorithm, using the 
 * current estimate in @bcv for initial values.  Return the residual sum of
 * squares (RSS) between the training set and the current SVD estimate
 * in $*rss.  Return 0 on success and non-zero on failure to compute the 
 * SVD of the current estimate of @x.
 */
bcv_error_t 
bcv_svd_wrep_impute_step (bcv_svd_wrep_t *bcv, bcv_index_t k, 
                          double *rss);


#endif /* _BCV_SVD_WOLD_REP_H */
