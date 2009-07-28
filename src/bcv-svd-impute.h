
#ifndef _BCV_SVD_IMPUTE_H
#define _BCV_SVD_IMPUTE_H

#include <stdlib.h>
#include "bcv-types.h"

/**
 * bcv_svd_impute_t:
 *
 * A #bcv_svd_impute_t is a workspace for imputing missing values in a
 * matrix using an SVD approximation estimated by the EM algorithm.
 */
typedef struct _bcv_svd_impute bcv_svd_impute_t;

/**
 * bcv_svd_impute_alloc:
 * @m: the number of rows in the matrix
 * @n: the number of columns in the matrix
 *
 * Allocate enough memory to hold an SVDImpute workspace for a matrix of the
 * given dimensions.  The memory should be freed by bcv_svd_impute_free().
 */
bcv_svd_impute_t *
bcv_svd_impute_alloc (bcv_index_t m, bcv_index_t n);

/**
 * bcv_svd_impute_free:
 * @impute: an SVDImpute workspace.
 *
 * Free the memory returned by bcv_svd_impute_alloc().
 */
void
bcv_svd_impute_free (bcv_svd_impute_t *impute);

/**
 * bcv_svd_impute_size:
 * @m: the number of rows in the matrix
 * @n: the number of columns in the matrix
 *
 * Return the size (in bytes) of an SVDImpute workspace for a matrix of the
 * given dimensions.
 */
size_t
bcv_svd_impute_size (bcv_index_t m, bcv_index_t n);

/**
 * bcv_svd_impute_align:
 *
 * Return the alignment of an SVDImpute workspace.
 */
size_t
bcv_svd_impute_align ();

/**
 * bcv_svd_impute_init:
 * @impute: memory for an SVDImpute workspace
 * @xhat: a matrix to hold imputed values
 * @x: a matrix with missing values
 * @indices: the column-major indices of the missing elements in @x
 * @num_indices: the number of elements in the @indices array
 *
 * Initialize the memory in @impute for an SVDImpute computation.  Set
 * @xhat to be equal to @x at the non-misisng values, and equal to the
 * column means at the missing values.
 *
 * This function does not allocate any memory.
 */
void
bcv_svd_impute_init (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices);

/**
 * bcv_svd_impute_step:
 * @impute: an initialized SVDImpute workspace
 * @xhat: a matrix to hold imputed values
 * @x: a matrix with missing values
 * @indices: the column-major indices of the missing elements in @x
 * @num_indices: the number of elements in the @indices array
 * @k: the rank of the SVD approximation
 * 
 * Peform one step of the EM algorithm for estimating a rank-@k SVD 
 * approximation of @x.  Use the initial values for the missing values in 
 * @xhat.  On returning, update @xhat with the new missing values and
 * update @impute with the new SVD and RSS.
 *
 * Return 0 on success and a non-zero value on failing to compute the SVD
 * of @x.
 */
bcv_error_t
bcv_svd_impute_step (bcv_svd_impute_t *impute, 
                     bcv_matrix_t *xhat, const bcv_matrix_t *x, 
                     const bcv_index_t *indices, bcv_index_t num_indices,
                     bcv_index_t k);

/**
 * bcv_svd_impute_get_rss:
 * @impute: an SVDImpute workspace
 *
 * Get the RSS stored in @impute.  This function can only be called after
 * bcv_svd_impute_step().
 */
double
bcv_svd_impute_get_rss (const bcv_svd_impute_t *impute);

/**
 * bcv_svd_impute_get_rss:
 * @impute: an SVDImpute workspace
 * @udvt: a matrix of the same dimensions as the imputed workspace.
 *
 * Replace @udvt with the estimate of the SVD approximation.  This function 
 * can only be called after bcv_svd_impute_step().
 */
void
bcv_svd_impute_get_svd (const bcv_svd_impute_t *impute, bcv_matrix_t *udvt);

#endif /* _BCV_SVD_IMPUTE_H */
