
#ifndef _BCV_SVD_WOLD_H
#define _BCV_SVD_WOLD_H

#include <assert.h>
#include <stdlib.h>
#include "bcv-partition.h"
#include "bcv-types.h"

/**
 * bcv_svd_wold_t:
 *
 * A #bcv_svd_wold_t provides a workspace for performing a Wold-style SVD 
 * cross-validation.
 */
typedef struct _bcv_svd_wold bcv_svd_wold_t;

/**
 * bcv_wold_holdout_t:
 *
 * A #bcv_wold_holdout_t holds the indices of the elements in held-out set.
 * They are assumed to be column-major indices in a matrix, and they must
 * be sorted in ascending order.
 */
typedef struct _bcv_wold_holdout {
    bcv_index_t *indices;
    bcv_index_t num_indices;
} bcv_wold_holdout_t;

#define _bcv_assert_valid_wold_holdout(x,M,N) \
    assert (x); \
    assert (0 <= (x)->num_indices && (x)->num_indices <= (M)*(N)); \
    assert ((x)->num_indices == 0 || (x)->indices)

/**
 * bcv_svd_wold_alloc:
 * @max_holdout: the size (in elements) of the bigest holdout set
 * @M: the number of rows in the matrix
 * @N: the number of columns in the matrix
 *
 * Allocate enough memory for a #bcv_svd_wold_t for a matrix of the given
 * dimensions and the given maximum holdout set size.  This memory should
 * be freed by bcv_svd_wold_free().
 */
bcv_svd_wold_t *
bcv_svd_wold_alloc (bcv_index_t max_holdout, bcv_index_t M, bcv_index_t N);

/**
 * bcv_svd_wold_alloc:
 * @max_holdout: the size (in elements) of the bigest holdout set
 * @M: the number of rows in the matrix
 * @N: the number of columns in the matrix
 *
 * Return the size (in bytes) for a #bcv_svd_wold_t for a matrix of the given
 * dimensions and the given maximum holdout set size.
 */
size_t
bcv_svd_wold_size (bcv_index_t max_holdout, bcv_index_t M, bcv_index_t N);

/**
 * bcv_svd_wold_align:
 *
 * Return the alignment of a #bcv_svd_wold_t.
 */
size_t
bcv_svd_wold_align ();

/**
 * bcv_svd_wold_init:
 * @bcv: a Wold CV workspace
 * @x: a matrix to cross-validate
 * @indices: a partition of the indices of @x into holdout sets
 *
 * Initialize the memory in @bcv for cross-validatign @x with the given
 * houldout sets.
 */
void
bcv_svd_wold_init (bcv_svd_wold_t *bcv, const bcv_matrix_t *x, 
                   const bcv_partition_t *indices);

/**
 * bcv_svd_wold_free:
 * @bcv: a Wold CV workspace
 *
 * Free the memory returned by bcv_svd_wold_alloc().
 */
void
bcv_svd_wold_free (bcv_svd_wold_t *bcv);

/**
 * bcv_svd_wold_get_press:
 * @bcv: an initialized Wold CV workspace
 * @i: the index of a hold-out set
 * @tol: the tolerance for determining when an SVD imputation has converged
 * @max_iter: the maximum number of imputation steps to take
 * @press: an array of length @max_rank+1
 * @max_rank: the rank to cross-validate up to
 *
 * Perform Wold-style cross-validation of the @i-th hold-out set.  Store
 * the prediction error sum of squares (PRESS) between the test set and its
 * estimate for ranks 0, 1, ..., @max_rank in @press[0], @press[1], ..., 
 * @press[@max_rank].
 *
 * Return 0 on success or nonzero on failure to compute an SVD in one of
 * the imputation steps.
 */
bcv_error_t
bcv_svd_wold_get_press (const bcv_svd_wold_t *bcv, bcv_index_t i,
                        double tol, bcv_index_t max_iter,
                        double *press, bcv_index_t max_rank);

/**
 * bcv_svd_wold_get_msep:
 * @bcv: an initialized Wold CV workspace
 * @i: the index of a hold-out set
 * @tol: the tolerance for determining when an SVD imputation has converged
 * @max_iter: the maximum number of imputation steps to take
 * @msep: an array of length @max_rank+1
 * @max_rank: the rank to cross-validate up to
 *
 * Perform Wold-style cross-validation of the @i-th hold-out set.  Store
 * the mean square error of prediction (MSEP) between the test set and its
 * estimate for ranks 0, 1, ..., @max_rank in @msep[0], @msep[1], ..., 
 * @msep[@max_rank].
 *
 * Return 0 on success or nonzero on failure to compute an SVD in one of
 * the imputation steps.
 */
bcv_error_t
bcv_svd_wold_get_msep (const bcv_svd_wold_t *bcv, bcv_index_t i,
                       double tol, bcv_index_t max_iter,
                       double *msep, bcv_index_t max_rank);

/**
 * bcv_svd_wold_get_max_rank:
 * @bcv: an initialized Wold CV workspace
 * @i: the index of a hold-out set
 *
 * Return the maximum rank allowable for cross-validating with the @i-th
 * holdout set.  This is equal to the minimum of the dimensions of the @x
 * that @bcv was initialized with.
 */
bcv_index_t
bcv_svd_wold_get_max_rank (const bcv_svd_wold_t *bcv, bcv_index_t i);

#endif /* _BCV_SVD_WOLD_H */
