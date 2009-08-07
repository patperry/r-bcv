 
#ifndef _BCV_SVD_GABRIEL_REP_H
#define _BCV_SVD_GABRIEL_REP_H

#include <stdlib.h>
#include "bcv-svd-gabriel.h"
#include "bcv-types.h"

/** 
 * bcv_svd_grep_t:
 *  
 * A #bcv_svd_grep_t is a workspace for a perfoming a single replicate of a 
 * Gabrial-style SVD cross-validation.
 *
 * After being initiliazed with either bcv_svd_grep_init() or 
 * bcv_svd_grep_init_with_perm(), this structure maintains the current
 * residual estimate of the held-out block along with the estimated factors. 
 * The residuals can we updated with bcv_svd_grep_update_resid() and queried
 * with bcv_svd_grep_get_rss().
 */
typedef struct _bcv_svd_grep bcv_svd_grep_t;

/** 
 * bcv_svd_grep_alloc:
 * @holdin: the size of the held-in matrix.
 * @M: the number of rows in the matrix
 * @N: the number of columns in the matrix
 * 
 * Allocate memory for a BCV computation large enough to be used by
 * bcv_svd_grep_init().  The workspace should be freed with 
 * bcv_svd_grep_free().
 *
 * This function is monotonic in all of its arguments.
 */
bcv_svd_grep_t *
bcv_svd_grep_alloc (bcv_gabriel_holdin_t holdin, bcv_index_t M, bcv_index_t N);

/**
 * bcv_svd_grep_size:
 * @holdin: the size of the held-in matrix.
 * @M: the number of rows in the matrix
 * @N: the number of columns in the matrix
 *
 * Returns the size (in bytes) of a BCV replicate workspace necessary to 
 * cross-validate a matrix with the given dimensions, or 0 if the workspace
 * is larger than %SIZE_MAX bytes.
 */
size_t
bcv_svd_grep_size (bcv_gabriel_holdin_t holdin, bcv_index_t M, bcv_index_t N);

/**
 * bcv_svd_grep_align:
 *
 * Get the alignment of a bcv_svd_grep_t structure.
 */
size_t
bcv_svd_grep_align ();

/**
 * bcv_svd_grep_free:
 * @bcv: the BCV workspace
 * 
 * Free the BCV workspace.
 */
void 
bcv_svd_grep_free (bcv_svd_grep_t *bcv);

/**
 * bcv_svd_grep_init:
 * @bcv: unitialized memory for a #bcv_svd_grep_t
 * @holdin: the size of the held-in matrix
 * @x: a matrix to cross-validate
 *
 * Initialize the @bcv structure for a BCV computation with the given holdin
 * dimensions.  Divide @x into blocks, copy it into the BCV workspace, and
 * take the SVD of the held-in block.  Return zero on success and a positive
 * number on failure to compute the SVD of the held-in block.
 *
 * The @bcv pointer must point to at  least bcv_svd_grep_size() free bytes.  
 * This function does not allocate any memory.
 */
bcv_error_t
bcv_svd_grep_init (bcv_svd_grep_t *bcv, bcv_gabriel_holdin_t holdin, 
                   const bcv_matrix_t *x);

/**
 * bcv_svd_grep_init_with_perm:
 * @bcv: uninitialized memory for a #bcv_svd_grep_t
 * @holdin: the size of the held-in matrix
 * @x: a matrix to cross-validate
 * @p: a row permutation
 * @q: a column permutation
 *
 * Initialize the @bcv structure for a BCV computation with the given holdin
 * dimensions.  Divide @x into blocks, copy it into the BCV workspace while
 * optionally permuting the rows or columns of the copy.  If either of
 * @p or @q is non-null, the matrix rows or columns get permuted.  
 * In this case that @p and @q are both non-null, @x[i,j] gets replaced by
 * @x[p[i], q[j]].  Finally, take the SVD of the held-in block.  Return zero 
 * on success and a positive number on failure to compute the SVD of the 
 * held-in block.
 *
 * The @bcv pointer must point to at  least bcv_svd_grep_size() free bytes.  
 * This function does not allocate any memory.
 */
bcv_error_t
bcv_svd_grep_init_with_perm (bcv_svd_grep_t *bcv, bcv_gabriel_holdin_t holdin,
                             const bcv_matrix_t *x, bcv_index_t *p,
                             bcv_index_t *q);

/**
 * bcv_svd_grep_get_resid:
 * @bcv: the BCV workspace
 * @resid: the residual matrix
 *
 * Initialize @resid@ with the dimensions and data pointer of the residual
 * from the held-out matrix.
 */
void 
bcv_svd_grep_get_resid (const bcv_svd_grep_t *bcv, bcv_matrix_t *resid);

/**
 * bcv_svd_get_press:
 * @bcv: the BCV workspace
 * 
 * Get the prediction error sum of squares of from the residual matrix.
 */
double 
bcv_svd_grep_get_press (const bcv_svd_grep_t *bcv);

/**
 * bcv_svd_grep_get_msep:
 * @bcv: the BCV workspace
 *
 * Get the mean square error of prediction from the residual matrix.
 */
double
bcv_svd_grep_get_msep (const bcv_svd_grep_t *bcv);

/**
 * bcv_svd_grep_get_max_rank:
 * @bcv: the BCV workspace
 *
 * Get the maximum SVD rank.  This is equal to the smallest dimension
 * of the held-in set.
 */
bcv_index_t
bcv_svd_grep_get_max_rank (const bcv_svd_grep_t *bcv);

/**
 * bcv_svd_grep_get_holdout_sizes:
 * @bcv: the BCV workspace
 * @m2: (output) the number of rows in the holdout set
 * @n2: (output) the number of columns in the holdout set
 *
 * Get the dimensions of the holdout matrix.
 */
void
bcv_svd_grep_get_holdout_sizes (const bcv_svd_grep_t *bcv, bcv_index_t *m2,
                                bcv_index_t *n2);

/**
 * bcv_svd_grep_update_resid
 * @bcv: the BCV workspace
 * @scale: a scalar to multiply the factor by before subtracting it
 *   from the residual matrix
 * @i: the factor index
 *
 * Replace the residual estimate x22 by x22 - scale d2[i] u2[i] v2[i]^T.
 * Normally, @scale will be 1.0.
 */
void 
bcv_svd_grep_update_resid (bcv_svd_grep_t *bcv, double scale, bcv_index_t k);

#endif /* _BCV_SVD_GABRIEL_REP_H */
