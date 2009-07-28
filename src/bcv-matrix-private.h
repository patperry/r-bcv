
#ifndef _BCV_MATRIX_PRIVATE_H
#define _BCV_MATRIX_PRIVATE_H

#include "bcv-types.h"

/**
 * bcv_matrix_transpose_t:
 *
 * A #bcv_matrix_transpose_t specifies a transpose operation to apply
 * to a matrix.  For matrix "A", it can specify one "A", "A^T", or "A^H"
 * by %BCV_MATRIX_NOTRANS, %BCV_MATRIX_TRANS, or %BCV_MATRIX_CONJTRANS,
 * respectively.
 */
typedef enum _bcv_matrix_transpose { 
    BCV_MATRIX_NOTRANS, 
    BCV_MATRIX_TRANS, 
    BCV_MATRIX_CONJTRANS 
} bcv_matrix_transpose_t;

/**
 * _bcv_matrix_set_constant
 * @a: a matrix
 * @value: a double value
 *
 * Set all elements of @a to the given @value.
 */
void
_bcv_matrix_set_constant (bcv_matrix_t *a, double value);

/**
 * _bcv_matrix_set_identity:
 * @a: a matrix
 *
 * Set the off-diagonal elements of @a to zero and set the diagonal elements
 * to one.
 */
void
_bcv_matrix_set_identity (bcv_matrix_t *a);

/**
 * _bcv_matrix_set_index
 * @a: a matrix
 * @value: a double value
 * @indices: an array of 0-based column-major indices into @a.
 * @num_indices: the length of the @indices array
 *
 * Set the specified elements in @a to the given @value.
 */
void
_bcv_matrix_set_indices (bcv_matrix_t *a, double value, 
                         const bcv_index_t *indices, bcv_index_t num_indices);

/**
 * _bcv_matrix_copy:
 * @dst: the destination matrix
 * @src: the source matrix
 *
 * Elementwise copy so that @dst[i,j] := src[i,j] for all i,j.  @dst
 * and @src must have the same dimensions.
 */
void
_bcv_matrix_copy (bcv_matrix_t *dst, const bcv_matrix_t *src);

/**
 * _bcv_matrix_permute_copy:
 * @dst: the destination matrix
 * @src: the source matrix
 * @p: a row permutation
 * @q: a column permutation
 *
 * Elementwise copy with a permutation so that @dst[p[i],q[j]] := src[i,j].
 * If p or q is NULL then the permutation is assumed to be the identity,
 * so that p[i] (respectively, q[i]) is equal to i.
 */
void
_bcv_matrix_permute_copy (bcv_matrix_t *dst, const bcv_matrix_t *src,
                          bcv_index_t *p, bcv_index_t *q);

/**
 * _bcv_matrix_norm_frob:
 * @a: a matrix
 *
 * Return the Frobenius norm of a matrix, equal to the square root of
 * the sum of squares of its elements.
 */
double
_bcv_matrix_norm_frob (const bcv_matrix_t *a);


/* low-level BLAS routines */

double
_bcv_blas_dnrm2 (const bcv_vector_t *x);

void
_bcv_blas_dscal (double alpha, bcv_vector_t *x);

void
_bcv_blas_dcopy (const bcv_vector_t *x, bcv_vector_t *y);

void
_bcv_blas_dgemv (bcv_matrix_transpose_t transA,
                 double alpha, const bcv_matrix_t *a, const bcv_vector_t *x,
                 double beta, bcv_vector_t *y);

void
_bcv_blas_dger (double alpha, const bcv_vector_t *x, const bcv_vector_t *y, 
                bcv_matrix_t *a);

void
_bcv_blas_dgemm (bcv_matrix_transpose_t transA,
                 bcv_matrix_transpose_t transB,
                 double alpha, const bcv_matrix_t *a,
                 const bcv_matrix_t *b, double beta,
                 bcv_matrix_t *c);

/* low-level LAPACK SVD routines */

typedef enum _bcv_matrix_uplo {
    BCV_MATRIX_UPPER,
    BCV_MATRIX_LOWER,
    BCV_MATRIX_UPLO_NA
} bcv_matrix_uplo_t;

typedef enum _bcv_matrix_side {
    BCV_MATRIX_LEFT,
    BCV_MATRIX_RIGHT
} bcv_matrix_side_t;
    
typedef enum _bcv_matrix_norm {
    BCV_MATRIX_NORM_MAXABS,
    BCV_MATRIX_NORM_ONE,
    BCV_MATRIX_NORM_INF,
    BCV_MATRIX_NORM_FROB
} bcv_matrix_norm_t;

typedef enum _bcv_matrix_svdjob {
    BCV_MATRIX_SVDJOB_ALL,
    BCV_MATRIX_SVDJOB_SOME,
    BCV_MATRIX_SVDJOB_OVERWRITE,
    BCV_MATRIX_SVDJOB_NONE
} bcv_matrix_svdjob_t;

typedef enum _bcv_matrix_vect {
    BCV_MATRIX_VECT_P,
    BCV_MATRIX_VECT_Q
} bcv_matrix_vect_t;

void
_bcv_lapack_dlacpy (bcv_matrix_uplo_t uplo, const bcv_matrix_t *a, 
                    bcv_matrix_t *b);
                    
double
_bcv_lapack_dlange (bcv_matrix_norm_t norm, const bcv_matrix_t *a, 
                    double *work);

bcv_index_t
_bcv_lapack_dlange_work_len (bcv_matrix_norm_t norm, 
                             bcv_index_t m, bcv_index_t n);

bcv_error_t
_bcv_lapack_dgesvd (bcv_matrix_svdjob_t jobu, bcv_matrix_svdjob_t jobvt,
                    bcv_matrix_t *a, double *s, bcv_matrix_t *u,
                    bcv_matrix_t *vt, double *work, bcv_index_t lwork);

bcv_index_t
_bcv_lapack_dgesvd_work_len (bcv_matrix_svdjob_t jobu, 
                             bcv_matrix_svdjob_t jobvt, bcv_index_t m,
                             bcv_index_t n);

bcv_error_t
_bcv_lapack_dgesdd (bcv_matrix_svdjob_t jobz,
                    bcv_matrix_t *a, double *s, bcv_matrix_t *u,
                    bcv_matrix_t *vt, double *work, bcv_index_t lwork,
                    bcv_index_t *iwork);

bcv_index_t
_bcv_lapack_dgesdd_work_len (bcv_matrix_svdjob_t jobz, 
                             bcv_index_t m, bcv_index_t n);

bcv_index_t
_bcv_lapack_dgesdd_iwork_len (bcv_index_t m, bcv_index_t n);


void
_bcv_lapack_dgebrd (bcv_matrix_t *a, double *d, double *e, double *tauq,
                    double *taup, double *work, bcv_index_t lwork);

bcv_index_t
_bcv_lapack_dgebrd_work_len (bcv_index_t m, bcv_index_t n);

void
_bcv_lapack_dormbr (bcv_matrix_vect_t vect, bcv_matrix_side_t side,
                    bcv_matrix_transpose_t trans, 
                    const bcv_matrix_t *a, double *tau,
                    bcv_matrix_t *c, double *work, bcv_index_t lwork);

bcv_index_t
_bcv_lapack_dormbr_work_len (bcv_matrix_vect_t vect, bcv_matrix_side_t side,
                             bcv_index_t ma, bcv_index_t na,
                             bcv_index_t mc, bcv_index_t nc);

bcv_error_t
_bcv_lapack_dbdsqr (bcv_matrix_uplo_t uplo, bcv_index_t n,
                    double *d, double *e, bcv_matrix_t *vt, bcv_matrix_t *u,
                    bcv_matrix_t *c, double *work);

bcv_index_t
_bcv_lapack_dbdsqr_work_len (bcv_index_t n, bcv_bool_t only_values);

#endif /* _BCV_MATRIX_PRIVATE_H */
