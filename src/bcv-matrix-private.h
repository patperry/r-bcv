
#ifndef _BCV_MATRIX_PRIVATE_H
#define _BCV_MATRIX_PRIVATE_H

#include "bcv-types.h"

typedef enum _bcv_matrix_transpose { 
    BCV_MATRIX_NOTRANS, 
    BCV_MATRIX_TRANS, 
    BCV_MATRIX_CONJTRANS 
} bcv_matrix_transpose_t;

void
_bcv_matrix_set_identity (bcv_matrix_t *a);

void
_bcv_matrix_copy (bcv_matrix_t *dst, const bcv_matrix_t *src);

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

typedef enum _bcv_matrix_vect {
    BCV_MATRIX_VECT_P,
    BCV_MATRIX_VECT_Q
} bcv_matrix_vect_t;

void
_bcv_lapack_dlacpy (bcv_matrix_uplo_t uplo, const bcv_matrix_t *a, 
                    bcv_matrix_t *b);
                    
double
_bcv_lapack_dlange (bcv_matrix_norm_t norm, const bcv_matrix_t *a);

void
_bcv_lapack_dgebrd (bcv_matrix_t *a, double *d, double *e, double *tauq,
                    double *taup, double *work, bcv_index_t lwork);

void
_bcv_lapack_dormbr (bcv_matrix_vect_t vect, bcv_matrix_side_t side,
                    bcv_matrix_transpose_t trans, 
                    const bcv_matrix_t *a, double *tau,
                    bcv_matrix_t *c, double *work, bcv_index_t lwork);

bcv_error_t
_bcv_lapack_dbdsqr (bcv_matrix_uplo_t uplo, bcv_index_t n,
                    double *d, double *e, bcv_matrix_t *vt, bcv_matrix_t *u,
                    bcv_matrix_t *c, double *work);

#endif /* _BCV_MATRIX_PRIVATE_H */
