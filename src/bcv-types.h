
#ifndef _BCV_TYPES_H
#define _BCV_TYPES_H

#include <limits.h>

/**
 * BCV_DBL_EPSILON:
 *
 * The smallest value x such that 1 + x != 1.
 */
#define BCV_DBL_EPSILON 2.2204460492503131e-16

/**
 * BCV_POSINF, BCV_NEGINF:
 *
 * Positive and negative infinities.
 */
#define BCV_DBL_POSINF ((double) 1.0 / 0.0)
#define BCV_DBL_NEGINF (-BCV_DBL_POSINF)

/** 
 * BCV_MIN:
 * @a: the first value
 * @b: the second value
 *
 * Return the minimum of @a and @b.
 */
#define BCV_MIN(a,b) ((a) < (b) ? (a) : (b))

/** 
 * BCV_MAX:
 * @a: the first value
 * @b: the second value
 *
 * Return the maximum of @a and @b.
 */
#define BCV_MAX(a,b) ((a) > (b) ? (a) : (b))

/**
 * bcv_bool_t:
 *
 * A boolean value.
 */
typedef enum _bcv_bool { 
    BCV_FALSE = 0, 
    BCV_TRUE 
} bcv_bool_t;

/**
 * bcv_index_t:
 *
 * A #bcv_index_t is the type of a matrix of vector index.  This type
 * is used for storing matrix and vector dimensions.
 */
typedef int bcv_index_t;


/**
 * BCV_MAX_INDEX:
 *
 * The maximum value of a #bcv_index_t.
 */
#define BCV_MAX_INDEX INT_MAX


/**
 * bcv_vector_t:
 * @n:    the number of elements in the vector
 * @data: a pointer to the vector data
 * @inc:  the stride of the vector
 * 
 * A #bcv_vector_t associates a pointer to vector-data with the vector
 * length and stride.
 */
typedef struct _bcv_vector {
    bcv_index_t n;
    double *data;
    bcv_index_t inc;
} bcv_vector_t;

#define _bcv_assert_valid_vector(x) \
    assert (x); \
    assert ((x)->n   >= 0); \
    assert ((x)->inc >= 1)


/**
 * bcv_matrix_t:
 * @m:    the number of rows in the matrix
 * @n:    the number of columns in the matrix
 * @data: a pointer to the matrix data
 * @lda:  the leading dimension of the matrix
 * 
 * A #bcv_matrix_t associates a pointer to matrix-data with the matrix
 * dimensions.  The matrix is assumed to be stored in column-major order.
 */
typedef struct _bcv_matrix {
    bcv_index_t m;
    bcv_index_t n;
    double *data;
    bcv_index_t lda;
} bcv_matrix_t;

#define _bcv_assert_valid_matrix(x) \
    assert (x); \
    assert ((x)->m >= 0); \
    assert ((x)->n >= 0); \
    assert ((x)->data); \
    assert ((x)->lda >= 1); \
    assert ((x)->lda >= (x)->m)


/**
 * bcv_error_t:
 *
 * A #bcv_error_t holds an error code.  This is 0 on success and nonzero
 * on error.
 */
typedef bcv_index_t bcv_error_t;

#endif
