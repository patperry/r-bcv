
#ifndef _BCV_VECTOR_PRIVATE_H
#define _BCV_VECTOR_PRIVATE_H

#include "bcv-types.h"

/**
 * _bcv_vector_set_constant
 * @x: a vector
 * @value: a double value
 *
 * Set all elements of @x to the given @value.
 */
void
_bcv_vector_set_constant (bcv_vector_t *x, double value);

#endif /* _BCV_VECTOR_PRIVATE_H */
