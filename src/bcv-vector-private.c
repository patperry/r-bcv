
#include <assert.h>
#include "bcv-vector-private.h"

void
_bcv_vector_set_constant (bcv_vector_t *x, double value)
{
    bcv_index_t n, inc;
    double *ptr, *start, *end;
    
    _bcv_assert_valid_vector (x);
    
    n     = x->n;
    inc   = x->inc;
    start = x->data;

    if (inc == 1)
    {
        end = start + n;

        for (ptr = start; ptr < end; ptr++)
        {
            *ptr = value;
        }
    }
    else
    {
        end = start + n*inc;
        
        for (ptr = start; ptr < end; ptr += inc)
        {
            *ptr = value;
        }
    }
}
