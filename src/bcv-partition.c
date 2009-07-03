
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "bcv-partition.h"

/* a _bcv_partition represents a partition of the integers
 * [0..n) into k sets. */
struct _bcv_partition
{
    /* the total number of elements in the partitions */
    bcv_index_t n;
    
    /* the number of sets in the partition */
    bcv_index_t k;

    /* an array of length n; 
     * sets[i], in the range [0..k), indicates which set
     * i belongs to. */
    bcv_index_t *sets; 
};


bcv_partition_t *
bcv_partition_alloc (bcv_index_t n, bcv_index_t k)
{
    bcv_partition_t *part = NULL;
    void *work;
    size_t size = bcv_partition_size (n, k);
    
    if (size > 0
        && (work = malloc (size)))
    {
        part       = work; work += sizeof (bcv_partition_t);
        part->sets = work;
    }
    
    return part;
}


size_t
bcv_partition_size (bcv_index_t n, bcv_index_t k)
{
    size_t result = 0;

    assert (n >= 0);
    
    if (n <= (SIZE_MAX - sizeof (bcv_partition_t)) / sizeof (bcv_index_t))
    {
        result = sizeof (bcv_partition_t) + n * sizeof (bcv_index_t);
    }
    
    return result;
}

void
bcv_partition_free (bcv_partition_t *part)
{
    if (part)
        free (part);
}


void 
bcv_partition_init (bcv_partition_t *part, bcv_index_t n, bcv_index_t k, 
                    const bcv_index_t *sets)
{
    assert (part);
    assert (n >= 0);
    assert (k >= 0);
    assert (sets || n == 0);
    
    part->n = n;
    part->k = k;
    memcpy (part->sets, sets, n * sizeof (bcv_index_t));
}


bcv_index_t
bcv_partition_get_perm (const bcv_partition_t *part, bcv_index_t part_index,
                        bcv_index_t *p)
{
    bcv_index_t i, m, mc, n, k, *sets;
    assert (part);
    
    n    = part->n;
    k    = part->k;
    sets = part->sets;
    m    = 0;
    mc   = n;

    assert (0 <= part_index && part_index < k);
    
    if (n > 0)
    {
        assert (sets);
        assert (p);
    
        /* put indices in the partition set at the end and indices in
         * the complement at the beginning */
        for (i = 0; i < n; i++)
        {
            if (sets[i] != part_index)
                p[i] = m++;
            else
                p[i] = --mc;
        }
    
        /* restore the relative ordering of the indices */
        for (i = 0; i < n; i++)
            if (p[i] >= m)
                p[i] = (n - 1) - (p[i] - m);
    }
    
    return m;
}
