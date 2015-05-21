
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "bcv-align-private.h"
#include "bcv-partition.h"


bcv_partition_t *
bcv_partition_alloc (bcv_index_t n)
{
    bcv_partition_t *result = NULL;
    size_t size = bcv_partition_size (n);
    
    if (size > 0)
        result = malloc (size);
    
    return result;
}


size_t
bcv_partition_size (bcv_index_t n)
{
    size_t result = 0;

    assert (n >= 0);
    
    if (n <= (SIZE_MAX - sizeof (bcv_partition_t)) / sizeof (bcv_index_t))
    {
        result = sizeof (bcv_partition_t) + n * sizeof (bcv_index_t);
    }
    
    return result;
}


size_t
bcv_partition_align ()
{
    return _bcv_alignof (bcv_partition_t);
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
    char *mem = (char *)part; mem += sizeof (bcv_partition_t);
    
    assert (mem);
    assert (n >= 0);
    assert (k >= 0);
    assert (sets || n == 0);
    
    part->n    = n;
    part->k    = k;
    part->sets = (void *)mem;
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


bcv_index_t
bcv_partition_get_set (const bcv_partition_t *part, bcv_index_t i,
                       bcv_index_t *set)
{
    bcv_index_t set_size = 0;
    bcv_index_t n, idx;
    bcv_index_t *sets;
    
    assert (part);
    n    = part->n;
    sets = part->sets;

    assert (sets || n == 0);
    assert (set);

    for (idx = 0; idx < n; idx++)
    {
        if (sets[idx] == i)
        {
            set[set_size++] = idx;
        }
    }
    
    return set_size;
}


bcv_index_t
bcv_partition_get_size (const bcv_partition_t *part, bcv_index_t i)
{
    bcv_index_t size = 0, n, k;
    bcv_index_t *sets, *sets_begin, *sets_end;
    
    assert (part);
    n = part->n;
    k = part->k;
    sets_begin = part->sets;
    sets_end   = sets_begin + n;

    assert (0 <= i && i < k);
    assert (sets_begin || n == 0);

    for (sets = sets_begin; sets < sets_end; sets++)
    {
        if (*sets == i)
        {
            size++;
        }
    }
    
    return size;
}


void
bcv_partition_get_sizes (const bcv_partition_t *part,
                         bcv_index_t *sizes)
{
    bcv_index_t k, n, i;
    bcv_index_t *sets, *sets_begin, *sets_end;
    
    assert (part);
    n = part->n;
    k = part->k;
    sets_begin = part->sets;
    sets_end   = sets_begin + n;

    assert (sets_begin || n == 0);
    assert (sizes      || k == 0);

    /* initialize sizes to 0 */
    memset (sizes, 0, k * sizeof (bcv_index_t));
    
    for (sets = sets_begin; sets < sets_end; sets++)
    {
        i = *sets;
        assert (0 <= i && i < k);
        
        sizes[i] = sizes[i] + 1;
    }
}
