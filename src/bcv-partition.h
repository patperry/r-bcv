
#ifndef _BCV_PARTITION_H
#define _BCV_PARTITION_H

#include <stdlib.h>
#include "bcv-types.h"

/**
 * bcv_partition_t:
 *
 * A #bcv_partition_t stores a partition of the integers [0..n) into
 * k different sets.
 */
typedef struct _bcv_partition
{
    /* the total number of elements in the partitions */
    bcv_index_t n;
    
    /* the number of sets in the partition */
    bcv_index_t k;

    /* an array of length n; 
     * sets[i], in the range [0..k), indicates which set
     * i belongs to. */
    bcv_index_t *sets; 
} bcv_partition_t;

/**
 * bcv_partition_alloc:
 * @n: the number of elements
 *
 * Allocate space large enough to hold a partition of @n
 * integers into k sets.  The memory should be freed with 
 * bcv_partition_free().
 */
bcv_partition_t *
bcv_partition_alloc (bcv_index_t n);

/**
 * bcv_partition_size:
 * @n: the number of elements
 *
 * Return the number of bytes needed to store a #bcv_partition_t of @n
 * integers into k sets.
 */
size_t
bcv_partition_size (bcv_index_t n);

/**
 * bcv_partition_align:
 *
 * Return the alignment of a #bcv_partition_t.
 */
size_t
bcv_partition_align ();

/**
 * bcv_partition_free:
 * @part: a partition
 *
 * Free a #bcv_partition_t object that was allocated by bcv_partition_alloc().
 */
void
bcv_partition_free (bcv_partition_t *part);

/**
 * bcv_partition_init:
 * @part: uninitialized memory
 * @n: the number of elements
 * @k: the number of sets in the partition
 * @sets: an array of length @n indicating which partition set the elements
 *     [0..n) belong to.
 *
 * Initialize the memory pointed to by @part with the given parameters.
 * The @part pointer must point to at least bcv_partition_size(n,k) free
 * bytes. 
 * 
 * This function does not allocate any memory.
 */
void
bcv_partition_init (bcv_partition_t *part, bcv_index_t n, bcv_index_t k, 
                    const bcv_index_t *sets);

/**
 * bcv_partition_get_perm:
 * @part: a partition of [0..n) into k sets
 * @i: the partition set index, a number in the range [0..k)
 * @p: an array of length n.
 *
 * Set @p to be a permutation of [0..n) that moves the elements in the 
 * @i-th partition set to the end and its complement to the beginning.
 * Return the size of the @i-th partition set.
 *
 * If m is the return value, then
 *     p[0], ..., p[n-m-1] are set to the indices of the complement, and
 *     p[n-m], ..., p[n-1] are set to the indices of the partition set.
 * The elements are sorted so that
 *     p[0] < ... < p[n-m-1], and
 *     p[n-m] < ... < p[n-1].
 */
bcv_index_t
bcv_partition_get_perm (const bcv_partition_t *part, bcv_index_t i,
                        bcv_index_t *p);

/**
 * bcv_partition_get_set:
 * @part: a partition of [0..n) into k sets
 * @i: the partition set index, a number in the range [0..k)
 * @set: an array of length equal to the size of the @i-th set.
 *
 * Set the elements of @set to contain the elements of the @i-th partition
 * set, sorted in ascending order.  Return the size of the set.
 */
bcv_index_t
bcv_partition_get_set (const bcv_partition_t *part, bcv_index_t i,
                       bcv_index_t *set);

/**
 * bcv_partition_get_size:
 * @part: a partition of [0..n) into k sets
 * @i: the partition set index, a number in the range [0..k)
 *
 * Return the size (in elements) of the @ith set in the partition.
 */
bcv_index_t
bcv_partition_get_size (const bcv_partition_t *part, bcv_index_t i);

/**
 * bcv_partition_get_sizes:
 * @part: a partition of [0..n) into k sets
 * @sizes: an array of length k.
 *
 * For i in [0..k), sets @sizes[i] to be the size of the ith set in
 * the partition, @part.
 */
void
bcv_partition_get_sizes (const bcv_partition_t *part,
                         bcv_index_t *sizes);


#endif /* _BCV_PARTITION_H */
