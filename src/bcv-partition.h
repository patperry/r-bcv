
#ifndef _BCV_PARTITION_H
#define _BCV_PARTITION_H

#include <stdlib.h>
#include "bcv-types.h"


typedef struct _bcv_partition bcv_partition_t;

bcv_partition_t *
bcv_partition_alloc (bcv_index_t n);

size_t
bcv_partition_size (bcv_index_t n);

void
bcv_partition_free (bcv_partition_t *part);

void 
bcv_partition_init (bcv_partition_t *part, bcv_index_t n, bcv_index_t k, 
                    const bcv_index_t *sets);

bcv_index_t
bcv_partition_get_perm (const bcv_partition_t *part, bcv_index_t part_index,
                        bcv_index_t *p);

#endif /* _BCV_PARTITION_H */
