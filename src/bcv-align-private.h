
#ifndef _BCV_ALIGN_PRIVATE_H
#define _BCV_ALIGN_PRIVATE_H

#include <stdint.h>

/**
 * _bcv_alignof:
 *
 * Get the alignment of a type.
 */
#define _bcv_alignof __alignof__

/**
 * BCV_ALIGN_PTR:
 * @p: a pointer
 * @a: an alignment (typically a small power of 2)
 *
 * Advance the pointer @p by the minimum number of bytes necessary so that
 * the result is divisible by @a.
 */
#define BCV_ALIGN_PTR(p,a) \
            ((void *) \
             (((char *) (p) + (size_t) (a) - 1) \
              - \
              ((uintptr_t) (((char *) (p) + (size_t) (a) - 1)) \
               & \
               ((size_t) (a)-1))))

#endif /* _BCV_ALIGN_PRIVATE_H */
