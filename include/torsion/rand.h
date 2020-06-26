/*!
 * rand.h - global RNG for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef _TORSION_RAND_H
#define _TORSION_RAND_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include "common.h"

/*
 * Global API
 */

TORSION_EXTERN int
torsion_is_reentrant(void);

TORSION_EXTERN int
torsion_has_tls(void);

TORSION_EXTERN int
torsion_getentropy(void *dst, size_t size);

TORSION_EXTERN int
torsion_getrandom(void *dst, size_t size);

TORSION_EXTERN int
torsion_random(uint32_t *out);

TORSION_EXTERN int
torsion_uniform(uint32_t *out, uint32_t max);

#ifdef __cplusplus
}
#endif

#endif /* _TORSION_RAND_H */
