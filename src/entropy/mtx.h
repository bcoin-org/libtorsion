/*!
 * mtx.h - mutexes for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef _TORSION_MTX_H
#define _TORSION_MTX_H

#include <stdint.h>

#if defined(__CloudABI__) \
 || defined(__wasi__) \
 || defined(__EMSCRIPTEN__) \
 || defined(__wasm__) \
 || defined(__asmjs__) \
 || defined(__vxworks) \
 || defined(__fuchsia__)
/* nothing */
#else
#  define TORSION_HAS_MUTEX
#endif

#if defined(_WIN32) && defined(TORSION_HAS_MUTEX)
#include <windows.h>
typedef struct torsion_mutex_s {
  CRITICAL_SECTION lock;
  LONG flag;
  int locked;
} torsion_mutex_t;
#elif defined(HAVE_PTHREAD) && defined(TORSION_HAS_MUTEX)
#include <pthread.h>
typedef pthread_mutex_t torsion_mutex_t;
#define TORSION_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER
#else
typedef struct torsion_mutex_s {
  int flag;
  int locked;
} torsion_mutex_t;
#endif

uint64_t
torsion_getpid(void);

void
torsion_mutex_lock(torsion_mutex_t *mtx);

void
torsion_mutex_unlock(torsion_mutex_t *mtx);

#endif /* _TORSION_MTX_H */
