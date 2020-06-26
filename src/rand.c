/*!
 * rand.c - global RNG for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

/**
 * Global RNG
 *
 * We expose a global fork-aware and thread-safe
 * RNG. However, instead of using locks or atomics
 * for thread safety, we use thread local storage
 * for the global context, meaning the context is
 * not actually global, rather, it is local to
 * each thread. This avoids us having to link to
 * pthreads and deal with other OS compat issues.
 *
 * If TLS is not supported, we try to fall back
 * to pthreads or atomics (if they are available).
 */

#include <stddef.h>
#include <stdint.h>
#include <torsion/rand.h>
#include <torsion/rng.h>
#include "entropy/entropy.h"
#include "internal.h"
#include "tls.h"

/*
 * Global Lock
 */

#undef HAVE_UNIX
#undef HAVE_TLS
#undef HAVE_PTHREAD
#undef HAVE_ATOMICS

#if (defined(__unix) || defined(__unix__))    \
 || (defined(__APPLE__) && defined(__MACH__))
#  define HAVE_UNIX
#endif

#if defined(TORSION_HAVE_TLS)
#  define HAVE_TLS
#elif defined(HAVE_UNIX) && defined(_REENTRANT)
#  include <pthread.h>
static pthread_mutex_t rng_lock = PTHREAD_MUTEX_INITIALIZER;
#  define HAVE_PTHREAD
#elif TORSION_GNUC_PREREQ(4, 1)
static volatile int rng_lock;
#  define HAVE_ATOMICS
#endif

static void
rng_global_lock(void) {
#if defined(HAVE_PTHREAD)
  if (pthread_mutex_lock(&rng_lock) != 0)
    abort();
#elif defined(HAVE_ATOMICS)
  while (__sync_val_compare_and_swap(&rng_lock, 0, 1) != 0) {
#if defined(__i386__) || defined(__amd64__) || defined(__x86_64__)
    __asm__ __volatile__("pause");
#elif defined(__aarch64__)
    __asm__ __volatile__("yield");
#endif
  }
#endif
}

static void
rng_global_unlock(void) {
#if defined(HAVE_PTHREAD)
  if (pthread_mutex_unlock(&rng_lock) != 0)
    abort();
#elif defined(HAVE_ATOMICS)
  __sync_synchronize();
  rng_lock = 0;
  __sync_synchronize();
#endif
}

/*
 * Global Context
 */

static TORSION_TLS struct {
  rng_t rng;
  int started;
  uint64_t pid;
} rng_state;

static int
rng_global_init(void) {
  uint64_t pid = torsion_getpid();

  if (!rng_state.started || rng_state.pid != pid) {
    if (!rng_init(&rng_state.rng))
      return 0;

    rng_state.started = 1;
    rng_state.pid = pid;
  }

  return 1;
}

#ifdef TORSION_TEST
TORSION_EXTERN uintptr_t
__torsion_global_rng_addr(void) {
  void *rng = (void *)&rng_state.rng;
  return (uintptr_t)rng;
}
#endif

/*
 * Global API
 */

int
torsion_is_reentrant(void) {
#if defined(HAVE_TLS) || defined(HAVE_PTHREAD) || defined(HAVE_ATOMICS)
  return 1;
#else
  return 0;
#endif
}

int
torsion_has_tls(void) {
#if defined(HAVE_TLS)
  return 1;
#else
  return 0;
#endif
}

int
torsion_getentropy(void *dst, size_t size) {
  return torsion_sysrand(dst, size);
}

int
torsion_getrandom(void *dst, size_t size) {
  rng_global_lock();

  if (!rng_global_init()) {
    rng_global_unlock();
    return 0;
  }

  rng_generate(&rng_state.rng, dst, size);
  rng_global_unlock();

  return 1;
}

int
torsion_random(uint32_t *out) {
  rng_global_lock();

  if (!rng_global_init()) {
    rng_global_unlock();
    return 0;
  }

  *out = rng_random(&rng_state.rng);

  rng_global_unlock();

  return 1;
}

int
torsion_uniform(uint32_t *out, uint32_t max) {
  rng_global_lock();

  if (!rng_global_init()) {
    rng_global_unlock();
    return 0;
  }

  *out = rng_uniform(&rng_state.rng, max);

  rng_global_unlock();

  return 1;
}
