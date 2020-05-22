/*!
 * mtx.c - mutexes for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include "mtx.h"

#ifndef __has_feature
#  define __has_feature(x) 0
#endif

#if defined(__GNUC__) && defined(__GNUC_MINOR__)
#  define TORSION_GNUC_PREREQ(maj, min) \
    ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
#  define TORSION_GNUC_PREREQ(maj, min) 0
#endif

#ifdef TORSION_HAS_MUTEX
#  if defined(_WIN32)
#    include <windows.h>
#    define HAVE_WINDOWS_MUTEX
#  elif defined(HAVE_PTHREAD)
#    include <pthread.h>
#    define HAVE_PTHREAD_MUTEX
#  elif defined(__clang__) && __has_feature(c_atomic)
#    define HAVE_UNIX_MUTEX
#    define HAVE_CLANG_ATOMICS
#  elif TORSION_GNUC_PREREQ(4, 7)
#    define HAVE_UNIX_MUTEX
#    define HAVE_GNUC_ATOMICS
#  elif defined(__GNUC__)
#    define HAVE_UNIX_MUTEX
#    define HAVE_LEGACY_ATOMICS
#  endif
#endif /* TORSION_HAS_MUTEX */

#ifdef HAVE_UNIX_MUTEX
#  include <unistd.h> /* getpid */
#endif

uint64_t
torsion_getpid(void) {
#ifdef TORSION_HAS_MUTEX
#ifdef _WIN32
  return GetCurrentProcessId();
#else
  return getpid();
#endif
#else /* TORSION_HAS_MUTEX */
  return 0;
#endif /* TORSION_HAS_MUTEX */
}

#if defined(HAVE_CLANG_ATOMICS)
static int
atomic_flag_test_and_set(int *flag) {
  int expected = 0;
  return __c11_atomic_compare_exchange_strong(flag, &expected, 1, 5, 5);
}
static void
atomic_flag_clear(int *flag) {
  __c11_atomic_store(flag, 0, 5);
}
#elif defined(HAVE_GNUC_ATOMICS)
static int
atomic_flag_test_and_set(int *flag) {
  int expected = 0;
  return __atomic_compare_exchange_n(flag, &expected, 1, 0, 5, 5);
}
static void
atomic_flag_clear(int *flag) {
  __atomic_store_n(flag, 0, 5);
}
#elif defined(HAVE_LEGACY_ATOMICS)
static int
atomic_flag_test_and_set(int *flag) {
  return __sync_val_compare_and_swap(flag, 0, 1) == 0;
}
static void
atomic_flag_clear(int *flag) {
  __sync_synchronize();
  *flag = 0;
  __sync_synchronize();
}
#endif

#if defined(HAVE_WINDOWS_MUTEX)
void
torsion_mutex_lock(torsion_mutex_t *mtx) {
  /* Borrowed from libsodium. */
  LONG status = 0;

  for (;;) {
    status = InterlockedCompareExchange(&mtx->flag, 1, 0);

    if (status != 1)
      break;

    Sleep(0);
  }

  switch (status) {
    case 0:
      InitializeCriticalSection(&mtx->lock);
      if (InterlockedExchange(&mtx->flag, 2) != 1)
        abort();
      break;
    case 2:
      break;
    default:
      abort();
      return;
  }

  EnterCriticalSection(&mtx->lock);

  if (mtx->locked != 0)
    abort();

  mtx->locked = 1;

  return 0;
}

void
torsion_mutex_unlock(torsion_mutex_t *mtx) {
  if (mtx->locked != 1)
    abort();

  mtx->locked = 0;

  LeaveCriticalSection(&mtx->lock);
}
#elif defined(HAVE_PTHREAD_MUTEX)
void
torsion_mutex_lock(torsion_mutex_t *mtx) {
  if (pthread_mutex_lock(mtx) != 0)
    abort();
}

void
torsion_mutex_unlock(torsion_mutex_t *mtx) {
  if (pthread_mutex_unlock(mtx) != 0)
    abort();
}
#elif defined(HAVE_UNIX_MUTEX)
void
torsion_mutex_lock(torsion_mutex_t *mtx) {
  /* Borrowed from libsodium. */
  while (!atomic_flag_test_and_set(&mtx->flag)) {
#if defined(__i386__) || defined(__x86_64__)
    __asm__ __volatile__("pause");
#elif defined(__aarch64__)
    __asm__ __volatile__("yield");
#endif
  }

  if (mtx->locked != 0)
    abort();

  mtx->locked = 1;
}

void
torsion_mutex_unlock(torsion_mutex_t *mtx) {
  if (mtx->locked != 1)
    abort();

  mtx->locked = 0;

  atomic_flag_clear(&mtx->flag);
}
#else
void
torsion_mutex_lock(torsion_mutex_t *mtx) {
  (void)mtx;
}

void
torsion_mutex_unlock(torsion_mutex_t *mtx) {
  (void)mtx;
}
#endif
