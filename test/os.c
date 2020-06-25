/*!
 * os.c - platform-specific functionality for libtorsion tests
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#if !defined(_WIN32) && !defined(_GNU_SOURCE)
/* For clock_gettime(2). */
#  define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdint.h>
#include "os.h"

#undef HAVE_GETTIMEOFDAY

#if defined(_WIN32)
#  include <windows.h> /* QueryPerformance{Counter,Frequency}, CreateThread */
#  pragma comment(lib, "kernel32.lib")
#else /* _WIN32 */
#  include <time.h> /* clock_gettime, time */
#  ifndef CLOCK_MONOTONIC
#    if defined(__unix) || defined(__unix__)     \
     || (defined(__APPLE__) && defined(__MACH__))
#      include <sys/time.h> /* gettimeofday */
#      define HAVE_GETTIMEOFDAY
#    endif
#  endif
#endif /* _WIN32 */

#ifdef TORSION_HAVE_PTHREAD
#  include <pthread.h>
#endif

uint64_t
torsion_hrtime(void) {
#if defined(_WIN32)
  static unsigned int scale = 1000000000;
  LARGE_INTEGER freq, ctr;
  double scaled, result;

  if (!QueryPerformanceFrequency(&freq))
    abort();

  if (!QueryPerformanceCounter(&ctr))
    abort();

  if (freq.QuadPart == 0)
    abort();

  /* We have no idea of the magnitude of `freq`,
   * so we must resort to double arithmetic[1].
   * Furthermore, we use some wacky arithmetic
   * to avoid a bug in Visual Studio 2019[2][3].
   *
   * [1] https://github.com/libuv/libuv/blob/7967448/src/win/util.c#L503
   * [2] https://github.com/libuv/libuv/issues/1633
   * [3] https://github.com/libuv/libuv/pull/2866
   */
  scaled = (double)freq.QuadPart / scale;
  result = (double)ctr.QuadPart / scaled;

  return (uint64_t)result;
#elif defined(CLOCK_MONOTONIC)
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
    abort();

  return (uint64_t)ts.tv_sec * 1000000000 + (uint64_t)ts.tv_nsec;
#elif defined(HAVE_GETTIMEOFDAY)
  struct timeval tv;

  if (gettimeofday(&tv, NULL) != 0)
    abort();

  return (uint64_t)tv.tv_sec * 1000000000 + (uint64_t)tv.tv_usec * 1000;
#else
  return (uint64_t)time(NULL) * 1000000000;
#endif
}

#ifdef TORSION_HAVE_THREADS
struct torsion_thread_s {
#ifdef _WIN32
  HANDLE handle;
#else
  pthread_t handle;
#endif
};

struct torsion_thread_s *
torsion_thread_alloc(void) {
  struct torsion_thread_s *thread = malloc(sizeof(struct torsion_thread_s));

  if (thread == NULL)
    abort();

  return thread;
}

void
torsion_thread_free(struct torsion_thread_s *thread) {
  free(thread);
}

#ifdef _WIN32
typedef struct torsion_thread_args_s {
  torsion_thread_start_f *start_routine;
  void *arg;
} torsion_thread_args_t;

static DWORD WINAPI
torsion_thread_run(void *ptr) {
  torsion_thread_args_t args;

  args = *((torsion_thread_args_t *)ptr);

  free(ptr);

  (void *)args.start_routine(args.arg);

  return ERROR_SUCCESS;
}
#endif

int
torsion_thread_create(struct torsion_thread_s *thread,
                      const torsion_thread_attr_t *attr,
                      torsion_thread_start_f *start_routine,
                      void *arg) {
#ifdef _WIN32
  torsion_thread_args_t *args;

  if (attr != NULL)
    return -1;

  args = malloc(sizeof(torsion_thread_args_t));

  if (args == NULL)
    return -1;

  args->start_routine = start_routine;
  args->arg = arg;

  thread->handle = CreateThread(NULL, 0, torsion_thread_run, args, 0, NULL);

  if (thread->handle == NULL) {
    free(args);
    return -1;
  }

  return 0;
#else
  if (attr != NULL)
    return -1;

  return pthread_create(&thread->handle, NULL, start_routine, arg);
#endif
}

int
torsion_thread_join(struct torsion_thread_s *thread, void **retval) {
#ifdef _WIN32
  if (retval != NULL)
    return -1;

  WaitForSingleObject(thread->handle, INFINITE);

  if (CloseHandle(thread->handle) == FALSE)
    return -1;

  return 0;
#else
  if (retval != NULL)
    return -1;

  return pthread_join(thread->handle, NULL);
#endif
}
#endif /* TORSION_HAVE_THREADS */
