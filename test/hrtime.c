/*!
 * hrtime.c - high-resolution time for libtorsion benchmarks
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifdef __linux__
#  define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdint.h>

#ifdef _WIN32
#  include <windows.h>
uint64_t
torsion_hrtime(void) {
  /* From libuv. See:
   *   https://github.com/libuv/libuv/blob/7967448/src/win/util.c#L77
   *   https://github.com/libuv/libuv/blob/7967448/src/win/util.c#L493
   *
   * Warning: not reentrant.
   */
  static double interval = 0;
  LARGE_INTEGER counter;

  if (interval == 0) {
    LARGE_INTEGER fequency;

    if (QueryPerformanceFrequency(&fequency))
      interval = 1.0 / fequency.QuadPart;
    else
      interval= 0;
  }

  if (!QueryPerformanceCounter(&counter))
    return 0;

  return (uint64_t)((double)counter.QuadPart * interval * 1000000000);
}
#else /* _WIN32 */
#  include <time.h>
#  ifdef CLOCK_MONOTONIC
uint64_t
torsion_hrtime(void) {
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0)
    abort();

  return (uint64_t)ts.tv_sec * 1000000000 + (uint64_t)ts.tv_nsec;
}
#  else /* CLOCK_MONOTONIC */
#    include <sys/time.h>
uint64_t
torsion_hrtime(void) {
  struct timeval tv;

  if (gettimeofday(&tv, NULL) != 0)
    abort();

  return (uint64_t)tv.tv_sec * 1000000000 + (uint64_t)tv.tv_usec * 1000;
}
#  endif /* CLOCK_MONOTONIC */
#endif /* _WIN32 */
