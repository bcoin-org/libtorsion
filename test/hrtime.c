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
#  pragma comment(lib, "kernel32.lib")
uint64_t
torsion_hrtime(void) {
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
