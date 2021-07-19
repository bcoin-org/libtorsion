/*!
 * hrtime.c - high-resolution time for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdint.h>
#include <stdlib.h>
#include "utils.h"

/*
 * Backend
 */

#if defined(_WIN32)
#  include <windows.h> /* QueryPerformanceCounter */
#  ifndef __MINGW32__
#    pragma comment(lib, "kernel32.lib")
#  endif
#elif defined(TORSION_HAVE_CLOCK_GETTIME)
#  include <time.h> /* clock_gettime */
#elif defined(TORSION_HAVE_GETTIMEOFDAY)
#  include <sys/time.h> /* gettimeofday */
#elif defined(TORSION_HAVE_TIME)
#  include <time.h> /* time */
#endif

/*
 * High-Resolution Time
 */

uint64_t
torsion_hrtime(void) {
#if defined(_WIN32)
  LARGE_INTEGER freq, ctr;
  double sec;

  if (!QueryPerformanceFrequency(&freq))
    abort(); /* LCOV_EXCL_LINE */

  if (freq.QuadPart == 0)
    abort(); /* LCOV_EXCL_LINE */

  if (!QueryPerformanceCounter(&ctr))
    abort(); /* LCOV_EXCL_LINE */

  sec = (double)ctr.QuadPart / (double)freq.QuadPart;

  return (uint64_t)(sec * 1000000000.0);
#elif defined(TORSION_HAVE_CLOCK_GETTIME)
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) {
    if (clock_gettime(CLOCK_REALTIME, &ts) != 0)
      abort(); /* LCOV_EXCL_LINE */
  }

  return (uint64_t)ts.tv_sec * 1000000000 + (uint64_t)ts.tv_nsec;
#elif defined(TORSION_HAVE_GETTIMEOFDAY)
  struct timeval tv;

  if (gettimeofday(&tv, NULL) != 0)
    abort(); /* LCOV_EXCL_LINE */

  return (uint64_t)tv.tv_sec * 1000000000 + (uint64_t)tv.tv_usec * 1000;
#elif defined(TORSION_HAVE_TIME)
  time_t ts = time(NULL);

  if (ts == (time_t)-1)
    abort(); /* LCOV_EXCL_LINE */

  return (uint64_t)ts * 1000000000;
#else
  abort();
  return 0;
#endif
}
