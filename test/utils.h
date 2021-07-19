/*!
 * utils.h - test utils for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef TORSION_UTILS_H
#define TORSION_UTILS_H

#include <stddef.h>
#include <stdint.h>

#undef ASSERT

#define ASSERT(expr) do {                                 \
  if (!(expr))                                            \
    torsion__test_assert_fail(__FILE__, __LINE__, #expr); \
} while (0)

#define ENTROPY_SIZE 32
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

void
torsion__test_assert_fail(const char *file, int line, const char *expr);

uint64_t
torsion_hrtime(void);

#ifndef TORSION_HAVE_CONFIG
#  undef TORSION_HAVE_CLOCK_GETTIME
#  undef TORSION_HAVE_FORK
#  undef TORSION_HAVE_GETTIMEOFDAY
#  undef TORSION_HAVE_TIME
#  if (defined(__linux__))                      \
   || (defined(__APPLE__) && defined(__MACH__)) \
   || (defined(__FreeBSD__))                    \
   || (defined(__OpenBSD__))                    \
   || (defined(__NetBSD__))                     \
   || (defined(__DragonFly__))                  \
   || (defined(__sun) && defined(__SVR4))       \
   || (defined(__CYGWIN__))
#    define TORSION_HAVE_GETTIMEOFDAY
#    define TORSION_HAVE_FORK
#    define TORSION_HAVE_TIME
#  endif
#endif

#endif /* TORSION_UTILS_H */
