/*!
 * utils.h - test utils for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef _TORSION_UTILS_H
#define _TORSION_UTILS_H

#include <stddef.h>

#undef ASSERT

#define ASSERT(expr) do {                                  \
  if (!(expr))                                             \
    __torsion_test_assert_fail(__FILE__, __LINE__, #expr); \
} while (0)

#define ENTROPY_SIZE 32
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
/* Avoid a GCC bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95189 */
#  define torsion_memcmp __torsion_test_memcmp
#else
#  include <string.h>
#  define torsion_memcmp memcmp
#endif

void
__torsion_test_assert_fail(const char *file, int line, const char *expr);

int
__torsion_test_memcmp(const void *s1, const void *s2, size_t n);

#endif /* _TORSION_UTILS_H */
