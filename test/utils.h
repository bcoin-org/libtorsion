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

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
/* Avoid a GCC bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95189 */
#  define torsion_memcmp torsion__test_memcmp
#else
#  include <string.h>
#  define torsion_memcmp memcmp
#endif

void
torsion__test_assert_fail(const char *file, int line, const char *expr);

int
torsion__test_memcmp(const void *s1, const void *s2, size_t n);

uint64_t
torsion_hrtime(void);

#endif /* TORSION_UTILS_H */
