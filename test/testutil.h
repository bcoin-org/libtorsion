/*!
 * testutil.h - test utils for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef _TORSION_TESTUTIL_H
#define _TORSION_TESTUTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef CHECK

#define CHECK(expr) do {                                   \
  if (!(expr))                                             \
    __torsion_test_assert_fail(__FILE__, __LINE__, #expr); \
} while (0)

static void
__torsion_test_assert_fail(const char *file, int line, const char *expr) {
  fprintf(stderr, "%s:%d: Assertion `%s' failed.\n", file, line, expr);
  fflush(stderr);
  abort();
}

#define ENTROPY_SIZE 32
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

#endif /* _TORSION_TESTUTIL_H */
