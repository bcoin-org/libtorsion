/*!
 * utils.c - test utils for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

void
__torsion_test_assert_fail(const char *file, int line, const char *expr) {
  fprintf(stderr, "%s:%d: Assertion `%s' failed.\n", file, line, expr);
  fflush(stderr);
  abort();
}

int
__torsion_test_memcmp(const void *s1, const void *s2, size_t n) {
  const unsigned char *x = (const unsigned char *)s1;
  const unsigned char *y = (const unsigned char *)s2;
  size_t i;

  for (i = 0; i < n; i++) {
    if (x[i] != y[i])
      return (int)x[i] - (int)y[i];
  }

  return 0;
}
