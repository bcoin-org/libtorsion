/*!
 * internal.c - internal utils for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
#  include <stdio.h>
#endif

#include <stdlib.h>
#include <string.h>
#include "internal.h"

void
__torsion_assert_fail(const char *file, int line, const char *expr) {
#if !defined(__EMSCRIPTEN__) && !defined(__wasm__)
  fprintf(stderr, "%s:%d: Assertion `%s' failed.\n", file, line, expr);
  fflush(stderr);
#else
  (void)file;
  (void)line;
  (void)expr;
#endif
  abort();
}

void
torsion_abort(void) {
  abort();
}
