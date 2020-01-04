#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>

#ifndef BCRYPTO_HAS_GMP
static unsigned long
mpz_zerobits(const mpz_t n) {
  /* Note: mpz_ptr is undocumented. */
  /* https://gmplib.org/list-archives/gmp-discuss/2009-May/003769.html */
  /* https://gmplib.org/list-archives/gmp-devel/2013-February/002775.html */
  int sgn = mpz_sgn(n);
  unsigned long bits;

  if (sgn == 0)
    return 0;

  if (sgn < 0)
    mpz_neg((mpz_ptr)n, n);

  bits = mpz_scan1(n, 0);

  if (sgn < 0)
    mpz_neg((mpz_ptr)n, n);

  return bits;
}

/* `mpz_jacobi` is not implemented in mini-gmp. */
/* https://github.com/golang/go/blob/aadaec5/src/math/big/int.go#L754 */
static int
mpz_jacobi(const mpz_t x, const mpz_t y) {
  mpz_t a, b, c;
  unsigned long s, bmod8;
  int j;

  assert(mpz_sgn(y) > 0);
  assert(mpz_odd_p(y));

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  /* a = x */
  mpz_set(a, x);

  /* b = y */
  mpz_set(b, y);

  j = 1;

  /* if b < 0 */
  if (mpz_sgn(b) < 0) {
    /* if a < 0 */
    if (mpz_sgn(a) < 0)
      j = -1;

    /* b = -b */
    mpz_neg(b, b);
  }

  for (;;) {
    /* if b == 1 */
    if (mpz_cmp_ui(b, 1) == 0)
      break;

    /* if a == 0 */
    if (mpz_sgn(a) == 0) {
      j = 0;
      break;
    }

    /* a = a mod b */
    mpz_mod(a, a, b);

    /* if a == 0 */
    if (mpz_sgn(a) == 0) {
      j = 0;
      break;
    }

    /* s = a factors of 2 */
    s = mpz_zerobits(a);

    if (s & 1) {
      /* bmod8 = b mod 8 */
      bmod8 = mpz_getlimbn(b, 0) & 7;

      if (bmod8 == 3 || bmod8 == 5)
        j = -j;
    }

    /* c = a >> s */
    mpz_tdiv_q_2exp(c, a, s);

    /* if b mod 4 == 3 and c mod 4 == 3 */
    if ((mpz_getlimbn(b, 0) & 3) == 3 && (mpz_getlimbn(c, 0) & 3) == 3)
      j = -j;

    /* a = b */
    mpz_set(a, b);

    /* b = c */
    mpz_set(b, c);
  }

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);

  return j;
}
#endif
