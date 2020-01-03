#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>

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

#ifndef BCRYPTO_HAS_GMP
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

static int
mpz_legendre(const mpz_t x, const mpz_t y) {
  mpz_t e, s;
  int a, b, c;

  assert(mpz_sgn(y) > 0);
  assert(mpz_odd_p(y));

  mpz_init(e);
  mpz_init(s);

  /* e = (y - 1) >> 1 */
  mpz_sub_ui(e, y, 1);
  mpz_tdiv_q_2exp(e, e, 1);

  /* s = x^e mod y */
  mpz_powm(s, x, e, y);

  a = mpz_sgn(s) == 0;
  b = mpz_cmp_ui(s, 1) == 0;
  mpz_add_ui(s, s, 1);
  c = mpz_cmp(s, y) == 0;

  assert(a + b + c == 1);

  mpz_clear(e);
  mpz_clear(s);

  /* 0, 1, or -1. */
  return b - c;
}
#endif

/* https://github.com/golang/go/blob/c86d464/src/math/big/int.go#L906 */
static int
mpz_sqrtm(mpz_t ret, const mpz_t num, const mpz_t p) {
  int r = 0;
  mpz_t x, e, t, a, s, n, y, b, g;
  unsigned long z, k;

  mpz_init(x);
  mpz_init(e);
  mpz_init(t);
  mpz_init(a);
  mpz_init(s);
  mpz_init(n);
  mpz_init(y);
  mpz_init(b);
  mpz_init(g);

  /* x = num */
  mpz_set(x, num);

  /* if p <= 0 or p mod 2 == 0 */
  if (mpz_sgn(p) <= 0 || mpz_even_p(p))
    goto fail;

  /* if x < 0 or x >= p */
  if (mpz_sgn(x) < 0 || mpz_cmp(x, p) >= 0) {
    /* x = x mod p */
    mpz_mod(x, x, p);
  }

  /* if p mod 4 == 3 */
  if ((mpz_getlimbn(p, 0) & 3) == 3) {
    /* b = x^((p + 1) / 4) mod p */
    mpz_add_ui(e, p, 1);
    mpz_tdiv_q_2exp(e, e, 2);
    mpz_powm(b, x, e, p);

    /* g = b^2 mod p */
    mpz_mul(g, b, b);
    mpz_mod(g, g, p);

    /* g != x */
    if (mpz_cmp(g, x) != 0)
      goto fail;

    /* ret = b */
    mpz_set(ret, b);

    goto succeed;
  }

  /* if p mod 8 == 5 */
  if ((mpz_getlimbn(p, 0) & 7) == 5) {
    /* t = x * 2 mod p */
    mpz_mul_2exp(t, x, 1);
    mpz_mod(t, t, p);

    /* a = t^((p - 5) / 8) mod p */
    mpz_tdiv_q_2exp(e, p, 3);
    mpz_powm(a, t, e, p);

    /* b = (a^2 * t - 1) * x * a mod p */
    mpz_mul(b, a, a);
    mpz_mod(b, b, p);
    mpz_mul(b, b, t);
    mpz_mod(b, b, p);
    mpz_sub_ui(b, b, 1);
    mpz_mod(b, b, p);
    mpz_mul(b, b, x);
    mpz_mod(b, b, p);
    mpz_mul(b, b, a);
    mpz_mod(b, b, p);

    /* g = b^2 mod p */
    mpz_mul(g, b, b);
    mpz_mod(g, g, p);

    /* g != x */
    if (mpz_cmp(g, x) != 0)
      goto fail;

    /* ret = b */
    mpz_set(ret, b);

    goto succeed;
  }

  /* if p == 1 */
  if (mpz_cmp_ui(p, 1) == 0)
    goto fail;

  switch (mpz_jacobi(x, p)) {
    case -1:
      goto fail;
    case 0:
      mpz_set_ui(ret, 0);
      goto succeed;
    case 1:
      break;
  }

  /* s = p - 1 */
  mpz_sub_ui(s, p, 1);

  /* z = s factors of 2 */
  z = mpz_zerobits(s);

  /* s = s >> z */
  mpz_tdiv_q_2exp(s, s, z);

  /* n = 2 */
  mpz_set_ui(n, 2);

  /* while n^((p - 1) / 2) != -1 mod p */
  while (mpz_jacobi(n, p) != -1) {
    /* n = n + 1 */
    mpz_add_ui(n, n, 1);
  }

  /* y = x^((s + 1) / 2) mod p */
  mpz_add_ui(y, s, 1);
  mpz_tdiv_q_2exp(y, y, 1);
  mpz_powm(y, x, y, p);

  /* b = x^s mod p */
  mpz_powm(b, x, s, p);

  /* g = n^s mod p */
  mpz_powm(g, n, s, p);

  /* k = z */
  k = z;

  for (;;) {
    unsigned long m = 0;

    /* t = b */
    mpz_set(t, b);

    /* while t != 1 */
    while (mpz_cmp_ui(t, 1) != 0) {
      /* t = t^2 mod p */
      mpz_mul(t, t, t);
      mpz_mod(t, t, p);
      m += 1;
    }

    /* if m == 0 */
    if (m == 0)
      break;

    /* if m >= k */
    if (m >= k)
      goto fail;

    /* t = g^(2^(k - m - 1)) mod p */
    mpz_set_ui(t, 1);
    mpz_mul_2exp(t, t, k - m - 1);
    mpz_powm(t, g, t, p);

    /* g = t^2 mod p */
    mpz_mul(g, t, t);
    mpz_mod(g, g, p);

    /* y = y * t mod p */
    mpz_mul(y, y, t);
    mpz_mod(y, y, p);

    /* b = b * g mod p */
    mpz_mul(b, b, g);
    mpz_mod(b, b, p);

    /* k = m */
    k = m;
  }

  /* ret = y */
  mpz_set(ret, y);
succeed:
  r = 1;
fail:
  mpz_clear(x);
  mpz_clear(e);
  mpz_clear(t);
  mpz_clear(a);
  mpz_clear(s);
  mpz_clear(n);
  mpz_clear(y);
  mpz_clear(b);
  mpz_clear(g);
  return r;
}
