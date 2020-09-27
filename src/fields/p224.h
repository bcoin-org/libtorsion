/*!
 * p224.h - p224 field element for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#if defined(TORSION_HAVE_INT128)
typedef uint64_t p224_fe_word_t;
#define P224_FIELD_WORDS 4
#include "p224_64.h"
#else
typedef uint32_t p224_fe_word_t;
#define P224_FIELD_WORDS 7
#include "p224_32.h"
#endif

typedef p224_fe_word_t p224_fe_t[P224_FIELD_WORDS];

#define p224_fe_add fiat_p224_add
#define p224_fe_sub fiat_p224_sub
#define p224_fe_neg fiat_p224_opp
#define p224_fe_mul fiat_p224_mul
#define p224_fe_sqr fiat_p224_square

#if defined(TORSION_HAVE_INT128)
static const p224_fe_t p224_zero = {0, 0, 0, 0};

static const p224_fe_t p224_one = {
  UINT64_C(0xffffffff00000000), UINT64_C(0xffffffffffffffff),
  UINT64_C(0x0000000000000000), UINT64_C(0x0000000000000000)
};

/* 11^(2^128 - 1) mod p */
/* mont: 0xa31b1da46d3e2af0dd915e4b7869be5d866c223b174131b85ee27c6c */
static const p224_fe_t p224_g = {
  UINT64_C(0x174131b85ee27c6c), UINT64_C(0x7869be5d866c223b),
  UINT64_C(0x6d3e2af0dd915e4b), UINT64_C(0x00000000a31b1da4)
};
#else
static const p224_fe_t p224_zero = {0, 0, 0, 0, 0, 0, 0};

static const p224_fe_t p224_one = {
  0xffffffff, 0xffffffff, 0xffffffff, 0x00000000,
  0x00000000, 0x00000000, 0x00000000
};

/* 11^(2^128 - 1) mod p */
/* mont: 0xa11d8394a31b1da46d3e2af0dd915e4ad74c3ac9866c223b174131b9 */
static const p224_fe_t p224_g = {
  0x174131b9, 0x866c223b, 0xd74c3ac9, 0xdd915e4a,
  0x6d3e2af0, 0xa31b1da4, 0xa11d8394
};
#endif

static void
p224_fe_set(p224_fe_t r, const p224_fe_t x) {
  r[0] = x[0];
  r[1] = x[1];
  r[2] = x[2];
  r[3] = x[3];
#if P224_FIELD_WORDS == 7
  r[4] = x[4];
  r[5] = x[5];
  r[6] = x[6];
#endif
}

static int
p224_fe_equal(const p224_fe_t x, const p224_fe_t y) {
  p224_fe_word_t z = 0;
  size_t i;

  for (i = 0; i < P224_FIELD_WORDS; i++)
    z |= x[i] ^ y[i];

  z = (z >> 1) | (z & 1);

  return (z - 1) >> (sizeof(z) * CHAR_BIT - 1);
}

static void
p224_fe_sqrn(p224_fe_t r, const p224_fe_t x, int rounds) {
  int i;

  /* Handle zero for the tonelli-shanks loop. */
  if (rounds == 0) {
    p224_fe_set(r, x);
    return;
  }

  p224_fe_sqr(r, x);

  for (i = 1; i < rounds; i++)
    p224_fe_sqr(r, r);
}

static void
p224_fe_pow_s(p224_fe_t r, const p224_fe_t x1) {
  /* Exponent: 2^128 - 1 */
  /* Bits: 128x1 */
  p224_fe_t t1, t2;

  /* x2 = x1^(2^1) * x1 */
  p224_fe_sqr(t1, x1);
  p224_fe_mul(t1, t1, x1);

  /* x4 = x2^(2^2) * x2 */
  p224_fe_sqrn(t2, t1, 2);
  p224_fe_mul(t2, t2, t1);

  /* x8 = x4^(2^4) * x4 */
  p224_fe_sqrn(t1, t2, 4);
  p224_fe_mul(t1, t1, t2);

  /* x16 = x8^(2^8) * x8 */
  p224_fe_sqrn(t2, t1, 8);
  p224_fe_mul(t2, t2, t1);

  /* x32 = x16^(2^16) * x16 */
  p224_fe_sqrn(t1, t2, 16);
  p224_fe_mul(t1, t1, t2);

  /* x64 = x32^(2^32) * x32 */
  p224_fe_sqrn(t2, t1, 32);
  p224_fe_mul(t2, t2, t1);

  /* x128 = x64^(2^64) * x64 */
  p224_fe_sqrn(r, t2, 64);
  p224_fe_mul(r, r, t2);
}

static void
p224_fe_pow_e(p224_fe_t r, const p224_fe_t x) {
  /* Exponent: 2^127 */
  /* Bits: 1x1 127x0 */

  /* r = x^(2^127) */
  p224_fe_sqrn(r, x, 127);
}

static void
p224_fe_invert(p224_fe_t r, const p224_fe_t x) {
  /* Exponent: p - 2 */
  /* Bits: 127x1 1x0 96x1 */
  p224_fe_t t0, t1, t2, t3, t4;

  /* x1 = x */
  p224_fe_set(t0, x);

  /* x2 = x1^(2^1) * x1 */
  p224_fe_sqr(t1, t0);
  p224_fe_mul(t1, t1, t0);

  /* x3 = x2^(2^1) * x1 */
  p224_fe_sqr(t1, t1);
  p224_fe_mul(t1, t1, t0);

  /* x6 = x3^(2^3) * x3 */
  p224_fe_sqrn(t2, t1, 3);
  p224_fe_mul(t2, t2, t1);

  /* x12 = x6^(2^6) * x6 */
  p224_fe_sqrn(t1, t2, 6);
  p224_fe_mul(t1, t1, t2);

  /* x24 = x12^(2^12) * x12 */
  p224_fe_sqrn(t3, t1, 12);
  p224_fe_mul(t3, t3, t1);

  /* x48 = x24^(2^24) * x24 */
  p224_fe_sqrn(t1, t3, 24);
  p224_fe_mul(t1, t1, t3);

  /* x96 = x48^(2^48) * x48 */
  p224_fe_sqrn(t4, t1, 48);
  p224_fe_mul(t4, t4, t1);

  /* x120 = x96^(2^24) * x24 */
  p224_fe_sqrn(r, t4, 24);
  p224_fe_mul(r, r, t3);

  /* x126 = x120^(2^6) * x6 */
  p224_fe_sqrn(r, r, 6);
  p224_fe_mul(r, r, t2);

  /* x127 = x126^(2^1) * x1 */
  p224_fe_sqr(r, r);
  p224_fe_mul(r, r, t0);

  /* r = r^(2^1) */
  p224_fe_sqr(r, r);

  /* r = r^(2^96) * x96 */
  p224_fe_sqrn(r, r, 96);
  p224_fe_mul(r, r, t4);
}

static int
p224_fe_sqrt_var(p224_fe_t r, const p224_fe_t x) {
  /* Tonelli-Shanks for P224.
   *
   * Algorithm:
   *
   *   s = 2^128 - 1 (0xffffffffffffffffffffffffffffffff)
   *   n = 11
   *   e = 2^127 (0x80000000000000000000000000000000)
   *   y = x^e mod p
   *   b = x^s mod p
   *   g = n^s mod p (0x6a0fec678598a7920c55b2d40b2d6ffbbea3d8cef3fb3632dc691b74)
   *   k = 96
   *
   *   loop:
   *     m = 0
   *     t = b
   *
   *     while t != 1:
   *       t = t^2 mod p
   *       m += 1
   *
   *     if m == 0:
   *       break
   *
   *     if m >= k:
   *       fail
   *
   *     t = g^(2^(k - m - 1)) mod p
   *     g = t^2 mod p
   *     y = y * t mod p
   *     b = b * g mod p
   *     k = m
   *
   *   ret = y
   */
  p224_fe_t y, b, g, t;
  int k, m;

  p224_fe_pow_e(y, x);
  p224_fe_pow_s(b, x);
  p224_fe_set(g, p224_g);
  p224_fe_set(r, p224_zero);

  /* Note that b happens to be the first
   * step of Euler's criterion. Squaring
   * it 95 times more gives us the Legendre
   * symbol.
   */
  p224_fe_sqrn(t, b, 95);

  /* Zero. */
  if (p224_fe_equal(t, p224_zero))
    return 1;

  /* Quadratic non-residue. */
  if (!p224_fe_equal(t, p224_one))
    return 0;

  /* Loop until we find a solution. */
  k = 96;

  for (;;) {
    m = 0;

    p224_fe_set(t, b);

    while (!p224_fe_equal(t, p224_one) && m < k) {
      p224_fe_sqr(t, t);
      m += 1;
    }

    if (m == 0)
      break;

    if (m >= k)
      return 0;

    p224_fe_sqrn(t, g, k - m - 1);
    p224_fe_sqr(g, t);
    p224_fe_mul(y, y, t);
    p224_fe_mul(b, b, g);

    k = m;
  }

  p224_fe_set(r, y);

  return 1;
}

TORSION_UNUSED static void
p224_fe_legendre(p224_fe_t r, const p224_fe_t x) {
  /* Exponent: (p - 1) / 2 */
  /* Bits: 128x1 95x0 */

  /* r = x^(2^128 - 1) */
  p224_fe_pow_s(r, x);

  /* r = r^(2^95) */
  p224_fe_sqrn(r, r, 95);
}
