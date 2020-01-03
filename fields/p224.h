#ifdef BCRYPTO_EC_64BIT
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
#define p224_fe_nonzero fiat_p224_nonzero

#ifdef BCRYPTO_EC_64BIT
static const p224_fe_t p224_zero = {0, 0, 0, 0};

/* 64 bit alignment */
static const p224_fe_t p224_one = {
  0xffffffff00000000, 0xffffffffffffffff,
  0x0000000000000000, 0x0000000000000000
};

/* 32 bit alignment as 64 bit */
static const p224_fe_t p224_one_32 = {
  0xffffffffffffffff, 0xffffffff00000000,
  0x0000000000000000, 0x0000000000000000
};

/* 11^(2^128 - 1) mod p */
/* mont: 0xa31b1da46d3e2af0dd915e4b7869be5d866c223b174131b85ee27c6c */
/* 64 bit alignment */
static const p224_fe_t p224_g = {
  0x174131b85ee27c6c, 0x7869be5d866c223b,
  0x6d3e2af0dd915e4b, 0x00000000a31b1da4
};

/* 11^(2^128 - 1) mod p */
/* mont: 0xa11d8394a31b1da46d3e2af0dd915e4ad74c3ac9866c223b174131b9 */
/* 32 bit alignment as 64 bit */
static const p224_fe_t p224_g_32 = {
  0x174131b9866c223b, 0xd74c3ac9dd915e4a,
  0x6d3e2af0a31b1da4, 0x00000000a11d8394
};
#else
static const p224_fe_t p224_zero = {0, 0, 0, 0, 0, 0, 0};

/* 32 bit alignment */
static const p224_fe_t p224_one = {
  0xffffffff, 0xffffffff, 0xffffffff, 0x00000000,
  0x00000000, 0x00000000, 0x00000000
};

/* 11^(2^128 - 1) mod p */
/* mont: 0xa11d8394a31b1da46d3e2af0dd915e4ad74c3ac9866c223b174131b9 */
/* 32 bit alignment */
static const p224_fe_t p224_g = {
  0x174131b9, 0x866c223b, 0xd74c3ac9, 0xdd915e4a,
  0x6d3e2af0, 0xa31b1da4, 0xa11d8394
};
#endif

static void
p224_fe_set(p224_fe_t out, const p224_fe_t in) {
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
  out[3] = in[3];
#if P224_FIELD_WORDS == 7
  out[4] = in[4];
  out[5] = in[5];
  out[6] = in[6];
#endif
}

static int
p224_fe_equal(const p224_fe_t a, const p224_fe_t b) {
  p224_fe_t c;
  p224_fe_word_t ret;

  p224_fe_sub(c, a, b);
  fiat_p224_nonzero(&ret, c);

  return ret == 0;
}

static void
p224_fe_sqrn(p224_fe_t out, const p224_fe_t in, int rounds) {
  int i;

  p224_fe_set(out, in);

  for (i = 0; i < rounds; i++)
    p224_fe_sqr(out, out);
}

/* https://github.com/openssl/openssl/blob/master/crypto/ec/ecp_nistp224.c#L701 */
/* TODO: optimize */
static void
p224_fe_invert(p224_fe_t out, const p224_fe_t in) {
  p224_fe_t ftmp, ftmp2, ftmp3, ftmp4;
  p224_fe_t tmp;
  unsigned int i;

  p224_fe_sqr(tmp, in);
  p224_fe_set(ftmp, tmp); /* 2 */
  p224_fe_mul(tmp, in, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^2 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^3 - 2 */
  p224_fe_mul(tmp, in, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^3 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp2, tmp); /* 2^4 - 2 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp2, tmp); /* 2^5 - 4 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp2, tmp); /* 2^6 - 8 */
  p224_fe_mul(tmp, ftmp2, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^6 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp2, tmp); /* 2^7 - 2 */

  for (i = 0; i < 5; ++i) { /* 2^12 - 2^6 */
    p224_fe_sqr(tmp, ftmp2);
    p224_fe_set(ftmp2, tmp);
  }

  p224_fe_mul(tmp, ftmp2, ftmp);
  p224_fe_set(ftmp2, tmp); /* 2^12 - 1 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp3, tmp); /* 2^13 - 2 */

  for (i = 0; i < 11; ++i) { /* 2^24 - 2^12 */
    p224_fe_sqr(tmp, ftmp3);
    p224_fe_set(ftmp3, tmp);
  }

  p224_fe_mul(tmp, ftmp3, ftmp2);
  p224_fe_set(ftmp2, tmp); /* 2^24 - 1 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp3, tmp); /* 2^25 - 2 */

  for (i = 0; i < 23; ++i) { /* 2^48 - 2^24 */
    p224_fe_sqr(tmp, ftmp3);
    p224_fe_set(ftmp3, tmp);
  }

  p224_fe_mul(tmp, ftmp3, ftmp2);
  p224_fe_set(ftmp3, tmp); /* 2^48 - 1 */
  p224_fe_sqr(tmp, ftmp3);
  p224_fe_set(ftmp4, tmp); /* 2^49 - 2 */

  for (i = 0; i < 47; ++i) { /* 2^96 - 2^48 */
    p224_fe_sqr(tmp, ftmp4);
    p224_fe_set(ftmp4, tmp);
  }

  p224_fe_mul(tmp, ftmp3, ftmp4);
  p224_fe_set(ftmp3, tmp); /* 2^96 - 1 */
  p224_fe_sqr(tmp, ftmp3);
  p224_fe_set(ftmp4, tmp); /* 2^97 - 2 */

  for (i = 0; i < 23; ++i) { /* 2^120 - 2^24 */
    p224_fe_sqr(tmp, ftmp4);
    p224_fe_set(ftmp4, tmp);
  }

  p224_fe_mul(tmp, ftmp2, ftmp4);
  p224_fe_set(ftmp2, tmp); /* 2^120 - 1 */

  for (i = 0; i < 6; ++i) { /* 2^126 - 2^6 */
    p224_fe_sqr(tmp, ftmp2);
    p224_fe_set(ftmp2, tmp);
  }

  p224_fe_mul(tmp, ftmp2, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^126 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^127 - 2 */
  p224_fe_mul(tmp, ftmp, in);
  p224_fe_set(ftmp, tmp); /* 2^127 - 1 */

  for (i = 0; i < 97; ++i) { /* 2^224 - 2^97 */
    p224_fe_sqr(tmp, ftmp);
    p224_fe_set(ftmp, tmp);
  }

  p224_fe_mul(tmp, ftmp, ftmp3);
  p224_fe_set(out, tmp); /* 2^224 - 2^96 - 1 */
}

static int
p224_fe_legendre(const p224_fe_t in) {
  p224_fe_t t;
  int i, a, b, c;

  p224_fe_set(t, in);

  /* TODO: optimize */
  for (i = 1; i < 128; i++) {
    p224_fe_sqr(t, t);
    p224_fe_mul(t, t, in);
  }

  p224_fe_sqrn(t, t, 95);

#if 0
  /* (p - 1) / 2 = 0x7fffffffffffffffffffffffffffffff800000000000000000000000 */
  /* 128 1s, 95 0s */
  static const uint32_t p224_legendre[7] = {
    0x00000000, 0x00000000, 0x80000000, 0xffffffff,
    0xffffffff, 0xffffffff, 0x7fffffff
  };

  uint32_t bit;

  p224_fe_set(t, p224_one);

  for (i = 224 - 1; i >= 0; i--) {
    bit = p224_legendre[i >> 5] >> (i & 31);

    p224_fe_sqr(t, t);

    if (bit & 1)
      p224_fe_mul(t, t, in);
  }
#endif

  a = p224_fe_equal(t, p224_zero);
  b = p224_fe_equal(t, p224_one);
  c = (a ^ 1) & (b ^ 1);

  assert((a | b | c) != 0);
  assert(a + b + c == 1);

  return b - c;
}

static void
p224_fe_pow_s(p224_fe_t out, const p224_fe_t in) {
  /* Compute x^(2^128 - 1) */
  p224_fe_t t;
  int i;

  p224_fe_set(t, in);

  for (i = 1; i < 128; i++) {
    p224_fe_sqr(t, t);
    p224_fe_mul(t, t, in);
  }

  p224_fe_set(out, t);
}

static void
p224_fe_pow_e(p224_fe_t out, const p224_fe_t in) {
  /* Compute x^(2^127) */
  p224_fe_sqrn(out, in, 127);
}

#if 0
const unsigned char p224_prime[28] = {
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x01
};

static int
p224_fe_jacobi(const p224_fe_t in) {
  unsigned char bytes[28];
  mp_limb_t xp[7];
  mp_limb_t yp[7];
  mpz_t x, y;

  fiat_p224_to_bytes(bytes, in);

  mpn_import_le(xp, 7, bytes, 28);
  mpn_import_be(yp, 7, p224_prime, 28);

  mpz_roinit_n(x, xp, 7);
  mpz_roinit_n(y, yp, 7);

  return mpz_jacobi(x, y);
}
#endif

static int
p224_fe_sqrt_var(p224_fe_t out, const p224_fe_t in) {
  /* Tonelli-Shanks for P224.
   *
   * Algorithm:
   *
   *   z = 96
   *   s = 2^128 - 1 (0xffffffffffffffffffffffffffffffff)
   *   n = 11
   *   e = 2^127 (0x80000000000000000000000000000000)
   *   y = x^e mod p
   *   b = x^s mod p
   *   g = n^s mod p (0x6a0fec678598a7920c55b2d40b2d6ffbbea3d8cef3fb3632dc691b74)
   *   k = z
   *
   *   loop:
   *     m = 0
   *     t = b
   *
   *     while t != 1:
   *       t = t^2 mod p
   *       m += 1;
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
  p224_fe_t y, b, g, t, l;
  int k, m;

#if 0
  switch (p224_fe_legendre(in)) {
    case 0:
      p224_fe_set(out, p224_zero);
      return 1;
    case -1:
      return 0;
  }
#endif

  p224_fe_pow_e(y, in);
  p224_fe_pow_s(b, in);
  p224_fe_set(g, p224_g);

  /* Note that b happens to be the first
   * step of Euler's criterion. Squaring
   * it 95 times more gives us the Legendre
   * symbol.
   */
  p224_fe_sqrn(l, b, 95);

  /* Zero. */
  if (p224_fe_equal(l, p224_zero)) {
    p224_fe_set(out, p224_zero);
    return 1;
  }

  /* Quadratic non-residue. */
  if (!p224_fe_equal(l, p224_one))
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

  p224_fe_set(out, y);

  return 1;
}
