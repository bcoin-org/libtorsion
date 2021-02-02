/*!
 * mpi_internal.h - mpi internal tests for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "data/jacobi_vectors.h"
#include "data/kronecker_vectors.h"
#include "data/mpz_vectors.h"
#include "data/prime_vectors.h"

#undef ASSERT
#define ASSERT(expr) ASSERT_ALWAYS(expr)

/*
 * Helpers
 */

static mp_limb_t
mp_random_limb(mp_rng_f *rng, void *arg) {
  mp_limb_t z;

  rng(&z, sizeof(z), arg);

  return z;
}

static mp_limb_t
mp_random_limb_nz(mp_rng_f *rng, void *arg) {
  mp_limb_t z;

  do {
    z = mp_random_limb(rng, arg);
  } while (z == 0);

  return z;
}

static mp_long_t
mp_random_long(mp_rng_f *rng, void *arg) {
  mp_limb_t z = mp_random_limb(rng, arg);
  mp_long_t s = z >> (MP_LIMB_BITS - 1);
  mp_long_t w = z & (MP_LIMB_HI - 1);

  return (-s | 1) * w;
}

static mp_long_t
mp_random_long_nz(mp_rng_f *rng, void *arg) {
  mp_long_t z;

  do {
    z = mp_random_long(rng, arg);
  } while (z == 0);

  return z;
}

static void
mpn_random_nz(mp_limb_t *zp, mp_size_t zn, mp_rng_f *rng, void *arg) {
  CHECK(zn != 0);

  do {
    mpn_random(zp, zn, rng, arg);
  } while (mpn_zero_p(zp, zn));
}

static void
mpz_random_nz(mpz_t z, mp_bits_t bits, mp_rng_f *rng, void *arg) {
  CHECK(bits != 0);

  do {
    mpz_urandomb(z, bits, rng, arg);
  } while (mpz_sgn(z) == 0);
}

static int
fake_puts(const char *s) {
  CHECK(s != NULL);
  return 0;
}

typedef struct arc4_s {
  uint8_t s[256];
  uint8_t i;
  uint8_t j;
} arc4_t;

static const arc4_t arc4_initial = {
  {
    /* Key: e7bf7a0072c8fa81ebf2af95deabe56d85842179ea0a0da70b7e0d4f31b70215 */
    0xa1, 0xa7, 0x23, 0x26, 0x36, 0x20, 0x05, 0x5a, 0xf8, 0xdf, 0x81, 0x79,
    0x22, 0xd4, 0xcd, 0x49, 0xde, 0xb0, 0xa6, 0x3b, 0x27, 0x75, 0xaa, 0x95,
    0x8a, 0xb2, 0x48, 0x06, 0xc8, 0x04, 0x4c, 0xfb, 0xf7, 0xd7, 0x62, 0x68,
    0x89, 0xe2, 0xa2, 0x87, 0x9a, 0xb5, 0xb7, 0x4e, 0x58, 0x92, 0xfe, 0x84,
    0x91, 0x46, 0x7a, 0xb9, 0x6b, 0xcf, 0xca, 0x3f, 0xad, 0x71, 0xb8, 0x24,
    0xaf, 0x1b, 0xf6, 0xa4, 0xa5, 0x2e, 0xff, 0x40, 0x00, 0xbf, 0xb1, 0x3c,
    0x73, 0x7b, 0x30, 0x45, 0xc3, 0x2c, 0x33, 0xd2, 0x7d, 0x52, 0xc4, 0x55,
    0x93, 0x29, 0xab, 0x4d, 0x72, 0xf5, 0xe6, 0x64, 0xc2, 0x85, 0x1c, 0x6c,
    0x2b, 0xf0, 0x5d, 0x0c, 0x65, 0x76, 0xe8, 0xd0, 0x5f, 0xb6, 0x44, 0x6a,
    0x86, 0x6e, 0xe9, 0xc5, 0x0b, 0xbe, 0xae, 0x67, 0xee, 0x4f, 0x0a, 0x78,
    0x94, 0xb3, 0x2a, 0xa3, 0x99, 0x07, 0x21, 0xf3, 0x8b, 0x28, 0x74, 0xcb,
    0x8f, 0xec, 0x25, 0x4a, 0x83, 0x02, 0x56, 0x51, 0x77, 0xa9, 0x3a, 0xb4,
    0xc9, 0xd3, 0x17, 0xfa, 0x1e, 0x14, 0x9d, 0x08, 0x88, 0x19, 0x98, 0x4b,
    0x32, 0x35, 0xc6, 0x7c, 0x8c, 0xd1, 0xed, 0xd9, 0x13, 0x5e, 0x18, 0xe3,
    0x01, 0xbd, 0x57, 0x53, 0xe0, 0x50, 0x09, 0xdc, 0x70, 0x8d, 0xea, 0xce,
    0xc1, 0x39, 0xd6, 0xc7, 0xfc, 0x16, 0xbb, 0x9f, 0xe1, 0x0d, 0xa0, 0xf4,
    0x63, 0x03, 0x6f, 0x41, 0x3e, 0x69, 0xf2, 0x3d, 0xd8, 0x80, 0x37, 0x12,
    0x34, 0x97, 0x38, 0xef, 0xdb, 0x0e, 0x7f, 0xd5, 0x9b, 0x54, 0x0f, 0x2f,
    0xbc, 0x1f, 0x61, 0xba, 0xac, 0x6d, 0x7e, 0x59, 0xcc, 0x47, 0x10, 0xf1,
    0xe7, 0x2d, 0x90, 0x96, 0xda, 0x82, 0x1a, 0x42, 0x8e, 0x60, 0xf9, 0x5b,
    0xc0, 0x11, 0x9c, 0xa8, 0x66, 0xeb, 0x15, 0x5c, 0xe4, 0xdd, 0x9e, 0x31,
    0xe5, 0x43, 0xfd, 0x1d
  },
  0,
  0
};

static void
arc4_rng(void *out, size_t size, void *arg) {
  arc4_t *ctx = arg;
  uint8_t *dst = out;
  uint8_t *s = ctx->s;
  uint8_t i = ctx->i;
  uint8_t j = ctx->j;
  uint8_t x, y;
  size_t k;

  for (k = 0; k < size; k++) {
    i = (i + 1) & 0xff;
    x = s[i];

    j = (j + x) & 0xff;
    y = s[j];

    s[i] = y;
    s[j] = x;

    dst[k] = s[(x + y) & 0xff];
  }

  ctx->i = i;
  ctx->j = j;
}

static TORSION_INLINE mp_bits_t
mp_popcount_simple(mp_limb_t x) {
  /* Kernighan's population count[1].
   *
   * Note that clang can optimize this
   * to a single popcnt instruction[2].
   *
   * Unless you're reading this many
   * years in the future, you probably
   * have to pass -march=native in order
   * for popcnt to be emitted on x86.
   *
   * See the LLVM code for more info[3].
   *
   * [1] https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan
   * [2] https://lemire.me/blog/2016/05/23/the-surprising-cleverness-of-modern-compilers/
   * [3] https://github.com/llvm/llvm-project/blob/51879a5
   *     /llvm/lib/Transforms/Scalar/LoopIdiomRecognize.cpp#L1285
   */
  mp_bits_t z = 0;

  while (x != 0) {
    x &= (x - 1);
    z += 1;
  }

  return z;
}

static TORSION_INLINE mp_bits_t
mp_clz_simple(mp_limb_t x) {
  /* Naive leading zero count.
   *
   * Note that clang will optimize this
   * to a single lzcnt/bsr instruction[1].
   *
   * [1] https://github.com/llvm/llvm-project/blob/51879a5
   *     /llvm/lib/Transforms/Scalar/LoopIdiomRecognize.cpp#L1417
   */
  mp_bits_t z = 0;

  while (x != 0) {
    x >>= 1;
    z += 1;
  }

  return MP_LIMB_BITS - z;
}

static TORSION_INLINE mp_bits_t
mp_ctz_simple(mp_limb_t x) {
  /* Naive trailing zero count.
   *
   * Note that clang will optimize this
   * to a single tzcnt/bsf instruction[1].
   *
   * [1] https://github.com/llvm/llvm-project/blob/51879a5
   *     /llvm/lib/Transforms/Scalar/LoopIdiomRecognize.cpp#L1417
   */
  mp_bits_t z = 0;

  while (x != 0) {
    x <<= 1;
    z += 1;
  }

  return MP_LIMB_BITS - z;
}

static TORSION_INLINE mp_bits_t
mp_bitlen_simple(mp_limb_t x) {
  /* Naive bit length calculation.
   *
   * Note that clang will optimize this
   * to a lzcnt/bsr instruction plus some
   * additional arithmetic[1].
   *
   * [1] https://github.com/llvm/llvm-project/blob/51879a5
   *     /llvm/lib/Transforms/Scalar/LoopIdiomRecognize.cpp#L1417
   */
  mp_bits_t z = 0;

  while (x != 0) {
    x >>= 1;
    z += 1;
  }

  return z;
}

/*
 * MPN/MPZ Functions
 */

static int
mpn_cmp3(const mp_limb_t *xp, mp_size_t xn,
         const mp_limb_t *yp, mp_size_t yn) {
  ASSERT(xn == 0 || xp[xn - 1] != 0);
  ASSERT(yn == 0 || yp[yn - 1] != 0);

  if (xn != yn)
    return xn < yn ? -1 : 1;

  return mpn_cmp(xp, yp, xn);
}

static void
mpz_mul_simple(mpz_t z, const mpz_t x, const mpz_t y) {
  mp_size_t xn, yn, zn;
  mp_limb_t *tp;

  if (x->size == 0 || y->size == 0) {
    z->size = 0;
    return;
  }

  xn = MP_ABS(x->size);
  yn = MP_ABS(y->size);
  zn = xn + yn;

  mpz_grow(z, zn);

  tp = mp_alloc_limbs(zn);

  mpn_mul(tp, x->limbs, xn, y->limbs, yn);

  zn -= (tp[zn - 1] == 0);

  mpn_copyi(z->limbs, tp, zn);

  mp_free_limbs(tp);

  z->size = (x->size ^ y->size) < 0 ? -zn : zn;
}

static mp_limb_t
mpn_divmod_1_simple(mp_limb_t *qp,
                    const mp_limb_t *np,
                    mp_size_t nn,
                    mp_limb_t d) {
  mp_bits_t s = mp_clz(d);
  mp_limb_t q, n1, n0;
  mp_limb_t r = 0;
  mp_size_t j;

  d <<= s;

  for (j = nn - 1; j >= 0; j--) {
    n1 = r;
    n0 = np[j];

    if (s != 0) {
      n1 = (n1 << s) | (n0 >> (MP_LIMB_BITS - s));
      n0 <<= s;
    }

    mp_div(&q, &r, n1, n0, d);

    r >>= s;

    if (qp != NULL)
      qp[j] = q;
  }

  return r;
}

static void
mpn_divmod_simple(mp_limb_t *qp, mp_limb_t *rp,
                  const mp_limb_t *np, mp_size_t nn,
                  const mp_limb_t *dp, mp_size_t dn) {
  mp_bits_t shift = mp_clz(dp[dn - 1]);
  mp_limb_t *up = mp_alloc_vla(nn + 1);
  mp_limb_t *vp = mp_alloc_vla(dn);
  mp_limb_t qhat, rhat, c;
  mp_size_t j;

  if (shift != 0) {
    mpn_lshift(vp, dp, dn, shift);
    up[nn] = mpn_lshift(up, np, nn, shift);
  } else {
    mpn_copyi(vp, dp, dn);
    mpn_copyi(up, np, nn);
    up[nn] = 0;
  }

  for (j = nn - dn; j >= 0; j--) {
    qhat = MP_LIMB_MAX;

    if (up[j + dn] != vp[dn - 1]) {
      mp_div(&qhat, &rhat, up[j + dn], up[j + dn - 1], vp[dn - 1]);

      while (mp_mul_gt_2(qhat, vp[dn - 2], rhat, up[j + dn - 2])) {
        qhat -= 1;
        rhat += vp[dn - 1];

        if (rhat < vp[dn - 1])
          break;
      }
    }

    c = mpn_submul_1(up + j, vp, dn, qhat);

    mp_sub(up[j + dn], c, up[j + dn], c);

    if (c != 0) {
      up[j + dn] += mpn_add_n(up + j, up + j, vp, dn);

      qhat -= 1;
    }

    ASSERT(up[j + dn] == 0);

    if (qp != NULL)
      qp[j] = qhat;
  }

  if (rp != NULL) {
    if (shift != 0)
      mpn_rshift(rp, up, dn, shift);
    else
      mpn_copyi(rp, up, dn);
  }

  mp_free_vla(up, nn + 1);
  mp_free_vla(vp, dn);
}

static mp_limb_t
mp_inv_mod_simple(mp_limb_t d) {
  mpz_t z, x, y;
  mp_limb_t m;

  mpz_init(z);
  mpz_init(x);
  mpz_init(y);

  mpz_set_ui(x, d);
  mpz_setbit(y, MP_LIMB_BITS);

  ASSERT(mpz_invert(z, x, y));
  ASSERT(mpz_fits_ui_p(z));

  m = mpz_get_ui(z);

  mpz_clear(z);
  mpz_clear(x);
  mpz_clear(y);

  return m;
}

static void
mp_div_2by1_simple(mp_limb_t *q, mp_limb_t *r,
                   mp_limb_t n1, mp_limb_t n0,
                   mp_limb_t d) {
  mp_limb_t np[2], qp[2];
  mp_limb_t r0;

  ASSERT(d != 0);

  np[0] = n0;
  np[1] = n1;

  r0 = mpn_divmod_1_simple(qp, np, 2, d);

  ASSERT(qp[1] == 0);

  *q = qp[0];
  *r = r0;
}

static void
mp_div_3by2_simple(mp_limb_t *q, mp_limb_t *r1, mp_limb_t *r0,
                   mp_limb_t n2, mp_limb_t n1, mp_limb_t n0,
                   mp_limb_t d1, mp_limb_t d0) {
  mp_limb_t np[3], dp[2], qp[3], rp[2];

  ASSERT((d0 | d1) != 0);

  np[0] = n0;
  np[1] = n1;
  np[2] = n2;

  dp[0] = d0;
  dp[1] = d1;

  if (d1 == 0) {
    rp[0] = mpn_divmod_1_simple(qp, np, 3, d0);

    ASSERT(qp[1] == 0);
    ASSERT(qp[2] == 0);

    *q = qp[0];
    *r0 = rp[0];
    *r1 = 0;
  } else {
    mpn_divmod_simple(qp, rp, np, 3, dp, 2);

    ASSERT(qp[1] == 0);

    *q = qp[0];
    *r0 = rp[0];
    *r1 = rp[1];
  }
}

static void
mpz_gcd_simple(mpz_t g, const mpz_t x, const mpz_t y) {
  mpz_t u, v, t;

  mpz_init(u);
  mpz_init(v);
  mpz_init(t);

  mpz_abs(u, x);
  mpz_abs(v, y);

  while (mpz_sgn(v) != 0) {
    mpz_rem(t, u, v);
    mpz_swap(u, v);
    mpz_swap(v, t);
  }

  mpz_swap(g, u);

  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(t);
}

static void
mpz_lcm_simple(mpz_t z, const mpz_t x, const mpz_t y) {
  mpz_t g, l;

  if (mpz_sgn(x) == 0 || mpz_sgn(y) == 0) {
    mpz_set_ui(z, 0);
    return;
  }

  mpz_init(g);
  mpz_init(l);

  mpz_gcd_simple(g, x, y);
  mpz_divexact(l, x, g);
  mpz_mul(z, l, y);
  mpz_abs(z, z);

  mpz_clear(g);
  mpz_clear(l);
}

static void
mpz_pow5(mpz_t z, const mpz_t x) {
  mpz_t t;

  mpz_init(t);
  mpz_set(t, x);

  mpz_set(z, x);
  mpz_sqr(z, z);
  mpz_sqr(z, z);
  mpz_mul(z, z, t);

  mpz_clear(t);
}

static void
mpz_pow_simple_ui(mpz_t z, const mpz_t x, mp_limb_t y) {
  mpz_t r, u;

  mpz_init(r);
  mpz_init(u);

  mpz_set_ui(r, 1);
  mpz_set(u, x);

  while (y > 0) {
    if (y & 1)
      mpz_mul_simple(r, r, u);

    mpz_mul_simple(u, u, u);

    y >>= 1;
  }

  mpz_swap(z, r);

  mpz_clear(r);
  mpz_clear(u);
}

static void
mpz_powm_simple(mpz_t z, const mpz_t x, const mpz_t y, const mpz_t m) {
  mp_bits_t i, bits;
  mpz_t r, u, v, t;

  mpz_init(r);
  mpz_init(u);
  mpz_init(v);
  mpz_init(t);

  mpz_set_ui(r, 1);

  mpz_mod(r, r, m);
  mpz_mod(u, x, m);
  mpz_abs(v, y);

  if (mpz_sgn(y) < 0)
    ASSERT(mpz_invert(u, u, m));

  bits = mpz_bitlen(v);

  for (i = 0; i < bits; i++) {
    if (mpz_tstbit(v, i)) {
      mpz_mul(t, r, u);
      mpz_mod(r, t, m);
    }

    mpz_sqr(t, u);
    mpz_mod(u, t, m);
  }

  mpz_swap(z, r);

  mpz_clear(r);
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(t);
}

static void
mpz_powm_simple_ui(mpz_t z, const mpz_t x, mp_limb_t y, const mpz_t m) {
  mpz_t r, u, t;

  mpz_init(r);
  mpz_init(u);
  mpz_init(t);

  mpz_set_ui(r, 1);

  mpz_mod(r, r, m);
  mpz_mod(u, x, m);

  while (y > 0) {
    if (y & 1) {
      mpz_mul(t, r, u);
      mpz_mod(r, t, m);
    }

    mpz_sqr(t, u);
    mpz_mod(u, t, m);

    y >>= 1;
  }

  mpz_swap(z, r);

  mpz_clear(r);
  mpz_clear(u);
  mpz_clear(t);
}

static void
mpn_powm_simple(mp_limb_t *zp,
                const mp_limb_t *xp, mp_size_t xn,
                const mp_limb_t *yp, mp_size_t yn,
                const mp_limb_t *mp, mp_size_t mn) {
  mpz_t z, x, y, m;

  mpz_init(z);

  mpz_roinit_n(x, xp, xn);
  mpz_roinit_n(y, yp, yn);
  mpz_roinit_n(m, mp, mn);

  mpz_powm_simple(z, x, y, m);

  mpn_copyi(zp, z->limbs, z->size);
  mpn_zero(zp + z->size, mn - z->size);

  mpz_clear(z);
}

static mp_size_t
mpn_gcd_simple(mp_limb_t *zp, const mp_limb_t *xp, mp_size_t xn,
                              const mp_limb_t *yp, mp_size_t yn) {
  mpz_t x, y, z;
  mp_size_t zn;

  CHECK(xn >= yn);
  CHECK(yn > 0);
  CHECK(yp[yn - 1] != 0);

  mpz_roinit_n(x, xp, xn);
  mpz_roinit_n(y, yp, yn);
  mpz_init(z);

  mpz_gcd_simple(z, x, y);

  zn = MP_ABS(z->size);

  CHECK(zn <= yn);

  mpn_copyi(zp, z->limbs, zn);
  mpn_zero(zp + zn, yn - zn);

  mpz_clear(z);

  return zn;
}

TORSION_UNUSED static int
mpz_cmp_str(const mpz_t x, const char *str) {
  mpz_t y;
  int cmp;

  mpz_init(y);

  CHECK(mpz_set_str(y, str, 0));

  cmp = mpz_cmp(x, y);

  mpz_clear(y);

  return cmp;
}

/*
 * To be implemented...
 */

static mp_size_t
mpn_sqrtrem(mp_limb_t *zp, mp_limb_t *rp, const mp_limb_t *xp, mp_size_t xn) {
  mp_size_t zn, rn;
  mpz_t x, z, r;

  CHECK(zp != xp);
  CHECK(xn > 0);
  CHECK(xp[xn - 1] != 0);

  mpz_roinit_n(x, xp, xn);
  mpz_init(z);
  mpz_init(r);

  mpz_sqrtrem(z, r, x);

  zn = MP_ABS(z->size);
  rn = MP_ABS(r->size);

  CHECK(zn <= (xn + 1) / 2);
  CHECK(rn <= xn);

  mpn_copyi(zp, z->limbs, zn);
  mpn_zero(zp + zn, (xn + 1) / 2 - zn);

  if (rp != NULL) {
    mpn_copyi(rp, r->limbs, rn);
    mpn_zero(rp + rn, xn - rn);
  }

  mpz_clear(z);
  mpz_clear(r);

  return rn;
}

static int
mpn_perfect_square_p(const mp_limb_t *xp, mp_size_t xn) {
  mpz_t x;

  CHECK(xn > 0);
  CHECK(xp[xn - 1] != 0);

  mpz_roinit_n(x, xp, xn);

  return mpz_perfect_square_p(x);
}

static mp_size_t
mpn_gcdext(mp_limb_t *gp,
           mp_limb_t *sp, mp_size_t *sn,
           mp_limb_t *xp, mp_size_t xn,
           mp_limb_t *yp, mp_size_t yn) {
  mpz_t x, y, g, s;
  mp_size_t gn, bn;

  CHECK(xn >= yn);
  CHECK(yn > 0);
  CHECK(yp[yn - 1] != 0);

  mpz_roinit_n(x, xp, xn);
  mpz_roinit_n(y, yp, yn);
  mpz_init(g);
  mpz_init(s);

  mpz_gcdext(g, s, NULL, x, y);
  mpz_rem(s, s, y);

  gn = MP_ABS(g->size);
  bn = MP_ABS(s->size);

  CHECK(gn <= yn);
  CHECK(bn <= yn + 1);

  mpn_copyi(gp, g->limbs, gn);
  mpn_zero(gp + gn, yn - gn);

  mpn_copyi(sp, s->limbs, bn);
  mpn_zero(sp + bn, yn + 1 - bn);

  *sn = s->size;

  mpz_clear(g);
  mpz_clear(s);

  return gn;
}

/*
 * P192
 */

#define MP_P192_LIMBS ((192 + MP_LIMB_BITS - 1) / MP_LIMB_BITS)
#define MP_P192_SHIFT (2 * MP_P192_LIMBS)

static void
mpn_p192_mod(mp_limb_t *zp) {
  /* z = 2^192 - 2^64 - 1 */
  mp_limb_t xp[MP_P192_LIMBS + 1];
  mp_limb_t yp[MP_P192_LIMBS + 1];

  mpn_zero(xp, MP_P192_LIMBS + 1);
  mpn_zero(yp, MP_P192_LIMBS + 1);

  mpn_setbit(xp, 192);
  mpn_setbit(yp, 64);

  mpn_sub_n(xp, xp, yp, MP_P192_LIMBS + 1);
  mpn_sub_1(zp, xp, MP_P192_LIMBS, 1);
}

static void
mpn_p192_exp(mp_limb_t *zp) {
  /* e = (p + 1) / 4 */
  mp_limb_t mp[MP_P192_LIMBS + 1];
  mp_limb_t ep[MP_P192_LIMBS + 1];

  mpn_zero(mp, MP_P192_LIMBS + 1);
  mpn_zero(ep, MP_P192_LIMBS + 1);

  mpn_p192_mod(mp);

  mpn_add_1(ep, mp, MP_P192_LIMBS + 1, 1);
  mpn_rshift(ep, ep, MP_P192_LIMBS + 1, 2);
  mpn_copyi(zp, ep, MP_P192_LIMBS);
}

static void
mpz_p192_mod(mpz_t z) {
  mpz_grow(z, MP_P192_LIMBS);

  mpn_p192_mod(z->limbs);

  z->size = MP_P192_LIMBS;
}

static void
mpz_p192_exp(mpz_t z) {
  mpz_grow(z, MP_P192_LIMBS);

  mpn_p192_exp(z->limbs);

  z->size = MP_P192_LIMBS;
}

/*
 * P224
 */

#define MP_P224_LIMBS ((224 + MP_LIMB_BITS - 1) / MP_LIMB_BITS)

TORSION_UNUSED static void
mpn_p224_mod(mp_limb_t *zp) {
  /* z = 2^224 - 2^96 + 1 */
  mp_limb_t xp[MP_P224_LIMBS + 1];
  mp_limb_t yp[MP_P224_LIMBS + 1];

  mpn_zero(xp, MP_P224_LIMBS + 1);
  mpn_zero(yp, MP_P224_LIMBS + 1);

  mpn_setbit(xp, 224);
  mpn_setbit(yp, 96);

  mpn_sub_n(xp, xp, yp, MP_P224_LIMBS + 1);
  mpn_add_1(zp, xp, MP_P224_LIMBS, 1);
}

/*
 * K256
 */

#define MP_K256_LIMBS ((256 + MP_LIMB_BITS - 1) / MP_LIMB_BITS)

TORSION_UNUSED static void
mpn_k256_mod(mp_limb_t *zp) {
  /* z = 2^256 - 2^32 - 977 */
  mp_limb_t xp[MP_K256_LIMBS + 1];
  mp_limb_t yp[MP_K256_LIMBS + 1];

  mpn_zero(xp, MP_K256_LIMBS + 1);
  mpn_zero(yp, MP_K256_LIMBS + 1);

  mpn_setbit(xp, 256);
  mpn_setbit(yp, 32);

  mpn_sub_n(xp, xp, yp, MP_K256_LIMBS + 1);
  mpn_sub_1(zp, xp, MP_K256_LIMBS, 977);
}

/*
 * P521
 */

#define MP_P521_LIMBS ((521 + MP_LIMB_BITS - 1) / MP_LIMB_BITS)

TORSION_UNUSED static void
mpn_p521_mod(mp_limb_t *zp) {
  /* z = 2^521 - 1 */
  mpn_zero(zp, MP_P521_LIMBS);
  mpn_setbit(zp, 521);
  mpn_sub_1(zp, zp, MP_P521_LIMBS, 1);
}

/*
 * P25519
 */

#define MP_P25519_LIMBS ((255 + MP_LIMB_BITS - 1) / MP_LIMB_BITS)

TORSION_UNUSED static void
mpn_p25519_mod(mp_limb_t *zp) {
  /* z = 2^255 - 19 */
  mpn_zero(zp, MP_P25519_LIMBS);
  mpn_setbit(zp, 255);
  mpn_sub_1(zp, zp, MP_P25519_LIMBS, 19);
}

/*
 * MP
 */

static void
test_mp_helpers(mp_rng_f *rng, void *arg) {
  mp_limb_t *xp;
  mp_limb_t x;
  int i;

  printf("  - MP helpers.\n");

  {
    xp = mp_alloc_limbs(2);

    xp[0] = 0;
    xp[1] = 0;

    mp_free_limbs(xp);
  }

  {
    xp = mp_alloc_vla(2);

    xp[0] = 0;
    xp[1] = 0;

    mp_free_vla(xp, 2);
  }

  {
    ASSERT(mp_popcount(0) == 0);
    ASSERT(mp_popcount(0x01010101) == 4);
    ASSERT(mp_popcount(0x0000ffff) == 16);
    ASSERT(mp_popcount(0xffffffff) == 32);
    ASSERT(mp_popcount(MP_LIMB_HI) == 1);
    ASSERT(mp_popcount(MP_LIMB_HI - 1) == MP_LIMB_BITS - 1);
    ASSERT(mp_popcount(MP_LIMB_MAX - 1) == MP_LIMB_BITS - 1);
    ASSERT(mp_popcount(MP_LIMB_MAX) == MP_LIMB_BITS);

    for (i = 0; i < MP_LIMB_BITS; i++)
      ASSERT(mp_popcount(MP_LIMB_MAX >> i) == MP_LIMB_BITS - i);

    for (i = 0; i < 100; i++) {
      x = mp_random_limb(rng, arg);

      ASSERT(mp_popcount(x) == mp_popcount_simple(x));
    }
  }

  {
    ASSERT(mp_clz(0) == MP_LIMB_BITS);
    ASSERT(mp_clz(1) == MP_LIMB_BITS - 1);
    ASSERT(mp_clz(2) == MP_LIMB_BITS - 2);
    ASSERT(mp_clz(MP_LIMB_MAX) == 0);

    for (i = 0; i < MP_LIMB_BITS; i++)
      ASSERT(mp_clz(MP_LIMB_HI >> i) == i);

    for (i = 0; i < 100; i++) {
      x = mp_random_limb(rng, arg);

      ASSERT(mp_clz(x) == mp_clz_simple(x));
    }
  }

  {
    ASSERT(mp_ctz(0) == MP_LIMB_BITS);
    ASSERT(mp_ctz(1) == 0);
    ASSERT(mp_ctz(2) == 1);
    ASSERT(mp_ctz(MP_LIMB_HI) == MP_LIMB_BITS - 1);
    ASSERT(mp_ctz(MP_LIMB_MAX) == 0);

    for (i = 0; i < MP_LIMB_BITS; i++)
      ASSERT(mp_ctz(MP_LIMB_HI >> i) == MP_LIMB_BITS - 1 - i);

    for (i = 0; i < 100; i++) {
      x = mp_random_limb(rng, arg);

      ASSERT(mp_ctz(x) == mp_ctz_simple(x));
    }
  }

  {
    ASSERT(mp_bitlen(0) == 0);
    ASSERT(mp_bitlen(1) == 1);
    ASSERT(mp_bitlen(2) == 2);
    ASSERT(mp_bitlen(MP_LIMB_HI) == MP_LIMB_BITS);
    ASSERT(mp_bitlen(MP_LIMB_MAX) == MP_LIMB_BITS);

    for (i = 0; i < MP_LIMB_BITS; i++)
      ASSERT(mp_bitlen(MP_LIMB_HI >> i) == MP_LIMB_BITS - i);

    for (i = 0; i < 100; i++) {
      x = mp_random_limb(rng, arg);

      ASSERT(mp_bitlen(x) == mp_bitlen_simple(x));
    }
  }

  ASSERT(!mp_mul_gt_2(0, 0, 0, 1));
  ASSERT(!mp_mul_gt_2(1, 0, 1, 1));
  ASSERT(mp_mul_gt_2(2, MP_LIMB_MAX, 1, 0));
  ASSERT(mp_mul_gt_2(1, 1, 0, 0));

  ASSERT(mp_limb_cast(MP_LONG_C(0)) == MP_LIMB_C(0));
  ASSERT(mp_long_cast(MP_LIMB_C(0), -1) == MP_LONG_C(0));

  ASSERT(mp_limb_cast(MP_LONG_C(100)) == MP_LIMB_C(100));
  ASSERT(mp_long_cast(MP_LIMB_C(100), 1) == MP_LONG_C(100));

  ASSERT(mp_limb_cast(MP_LONG_C(-100)) == MP_LIMB_C(100));
  ASSERT(mp_long_cast(MP_LIMB_C(100), -1) == MP_LONG_C(-100));

  ASSERT(mp_long_cast(MP_LIMB_HI, 1) == 0);
  ASSERT(mp_long_cast(MP_LIMB_HI + 1, 1) == 1);
  ASSERT(mp_long_cast(MP_LIMB_HI + 1, -1) == -1);

  ASSERT(mp_limb_cast(MP_LONG_MIN) == MP_LIMB_HI);
  ASSERT(mp_long_cast(MP_LIMB_HI, -1) == MP_LONG_MIN);

  ASSERT(mp_limb_cast(MP_LONG_MIN + 1) == MP_LIMB_HI - 1);
  ASSERT(mp_long_cast(MP_LIMB_HI - 1, -1) == MP_LONG_MIN + 1);
}

static void
test_mp_div(mp_rng_f *rng, void *arg) {
  mp_limb_t z1, z0, n1, n0, d, q, r;
  int i = 0;

  printf("  - MP div.\n");

  for (i = 0; i < 100; i++) {
    do {
      n1 = mp_random_limb(rng, arg);
      n0 = mp_random_limb(rng, arg);
      d = mp_random_limb(rng, arg);
    } while (n1 >= d);

    ASSERT(d != 0);

    mp_div(&q, &r, n1, n0, d);

    z0 = r;
    z1 = mpn_addmul_1(&z0, &d, 1, q);

    ASSERT(z0 == n0);
    ASSERT(z1 == n1);
  }

  for (i = 0; i < 100; i++) {
    n1 = 0;
    n0 = mp_random_limb(rng, arg);
    d = mp_random_limb_nz(rng, arg);

    mp_div(&q, &r, n1, n0, d);

    ASSERT(q == n0 / d);
    ASSERT(r == n0 % d);
  }
}

static void
test_mp_inv_2by1(mp_rng_f *rng, void *arg) {
  mp_limb_t bp[2], np[2], dp[2], qp[2];
  mp_limb_t d, v, q;
  int i;

  printf("  - MP inv (2-by-1).\n");

  bp[0] = 0;
  bp[1] = 1;

  np[0] = MP_LIMB_MAX;
  np[1] = MP_LIMB_MAX;

  for (i = 0; i < 100; i++) {
    /* v = (B^2 - 1) / d - B */
    /* d = (B^2 - 1) / (B + v) */
    d = mp_random_limb(rng, arg) | MP_LIMB_HI;
    v = mp_inv_2by1(d);

    mpn_divmod_1_simple(qp, np, 2, d);
    mpn_sub_n(qp, qp, bp, 2);

    ASSERT(qp[0] == v);
    ASSERT(qp[1] == 0);
  }

  for (i = 0; i < 100; i++) {
    /* v = (B^2 - 1) / d - B */
    /* d = (B^2 - 1) / (B + v) */
    d = mp_random_limb(rng, arg) | MP_LIMB_HI;
    v = mp_inv_2by1(d);

    dp[0] = v;
    dp[1] = 1;

    mpn_divmod_simple(&q, NULL, np, 2, dp, 2);

    ASSERT(q == d);
  }
}

static void
test_mp_div_2by1(mp_rng_f *rng, void *arg) {
  mp_limb_t z1, z0, n1, n0, d, m, q, r;
  mp_bits_t s;
  int i = 0;

  printf("  - MP div (2-by-1).\n");

  for (i = 0; i < 100; i++) {
    do {
      n1 = mp_random_limb(rng, arg);
      n0 = mp_random_limb(rng, arg);
      d = mp_random_limb(rng, arg);
    } while (n1 >= d);

    ASSERT(d != 0);

    s = mp_clz(d);
    d <<= s;
    m = mp_inv_2by1(d);

    if (s != 0) {
      n1 = (n1 << s) | (n0 >> (MP_LIMB_BITS - s));
      n0 <<= s;
    }

    mp_div_2by1(&q, &r, n1, n0, d, m);

    r >>= s;
    d >>= s;

    if (s != 0) {
      n0 = (n0 >> s) | (n1 << (MP_LIMB_BITS - s));
      n1 >>= s;
    }

    z0 = r;
    z1 = mpn_addmul_1(&z0, &d, 1, q);

    ASSERT(z0 == n0);
    ASSERT(z1 == n1);

    mp_div(&z0, &z1, n1, n0, d);

    ASSERT(z0 == q);
    ASSERT(z1 == r);

    mp_div_2by1_simple(&z0, &z1, n1, n0, d);

    ASSERT(z0 == q);
    ASSERT(z1 == r);
  }
}

static void
test_mp_inv_3by2(mp_rng_f *rng, void *arg) {
  mp_limb_t bp[2], np[3], dp[2], qp[2];
  mp_limb_t d1, d0, v;
  int i;

  printf("  - MP inv (3-by-2).\n");

  bp[0] = 0;
  bp[1] = 1;

  np[0] = MP_LIMB_MAX;
  np[1] = MP_LIMB_MAX;
  np[2] = MP_LIMB_MAX;

  for (i = 0; i < 100; i++) {
    /* v = (B^3 - 1) / (d1, d0) - B */
    /* (d1, d0) = (B^3 - 1) / (B + v) */
    d1 = mp_random_limb(rng, arg) | MP_LIMB_HI;
    d0 = mp_random_limb(rng, arg);
    v = mp_inv_3by2(d1, d0);

    dp[0] = d0;
    dp[1] = d1;

    mpn_divmod_simple(qp, NULL, np, 3, dp, 2);
    mpn_sub_n(qp, qp, bp, 2);

    ASSERT(qp[0] == v);
    ASSERT(qp[1] == 0);
  }

  for (i = 0; i < 100; i++) {
    /* v = (B^3 - 1) / (d1, d0) - B */
    /* (d1, d0) = (B^3 - 1) / (B + v) */
    d1 = mp_random_limb(rng, arg) | MP_LIMB_HI;
    d0 = mp_random_limb(rng, arg);
    v = mp_inv_3by2(d1, d0);

    dp[0] = v;
    dp[1] = 1;

    mpn_divmod_simple(qp, NULL, np, 3, dp, 2);

    ASSERT(qp[1] == d1 - 1
        || qp[1] == d1 + 0
        || qp[1] == d1 + 1);
  }
}

static void
test_mp_div_3by2(mp_rng_f *rng, void *arg) {
  mp_limb_t n2, n1, n0, d1, d0, m, q, r1, r0;
  mp_limb_t zp[3], dp[2];
  mp_bits_t s;
  int i = 0;

  printf("  - MP div (3-by-2).\n");

  for (i = 0; i < 100; i++) {
    do {
      n2 = mp_random_limb(rng, arg);
      n1 = mp_random_limb(rng, arg);
      n0 = mp_random_limb(rng, arg);
      d1 = mp_random_limb(rng, arg);
      d0 = mp_random_limb(rng, arg);
    } while (n2 > d1 || (n2 == d1 && n1 >= d0));

    ASSERT((d0 | d1) != 0);

    s = mp_clz(d1);

    if (s != 0) {
      d1 = (d1 << s) | (d0 >> (MP_LIMB_BITS - s));
      d0 <<= s;
    }

    m = mp_inv_3by2(d1, d0);

    if (s != 0) {
      n2 = (n2 << s) | (n1 >> (MP_LIMB_BITS - s));
      n1 = (n1 << s) | (n0 >> (MP_LIMB_BITS - s));
      n0 <<= s;
    }

    mp_div_3by2(&q, &r1, &r0, n2, n1, n0, d1, d0, m);

    if (s != 0) {
      r0 = (r0 >> s) | (r1 << (MP_LIMB_BITS - s));
      r1 >>= s;

      d0 = (d0 >> s) | (d1 << (MP_LIMB_BITS - s));
      d1 >>= s;

      n0 = (n0 >> s) | (n1 << (MP_LIMB_BITS - s));
      n1 = (n1 >> s) | (n2 << (MP_LIMB_BITS - s));
      n2 >>= s;
    }

    dp[0] = d0;
    dp[1] = d1;

    zp[0] = r0;
    zp[1] = r1;
    zp[2] = mpn_addmul_1(zp, dp, 2, q);

    ASSERT(zp[0] == n0);
    ASSERT(zp[1] == n1);
    ASSERT(zp[2] == n2);

    mp_div_3by2_simple(&zp[0], &zp[1], &zp[2],
                       n2, n1, n0, d1, d0);

    ASSERT(zp[0] == q);
    ASSERT(zp[1] == r1);
    ASSERT(zp[2] == r0);
  }
}

static void
test_mp_inv_mod(mp_rng_f *rng, void *arg) {
  {
    mp_limb_t d, m;
    int i;

    printf("  - MP inv mod.\n");

    for (i = 0; i < 100; i++) {
      d = mp_random_limb(rng, arg) | 1;
      m = mp_inv_mod(d);

      ASSERT(m == mp_inv_mod_simple(d));
      ASSERT(d * m == 1);
      ASSERT(3 * d * m == 3);
    }
  }

  {
    mp_limb_t d, m, q, c, hi, lo;
    mp_limb_t np[3];
    mp_limb_t qp[3];

    d = mp_random_limb(rng, arg) | MP_LIMB_HI | 1;
    m = mp_inv_mod(d);
    q = mp_random_limb(rng, arg);

    /* Multiply. */
    np[0] = q;
    np[1] = mpn_mul_1(np, np, 1, d);
    np[2] = mpn_mul_1(np, np, 2, d);

    /* Divide. */
    qp[0] = np[0] * m;

    /* qp[1] = (np[1] - (qp[0] * d)) * m; */
    mp_mul(hi, lo, qp[0], d);
    mp_sub(qp[1], c, np[1], hi);
    qp[1] *= m;

    /* qp[2] = (np[2] - (qp[1] * d) - c) * m; */
    mp_mul(hi, lo, qp[1], d);
    c += hi;
    mp_sub(qp[2], c, np[2], c);
    qp[2] *= m;

    /* New numerator. */
    np[0] = qp[0];
    np[1] = qp[1];

    ASSERT(qp[2] == 0);

    /* Divide again. */
    qp[0] = np[0] * m;

    /* qp[1] = (np[1] - (qp[0] * d)) * m; */
    mp_mul(hi, lo, qp[0], d);
    mp_sub(qp[1], c, np[1], hi);
    qp[1] *= m;

    (void)lo;

    ASSERT(qp[0] == q);
    ASSERT(qp[1] == 0);
  }
}

static void
test_mp_eratosthenes(void) {
  mp_size_t len = ARRAY_SIZE(mpz_test_primes);
  mp_limb_t n = mpz_test_primes[len - 1];
  mp_limb_t *sieve;
  mp_size_t i = 0;
  mp_limb_t p;

  printf("  - Sieve of Eratosthenes.\n");

  sieve = mp_eratosthenes(n);

  for (p = 0; p <= n; p++) {
    if (mpn_tstbit(sieve, p)) {
      ASSERT(i < len);
      ASSERT(p == (mp_limb_t)mpz_test_primes[i++]);
    }
  }

  ASSERT(i == len);

  mp_free_limbs(sieve);
}

/*
 * MPN
 */

static void
test_mpn_init(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4];

  printf("  - MPN init.\n");

  mpn_random_nz(xp, 4, rng, arg);

  ASSERT(!mpn_zero_p(xp, 4));

  mpn_zero(xp, 4);

  ASSERT(mpn_zero_p(xp, 4));

  mpn_cleanse(xp, 4);
}

static void
test_mpn_assign(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4];

  printf("  - MPN assignment.\n");

  mpn_random_nz(xp, 4, rng, arg);
  mpn_zero(yp, 4);

  ASSERT(mpn_cmp(xp, yp, 4) != 0);

  mpn_copyi(yp, xp, 4);

  ASSERT(mpn_cmp(xp, yp, 4) == 0);

  mpn_zero(yp, 4);

  ASSERT(mpn_cmp(xp, yp, 4) != 0);

  mpn_copyd(yp, xp, 4);

  ASSERT(mpn_cmp(xp, yp, 4) == 0);

  mpn_set_1(yp, 4, 1);

  ASSERT(mpn_strip(yp, 4) == 1);
  ASSERT(yp[0] == 1);
}

static void
test_mpn_cmp(void) {
  static const mp_limb_t zero[4] = {0, 0, 0, 0};
  static const mp_limb_t mod[4] = {4, 3, 2, 1};
  static const mp_limb_t minus1[4] = {3, 3, 2, 1};
  static const mp_limb_t plus1[4] = {5, 3, 2, 1};
  static const mp_limb_t full[4] = {0xff, 0xff, 0xff, 0xff};

  printf("  - MPN comparison.\n");

  ASSERT(mpn_zero_p(zero, 4));
  ASSERT(!mpn_zero_p(mod, 4));
  ASSERT(mpn_cmp(minus1, mod, 4) == -1);
  ASSERT(mpn_cmp(mod, mod, 4) == 0);
  ASSERT(mpn_cmp(plus1, mod, 4) == 1);
  ASSERT(mpn_cmp(mod, full, 4) == -1);
  ASSERT(mpn_cmp(full, mod, 4) == 1);
}

static void
test_mpn_addsub(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[6], yp[6], zp[6];
  mp_size_t xn, yn, zn;
  int i;

  printf("  - MPN add/sub.\n");

  for (i = 0; i < 100; i++) {
    mpn_zero(xp, 6);
    mpn_zero(yp, 6);
    mpn_zero(zp, 6);

    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(yp, 4, rng, arg);

    xp[3] |= (xp[3] == 0);
    yp[3] |= (yp[3] == 0);

    zp[4] = mpn_add_n(zp, xp, yp, 4);

    ASSERT(mpn_cmp(zp, xp, 5) > 0);
    ASSERT(mpn_cmp(zp, yp, 5) > 0);

    zp[5] = mpn_add_n(zp, zp, yp, 5);

    ASSERT(mpn_cmp(zp, xp, 6) > 0);
    ASSERT(mpn_cmp(zp, yp, 6) > 0);

    ASSERT(mpn_sub_n(zp, zp, yp, 6) == 0);

    ASSERT(zp[5] == 0);
    ASSERT(mpn_cmp(zp, xp, 5) > 0);
    ASSERT(mpn_cmp(zp, yp, 5) > 0);

    ASSERT(mpn_sub_n(zp, zp, yp, 5) == 0);

    ASSERT(zp[4] == 0);
    ASSERT(mpn_cmp(zp, xp, 4) == 0);

    ASSERT(mpn_sub_n(zp, zp, xp, 4) == 0);

    ASSERT(mpn_strip(zp, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(yp, 4, rng, arg);

    xn = mpn_strip(xp, 4);
    yn = mpn_strip(yp, 4);
    zn = MP_MAX(xn, yn);

    if (xn >= yn)
      zp[zn] = mpn_add(zp, xp, xn, yp, yn);
    else
      zp[zn] = mpn_add(zp, yp, yn, xp, xn);

    zn += (zp[zn] != 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(mpn_cmp3(zp, zn, yp, yn) > 0);

    zp[zn] = mpn_add(zp, zp, zn, yp, yn);

    zn += (zp[zn] != 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(mpn_cmp3(zp, zn, yp, yn) > 0);

    ASSERT(mpn_sub(zp, zp, zn, yp, yn) == 0);

    zn = mpn_strip(zp, zn);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(mpn_cmp3(zp, zn, yp, yn) > 0);

    ASSERT(mpn_sub(zp, zp, zn, yp, yn) == 0);

    zn = mpn_strip(zp, zn);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) == 0);

    ASSERT(mpn_sub(zp, zp, zn, xp, xn) == 0);

    zn = mpn_strip(zp, zn);

    ASSERT(zn == 0);
  }
}

static void
test_mpn_addsub_1(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[6], zp[6];
  mp_size_t xn, zn;
  mp_limb_t y;
  int i;

  printf("  - MPN add/sub (1 limb).\n");

  for (i = 0; i < 100; i++) {
    mpn_zero(xp, 6);
    mpn_zero(zp, 6);

    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    xp[3] |= (xp[3] == 0);

    zp[4] = mpn_add_1(zp, xp, 4, y);

    ASSERT(mpn_cmp(zp, xp, 5) > 0);

    zp[5] = mpn_add_1(zp, zp, 5, y);

    ASSERT(mpn_cmp(zp, xp, 6) > 0);

    ASSERT(mpn_sub_1(zp, zp, 6, y) == 0);

    ASSERT(zp[5] == 0);
    ASSERT(mpn_cmp(zp, xp, 5) > 0);

    ASSERT(mpn_sub_1(zp, zp, 5, y) == 0);

    ASSERT(zp[4] == 0);
    ASSERT(mpn_cmp(zp, xp, 4) == 0);

    ASSERT(mpn_sub_n(zp, zp, xp, 4) == 0);

    ASSERT(mpn_strip(zp, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    xn = mpn_strip(xp, 4);
    zn = xn;

    zp[zn] = mpn_add_1(zp, xp, xn, y);

    zn += (zp[zn] != 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(zn > 1 || zp[0] > y);

    zp[zn] = mpn_add_1(zp, zp, zn, y);

    zn += (zp[zn] != 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(zn > 1 || zp[0] > y);

    ASSERT(mpn_sub_1(zp, zp, zn, y) == 0);

    zn = mpn_strip(zp, zn);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(zn > 1 || zp[0] > y);

    ASSERT(mpn_sub_1(zp, zp, zn, y) == 0);

    zn = mpn_strip(zp, zn);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) == 0);

    ASSERT(mpn_sub(zp, zp, zn, xp, xn) == 0);

    zn = mpn_strip(zp, zn);

    ASSERT(zn == 0);
  }
}

static void
test_mpn_add_1(void) {
  printf("  - MPN add (1 limb).\n");

  {
    mp_limb_t xp[4] = {1, 2, 3, 4};
    mp_limb_t ep[4] = {2, 2, 3, 4};
    mp_limb_t zp[4];

    mpn_zero(zp, 4);

    ASSERT(mpn_add_1(zp, xp, 4, 1) == 0);
    ASSERT(mpn_cmp(zp, ep, 4) == 0);
  }

  {
    mp_limb_t xp[4] = {MP_LIMB_MAX, 2, 3, 4};
    mp_limb_t ep[4] = {0, 3, 3, 4};
    mp_limb_t zp[4];

    mpn_zero(zp, 4);

    ASSERT(mpn_add_1(zp, xp, 4, 1) == 0);
    ASSERT(mpn_cmp(zp, ep, 4) == 0);
  }

  {
    mp_limb_t xp[4] = {MP_LIMB_MAX, MP_LIMB_MAX, MP_LIMB_MAX, MP_LIMB_MAX};
    mp_limb_t ep[4] = {0, 0, 0, 0};
    mp_limb_t zp[4];

    mpn_zero(zp, 4);

    ASSERT(mpn_add_1(zp, xp, 4, 1) == 1);
    ASSERT(mpn_cmp(zp, ep, 4) == 0);
  }

  {
    mp_limb_t xp[4] = {MP_LIMB_MAX, MP_LIMB_MAX, MP_LIMB_MAX, MP_LIMB_MAX};
    mp_limb_t ep[4] = {MP_LIMB_MAX - 1, 0, 0, 0};
    mp_limb_t zp[4];

    mpn_zero(zp, 4);

    ASSERT(mpn_add_1(zp, xp, 4, MP_LIMB_MAX) == 1);
    ASSERT(mpn_cmp(zp, ep, 4) == 0);
  }

  {
    mp_limb_t xp[4] = {MP_LIMB_MAX, 1, 0, 0};
    mp_limb_t ep[4] = {MP_LIMB_MAX - 1, 2, 0, 0};
    mp_limb_t zp[4];

    mpn_zero(zp, 4);

    ASSERT(mpn_add_1(zp, xp, 4, MP_LIMB_MAX) == 0);
    ASSERT(mpn_cmp(zp, ep, 4) == 0);
  }
}

static void
test_mpn_muldiv(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[12], yp[12], up[12], vp[12];
  mp_size_t xn, yn, un, vn;
  int i;

  printf("  - MPN mul/div.\n");

  for (i = 0; i < 100; i++) {
    mpn_zero(xp, 12);
    mpn_zero(yp, 12);
    mpn_zero(up, 12);
    mpn_zero(vp, 12);

    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(yp, 4, rng, arg);

    xp[3] |= (xp[3] == 0);
    yp[3] |= (yp[3] == 0);

    mpn_mul_n(up, xp, yp, 4);

    ASSERT(mpn_cmp(up, xp, 8) > 0);
    ASSERT(mpn_cmp(up, yp, 8) > 0);

    mpn_mul(vp, up, 8, yp, 4);

    ASSERT(mpn_cmp(vp, up, 12) > 0);
    ASSERT(mpn_cmp(vp, yp, 12) > 0);

    mpn_div(up, vp, 12, yp, 4);

    ASSERT(mpn_cmp(up, xp, 9) > 0);
    ASSERT(mpn_cmp(up, yp, 9) > 0);

    mpn_div(vp, up, 9, yp, 4);

    ASSERT(mpn_cmp(vp, xp, 6) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(yp, 4, rng, arg);

    xn = mpn_strip(xp, 4);
    yn = mpn_strip(yp, 4);
    un = xn + yn;

    xp[0] += (xn == 1 && xp[0] == 1);
    yp[0] += (yn == 1 && yp[0] == 1);

    mpn_mul(up, xp, xn, yp, yn);

    un -= (up[un - 1] == 0);

    ASSERT(mpn_cmp3(up, un, xp, xn) > 0);
    ASSERT(mpn_cmp3(up, un, yp, yn) > 0);

    vn = un + yn;

    mpn_mul(vp, up, un, yp, yn);

    vn -= (vp[vn - 1] == 0);

    ASSERT(mpn_cmp3(vp, vn, up, un) > 0);
    ASSERT(mpn_cmp3(vp, vn, yp, yn) > 0);

    un = vn - yn + 1;

    mpn_div(up, vp, vn, yp, yn);

    un -= (up[un - 1] == 0);

    ASSERT(mpn_cmp3(up, un, xp, xn) > 0);
    ASSERT(mpn_cmp3(up, un, yp, yn) > 0);

    vn = un - yn + 1;

    mpn_div(vp, up, un, yp, yn);

    vn -= (vp[vn - 1] == 0);

    ASSERT(mpn_cmp3(vp, vn, xp, xn) == 0);
  }
}

static void
test_mpn_muldiv_1(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[6], zp[6];
  mp_size_t xn, zn;
  mp_limb_t y;
  int i;

  printf("  - MPN mul/div (1 limb).\n");

  for (i = 0; i < 100; i++) {
    mpn_zero(xp, 6);
    mpn_zero(zp, 6);

    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    xp[3] |= (xp[3] == 0);

    zp[4] = mpn_mul_1(zp, xp, 4, y);

    ASSERT(mpn_cmp(zp, xp, 5) > 0);

    zp[5] = mpn_mul_1(zp, zp, 5, y);

    ASSERT(mpn_cmp(zp, xp, 6) > 0);

    mpn_div_1(zp, zp, 6, y);

    ASSERT(zp[5] == 0);

    ASSERT(mpn_cmp(zp, xp, 5) > 0);

    mpn_div_1(zp, zp, 5, y);

    ASSERT(zp[4] == 0);

    ASSERT(mpn_cmp(zp, xp, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    xn = mpn_strip(xp, 4);
    zn = xn;

    xp[0] += (xn == 1 && xp[0] == 1);

    y += (y == 1);

    zp[zn] = mpn_mul_1(zp, xp, xn, y);

    zn += (zp[zn] != 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(zn > 1 || zp[0] > y);

    zp[zn] = mpn_mul_1(zp, zp, zn, y);

    zn += (zp[zn] != 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(zn > 1 || zp[0] > y);

    mpn_div_1(zp, zp, zn, y);

    zn -= (zp[zn - 1] == 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) > 0);
    ASSERT(zn > 1 || zp[0] > y);

    mpn_div_1(zp, zp, zn, y);

    zn -= (zp[zn - 1] == 0);

    ASSERT(mpn_cmp3(zp, zn, xp, xn) == 0);
  }
}

static void
test_mpn_muldivexact(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[2], tp[6], zp[8];
  int i;

  printf("  - MPN mul/divexact.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(yp, 2, rng, arg);

    xp[3] |= (xp[3] == 0);
    yp[1] |= (yp[1] == 0);

    mpn_mul(tp, xp, 4, yp, 2);
    mpn_mul(zp, tp, 6, yp, 2);

    mpn_divexact(zp, zp, 8, yp, 2);

    ASSERT(zp[6] == 0);

    mpn_divexact(zp, zp, 6, yp, 2);

    ASSERT(zp[4] == 0);

    ASSERT(mpn_cmp(zp, xp, 4) == 0);

    mpn_divexact(zp, zp, 4, xp, 4);

    ASSERT(zp[0] == 1);
  }
}

static void
test_mpn_muldivexact_1(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], zp[6];
  mp_limb_t y;
  int i;

  printf("  - MPN mul/divexact.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    xp[3] |= (xp[3] == 0);

    y += (y == 1);

    zp[4] = mpn_mul_1(zp, xp, 4, y);
    zp[5] = mpn_mul_1(zp, zp, 5, y);

    mpn_divexact_1(zp, zp, 6, y);

    ASSERT(zp[5] == 0);

    mpn_divexact_1(zp, zp, 5, y);

    ASSERT(zp[4] == 0);

    ASSERT(mpn_cmp(zp, xp, 4) == 0);

    mpn_divexact(zp, zp, 4, xp, 4);

    ASSERT(zp[0] == 1);
  }
}

static void
test_mpn_addmul_1(mp_rng_f *rng, void *arg) {
  mp_limb_t zp[5], xp[5], ep[5], tp[5];
  mp_limb_t y;
  int i;

  printf("  - MPN add+mul.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(zp, 4, rng, arg);
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    zp[3] |= (zp[3] == 0);
    zp[4] = 0;

    xp[3] |= (xp[3] == 0);
    xp[4] = 0;

    y += (y == 1);

    mpn_copyi(ep, zp, 5);

    tp[4] = mpn_mul_1(tp, xp, 4, y);

    ASSERT(mpn_add_n(ep, ep, tp, 5) == 0);

    ASSERT(mpn_addmul_1(zp, xp, 5, y) == 0);

    ASSERT(mpn_cmp(zp, ep, 5) == 0);
  }
}

static void
test_mpn_submul_1(mp_rng_f *rng, void *arg) {
  mp_limb_t zp[6], xp[6], ep[6], tp[6];
  mp_limb_t y;
  int i;

  printf("  - MPN sub+mul.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(zp, 6, rng, arg);
    mpn_random_nz(xp, 4, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    zp[5] |= (zp[5] == 0);

    xp[3] |= (xp[3] == 0);
    xp[4] = 0;
    xp[5] = 0;

    y += (y == 1);

    mpn_copyi(ep, zp, 6);

    ASSERT(mpn_mul_1(tp, xp, 6, y) == 0);
    ASSERT(mpn_sub_n(ep, ep, tp, 6) == 0);

    ASSERT(mpn_submul_1(zp, xp, 6, y) == 0);

    ASSERT(mpn_cmp(zp, ep, 6) == 0);
  }
}

static void
test_mpn_sqr(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], zp[8], ep[8];
  int i;

  printf("  - MPN sqr.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    mpn_sqr(zp, xp, 4, ep);
    mpn_mul_n(ep, xp, xp, 4);

    ASSERT(mpn_cmp(zp, ep, 8) == 0);
  }
}

static void
test_mpn_mod(mp_rng_f *rng, void *arg) {
  mp_limb_t np[8], dp[4], qp[5], rp[4], tp[9];
  mp_size_t nn, dn, qn, rn, tn;
  int i;

  printf("  - MPN mod.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(np, 8, rng, arg);
    mpn_random_nz(dp, 4, rng, arg);

    np[7] |= (np[7] == 0);
    dp[3] |= (dp[3] == 0);

    mpn_divmod(qp, rp, np, 8, dp, 4);

    mpn_mul(tp, qp, 5, dp, 4);

    ASSERT(tp[8] == 0);

    ASSERT(mpn_sub_n(tp, np, tp, 8) == 0);

    ASSERT(mpn_strip(tp, 8) <= 4);
    ASSERT(mpn_cmp(rp, tp, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    do {
      mpn_random_nz(np, 8, rng, arg);
      mpn_random_nz(dp, 4, rng, arg);
    } while (mpn_strip(np, 8) < 5);

    nn = mpn_strip(np, 8);
    dn = mpn_strip(dp, 4);
    qn = nn - dn + 1;

    mpn_divmod(qp, rp, np, nn, dp, dn);

    qn -= (qp[qn - 1] == 0);
    rn = mpn_strip(rp, dn);
    tn = qn + dn;

    mpn_mul(tp, qp, qn, dp, dn);

    tn -= (tp[tn - 1] == 0);

    ASSERT(mpn_sub(tp, np, nn, tp, tn) == 0);

    tn = mpn_strip(tp, tn);

    ASSERT(rn == tn);
    ASSERT(mpn_cmp(rp, tp, tn) == 0);
  }
}

static void
test_mpn_mod_1(mp_rng_f *rng, void *arg) {
  mp_limb_t np[8], qp[8], tp[9];
  mp_size_t nn, qn, tn;
  mp_limb_t d, r;
  int i;

  printf("  - MPN mod (1 limb).\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(np, 8, rng, arg);
    mpn_random_nz(&d, 1, rng, arg);

    np[7] |= (np[7] == 0);

    d += (d == 1);

    r = mpn_divmod_1(qp, np, 8, d);

    ASSERT(mpn_submul_1(np, qp, 8, d) == 0);

    ASSERT(mpn_strip(np, 8) <= 1);
    ASSERT(np[0] == r);
  }

  for (i = 0; i < 100; i++) {
    do {
      mpn_random_nz(np, 8, rng, arg);
      mpn_random_nz(&d, 1, rng, arg);
    } while (mpn_strip(np, 8) < 2);

    d += (d == 1);

    nn = mpn_strip(np, 8);
    qn = nn;

    r = mpn_divmod_1(qp, np, nn, d);

    qn -= (qp[qn - 1] == 0);
    tn = qn;

    tp[tn] = mpn_mul_1(tp, qp, qn, d);

    tn += (tp[tn] != 0);

    ASSERT(mpn_sub(tp, np, nn, tp, tn) == 0);

    tn = mpn_strip(tp, tn);

    ASSERT(tn <= 1 && tp[0] == r);
  }
}

static void
test_mpn_divround(mp_rng_f *rng, void *arg) {
  mp_limb_t zp[7], np[7], qp[5], dp[4], hp[2];
  int i;

  printf("  - MPN divround.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(qp, 4, rng, arg);
    mpn_random_nz(dp, 2, rng, arg);

    dp[0] &= ~MP_LIMB_HI;
    dp[1] |= 1;

    mpn_rshift(hp, dp, 2, 1);

    mpn_mul(np, qp, 4, dp, 2);

    if (i & 1) {
      hp[0] += 1;
      np[6] = mpn_add(np, np, 6, hp, 2);
      qp[4] = mpn_add_1(qp, qp, 4, 1);
    } else {
      hp[0] -= 1;
      np[6] = mpn_add(np, np, 6, hp, 2);
      qp[4] = 0;
    }

    mpn_divround(zp, np, 7, dp, 2);

    ASSERT(mpn_cmp(zp, qp, 5) == 0);
    ASSERT(zp[5] == 0);
    ASSERT(zp[6] == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(dp, 4, rng, arg);

    dp[3] |= 1;

    mpn_sub_1(np, dp, 4, 1);

    mpn_divround(qp, np, 4, dp, 4);

    ASSERT(qp[0] == 1);
    ASSERT(qp[1] == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(dp, 1, rng, arg);

    dp[0] |= MP_LIMB_HI;
    np[0] = dp[0] - 1;

    mpn_divround(qp, np, 1, dp, 1);

    ASSERT(qp[0] == 1);
    ASSERT(qp[1] == 0);

    mpn_divround(qp, np, 0, dp, 1);

    ASSERT(qp[0] == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(dp, 3, rng, arg);

    if (i & 1)
      dp[0] &= ~MP_LIMB_C(1);
    else
      dp[0] |= 1;

    dp[3] = 1;

    mpn_rshift(np, dp, 4, 1);

    mpn_divround(qp, np, 3, dp, 4);

    ASSERT(qp[0] == (mp_limb_t)(i & 1));
  }
}

static void
test_mpn_divround_1(mp_rng_f *rng, void *arg) {
  mp_limb_t zp[7], np[6], qp[5];
  mp_limb_t d;
  int i;

  printf("  - MPN divround (1 limb).\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(qp, 4, rng, arg);
    mpn_random_nz(&d, 1, rng, arg);

    d += (d == 1);

    np[4] = mpn_mul_1(np, qp, 4, d);

    if (i & 1) {
      np[5] = mpn_add_1(np, np, 5, (d >> 1) + 1);
      qp[4] = mpn_add_1(qp, qp, 4, 1);
    } else {
      np[5] = mpn_add_1(np, np, 5, (d >> 1) - 1);
      qp[4] = 0;
    }

    mpn_divround_1(zp, np, 6, d);

    ASSERT(mpn_cmp(zp, qp, 5) == 0);
    ASSERT(zp[5] == 0);
    ASSERT(zp[6] == 0);
  }
}

static void
test_mpn_roots(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], zp[4], rp[4], tp[4], sp[4];
  mp_size_t xn = 4;
  mp_size_t zn = 2;
  mp_size_t rn;
  int i;

  printf("  - MPN roots.\n");

  for (i = 0; i < 50; i++) {
    mpn_random_nz(tp, 2, rng, arg);

    mpn_sqr(xp, tp, 2, sp);

    ASSERT(mpn_perfect_square_p(xp, xn));

    rn = mpn_sqrtrem(zp, rp, xp, xn);

    ASSERT(mpn_cmp(zp, tp, zn) == 0);
    ASSERT(rn == 0);
    ASSERT(mpn_strip(rp, xn) == 0);
  }
}

static void
test_mpn_and(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4], tp[4], ep[4];
  int i, j;

  printf("  - MPN AND.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = xp[j] & yp[j];

    mpn_and_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = xp[j] & ~yp[j];

    mpn_andn_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ~(xp[j] & yp[j]);

    mpn_nand_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_ior(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4], tp[4], ep[4];
  int i, j;

  printf("  - MPN OR.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = xp[j] | yp[j];

    mpn_ior_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = xp[j] | ~yp[j];

    mpn_iorn_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ~(xp[j] | yp[j]);

    mpn_nior_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_xor(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4], tp[4], ep[4];
  int i, j;

  printf("  - MPN XOR.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = xp[j] ^ yp[j];

    mpn_xor_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ~(xp[j] ^ yp[j]);

    mpn_nxor_n(tp, xp, yp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_com(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], tp[4], ep[4];
  int i, j;

  printf("  - MPN NOT.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ~xp[j];

    mpn_com(tp, xp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_lshift(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[5], yp[5];
  int i;

  printf("  - MPN lshift.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    yp[4] = mpn_mul_1(yp, xp, 4, 8);
    xp[4] = mpn_lshift(xp, xp, 4, 3);

    ASSERT(mpn_cmp(xp, yp, 5) == 0);
  }
}

static void
test_mpn_rshift(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4];
  int i;

  printf("  - MPN rshift.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    mpn_div_1(yp, xp, 4, 8);
    mpn_rshift(xp, xp, 4, 3);

    ASSERT(mpn_cmp(xp, yp, 4) == 0);
  }
}

static void
test_mpn_shift(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[5], yp[5];
  int i;

  printf("  - MPN shift.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    xp[4] = 0;
    yp[4] = mpn_lshift(yp, xp, 4, 15);

    ASSERT(mpn_cmp(yp, xp, 5) > 0);

    mpn_rshift(yp, yp, 5, 15);

    ASSERT(mpn_cmp(yp, xp, 5) == 0);
  }
}

static void
test_mpn_bits(mp_rng_f *rng, void *arg) {
  mp_bits_t bits = 4 * MP_LIMB_BITS;
  mp_bits_t steps = (bits + 5 - 1) / 5;
  mp_limb_t xp[4], yp[4];
  mp_limb_t b, b1, b2;
  mp_bits_t j, r;
  mp_size_t s;
  int i;

  printf("  - MPN bits.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_copyi(yp, xp, 4);

    for (j = 0; j < bits; j++) {
      ASSERT(mpn_getbit(xp, 4, j) == (yp[0] & 1));
      ASSERT(mpn_tstbit(xp, j) == (int)(yp[0] & 1));

      mpn_rshift(yp, yp, 4, 1);
    }

    ASSERT(mpn_getbit(xp, 4, j) == 0);

    for (j = steps - 1; j >= 0; j--) {
      b1 = mpn_getbits(xp, 4, j * 5, 5);

      /* y = x >> (j * 5) */
      s = (j * 5) / MP_LIMB_BITS;
      r = (j * 5) % MP_LIMB_BITS;

      if (s >= 4) {
        mpn_zero(yp, 4);
      } else {
        if (r != 0)
          mpn_rshift(yp, xp + s, 4 - s, r);
        else
          mpn_copyi(yp, xp + s, 4 - s);
      }

      b2 = yp[0] & MP_MASK(5);

      ASSERT(b1 == b2);
    }

    ASSERT(mpn_getbits(xp, 4, steps * 5, 5) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_copyi(yp, xp, 4);

    for (j = 0; j < bits; j++) {
      b = mpn_getbit(xp, 4, j);

      if (b)
        mpn_clrbit(yp, j);
      else
        mpn_setbit(yp, j);

      ASSERT(mpn_cmp(yp, xp, 4) != 0);
      ASSERT(mpn_getbit(yp, 4, j) == (b ^ 1));

      if (b)
        mpn_setbit(yp, j);
      else
        mpn_clrbit(yp, j);

      ASSERT(mpn_cmp(yp, xp, 4) == 0);
      ASSERT(mpn_getbit(yp, 4, j) == b);

      mpn_combit(yp, j);

      ASSERT(mpn_cmp(yp, xp, 4) != 0);
      ASSERT(mpn_getbit(yp, 4, j) == (b ^ 1));

      mpn_combit(yp, j);

      ASSERT(mpn_cmp(yp, xp, 4) == 0);
      ASSERT(mpn_getbit(yp, 4, j) == b);
    }
  }
}

static void
test_mpn_scan(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], tp[4];
  mp_size_t xn = 1;
  int i;

  printf("  - MPN scan.\n");

  ASSERT(mpn_set_str(xp, xn, "10101000", 2));

  ASSERT(mpn_scan1(xp, xn, 0) == 3);
  ASSERT(mpn_scan1(xp, xn, 1) == 3);
  ASSERT(mpn_scan1(xp, xn, 2) == 3);
  ASSERT(mpn_scan1(xp, xn, 3) == 3);
  ASSERT(mpn_scan1(xp, xn, 4) == 5);
  ASSERT(mpn_scan1(xp, xn, 5) == 5);
  ASSERT(mpn_scan1(xp, xn, 6) == 7);
  ASSERT(mpn_scan1(xp, xn, 7) == 7);
  ASSERT(mpn_scan1(xp, xn, 8) == MP_BITS_MAX);

  ASSERT(mpn_scan0(xp, xn, 0) == 0);
  ASSERT(mpn_scan0(xp, xn, 1) == 1);
  ASSERT(mpn_scan0(xp, xn, 2) == 2);
  ASSERT(mpn_scan0(xp, xn, 3) == 4);
  ASSERT(mpn_scan0(xp, xn, 4) == 4);
  ASSERT(mpn_scan0(xp, xn, 5) == 6);
  ASSERT(mpn_scan0(xp, xn, 6) == 6);
  ASSERT(mpn_scan0(xp, xn, 7) == 8);
  ASSERT(mpn_scan0(xp, xn, 8) == 8);

  /* Shift left 71 bits. */
#if MP_LIMB_BITS == 64
  xp[1] = xp[0] << 7;
  xp[0] = 0;
  xn = 2;
#else
  xp[2] = xp[0] << 7;
  xp[1] = 0;
  xp[0] = 0;
  xn = 3;
#endif

  ASSERT(mpn_scan1(xp, xn, 0) == 71 + 3);
  ASSERT(mpn_scan1(xp, xn, 71 + 0) == 71 + 3);
  ASSERT(mpn_scan1(xp, xn, 71 + 1) == 71 + 3);
  ASSERT(mpn_scan1(xp, xn, 71 + 2) == 71 + 3);
  ASSERT(mpn_scan1(xp, xn, 71 + 3) == 71 + 3);
  ASSERT(mpn_scan1(xp, xn, 71 + 4) == 71 + 5);
  ASSERT(mpn_scan1(xp, xn, 71 + 5) == 71 + 5);
  ASSERT(mpn_scan1(xp, xn, 71 + 6) == 71 + 7);
  ASSERT(mpn_scan1(xp, xn, 71 + 7) == 71 + 7);
  ASSERT(mpn_scan1(xp, xn, 71 + 8) == MP_BITS_MAX);

  ASSERT(mpn_scan0(xp, xn, 0) == 0);
  ASSERT(mpn_scan0(xp, xn, 71 + 0) == 71 + 0);
  ASSERT(mpn_scan0(xp, xn, 71 + 1) == 71 + 1);
  ASSERT(mpn_scan0(xp, xn, 71 + 2) == 71 + 2);
  ASSERT(mpn_scan0(xp, xn, 71 + 3) == 71 + 4);
  ASSERT(mpn_scan0(xp, xn, 71 + 4) == 71 + 4);
  ASSERT(mpn_scan0(xp, xn, 71 + 5) == 71 + 6);
  ASSERT(mpn_scan0(xp, xn, 71 + 6) == 71 + 6);
  ASSERT(mpn_scan0(xp, xn, 71 + 7) == 71 + 8);
  ASSERT(mpn_scan0(xp, xn, 71 + 8) == 71 + 8);

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, 4, rng, arg);

    ASSERT(mpn_scan1(xp, 4, 0) == mpn_ctz(xp, 4));

    mpn_com(tp, xp, 4);

    ASSERT(mpn_scan0(tp, 4, 0) == mpn_ctz(xp, 4));

    mpn_random_nz(xp, 4, rng, arg);
    mpn_zero(xp, 2);

    ASSERT(mpn_scan1(xp, 4, 0) == mpn_ctz(xp, 4));

    mpn_com(tp, xp, 4);

    ASSERT(mpn_scan0(tp, 4, 0) == mpn_ctz(xp, 4));
  }
}

static void
test_mpn_popcount(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4];
  mp_bits_t c;
  int i, j;

  printf("  - MPN hamming weight.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    c = 0;

    for (j = 0; j < 4; j++)
      c += mp_popcount_simple(xp[j]);

    ASSERT(mpn_popcount(xp, 4) == c);
  }
}

static void
test_mpn_hamdist(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4];
  mp_bits_t c;
  int i, j;

  printf("  - MPN hamming distance.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);

    c = 0;

    for (j = 0; j < 4; j++)
      c += mp_popcount_simple(xp[j] ^ yp[j]);

    ASSERT(mpn_hamdist(xp, yp, 4) == c);
  }
}

static void
test_mpn_mask(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], mp[4], tp[4], ep[4];
  int i;

  printf("  - MPN mask.\n");

  mpn_zero(mp, 4);
  mpn_setbit(mp, 100);
  mpn_sub_1(mp, mp, 4, 1);

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    mpn_and_n(ep, xp, mp, 4);
    mpn_mask(tp, xp, 4, 100);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_negate(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4], cp[4];
  mp_limb_t c1, c2;
  int i;

  printf("  - MPN negation.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_copyi(yp, xp, 4);
    mpn_zero(cp, 4);

    c1 = mpn_sub_n(xp, cp, xp, 4);
    c2 = mpn_neg(yp, yp, 4);

    ASSERT(c1 == c2);
    ASSERT(mpn_cmp(xp, yp, 4) == 0);
  }
}

static void
test_mpn_mulshift(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch[MPN_MULSHIFT_ITCH(MP_P192_LIMBS)];
  mp_limb_t xp[MP_P192_LIMBS * 2];
  mp_limb_t yp[MP_P192_LIMBS * 2];
  mp_limb_t zp[MP_P192_LIMBS];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_limb_t qp[MP_P192_LIMBS * 2];
  mp_limb_t dp[MP_P192_LIMBS * 2];
  mp_size_t mn = MP_P192_LIMBS;
  mp_size_t sn = mn * 2;
  mp_size_t zn, dn, qn;
  int i;

  printf("  - MPN mulshift.\n");

  mpn_p192_mod(mp);

  mpn_zero(dp, mn * 2);
  mpn_setbit(dp, 196);

  dn = mpn_strip(dp, mn * 2);

  for (i = 0; i < 100; i++) {
    mpn_random(xp, mn * 2, rng, arg);
    mpn_random(yp, mn * 2, rng, arg);

    mpn_mod(xp, xp, mn * 2, mp, mn);
    mpn_mod(yp, yp, mn * 2, mp, mn);

    ASSERT(mpn_mulshift(zp, xp, yp, mn, 196, scratch) == 0);

    zn = mpn_strip(zp, mn);

    mpn_mul_n(sp, xp, yp, mn);
    mpn_divround(qp, sp, sn, dp, dn);

    qn = mpn_strip(qp, sn - dn + 2);

    ASSERT(zn == qn);
    ASSERT(mpn_cmp(zp, qp, qn) == 0);
  }
}

static void
test_mpn_reduce_weak(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch[MPN_REDUCE_WEAK_ITCH(4)];
  mp_limb_t xp[4], yp[4], zp[4];
  mp_limb_t c;
  int i, j;

  printf("  - MPN reduction (weak).\n");

  for (i = 0; i < 100; i++) {
    do {
      mpn_random(xp, 4, rng, arg);
      mpn_random(yp, 4, rng, arg);

      xp[0] |= 1;
      yp[0] |= 1;
      xp[3] |= MP_LIMB_HI;
      yp[3] |= MP_LIMB_HI;
    } while (mpn_cmp(xp, yp, 4) < 0);

    ASSERT(mpn_reduce_weak(zp, xp, yp, 4, 0, scratch));
    ASSERT(mpn_cmp(zp, yp, 4) < 0);

    ASSERT(!mpn_reduce_weak(zp, zp, yp, 4, 0, scratch));

    do {
      mpn_random(xp, 4, rng, arg);
    } while (mpn_cmp(xp, yp, 4) >= 0);

    mpn_copyi(zp, xp, 4);

    j = 0;

    do {
      j += 1;
      c = mpn_add_n(zp, zp, yp, 4);
    } while (c == 0);

    ASSERT(mpn_reduce_weak(zp, zp, yp, 4, c, scratch));

    j -= 1;

    while (j--)
      ASSERT(mpn_reduce_weak(zp, zp, yp, 4, 0, scratch));

    ASSERT(mpn_cmp(zp, yp, 4) < 0);
    ASSERT(mpn_cmp(zp, xp, 4) == 0);

    ASSERT(!mpn_reduce_weak(zp, zp, yp, 4, 0, scratch));

    ASSERT(mpn_cmp(zp, xp, 4) == 0);
  }
}

static void
test_mpn_barrett(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch1[MPN_BARRETT_ITCH(MP_P192_SHIFT)];
  mp_limb_t scratch2[MPN_REDUCE_ITCH(MP_P192_LIMBS, MP_P192_SHIFT)];
  mp_limb_t mp[MP_P192_SHIFT - MP_P192_LIMBS + 1]; /* MP_P192_LIMBS + 1 */
  mp_limb_t xp[MP_P192_SHIFT];
  mp_limb_t ep[MP_P192_LIMBS];
  mp_limb_t zp[MP_P192_LIMBS];
  mp_limb_t np[MP_P192_LIMBS];
  mp_size_t shift = MP_P192_SHIFT;
  mp_size_t n = MP_P192_LIMBS;
  int i;

  printf("  - MPN reduction (barrett).\n");

  mpn_p192_mod(np);

  mpn_barrett(mp, np, n, shift, scratch1);

  for (i = 0; i < 100; i++) {
    mpn_random(xp, n * 2, rng, arg);
    mpn_mod(ep, xp, n * 2, np, n);

    mpn_reduce(zp, xp, mp, np, n, shift, scratch2);

    ASSERT(mpn_cmp(zp, ep, n) == 0);
  }
}

static void
test_mpn_mont(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch1[MPN_MONT_ITCH(MP_P192_LIMBS)];
  mp_limb_t scratch2[MPN_MONTMUL_ITCH(MP_P192_LIMBS)];
  mp_limb_t up[MP_P192_LIMBS * 2];
  mp_limb_t vp[MP_P192_LIMBS * 2];
  mp_limb_t xp[MP_P192_LIMBS];
  mp_limb_t yp[MP_P192_LIMBS];
  mp_limb_t rp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_limb_t ep[MP_P192_LIMBS];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_size_t mn = MP_P192_LIMBS;
  mp_limb_t k;
  int i;

  printf("  - MPN reduction (montgomery).\n");

  mpn_p192_mod(mp);

  mpn_mont(&k, rp, mp, mn, scratch1);

  for (i = 0; i < 100; i++) {
    mpn_random(up, mn * 2, rng, arg);
    mpn_random(vp, mn * 2, rng, arg);

    mpn_mod(up, up, mn * 2, mp, mn);
    mpn_mod(vp, vp, mn * 2, mp, mn);

    mpn_mul_n(sp, up, vp, mn);
    mpn_mod(ep, sp, mn * 2, mp, mn);

    mpn_montmul(xp, up, rp, mp, mn, k, scratch2);
    mpn_montmul(yp, vp, rp, mp, mn, k, scratch2);
    mpn_montmul(sp, xp, yp, mp, mn, k, scratch2);
    mpn_set_1(yp, mn, 1);
    mpn_montmul(sp, sp, yp, mp, mn, k, scratch2);

    ASSERT(mpn_cmp(sp, ep, mn) == 0);

    mpn_montmul_var(xp, up, rp, mp, mn, k, scratch2);
    mpn_montmul_var(yp, vp, rp, mp, mn, k, scratch2);
    mpn_montmul_var(sp, xp, yp, mp, mn, k, scratch2);
    mpn_set_1(yp, mn, 1);
    mpn_montmul_var(sp, sp, yp, mp, mn, k, scratch2);
    mpn_mod(sp, sp, mn, mp, mn);

    ASSERT(mpn_cmp(sp, ep, mn) == 0);
  }
}

static void
test_mpn_gcd(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch[MPN_GCD_ITCH(4, 4)];
  mp_limb_t xp[4], yp[4], gp[4], zp[4];
  mp_size_t n = 4;
  mp_limb_t y, g;
  int i;

  printf("  - MPN gcd.\n");

  for (i = 0; i < 50; i++) {
    mpn_random_nz(xp, n, rng, arg);
    mpn_random_nz(yp, n, rng, arg);

    xp[n - 1] |= (xp[n - 1] == 0);
    yp[n - 1] |= (yp[n - 1] == 0);

    mpn_gcd(gp, xp, n, yp, n, scratch);
    mpn_gcd_simple(zp, xp, n, yp, n);

    ASSERT(mpn_cmp(gp, zp, n) == 0);
    ASSERT(mpn_cmp(gp, zp, n) == 0);
  }

  for (i = 0; i < 50; i++) {
    mpn_random_nz(xp, n, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    xp[n - 1] |= (xp[n - 1] == 0);

    g = mpn_gcd_1(xp, n, y, scratch);

    mpn_gcd_simple(zp, xp, n, &y, 1);

    ASSERT(g == zp[0]);
  }
}

static void
test_mpn_gcdext(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[MP_P192_LIMBS + 1];
  mp_limb_t yp[MP_P192_LIMBS + 1];
  mp_limb_t gp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS + 1];
  mp_limb_t tp[MP_P192_LIMBS * 2 + 2];
  mp_limb_t up[MP_P192_LIMBS * 2];
  mp_limb_t vp[MP_P192_LIMBS * 2];
  mp_limb_t q1[MP_P192_LIMBS + 1];
  mp_limb_t q2[MP_P192_LIMBS + 1];
  mp_size_t gn, sn, qn, q1n, q2n;
  mp_size_t n = MP_P192_LIMBS;
  int i;

  printf("  - MPN gcdext.\n");

  mpn_p192_mod(yp);

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, n, rng, arg);

    xp[n - 1] |= (xp[n - 1] == 0);

    gn = mpn_gcdext(gp, sp, &sn, xp, n, yp, n);

    ASSERT(gn == 1);
    ASSERT(gp[0] == 1);

    if (sn < 0)
      mpn_sub(sp, yp, n, sp, -sn);

    mpn_mul_n(up, xp, sp, n);
    mpn_mod(up, up, n * 2, yp, n);

    ASSERT(mpn_strip(up, n) == 1);
    ASSERT(up[0] == 1);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(xp, n, rng, arg);
    mpn_random_nz(yp, n, rng, arg);

    xp[n - 1] |= (xp[n - 1] == 0);
    yp[n - 1] |= (yp[n - 1] == 0);

    gn = mpn_gcdext(gp, sp, &sn, xp, n, yp, n);

    ASSERT(gn > 0);
    ASSERT(mpn_strip(gp, n) == gn);

    if (sn < 0)
      mpn_sub(sp, yp, n, sp, -sn);

    /* (g - x * s) / y == t */
    mpn_mul_n(tp, xp, sp, n);
    mpn_sub(tp, tp, n * 2, gp, n);
    mpn_div(tp, tp, n * 2, yp, n);

    /* x * s + y * t == g */
    mpn_mul_n(up, xp, sp, n);
    mpn_mul_n(vp, yp, tp, n);

    if (mpn_cmp(up, vp, n * 2) >= 0)
      mpn_sub_n(up, up, vp, n * 2);
    else
      mpn_sub_n(up, vp, up, n * 2);

    ASSERT(mpn_cmp(up, gp, n) == 0);

    mpn_div(q1, xp, n, gp, gn);
    mpn_div(q2, yp, n, gp, gn);

    qn = n - gn + 1;
    q1n = qn - (q1[qn - 1] == 0);
    q2n = qn - (q2[qn - 1] == 0);

    mpn_zero(q1 + qn, n - qn);
    mpn_zero(q2 + qn, n - qn);

    gn = mpn_gcdext(gp, sp, &sn, q1, q1n, q2, q2n);

    ASSERT(gn == 1);
    ASSERT(gp[0] == 1);

    if (sn < 0)
      mpn_sub(sp, yp, n, sp, -sn);

    /* (g - x * s) / y == t */
    mpn_mul_n(tp, xp, sp, n);
    mpn_sub(tp, tp, n * 2, gp, n);
    mpn_div(tp, tp, n * 2, yp, n);

    /* q1 * s + q2 * t == 1 */
    mpn_mul_n(up, q1, sp, n);
    mpn_mul_n(vp, q2, tp, n);

    if (mpn_cmp(up, vp, n * 2) >= 0)
      mpn_sub_n(up, up, vp, n * 2);
    else
      mpn_sub_n(up, vp, up, n * 2);

    ASSERT(mpn_strip(up, n * 2) == 1);
    ASSERT(up[0] == 1);
  }
}

static void
test_mpn_invert(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch[MPN_INVERT_ITCH(MP_P192_LIMBS)];
  mp_limb_t xp[MP_P192_LIMBS * 2];
  mp_limb_t zp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_size_t mn = MP_P192_LIMBS;
  int i = 0;

  printf("  - MPN invert.\n");

  mpn_p192_mod(mp);

  mpn_zero(xp, mn);

  ASSERT(mpn_invert_n(zp, xp, mp, mn, scratch) == 0);

  mpn_copyi(xp, mp, mn);

  ASSERT(mpn_invert_n(zp, xp, mp, mn, scratch) == 0);

  for (i = 0; i < 100; i++) {
    do {
      mpn_random(xp, mn * 2, rng, arg);
      mpn_mod(xp, xp, mn * 2, mp, mn);
    } while (mpn_zero_p(xp, mn));

    ASSERT(mpn_invert_n(zp, xp, mp, mn, scratch) == 1);

    mpn_mul_n(sp, xp, zp, mn);
    mpn_mod(sp, sp, mn * 2, mp, mn);

    ASSERT(mpn_strip(sp, mn) == 1);
    ASSERT(sp[0] == 1);
  }
}

static void
test_mpn_jacobi(void) {
  mp_limb_t scratch[MPN_JACOBI_ITCH(1)];
  mp_limb_t x, y;
  const int *v;
  size_t i;
  int j;

  printf("  - MPN jacobi.\n");

  for (i = 0; i < ARRAY_SIZE(jacobi_vectors); i++) {
    v = jacobi_vectors[i];
    x = MP_ABS(v[0]);
    y = MP_ABS(v[1]);

    if (x >= y)
      x %= y;

    if (v[0] < 0)
      x = y - x;

    j = mpn_jacobi_n(&x, &y, 1, scratch);

    if ((v[0] & v[1]) < 0)
      j = -j;

    ASSERT(j == v[2]);
  }
}

static void
test_mpn_powm(mp_rng_f *rng, void *arg) {
  mp_limb_t scratch1[MPN_JACOBI_ITCH(MP_P192_LIMBS)];
  mp_limb_t scratch2[MPN_POWM_ITCH(MP_P192_LIMBS, MP_P192_LIMBS)];
  mp_limb_t scratch3[MPN_SEC_POWM_ITCH(MP_P192_LIMBS)];
  mp_limb_t xp[MP_P192_LIMBS * 2];
  mp_limb_t yp[MP_P192_LIMBS];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_limb_t zp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_size_t mn = MP_P192_LIMBS;
  mp_size_t xn, yn;
  int i = 0;
  int k = 0;
  int j;

  printf("  - MPN powm.\n");

  mpn_p192_mod(mp);
  mpn_p192_exp(yp);

  mpn_zero(xp, mn);

  ASSERT(mpn_jacobi_n(xp, mp, mn, scratch1) == 0);

  mpn_copyi(xp, mp, mn);

  ASSERT(mpn_jacobi_n(xp, mp, mn, scratch1) == 0);

  while (i < 100 || k < 100) {
    mpn_random(xp, mn * 2, rng, arg);
    mpn_mod(xp, xp, mn * 2, mp, mn);

    xn = mpn_strip(xp, mn);

    j = mpn_jacobi_n(xp, mp, mn, scratch1);

    if (j == 0) {
      ASSERT(mpn_zero_p(xp, mn));
      continue;
    }

    ASSERT(j == 1 || j == -1);

    mpn_powm(zp, xp, xn, yp, mn, mp, mn, scratch2);
    mpn_sqr(sp, zp, mn, scratch1);
    mpn_mod(zp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(zp, xp, mn) == 0) == (j == 1));

    mpn_sec_powm(zp, xp, xn, yp, mn, mp, mn, scratch3);
    mpn_sqr(sp, zp, mn, scratch1);
    mpn_mod(zp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(zp, xp, mn) == 0) == (j == 1));

    mpn_mont_powm(zp, xp, xn, yp, mn, mp, mn, scratch2);
    mpn_sqr(sp, zp, mn, scratch1);
    mpn_mod(zp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(zp, xp, mn) == 0) == (j == 1));

    mpn_div_powm(zp, xp, xn, yp, mn, mp, mn, scratch2);
    mpn_sqr(sp, zp, mn, scratch1);
    mpn_mod(zp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(zp, xp, mn) == 0) == (j == 1));

    i += (j == 1);
    k += (j == -1);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, MP_P192_LIMBS, rng, arg);
    mpn_random_nz(yp, MP_P192_LIMBS, rng, arg);
    mpn_random_nz(mp, MP_P192_LIMBS, rng, arg);

    xn = mpn_strip(xp, MP_P192_LIMBS);
    yn = mpn_strip(yp, MP_P192_LIMBS);
    mn = mpn_strip(mp, MP_P192_LIMBS);

    mpn_powm(zp, xp, xn, yp, yn, mp, mn, scratch2);
    mpn_powm_simple(sp, xp, xn, yp, yn, mp, mn);

    ASSERT(mpn_cmp(zp, sp, mn) == 0);

    if (mp[0] & 1) {
      mpn_sec_powm(zp, xp, xn, yp, yn, mp, mn, scratch3);
      mpn_powm_simple(sp, xp, xn, yp, yn, mp, mn);

      ASSERT(mpn_cmp(zp, sp, mn) == 0);

      mpn_mont_powm(zp, xp, xn, yp, yn, mp, mn, scratch2);
      mpn_powm_simple(sp, xp, xn, yp, yn, mp, mn);

      ASSERT(mpn_cmp(zp, sp, mn) == 0);
    }

    mpn_div_powm(zp, xp, xn, yp, yn, mp, mn, scratch2);
    mpn_powm_simple(zp, xp, xn, yp, yn, mp, mn);

    ASSERT(mpn_cmp(zp, sp, mn) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random(xp, MP_P192_LIMBS, rng, arg);
    mpn_random_nz(yp, 1, rng, arg);
    mpn_random_nz(mp, MP_P192_LIMBS, rng, arg);

    xn = mpn_strip(xp, MP_P192_LIMBS);
    yn = mpn_strip(yp, 1);
    mn = mpn_strip(mp, MP_P192_LIMBS);

    mpn_powm(zp, xp, xn, yp, yn, mp, mn, scratch2);
    mpn_powm_simple(sp, xp, xn, yp, yn, mp, mn);

    ASSERT(mpn_cmp(zp, sp, mn) == 0);

    if (mp[0] & 1) {
      mpn_sec_powm(zp, xp, xn, yp, yn, mp, mn, scratch3);
      mpn_powm_simple(sp, xp, xn, yp, yn, mp, mn);

      ASSERT(mpn_cmp(zp, sp, mn) == 0);

      mpn_mont_powm(zp, xp, xn, yp, yn, mp, mn, scratch2);
      mpn_powm_simple(sp, xp, xn, yp, yn, mp, mn);

      ASSERT(mpn_cmp(zp, sp, mn) == 0);
    }

    mpn_div_powm(zp, xp, xn, yp, yn, mp, mn, scratch2);
    mpn_powm_simple(zp, xp, xn, yp, yn, mp, mn);

    ASSERT(mpn_cmp(zp, sp, mn) == 0);
  }
}

static void
test_mpn_helpers(void) {
  static const mp_limb_t trail[4] = {4, 3, 2, 0};
  static const mp_limb_t odd[4] = {3, 3, 2, 1};
  static const mp_limb_t even[4] = {4, 3, 2, 1};
  static const mp_limb_t tz[4] = {0, 3, 2, 1};

  printf("  - MPN helpers.\n");

  ASSERT(mpn_strip(trail, 4) == 3);
  ASSERT(mpn_odd_p(odd, 4));
  ASSERT(!mpn_odd_p(even, 4));
  ASSERT(mpn_even_p(even, 4));
  ASSERT(!mpn_even_p(odd, 4));
  ASSERT(mpn_ctz(tz, 4) == MP_LIMB_BITS);
  ASSERT(mpn_bitlen(trail, 4) == 2 * MP_LIMB_BITS + 2);
  ASSERT(mpn_bytelen(trail, 4) == ((size_t)mpn_bitlen(trail, 4) + 7) / 8);
  ASSERT(mpn_sizeinbase(trail, 4, 256) == (size_t)mpn_bytelen(trail, 4));
}

static void
test_mpn_select(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4], zp[4];

  printf("  - MPN select.\n");

  mpn_random_nz(xp, 4, rng, arg);
  mpn_random_nz(yp, 4, rng, arg);

  mpn_select(zp, xp, yp, 4, 0);

  ASSERT(mpn_cmp(zp, xp, 4) == 0);
  ASSERT(mpn_cmp(zp, yp, 4) != 0);

  mpn_select(zp, xp, yp, 4, 1);

  ASSERT(mpn_cmp(zp, yp, 4) == 0);
  ASSERT(mpn_cmp(zp, xp, 4) != 0);

  mpn_select_zero(zp, xp, 4, 0);

  ASSERT(mpn_cmp(zp, xp, 4) == 0);
  ASSERT(!mpn_zero_p(zp, 4));

  mpn_select_zero(zp, xp, 4, 1);

  ASSERT(mpn_cmp(zp, xp, 4) != 0);
  ASSERT(mpn_zero_p(zp, 4));
}

static void
test_mpn_sec_cmp(void) {
  static const mp_limb_t zero[4] = {0, 0, 0, 0};
  static const mp_limb_t mod[4] = {4, 3, 2, 1};
  static const mp_limb_t minus1[4] = {3, 3, 2, 1};
  static const mp_limb_t plus1[4] = {5, 3, 2, 1};
  static const mp_limb_t full[4] = {0xff, 0xff, 0xff, 0xff};

  printf("  - MPN secure comparison.\n");

  ASSERT(mpn_sec_zero_p(zero, 4));
  ASSERT(!mpn_sec_zero_p(mod, 4));

  ASSERT(!mpn_sec_equal_p(mod, minus1, 4));
  ASSERT(mpn_sec_equal_p(mod, mod, 4));

  ASSERT(mpn_sec_lt_p(minus1, mod, 4));
  ASSERT(!mpn_sec_lt_p(mod, mod, 4));
  ASSERT(!mpn_sec_lt_p(plus1, mod, 4));
  ASSERT(mpn_sec_lt_p(mod, full, 4));
  ASSERT(!mpn_sec_lt_p(full, mod, 4));

  ASSERT(mpn_sec_lte_p(minus1, mod, 4));
  ASSERT(mpn_sec_lte_p(mod, mod, 4));
  ASSERT(!mpn_sec_lte_p(plus1, mod, 4));
  ASSERT(mpn_sec_lte_p(mod, full, 4));
  ASSERT(!mpn_sec_lte_p(full, mod, 4));

  ASSERT(!mpn_sec_gt_p(minus1, mod, 4));
  ASSERT(!mpn_sec_gt_p(mod, mod, 4));
  ASSERT(mpn_sec_gt_p(plus1, mod, 4));
  ASSERT(!mpn_sec_gt_p(mod, full, 4));
  ASSERT(mpn_sec_gt_p(full, mod, 4));

  ASSERT(!mpn_sec_gte_p(minus1, mod, 4));
  ASSERT(mpn_sec_gte_p(mod, mod, 4));
  ASSERT(mpn_sec_gte_p(plus1, mod, 4));
  ASSERT(!mpn_sec_gte_p(mod, full, 4));
  ASSERT(mpn_sec_gte_p(full, mod, 4));

  ASSERT(mpn_sec_cmp(minus1, mod, 4) == -1);
  ASSERT(mpn_sec_cmp(mod, mod, 4) == 0);
  ASSERT(mpn_sec_cmp(plus1, mod, 4) == 1);
  ASSERT(mpn_sec_cmp(mod, full, 4) == -1);
  ASSERT(mpn_sec_cmp(full, mod, 4) == 1);
}

static void
test_mpn_sec_cmp_rand(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[4], yp[4];
  int i;

  printf("  - MPN secure comparison (random).\n");

  for (i = 0; i < 1000; i++) {
    mpn_random(xp, 4, rng, arg);
    mpn_random(yp, 4, rng, arg);

    ASSERT(mpn_sec_zero_p(xp, 4) == mpn_zero_p(xp, 4));
    ASSERT(mpn_sec_zero_p(yp, 4) == mpn_zero_p(yp, 4));

    ASSERT(mpn_sec_equal_p(xp, xp, 4));
    ASSERT(mpn_sec_equal_p(xp, yp, 4) == (mpn_cmp(xp, yp, 4) == 0));

    ASSERT(mpn_sec_lt_p(xp, xp, 4) == 0);
    ASSERT(mpn_sec_lt_p(yp, yp, 4) == 0);
    ASSERT(mpn_sec_lt_p(xp, yp, 4) == (mpn_cmp(xp, yp, 4) < 0));
    ASSERT(mpn_sec_lt_p(yp, xp, 4) == (mpn_cmp(yp, xp, 4) < 0));

    ASSERT(mpn_sec_lte_p(xp, xp, 4) == 1);
    ASSERT(mpn_sec_lte_p(yp, yp, 4) == 1);
    ASSERT(mpn_sec_lte_p(xp, yp, 4) == (mpn_cmp(xp, yp, 4) <= 0));
    ASSERT(mpn_sec_lte_p(yp, xp, 4) == (mpn_cmp(yp, xp, 4) <= 0));

    ASSERT(mpn_sec_gt_p(xp, xp, 4) == 0);
    ASSERT(mpn_sec_gt_p(yp, yp, 4) == 0);
    ASSERT(mpn_sec_gt_p(xp, yp, 4) == (mpn_cmp(xp, yp, 4) > 0));
    ASSERT(mpn_sec_gt_p(yp, xp, 4) == (mpn_cmp(yp, xp, 4) > 0));

    ASSERT(mpn_sec_gte_p(xp, xp, 4) == 1);
    ASSERT(mpn_sec_gte_p(yp, yp, 4) == 1);
    ASSERT(mpn_sec_gte_p(xp, yp, 4) == (mpn_cmp(xp, yp, 4) >= 0));
    ASSERT(mpn_sec_gte_p(yp, xp, 4) == (mpn_cmp(yp, xp, 4) >= 0));

    ASSERT(mpn_sec_cmp(xp, xp, 4) == 0);
    ASSERT(mpn_sec_cmp(yp, yp, 4) == 0);
    ASSERT(mpn_sec_cmp(xp, yp, 4) == mpn_cmp(xp, yp, 4));
    ASSERT(mpn_sec_cmp(yp, xp, 4) == mpn_cmp(yp, xp, 4));
  }
}

static void
test_mpn_io(mp_rng_f *rng, void *arg) {
  unsigned char raw[4 * MP_LIMB_BYTES];
  mp_limb_t xp[4], yp[4];
  int i;

  printf("  - MPN I/O.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    mpn_zero(yp, 4);
    mpn_export(raw, sizeof(raw), xp, 4, 1);
    mpn_import(yp, 4, raw, sizeof(raw), 1);

    ASSERT(mpn_cmp(xp, yp, 4) == 0);

    mpn_zero(yp, 4);
    mpn_export(raw, sizeof(raw), xp, 4, -1);
    mpn_import(yp, 4, raw, sizeof(raw), -1);

    ASSERT(mpn_cmp(xp, yp, 4) == 0);
  }
}

static void
test_mpn_io_str(mp_rng_f *rng, void *arg) {
  /* Base-2 size = 256 + 1 */
  /* Base-10 size = 78 + 1 */
  /* Base-16 size = 64 + 1 */
  mp_limb_t xp[4], yp[4];
  char str[258];
  int i;

  printf("  - MPN string I/O.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 4, rng, arg);

    mpn_zero(yp, 4);
    mpn_get_str(str, xp, 4, 10);

    ASSERT(mpn_set_str(yp, 4, str, 10));
    ASSERT(mpn_cmp(xp, yp, 4) == 0);

    mpn_zero(yp, 4);
    mpn_get_str(str, xp, 4, 16);

    ASSERT(mpn_set_str(yp, 4, str, 16));
    ASSERT(mpn_cmp(xp, yp, 4) == 0);

    mpn_print(xp, 4, 10, fake_puts);
    mpn_print(xp, 4, 16, fake_puts);
  }

  for (i = 2; i <= 62; i++) {
    mpn_random(xp, 4, rng, arg);

    mpn_zero(yp, 4);
    mpn_get_str(str, xp, 4, i);

    ASSERT(mpn_set_str(yp, 4, str, i));
    ASSERT(mpn_cmp(xp, yp, 4) == 0);

    mpn_print(xp, 4, i, fake_puts);
  }
}

static void
test_mpn_random(mp_rng_f *rng, void *arg) {
  mp_limb_t xp[8];
  int i;

  printf("  - MPN RNG.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(xp, 8, rng, arg);

    ASSERT(!mpn_zero_p(xp, 8));
  }
}

/*
 * MPZ
 */

static void
test_mpz_init(mp_rng_f *rng, void *arg) {
  static const mp_limb_t yp[1] = {1};
  static const mpz_t y = MPZ_ROINIT_N(yp, 1);
  mpz_t x;

  printf("  - MPZ init.\n");

  mpz_init(x);

  ASSERT(mpz_sgn(x) == 0);

  mpz_random_nz(x, 256, rng, arg);

  ASSERT(mpz_sgn(x) != 0);

  mpz_set_ui(x, 0);

  ASSERT(mpz_sgn(x) == 0);

  mpz_clear(x);

  mpz_init(x);
  mpz_cleanse(x);

  mpz_init_set(x, y);
  ASSERT(mpz_cmp_ui(x, 1) == 0);
  mpz_clear(x);

  mpz_init_set_ui(x, 2);
  ASSERT(mpz_cmp_ui(x, 2) == 0);
  mpz_clear(x);

  mpz_init_set_si(x, -3);
  ASSERT(mpz_cmp_si(x, -3) == 0);
  mpz_clear(x);

  ASSERT(mpz_init_set_str(x, "4", 10));
  ASSERT(mpz_cmp_ui(x, 4) == 0);
  mpz_clear(x);

  mpz_init2(x, 10 * MP_LIMB_BITS);

  ASSERT(x->alloc == 10);

  mpz_clear(x);
}

static void
test_mpz_assign(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;

  printf("  - MPZ assignment.\n");

  mpz_init(x);
  mpz_init(y);

  mpz_random_nz(x, 256, rng, arg);
  mpz_set_ui(y, 0);

  ASSERT(mpz_cmp(x, y) != 0);

  mpz_set(y, x);

  ASSERT(mpz_cmp(x, y) == 0);

  mpz_set_ui(y, 1);

  ASSERT(mpz_cmp_ui(y, 0) != 0);
  ASSERT(mpz_cmp_ui(y, 1) == 0);
  ASSERT(mpz_cmpabs_ui(y, 0) != 0);
  ASSERT(mpz_cmpabs_ui(y, 1) == 0);

  mpz_set_si(y, -1);

  ASSERT(mpz_cmp_si(y, 0) != 0);
  ASSERT(mpz_cmp_si(y, 1) != 0);
  ASSERT(mpz_cmp_si(y, -1) == 0);
  ASSERT(mpz_cmpabs_si(y, 1) == 0);
  ASSERT(mpz_cmpabs_si(y, -1) == 0);

  mpz_roset(z, x);

  ASSERT(mpz_cmp(x, z) == 0);

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_convert(void) {
  mpz_t x;

  printf("  - MPZ conversion.\n");

  mpz_init(x);
  mpz_set_ui(x, 0x01020304);

  ASSERT(mpz_get_ui(x) == 0x01020304);

  mpz_set_si(x, -0x01020304);

  ASSERT(mpz_get_si(x) == -0x01020304);

  mpz_clear(x);
}

static void
test_mpz_cmp(void) {
  mpz_t zero, mod, minus1, plus1, full;

  printf("  - MPZ comparison.\n");

  mpz_init(zero);
  mpz_init(mod);
  mpz_init(minus1);
  mpz_init(plus1);
  mpz_init(full);

  mpz_set_ui(zero, 0x00000000);
  mpz_set_ui(mod, 0x01020304);
  mpz_set_ui(minus1, 0x01020303);
  mpz_set_ui(plus1, 0x01020305);
  mpz_set_ui(full, 0xffffffff);

  ASSERT(mpz_sgn(zero) == 0);
  ASSERT(mpz_sgn(mod) != 0);

  ASSERT(mpz_cmp(minus1, mod) == -1);
  ASSERT(mpz_cmp(mod, mod) == 0);
  ASSERT(mpz_cmp(plus1, mod) == 1);
  ASSERT(mpz_cmp(mod, full) == -1);
  ASSERT(mpz_cmp(full, mod) == 1);

  ASSERT(mpz_cmp_ui(minus1, 0x01020304) == -1);
  ASSERT(mpz_cmp_ui(mod, 0x01020304) == 0);
  ASSERT(mpz_cmp_ui(plus1, 0x01020304) == 1);
  ASSERT(mpz_cmp_ui(mod, 0xffffffff) == -1);
  ASSERT(mpz_cmp_ui(full, 0x01020304) == 1);

  ASSERT(mpz_cmpabs(minus1, mod) == -1);
  ASSERT(mpz_cmpabs(mod, mod) == 0);
  ASSERT(mpz_cmpabs(plus1, mod) == 1);
  ASSERT(mpz_cmpabs(mod, full) == -1);
  ASSERT(mpz_cmpabs(full, mod) == 1);

  ASSERT(mpz_cmpabs_ui(minus1, 0x01020304) == -1);
  ASSERT(mpz_cmpabs_ui(mod, 0x01020304) == 0);
  ASSERT(mpz_cmpabs_ui(plus1, 0x01020304) == 1);
  ASSERT(mpz_cmpabs_ui(mod, 0xffffffff) == -1);
  ASSERT(mpz_cmpabs_ui(full, 0x01020304) == 1);

  mpz_clear(zero);
  mpz_clear(mod);
  mpz_clear(minus1);
  mpz_clear(plus1);
  mpz_clear(full);
}

static void
test_mpz_cmp_signed(void) {
  mpz_t zero, mod, minus1, plus1, full;

  printf("  - MPZ comparison (signed).\n");

  mpz_init(zero);
  mpz_init(mod);
  mpz_init(minus1);
  mpz_init(plus1);
  mpz_init(full);

  mpz_set_si(zero, 0x00000000);
  mpz_set_si(mod, -0x01020304);
  mpz_set_si(minus1, -0x01020303);
  mpz_set_si(plus1, -0x01020305);
  mpz_set_ui(full, 0xffffffff);

  ASSERT(mpz_sgn(zero) == 0);
  ASSERT(mpz_sgn(mod) != 0);

  ASSERT(mpz_cmp(minus1, mod) == 1);
  ASSERT(mpz_cmp(mod, mod) == 0);
  ASSERT(mpz_cmp(plus1, mod) == -1);
  ASSERT(mpz_cmp(mod, full) == -1);
  ASSERT(mpz_cmp(full, mod) == 1);

  ASSERT(mpz_cmp_si(minus1, -0x01020304) == 1);
  ASSERT(mpz_cmp_si(mod, -0x01020304) == 0);
  ASSERT(mpz_cmp_si(plus1, -0x01020304) == -1);
  ASSERT(mpz_cmp_ui(mod, 0xffffffff) == -1);
  ASSERT(mpz_cmp_si(full, -0x01020304) == 1);

  ASSERT(mpz_cmpabs(minus1, mod) == -1);
  ASSERT(mpz_cmpabs(mod, mod) == 0);
  ASSERT(mpz_cmpabs(plus1, mod) == 1);
  ASSERT(mpz_cmpabs(mod, full) == -1);
  ASSERT(mpz_cmpabs(full, mod) == 1);

  ASSERT(mpz_cmpabs_si(minus1, 0x01020304) == -1);
  ASSERT(mpz_cmpabs_si(mod, 0x01020304) == 0);
  ASSERT(mpz_cmpabs_si(plus1, 0x01020304) == 1);
  ASSERT(mpz_cmpabs_ui(mod, 0xffffffff) == -1);
  ASSERT(mpz_cmpabs_si(full, 0x01020304) == 1);

  mpz_clear(zero);
  mpz_clear(mod);
  mpz_clear(minus1);
  mpz_clear(plus1);
  mpz_clear(full);
}

static void
test_mpz_addsub(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ add/sub.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_add(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_add(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_sub(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_sub(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_sub(z, z, x);

    ASSERT(mpz_sgn(z) == 0);

    mpz_sub(z, z, x);

    ASSERT(mpz_sgn(z) == -1);
    ASSERT(mpz_cmp(z, x) == -1);
    ASSERT(mpz_cmpabs(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_add(z, x, y);
    mpz_add(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp(z, y) != 0);

    mpz_sub(z, z, y);
    mpz_sub(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_sub(z, z, x);

    ASSERT(mpz_sgn(z) == 0);

    mpz_sub(z, z, x);

    ASSERT(mpz_cmpabs(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_addsub_ui(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_limb_t y;
  int i;

  printf("  - MPZ add/sub (1 limb).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    mpz_add_ui(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_add_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_sub_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_sub_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_add_ui(z, x, y);
    mpz_add_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_ui(z, y) != 0);

    mpz_sub_ui(z, z, y);
    mpz_sub_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_addsub_si(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_long_t y;
  int i;

  printf("  - MPZ add/sub (1 limb, signed).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long_nz(rng, arg);

    mpz_add_si(z, x, y);
    mpz_add_si(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_si(z, y) != 0);

    mpz_sub_si(z, z, y);
    mpz_sub_si(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_ui_sub(mp_rng_f *rng, void *arg) {
  mpz_t y, z;
  mp_limb_t x;
  int i;

  printf("  - MPZ sub (1 limb).\n");

  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpn_random_nz(&x, 1, rng, arg);
    mpz_random_nz(y, 250, rng, arg);

    mpz_ui_sub(z, x, y);

    ASSERT(mpz_sgn(z) < 0);
    ASSERT(mpz_cmp(z, y) < 0);

    mpz_add(z, z, y);

    ASSERT(mpz_cmp_ui(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpn_random_nz(&x, 1, rng, arg);
    mpz_random_nz(y, 250, rng, arg);

    if (i & 1)
      mpz_neg(y, y);

    mpz_ui_sub(z, x, y);

    ASSERT(mpz_cmp(z, y) != 0);

    mpz_add(z, z, y);

    ASSERT(mpz_cmp_ui(z, x) == 0);
  }

  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_si_sub(mp_rng_f *rng, void *arg) {
  mpz_t y, z;
  mp_long_t x;
  int i;

  printf("  - MPZ sub (1 limb, signed).\n");

  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    x = mp_random_long_nz(rng, arg);

    mpz_random_nz(y, 250, rng, arg);

    if (i & 1)
      mpz_neg(y, y);

    mpz_si_sub(z, x, y);

    ASSERT(mpz_cmp(z, y) != 0);

    mpz_add(z, z, y);

    ASSERT(mpz_cmp_si(z, x) == 0);
  }

  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_mulquorem(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t;
  int i;

  printf("  - MPZ mul/quorem.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_mul(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_quorem(z, t, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_quorem(z, t, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_quorem(z, t, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_mul(z, x, y);
    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp(z, y) != 0);

    mpz_quorem(z, t, z, y);
    mpz_quorem(z, t, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_quorem(z, t, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
test_mpz_mulquo(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ mul/quo.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_mul(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_quo(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_quo(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_quo(z, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_mul(z, x, y);
    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp(z, y) != 0);

    mpz_quo(z, z, y);
    mpz_quo(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_quo(z, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_mulquo_ui(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_limb_t y;
  int i;

  printf("  - MPZ mul/quo (1 limb).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    mpz_mul_ui(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_mul_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_quo_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_quo_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul_ui(z, x, y);
    mpz_mul_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_ui(z, y) != 0);

    mpz_quo_ui(z, z, y);
    mpz_quo_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_mulquo_si(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_long_t y;
  int i;

  printf("  - MPZ mul/quo (1 limb, signed).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long_nz(rng, arg);

    mpz_mul_si(z, x, y);
    mpz_mul_si(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_si(z, y) != 0);

    mpz_quo_si(z, z, y);
    mpz_quo_si(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_muldivmod(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t;
  int i;

  printf("  - MPZ mul/divmod.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_mul(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_divmod(z, t, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_divmod(z, t, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_divmod(z, t, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_mul(z, x, y);
    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp(z, y) != 0);

    mpz_divmod(z, t, z, y);
    mpz_divmod(z, t, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_divmod(z, t, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
test_mpz_muldiv(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ mul/div.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_mul(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_div(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_div(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_div(z, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_mul(z, x, y);
    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp(z, y) != 0);

    mpz_div(z, z, y);
    mpz_div(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_div(z, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_muldiv_ui(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_limb_t y;
  int i;

  printf("  - MPZ mul/div (1 limb).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    mpz_mul_ui(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_mul_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_div_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_div_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_muldiv_si(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_long_t y;
  int i;

  printf("  - MPZ mul/div (1 limb, signed).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long_nz(rng, arg);

    mpz_mul_si(z, x, y);
    mpz_mul_si(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_si(z, y) != 0);

    mpz_div_si(z, z, y);
    mpz_div_si(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_muldivexact(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ mul/divexact.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_mul(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_divexact(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp(z, y) > 0);

    mpz_divexact(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_divexact(z, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_mul(z, x, y);
    mpz_mul(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp(z, y) != 0);

    mpz_divexact(z, z, y);
    mpz_divexact(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_divexact(z, z, x);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_muldivexact_ui(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  mp_limb_t y;
  int i;

  printf("  - MPZ mul/divexact (1 limb).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    mpz_mul_ui(z, x, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_mul_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_divexact_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) > 0);
    ASSERT(mpz_cmp_ui(z, y) > 0);

    mpz_divexact_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpn_random_nz(&y, 1, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul_ui(z, x, y);
    mpz_mul_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_ui(z, y) != 0);

    mpz_divexact_ui(z, z, y);
    mpz_divexact_ui(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_muldivexact_si(mp_rng_f *rng, void *arg) {
  mp_long_t y;
  mpz_t x, z;
  int i;

  printf("  - MPZ mul/divexact (1 limb, signed).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long_nz(rng, arg);

    mpz_mul_si(z, x, y);
    mpz_mul_si(z, z, y);

    ASSERT(mpz_cmp(z, x) != 0);
    ASSERT(mpz_cmp_si(z, y) != 0);

    mpz_divexact_si(z, z, y);
    mpz_divexact_si(z, z, y);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_divisibility(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;

  printf("  - MPZ divisibility.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  mpz_set_ui(x, 0);
  mpz_setbit(x, 100);

  ASSERT(mpz_divisible_2exp_p(x, 100));

  mpz_sub_ui(x, x, 1);

  ASSERT(!mpz_divisible_2exp_p(x, 100));

  mpz_add_ui(x, x, 2);

  ASSERT(!mpz_divisible_2exp_p(x, 100));

  mpz_random_nz(x, 250, rng, arg);
  mpz_random_nz(y, 250, rng, arg);

  mpz_mul(z, x, y);

  ASSERT(mpz_divisible_p(z, x));
  ASSERT(mpz_divisible_p(z, y));

  mpz_sub_ui(x, x, 1);
  mpz_sub_ui(y, y, 1);

  ASSERT(!mpz_divisible_p(z, x));
  ASSERT(!mpz_divisible_p(z, y));

  mpz_random_nz(x, 250, rng, arg);
  mpz_random_nz(y, 32, rng, arg);

  mpz_mul(z, x, y);

  ASSERT(mpz_divisible_ui_p(z, mpz_get_ui(y)));

  mpz_sub_ui(y, y, 1);

  ASSERT(!mpz_divisible_ui_p(z, mpz_get_ui(y)));

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_congruence(mp_rng_f *rng, void *arg) {
  mpz_t x, y, r1, r2, m;
  mp_bits_t bits;
  mp_limb_t d;
  int i;

  printf("  - MPZ congruence.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(r1);
  mpz_init(r2);
  mpz_init(m);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(r1, 250, rng, arg);
    mpz_random_nz(r2, 250, rng, arg);
    mpz_random_nz(m, 250, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(r1, r1);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(r2, r2);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(m, m);

    mpz_set(y, x);

    mpz_mul(r1, r1, m);
    mpz_mul(r2, r2, m);

    mpz_add(x, x, r1);
    mpz_add(y, y, r2);

    ASSERT(mpz_congruent_p(x, y, m));

    mpz_add_ui(x, x, 1);

    ASSERT(!mpz_congruent_p(x, y, m));
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(r1, 250, rng, arg);
    mpz_random_nz(r2, 250, rng, arg);

    d = mp_random_limb_nz(rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(r1, r1);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(r2, r2);

    mpz_set(y, x);

    mpz_mul_ui(r1, r1, d);
    mpz_mul_ui(r2, r2, d);

    mpz_add(x, x, r1);
    mpz_add(y, y, r2);

    ASSERT(mpz_congruent_ui_p(x, y, d));

    mpz_add_ui(x, x, 1);

    ASSERT(!mpz_congruent_ui_p(x, y, d));
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(r1, 250, rng, arg);
    mpz_random_nz(r2, 250, rng, arg);

    do {
      bits = mp_random_limb(rng, arg) & 0xff;
    } while (bits == 0);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(r1, r1);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(r2, r2);

    mpz_set(y, x);

    mpz_mul_2exp(r1, r1, bits);
    mpz_mul_2exp(r2, r2, bits);

    mpz_add(x, x, r1);
    mpz_add(y, y, r2);

    ASSERT(mpz_congruent_2exp_p(x, y, bits));

    mpz_add_ui(x, x, 1);

    ASSERT(!mpz_congruent_2exp_p(x, y, bits));
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(r1);
  mpz_clear(r2);
  mpz_clear(m);
}

static void
test_mpz_sqr(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ sqr.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul_simple(y, x, x);
    mpz_sqr(z, x);

    ASSERT(mpz_sgn(z) != -1);
    ASSERT(mpz_cmp(z, y) == 0);

    mpz_div(z, z, x);

    ASSERT(mpz_cmpabs(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_addmul(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t;
  int i;

  printf("  - MPZ add+mul.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);
    mpz_random_nz(z, 100, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(z, z);

    mpz_mul(t, x, y);
    mpz_add(t, z, t);

    mpz_addmul(z, x, y);

    ASSERT(mpz_cmp(z, t) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS, rng, arg);
    mpz_random_nz(z, 100, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(z, z);

    mpz_mul(t, x, y);
    mpz_add(t, z, t);

    mpz_addmul_ui(z, x, mpz_get_ui(y));

    ASSERT(mpz_cmp(z, t) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS - 1, rng, arg);
    mpz_random_nz(z, 100, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(z, z);

    mpz_mul(t, x, y);
    mpz_add(t, z, t);

    mpz_addmul_si(z, x, mpz_get_si(y));

    ASSERT(mpz_cmp(z, t) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
test_mpz_submul(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t;
  int i;

  printf("  - MPZ sub+mul.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);
    mpz_random_nz(z, 100, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(z, z);

    mpz_mul(t, x, y);
    mpz_sub(t, z, t);

    mpz_submul(z, x, y);

    ASSERT(mpz_cmp(z, t) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS, rng, arg);
    mpz_random_nz(z, 100, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(z, z);

    mpz_mul(t, x, y);
    mpz_sub(t, z, t);

    mpz_submul_ui(z, x, mpz_get_ui(y));

    ASSERT(mpz_cmp(z, t) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS - 1, rng, arg);
    mpz_random_nz(z, 100, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(z, z);

    mpz_mul(t, x, y);
    mpz_sub(t, z, t);

    mpz_submul_si(z, x, mpz_get_si(y));

    ASSERT(mpz_cmp(z, t) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
test_mpz_rem(mp_rng_f *rng, void *arg) {
  mpz_t n, d, q, r, t;
  int i;

  printf("  - MPZ rem.\n");

  mpz_init(n);
  mpz_init(d);
  mpz_init(q);
  mpz_init(r);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(n, 250, rng, arg);
    mpz_random_nz(d, 125, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(n, n);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(d, d);

    mpz_quo(q, n, d);
    mpz_rem(r, n, d);

    mpz_mul(t, q, d);
    mpz_sub(t, n, t);

    ASSERT(mpz_cmp(r, t) == 0);
  }

  mpz_clear(n);
  mpz_clear(d);
  mpz_clear(q);
  mpz_clear(r);
  mpz_clear(t);
}

static void
test_mpz_rem_ui(mp_rng_f *rng, void *arg) {
  mpz_t n, q, t;
  mp_limb_t d, r;
  int i;

  printf("  - MPZ rem (1 limb).\n");

  mpz_init(n);
  mpz_init(q);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(n, 250, rng, arg);
    mpn_random_nz(&d, 1, rng, arg);

    mpz_quo_ui(q, n, d);

    r = mpz_rem_ui(n, d);

    mpz_mul_ui(t, q, d);
    mpz_sub(t, n, t);

    ASSERT(mpz_cmp_ui(t, r) == 0);
  }

  mpz_clear(n);
  mpz_clear(q);
  mpz_clear(t);
}

static void
test_mpz_rem_si(mp_rng_f *rng, void *arg) {
  mpz_t n, q, t;
  mp_long_t d, r;
  int i;

  printf("  - MPZ rem (1 limb, signed).\n");

  mpz_init(n);
  mpz_init(q);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(n, 250, rng, arg);

    if (i & 1)
      mpz_neg(n, n);

    d = mp_random_limb_nz(rng, arg);

    mpz_quo_si(q, n, d);

    r = mpz_rem_si(n, d);

    mpz_mul_si(t, q, d);
    mpz_sub(t, n, t);

    ASSERT(mpz_cmp_si(t, r) == 0);
  }

  mpz_clear(n);
  mpz_clear(q);
  mpz_clear(t);
}

static void
test_mpz_mod(mp_rng_f *rng, void *arg) {
  mpz_t n, d, q, r, t;
  int i;

  printf("  - MPZ mod.\n");

  mpz_init(n);
  mpz_init(d);
  mpz_init(q);
  mpz_init(r);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(n, 250, rng, arg);
    mpz_random_nz(d, 125, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(n, n);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(d, d);

    mpz_div(q, n, d);
    mpz_mod(r, n, d);

    mpz_mul(t, q, d);
    mpz_sub(t, n, t);

    ASSERT(mpz_cmp(r, t) == 0);
  }

  mpz_clear(n);
  mpz_clear(d);
  mpz_clear(q);
  mpz_clear(r);
  mpz_clear(t);
}

static void
test_mpz_mod_ui(mp_rng_f *rng, void *arg) {
  mpz_t n, q, t;
  mp_limb_t d, r;
  int i;

  printf("  - MPZ mod (1 limb).\n");

  mpz_init(n);
  mpz_init(q);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(n, 250, rng, arg);
    mpn_random_nz(&d, 1, rng, arg);

    mpz_div_ui(q, n, d);

    r = mpz_mod_ui(n, d);

    mpz_mul_ui(t, q, d);
    mpz_sub(t, n, t);

    ASSERT(mpz_cmp_ui(t, r) == 0);
  }

  mpz_clear(n);
  mpz_clear(q);
  mpz_clear(t);
}

static void
test_mpz_mod_si(mp_rng_f *rng, void *arg) {
  mpz_t n, q, t;
  mp_long_t d, r;
  int i;

  printf("  - MPZ mod (1 limb, signed).\n");

  mpz_init(n);
  mpz_init(q);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(n, 250, rng, arg);

    if (i & 1)
      mpz_neg(n, n);

    d = mp_random_long_nz(rng, arg);

    mpz_div_si(q, n, d);

    r = mpz_mod_si(n, d);

    mpz_mul_si(t, q, d);
    mpz_sub(t, n, t);

    ASSERT(mpz_cmp_si(t, r) == 0);
  }

  mpz_clear(n);
  mpz_clear(q);
  mpz_clear(t);
}

static void
test_mpz_divround(mp_rng_f *rng, void *arg) {
  mpz_t z, n, q, d, h;
  int i;

  printf("  - MPZ divround.\n");

  mpz_init(z);
  mpz_init(n);
  mpz_init(q);
  mpz_init(d);
  mpz_init(h);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(q, 250, rng, arg);
    mpz_random_nz(d, 125, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(q, q);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(d, d);

    mpz_quo_2exp(h, d, 1);

    mpz_mul(n, q, d);

    if (i & 1) {
      if (mpz_sgn(d) < 0)
        mpz_sub_ui(h, h, 1);
      else
        mpz_add_ui(h, h, 1);

      if (mpz_sgn(q) < 0)
        mpz_sub_ui(q, q, 1);
      else
        mpz_add_ui(q, q, 1);
    } else {
      if (mpz_sgn(d) < 0)
        mpz_add_ui(h, h, 1);
      else
        mpz_sub_ui(h, h, 1);
    }

    if ((mpz_sgn(n) ^ mpz_sgn(h)) < 0)
      mpz_sub(n, n, h);
    else
      mpz_add(n, n, h);

    mpz_divround(z, n, d);

    ASSERT(mpz_cmp(z, q) == 0);
  }

  mpz_clear(z);
  mpz_clear(n);
  mpz_clear(q);
  mpz_clear(d);
  mpz_clear(h);
}

static void
test_mpz_divround_ui(mp_rng_f *rng, void *arg) {
  mp_limb_t d, h;
  mpz_t z, n, q;
  int i;

  printf("  - MPZ divround (1 limb).\n");

  mpz_init(z);
  mpz_init(n);
  mpz_init(q);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(q, 250, rng, arg);
    mpn_random_nz(&d, 1, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(q, q);

    mpz_mul_ui(n, q, d);

    if (i & 1) {
      h = (d >> 1) + 1;

      if (mpz_sgn(q) < 0)
        mpz_sub_ui(q, q, 1);
      else
        mpz_add_ui(q, q, 1);
    } else {
      h = (d >> 1) - 1;
    }

    if (mpz_sgn(n) < 0)
      mpz_sub_ui(n, n, h);
    else
      mpz_add_ui(n, n, h);

    mpz_divround_ui(z, n, d);

    ASSERT(mpz_cmp(z, q) == 0);
  }

  mpz_clear(z);
  mpz_clear(n);
  mpz_clear(q);
}

static void
test_mpz_divround_si(mp_rng_f *rng, void *arg) {
  mp_long_t d, h;
  mpz_t z, n, q;
  int i;

  printf("  - MPZ divround (1 limb, signed).\n");

  mpz_init(z);
  mpz_init(n);
  mpz_init(q);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(q, 250, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(q, q);

    d = mp_random_long_nz(rng, arg);

    mpz_mul_si(n, q, d);

    if (i & 1) {
      if (d < 0)
        h = (d / 2) - 1;
      else
        h = (d / 2) + 1;

      if (mpz_sgn(q) < 0)
        mpz_sub_ui(q, q, 1);
      else
        mpz_add_ui(q, q, 1);
    } else {
      if (d < 0)
        h = (d / 2) + 1;
      else
        h = (d / 2) - 1;
    }

    if ((mpz_sgn(n) ^ h) < 0)
      mpz_sub_si(n, n, h);
    else
      mpz_add_si(n, n, h);

    mpz_divround_si(z, n, d);

    ASSERT(mpz_cmp(z, q) == 0);
  }

  mpz_clear(z);
  mpz_clear(n);
  mpz_clear(q);
}

static void
test_mpz_pow(mp_rng_f *rng, void *arg) {
  mpz_t x, z1, z2;
  mp_limb_t y;
  int i;

  printf("  - MPZ pow.\n");

  mpz_init(x);
  mpz_init(z1);
  mpz_init(z2);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 32, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_pow_ui(z1, x, 5);
    mpz_pow5(z2, x);

    ASSERT(mpz_cmp(z1, z2) == 0);

    mpz_pow_simple_ui(z2, x, 5);

    ASSERT(mpz_cmp(z1, z2) == 0);
  }

  for (i = 0; i < 100; i++) {
    y = mp_random_limb(rng, arg);

    mpz_ui_pow_ui(z1, y, 5);

    mpz_set_ui(z2, y);
    mpz_pow5(z2, z2);

    ASSERT(mpz_cmp(z1, z2) == 0);

    mpz_set_ui(z2, y);
    mpz_pow_simple_ui(z2, z2, 5);

    ASSERT(mpz_cmp(z1, z2) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 32, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_pow_ui(z1, x, 16);
    mpz_pow_simple_ui(z2, x, 16);

    ASSERT(mpz_cmp(z1, z2) == 0);
  }

  mpz_clear(x);
  mpz_clear(z1);
  mpz_clear(z2);
}

static void
test_mpz_roots(mp_rng_f *rng, void *arg) {
  static const int roots[][3] = {
    { 0, 0, 1 },
    { 1, 1, 1 },
    { 4, 2, 1 },
    { 5, 2, 0 },
    { 8, 2, 0 },
    { 9, 3, 1 },
    { 121, 11, 1 },
    { 122, 11, 0 },
    { 1024, 32, 1 },
    { 1025, 32, 0 }
  };

  mpz_t x, z, r;
  const int *v;
  size_t i;

  printf("  - MPZ roots.\n");

  mpz_init(x);
  mpz_init(z);
  mpz_init(r);

  for (i = 0; i < ARRAY_SIZE(roots); i++) {
    v = roots[i];

    mpz_set_ui(x, v[0]);
    mpz_sqrt(x, x);

    ASSERT(mpz_cmp_ui(x, v[1]) == 0);

    mpz_set_ui(x, v[0]);

    ASSERT(mpz_perfect_square_p(x) == v[2]);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_sqr(z, x);

    ASSERT(mpz_perfect_square_p(z));

    mpz_sqrt(z, z);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    mpz_sqrtrem(z, r, x);

    if (mpz_sgn(r) == 0)
      ASSERT(mpz_perfect_square_p(x));

    mpz_sqr(z, z);
    mpz_add(z, z, r);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 64, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_pow_ui(z, x, 5);

    ASSERT(mpz_root(z, z, 5));
    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_rootrem(z, r, x, 5);

    if (mpz_sgn(r) == 0)
      ASSERT(mpz_root(NULL, x, 5));

    mpz_pow_ui(z, z, 5);
    mpz_add(z, z, r);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  for (i = 0; i < 4; i++) {
    mpz_randprime(x, 125, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    ASSERT(!mpz_perfect_power_p(x));

    mpz_pow_ui(x, x, i + 2);

    ASSERT(mpz_perfect_power_p(x));
  }

  mpz_clear(x);
  mpz_clear(z);
  mpz_clear(r);
}

static void
test_mpz_and(mp_rng_f *rng, void *arg) {
  mpz_t x, m, z, e, r;
  mp_bits_t j;
  int i;

  printf("  - MPZ AND.\n");

  mpz_init(x);
  mpz_init(m);
  mpz_init(z);
  mpz_init(e);
  mpz_init(r);

  mpz_setbit(m, 100);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    mpz_set_ui(e, 0);

    for (j = 0; j < 100; j++) {
      if (mpz_tstbit(x, j))
        mpz_setbit(e, j);
    }

    mpz_and(z, x, m);

    ASSERT(mpz_cmp(z, e) == 0);

    ASSERT(mpz_and_ui(x, 13) == (mpz_getlimbn(x, 0) & 13));

    mpz_and_si(r, x, 13);

    ASSERT(mpz_get_si(r) == (mp_long_t)(mpz_getlimbn(x, 0) & 13));
  }

  mpz_clear(x);
  mpz_clear(m);
  mpz_clear(z);
  mpz_clear(e);
  mpz_clear(r);
}

static void
test_mpz_ior(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, r;
  int i;

  printf("  - MPZ OR.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(r);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 31, rng, arg);

    mpz_mul_2exp(z, x, 32);
    mpz_ior(z, z, y);

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_quo_2exp(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_mul_2exp(z, x, 32);
    mpz_ior_ui(z, z, mpz_get_ui(y));

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_div_2exp(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_mul_2exp(z, x, 32);
    mpz_ior_si(z, z, mpz_get_si(y));

    mpz_and_si(r, z, INT32_MAX);

    ASSERT(mpz_cmp(r, y) == 0);

    mpz_div_2exp(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(r);
}

static void
test_mpz_xor(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, r;
  int i;

  printf("  - MPZ XOR.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(r);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 31, rng, arg);

    mpz_mul_2exp(z, x, 32);
    mpz_xor(z, z, y);

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_div_2exp(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_mul_2exp(z, x, 32);
    mpz_xor_ui(z, z, mpz_get_ui(y));

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_div_2exp(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_mul_2exp(z, x, 32);
    mpz_xor_si(z, z, mpz_get_si(y));

    mpz_and_si(r, z, INT32_MAX);

    ASSERT(mpz_cmp(r, y) == 0);

    mpz_div_2exp(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(r);
}

static void
test_mpz_com(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  int i;

  printf("  - MPZ NOT.\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_com(z, x);
    mpz_com(z, z);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_lshift(mp_rng_f *rng, void *arg) {
  mpz_t x, y;
  int i;

  printf("  - MPZ lshift.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul_ui(y, x, 8);
    mpz_mul_2exp(x, x, 3);

    ASSERT(mpz_cmp(x, y) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_rshift(mp_rng_f *rng, void *arg) {
  mpz_t x, y;
  int i;

  printf("  - MPZ rshift.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1) {
      mpz_neg(x, x);
      mpz_div_ui(y, x, 8);
    } else {
      mpz_quo_ui(y, x, 8);
    }

    mpz_div_2exp(x, x, 3);

    ASSERT(mpz_cmp(x, y) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_shift(mp_rng_f *rng, void *arg) {
  mpz_t x, y;
  int i;

  printf("  - MPZ shift.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul_2exp(y, x, 173);

    ASSERT(mpz_cmpabs(y, x) > 0);

    mpz_div_2exp(y, y, 173);

    ASSERT(mpz_cmp(x, y) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_bits(mp_rng_f *rng, void *arg) {
  mp_bits_t bits = 4 * MP_LIMB_BITS;
  mp_bits_t j;
  mp_limb_t w;
  mpz_t x, y;
  int i, b;

  printf("  - MPZ bits.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, bits, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    for (j = 0; j < bits; j++) {
      if (mpz_sgn(x) < 0) {
        mpz_com(y, x);
        mpz_div_2exp(y, y, j);
      } else {
        mpz_div_2exp(y, x, j);
      }

      w = mpz_getlimbn(y, 0);

      if (mpz_sgn(x) < 0)
        w = ~w;

      ASSERT((mp_limb_t)mpz_tstbit(x, j) == (w & 1));
    }

    ASSERT(mpz_tstbit(x, j) == (mpz_sgn(x) < 0));
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, bits, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_set(y, x);

    for (j = 0; j < bits; j++) {
      b = mpz_tstbit(x, j);

      if (b)
        mpz_clrbit(y, j);
      else
        mpz_setbit(y, j);

      ASSERT(mpz_cmp(y, x) != 0);
      ASSERT(mpz_tstbit(y, j) == (b ^ 1));

      if (b)
        mpz_setbit(y, j);
      else
        mpz_clrbit(y, j);

      ASSERT(mpz_cmp(y, x) == 0);
      ASSERT(mpz_tstbit(y, j) == b);

      mpz_combit(y, j);

      ASSERT(mpz_cmp(y, x) != 0);
      ASSERT(mpz_tstbit(y, j) == (b ^ 1));

      mpz_combit(y, j);

      ASSERT(mpz_cmp(y, x) == 0);
      ASSERT(mpz_tstbit(y, j) == b);
    }
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_scan(void) {
  mpz_t x;

  printf("  - MPZ scan.\n");

  mpz_init(x);

  ASSERT(mpz_set_str(x, "10101000", 2));

  ASSERT(mpz_scan1(x, 0) == 3);
  ASSERT(mpz_scan1(x, 1) == 3);
  ASSERT(mpz_scan1(x, 2) == 3);
  ASSERT(mpz_scan1(x, 3) == 3);
  ASSERT(mpz_scan1(x, 4) == 5);
  ASSERT(mpz_scan1(x, 5) == 5);
  ASSERT(mpz_scan1(x, 6) == 7);
  ASSERT(mpz_scan1(x, 7) == 7);
  ASSERT(mpz_scan1(x, 8) == MP_BITS_MAX);

  ASSERT(mpz_scan0(x, 0) == 0);
  ASSERT(mpz_scan0(x, 1) == 1);
  ASSERT(mpz_scan0(x, 2) == 2);
  ASSERT(mpz_scan0(x, 3) == 4);
  ASSERT(mpz_scan0(x, 4) == 4);
  ASSERT(mpz_scan0(x, 5) == 6);
  ASSERT(mpz_scan0(x, 6) == 6);
  ASSERT(mpz_scan0(x, 7) == 8);
  ASSERT(mpz_scan0(x, 8) == 8);

  mpz_mul_2exp(x, x, 71);

  ASSERT(mpz_scan1(x, 0) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 0) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 1) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 2) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 3) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 4) == 71 + 5);
  ASSERT(mpz_scan1(x, 71 + 5) == 71 + 5);
  ASSERT(mpz_scan1(x, 71 + 6) == 71 + 7);
  ASSERT(mpz_scan1(x, 71 + 7) == 71 + 7);
  ASSERT(mpz_scan1(x, 71 + 8) == MP_BITS_MAX);

  ASSERT(mpz_scan0(x, 0) == 0);
  ASSERT(mpz_scan0(x, 71 + 0) == 71 + 0);
  ASSERT(mpz_scan0(x, 71 + 1) == 71 + 1);
  ASSERT(mpz_scan0(x, 71 + 2) == 71 + 2);
  ASSERT(mpz_scan0(x, 71 + 3) == 71 + 4);
  ASSERT(mpz_scan0(x, 71 + 4) == 71 + 4);
  ASSERT(mpz_scan0(x, 71 + 5) == 71 + 6);
  ASSERT(mpz_scan0(x, 71 + 6) == 71 + 6);
  ASSERT(mpz_scan0(x, 71 + 7) == 71 + 8);
  ASSERT(mpz_scan0(x, 71 + 8) == 71 + 8);

  ASSERT(mpz_set_str(x, "-10101000", 2)); /* 01011000 */

  ASSERT(mpz_scan1(x, 0) == 3);
  ASSERT(mpz_scan1(x, 1) == 3);
  ASSERT(mpz_scan1(x, 2) == 3);
  ASSERT(mpz_scan1(x, 3) == 3);
  ASSERT(mpz_scan1(x, 4) == 4);
  ASSERT(mpz_scan1(x, 5) == 6);
  ASSERT(mpz_scan1(x, 6) == 6);
  ASSERT(mpz_scan1(x, 7) == 8);
  ASSERT(mpz_scan1(x, 8) == 8);

  ASSERT(mpz_scan0(x, 0) == 0);
  ASSERT(mpz_scan0(x, 1) == 1);
  ASSERT(mpz_scan0(x, 2) == 2);
  ASSERT(mpz_scan0(x, 3) == 5);
  ASSERT(mpz_scan0(x, 4) == 5);
  ASSERT(mpz_scan0(x, 5) == 5);
  ASSERT(mpz_scan0(x, 6) == 7);
  ASSERT(mpz_scan0(x, 7) == 7);
  ASSERT(mpz_scan0(x, 8) == MP_BITS_MAX);

  mpz_mul_2exp(x, x, 71);

  ASSERT(mpz_scan1(x, 0) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 0) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 1) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 2) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 3) == 71 + 3);
  ASSERT(mpz_scan1(x, 71 + 4) == 71 + 4);
  ASSERT(mpz_scan1(x, 71 + 5) == 71 + 6);
  ASSERT(mpz_scan1(x, 71 + 6) == 71 + 6);
  ASSERT(mpz_scan1(x, 71 + 7) == 71 + 8);
  ASSERT(mpz_scan1(x, 71 + 8) == 71 + 8);

  ASSERT(mpz_scan0(x, 0) == 0);
  ASSERT(mpz_scan0(x, 71 + 0) == 71 + 0);
  ASSERT(mpz_scan0(x, 71 + 1) == 71 + 1);
  ASSERT(mpz_scan0(x, 71 + 2) == 71 + 2);
  ASSERT(mpz_scan0(x, 71 + 3) == 71 + 5);
  ASSERT(mpz_scan0(x, 71 + 4) == 71 + 5);
  ASSERT(mpz_scan0(x, 71 + 5) == 71 + 5);
  ASSERT(mpz_scan0(x, 71 + 6) == 71 + 7);
  ASSERT(mpz_scan0(x, 71 + 7) == 71 + 7);
  ASSERT(mpz_scan0(x, 71 + 8) == MP_BITS_MAX);

  mpz_clear(x);
}

static void
test_mpz_popcount(mp_rng_f *rng, void *arg) {
  mp_size_t j, n;
  mp_bits_t c;
  mpz_t x, y;
  int i;

  printf("  - MPZ hamming weight.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    n = mpz_size(x);
    c = 0;

    if (mpz_sgn(x) < 0) {
      c = MP_BITS_MAX;
    } else {
      for (j = 0; j < n; j++)
        c += mp_popcount_simple(mpz_getlimbn(x, j));
    }

    ASSERT(mpz_popcount(x) == c);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_hamdist(mp_rng_f *rng, void *arg) {
  mp_size_t xn, yn, n, j;
  mpz_t x, y, xc, yc;
  mp_limb_t xw, yw;
  mp_bits_t c;
  int i;

  printf("  - MPZ hamming distance.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(xc);
  mpz_init(yc);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 150, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    xn = mpz_size(x);
    yn = mpz_size(y);
    n = MP_MAX(xn, yn);
    c = 0;

    if ((mpz_sgn(x) & mpz_sgn(y)) < 0) {
      mpz_com(xc, x);
      mpz_com(yc, y);

      for (j = 0; j < n; j++) {
        xw = ~mpz_getlimbn(xc, j);
        yw = ~mpz_getlimbn(yc, j);
        c += mp_popcount_simple(xw ^ yw);
      }
    } else if ((mpz_sgn(x) | mpz_sgn(y)) < 0) {
      c = MP_BITS_MAX;
    } else {
      for (j = 0; j < n; j++) {
        xw = mpz_getlimbn(x, j);
        yw = mpz_getlimbn(y, j);
        c += mp_popcount_simple(xw ^ yw);
      }
    }

    ASSERT(mpz_hamdist(x, y) == c);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(xc);
  mpz_clear(yc);
}

static void
test_mpz_mask(mp_rng_f *rng, void *arg) {
  mpz_t x, m, z, e;
  int i;

  printf("  - MPZ mask.\n");

  mpz_init(x);
  mpz_init(m);
  mpz_init(z);
  mpz_init(e);

  mpz_set_ui(m, 0);
  mpz_setbit(m, 100);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    mpz_and(e, x, m);
    mpz_rem_2exp(z, x, 100);

    ASSERT(mpz_cmp(z, e) == 0);
  }

  mpz_set_ui(m, 0);
  mpz_setbit(m, 500);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    mpz_and(e, x, m);
    mpz_rem_2exp(z, x, 500);

    ASSERT(mpz_cmp(z, e) == 0);
  }

  mpz_set_ui(m, 0);
  mpz_setbit(m, 100);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_and(e, x, m);
    mpz_mod_2exp(z, x, 100);

    ASSERT(mpz_cmp(z, e) == 0);
  }

  mpz_set_ui(m, 0);
  mpz_setbit(m, 500);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_and(e, x, m);
    mpz_mod_2exp(z, x, 500);

    ASSERT(mpz_cmp(z, e) == 0);
  }

  mpz_clear(x);
  mpz_clear(m);
  mpz_clear(z);
  mpz_clear(e);
}

static void
test_mpz_negate(mp_rng_f *rng, void *arg) {
  mpz_t x, y;
  int i;

  printf("  - MPZ negation.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_abs(y, x);

    ASSERT(mpz_sgn(y) != -1);
    ASSERT(mpz_cmp(x, y) == 0);
    ASSERT(mpz_cmpabs(x, y) == 0);

    mpz_neg(y, y);

    ASSERT(mpz_sgn(y) == -1);
    ASSERT(mpz_cmp(x, y) != 0);
    ASSERT(mpz_cmpabs(x, y) == 0);

    mpz_abs(y, y);

    ASSERT(mpz_sgn(y) != -1);
    ASSERT(mpz_cmp(x, y) == 0);
    ASSERT(mpz_cmpabs(x, y) == 0);

    mpz_neg(y, y);

    ASSERT(mpz_sgn(y) == -1);

    mpz_neg(y, y);

    ASSERT(mpz_sgn(y) != -1);
    ASSERT(mpz_cmp(x, y) == 0);
    ASSERT(mpz_cmpabs(x, y) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_mulshift(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, m, q, d;
  int i;

  printf("  - MPZ mulshift.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(m);
  mpz_init(q);
  mpz_init(d);

  mpz_p192_mod(m);

  mpz_setbit(d, 196);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 384, rng, arg);
    mpz_random_nz(y, 384, rng, arg);

    mpz_mod(x, x, m);
    mpz_mod(y, y, m);

    mpz_mulshift(z, x, y, 196);

    mpz_mul(q, x, y);
    mpz_divround(q, q, d);

    ASSERT(mpz_cmp(z, q) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(m);
  mpz_clear(q);
  mpz_clear(d);
}

static void
test_mpz_gcd(mp_rng_f *rng, void *arg) {
  mpz_t g1, g2, g3, x, y;
  mp_limb_t g;
  int i;

  printf("  - MPZ gcd.\n");

  mpz_init(g1);
  mpz_init(g2);
  mpz_init(g3);
  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 50; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_gcd(g1, x, y);
    mpz_gcdext(g2, NULL, NULL, x, y);
    mpz_gcd_simple(g3, x, y);

    ASSERT(mpz_cmp(g1, g2) == 0);
    ASSERT(mpz_cmp(g1, g3) == 0);
  }

  for (i = 0; i < 50; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS, rng, arg);

    g = mpz_gcd_ui(g1, x, mpz_get_ui(y));

    mpz_gcdext(g2, NULL, NULL, x, y);
    mpz_gcd_simple(g3, x, y);

    ASSERT(mpz_cmp_ui(g1, g) == 0);
    ASSERT(mpz_cmp(g1, g2) == 0);
    ASSERT(mpz_cmp(g1, g3) == 0);
  }

  mpz_clear(g1);
  mpz_clear(g2);
  mpz_clear(g3);
  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_lcm(mp_rng_f *rng, void *arg) {
  mpz_t l1, l2, l3, g, x, y;
  int i;

  printf("  - MPZ lcm.\n");

  mpz_init(l1);
  mpz_init(l2);
  mpz_init(l3);
  mpz_init(g);
  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 50; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_lcm(l1, x, y);

    mpz_gcd(g, x, y);
    mpz_divexact(l2, x, g);
    mpz_mul(l2, l2, y);

    mpz_lcm_simple(l3, x, y);

    ASSERT(mpz_cmp(l1, l2) == 0);
    ASSERT(mpz_cmp(l1, l3) == 0);
  }

  for (i = 0; i < 50; i++) {
    mpz_random_nz(x, 250, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS, rng, arg);

    mpz_lcm_ui(l1, x, mpz_get_ui(y));

    mpz_gcd(g, x, y);
    mpz_divexact(l2, x, g);
    mpz_mul(l2, l2, y);

    mpz_lcm_simple(l3, x, y);

    ASSERT(mpz_cmp(l1, l2) == 0);
    ASSERT(mpz_cmp(l1, l3) == 0);
  }

  mpz_clear(l1);
  mpz_clear(l2);
  mpz_clear(l3);
  mpz_clear(g);
  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_gcdext(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, m;
  mpz_t g1, s1, t1;
  mpz_t g2, s2, t2;
  mpz_t q1, q2;
  int i, j, mode;

  printf("  - MPZ gcdext.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(m);
  mpz_init(g1);
  mpz_init(s1);
  mpz_init(t1);
  mpz_init(g2);
  mpz_init(s2);
  mpz_init(t2);
  mpz_init(q1);
  mpz_init(q2);

  mpz_p192_mod(m);

  for (i = 0; i < 100; i++) {
    do {
      mpz_random_nz(x, 384, rng, arg);
      mpz_mod(x, x, m);
    } while (mpz_sgn(x) == 0);

    mpz_gcdext(g1, z, NULL, x, m);

    ASSERT(mpz_cmp_ui(g1, 1) == 0);

    mpz_mod(z, z, m);
    mpz_mul(z, z, x);
    mpz_mod(z, z, m);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 10; i++) {
    mpz_set_ui(x, 0);
    mpz_random_nz(y, 192, rng, arg);

    if (i & 1)
      mpz_neg(y, y);

    mpz_gcdext(g1, s1, t1, x, y);

    ASSERT(mpz_cmpabs(g1, y) == 0 && mpz_sgn(g1) >= 0);
    ASSERT(mpz_sgn(s1) == 0);
    ASSERT(mpz_cmp_si(t1, mpz_sgn(y)) == 0);
  }

  for (i = 0; i < 10; i++) {
    mpz_random_nz(x, 192, rng, arg);
    mpz_set_ui(y, 0);

    if (i & 1)
      mpz_neg(x, x);

    mpz_gcdext(g1, s1, t1, x, y);

    ASSERT(mpz_cmpabs(g1, x) == 0 && mpz_sgn(g1) >= 0);
    ASSERT(mpz_cmp_si(s1, mpz_sgn(x)) == 0);
    ASSERT(mpz_sgn(t1) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 256, rng, arg);
    mpz_random_nz(y, 256, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(y, y);

    mpz_gcdext(g1, s1, t1, x, y);

    /* (g - x * s) / y == t */
    mpz_set(z, g1);
    mpz_submul(z, x, s1);
    mpz_divexact(z, z, y);

    ASSERT(mpz_cmp(z, t1) == 0);

    mpz_quo(q1, x, g1);
    mpz_quo(q2, y, g1);

    mpz_gcdext(g2, s2, t2, q1, q2);

    /* x * s + y * t == g */
    mpz_mul(x, x, s1);
    mpz_mul(y, y, t1);
    mpz_add(x, x, y);

    /* q1 * s2 + q2 * t2 == 1 */
    mpz_mul(q1, q1, s2);
    mpz_mul(q2, q2, t2);
    mpz_add(q1, q1, q2);

    ASSERT(mpz_cmp(g1, x) == 0);
    ASSERT(mpz_cmp_ui(g2, 1) == 0);
    ASSERT(mpz_cmp_ui(q1, 1) == 0);
  }

  for (mode = 0; mode < 3; mode++) {
    i = 0;
    j = 0;

    while (i < 20 || j < 20) {
      mpz_random_nz(x, 256, rng, arg);
      mpz_random_nz(y, 256, rng, arg);

      if (mp_random_limb(rng, arg) & 1)
        mpz_neg(x, x);

      if (mp_random_limb(rng, arg) & 1)
        mpz_neg(y, y);

      if (mode == 1)
        mpz_clrbit(x, 0);
      else if (mode == 2)
        mpz_clrbit(y, 0);

      mpz_gcdext(g1, s1, t1, x, y);

      if (mpz_cmp_ui(g1, 1) != 0) {
        j += 1;
        continue;
      }

      mpz_mul(x, x, s1);
      mpz_mod(x, x, y);

      ASSERT(mpz_cmp_ui(x, 1) == 0);

      i += 1;
    }
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(m);
  mpz_clear(g1);
  mpz_clear(s1);
  mpz_clear(t1);
  mpz_clear(g2);
  mpz_clear(s2);
  mpz_clear(t2);
  mpz_clear(q1);
  mpz_clear(q2);
}

static void
test_mpz_invert(mp_rng_f *rng, void *arg) {
  mpz_t x, z, m;
  int i = 0;

  printf("  - MPZ invert.\n");

  mpz_init(x);
  mpz_init(z);
  mpz_init(m);

  mpz_p192_mod(m);

  mpz_set_ui(x, 0);

  ASSERT(mpz_invert(z, x, m) == 0);

  for (i = 0; i < 100; i++) {
    do {
      mpz_random_nz(x, 384, rng, arg);
      mpz_mod(x, x, m);
    } while (mpz_sgn(x) == 0);

    if (i & 1)
      mpz_neg(x, x);

    ASSERT(mpz_invert(z, x, m) == 1);

    mpz_mul(z, z, x);
    mpz_mod(z, z, m);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  for (i = 0; i < 100; i++) {
    do {
      mpz_random_nz(x, 384, rng, arg);
      mpz_random_nz(m, 256, rng, arg);
      mpz_mod(x, x, m);

      if (i & 1)
        mpz_neg(x, x);
    } while (!mpz_invert(z, x, m));

    mpz_mul(z, z, x);
    mpz_mod(z, z, m);

    ASSERT(mpz_cmp_ui(z, 1) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
  mpz_clear(m);
}

static void
test_mpz_jacobi(void) {
  const int *v;
  mpz_t x, y;
  size_t i;

  printf("  - MPZ jacobi.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < ARRAY_SIZE(jacobi_vectors); i++) {
    v = jacobi_vectors[i];

    mpz_set_si(x, v[0]);
    mpz_set_si(y, v[1]);

    ASSERT(mpz_jacobi(x, y) == v[2]);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_kronecker(void) {
  const int *v;
  mpz_t x, y;
  size_t i;

  printf("  - MPZ kronecker.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < ARRAY_SIZE(kronecker_vectors); i++) {
    v = kronecker_vectors[i];

    mpz_set_si(x, v[0]);
    mpz_set_si(y, v[1]);

    ASSERT(mpz_kronecker(x, y) == v[2]);
    ASSERT(mpz_si_kronecker(v[0], y) == v[2]);
    ASSERT(mpz_kronecker_si(x, v[1]) == v[2]);

    if (v[0] >= 0)
      ASSERT(mpz_ui_kronecker(v[0], y) == v[2]);

    if (v[1] >= 0)
      ASSERT(mpz_kronecker_ui(x, v[1]) == v[2]);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_powm(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t, m;
  int i = 0;
  int k = 0;
  int j;

  printf("  - MPZ powm.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);
  mpz_init(m);

  mpz_p192_exp(y);
  mpz_p192_mod(m);

  mpz_set_ui(x, 0);

  ASSERT(mpz_jacobi(x, m) == 0);

  mpz_set(x, m);

  ASSERT(mpz_jacobi(x, m) == 0);

  while (i < 100 || k < 100) {
    mpz_random_nz(x, 384, rng, arg);
    mpz_mod(x, x, m);

    if (i & 1) {
      mpz_neg(x, x);
      mpz_mod(t, x, m);
    } else {
      mpz_set(t, x);
    }

    j = mpz_jacobi(x, m);

    ASSERT(j == mpz_legendre(x, m));

    if (j == 0) {
      ASSERT(mpz_sgn(x) == 0);
      continue;
    }

    ASSERT(j == 1 || j == -1);

    mpz_powm(z, x, y, m);
    mpz_sqr(z, z);
    mpz_mod(z, z, m);

    ASSERT((mpz_cmp(z, t) == 0) == (j == 1));

    mpz_powm_sec(z, x, y, m);
    mpz_sqr(z, z);
    mpz_mod(z, z, m);

    ASSERT((mpz_cmp(z, t) == 0) == (j == 1));

    i += (j == 1);
    k += (j == -1);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 192, rng, arg);
    mpz_random_nz(y, 128, rng, arg);
    mpz_random_nz(m, 192, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1) {
      if (mpz_invert(t, x, m))
        mpz_neg(y, y);
    }

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(m, m);

    mpz_powm(z, x, y, m);
    mpz_powm_simple(t, x, y, m);

    ASSERT(mpz_cmp(z, t) == 0);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 192, rng, arg);
    mpz_random_nz(y, MP_LIMB_BITS, rng, arg);
    mpz_random_nz(m, 192, rng, arg);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(x, x);

    if (mp_random_limb(rng, arg) & 1)
      mpz_neg(m, m);

    mpz_powm_ui(z, x, mpz_get_ui(y), m);
    mpz_powm_simple_ui(t, x, mpz_get_ui(y), m);

    ASSERT(mpz_cmp(z, t) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
  mpz_clear(m);
}

static void
test_mpz_sqrtm(mp_rng_f *rng, void *arg) {
  mpz_t x, z, t, p;
  int i, j, mode;

  printf("  - MPZ sqrtm.\n");

  mpz_init(x);
  mpz_init(z);
  mpz_init(t);
  mpz_init(p);

  for (mode = 0; mode < 3; mode++) {
    for (;;) {
      mpz_randprime(p, 128, rng, arg);

      if (mode == 0 && mpz_and_ui(p, 3) == 3)
        break;

      if (mode == 1 && mpz_and_ui(p, 7) == 5)
        break;

      if (mode == 2 && mpz_and_ui(p, 15) == 1)
        break;
    }

    i = 0;
    j = 0;

    while (i < 20 || j < 20) {
      mpz_random_nz(x, 128, rng, arg);

      if (i & 1)
        mpz_neg(x, x);

      mpz_mod(t, x, p);

      if (!mpz_sqrtm(z, x, p)) {
        ASSERT(mpz_jacobi(x, p) == -1);
        j += 1;
        continue;
      }

      ASSERT(mpz_jacobi(x, p) >= 0);

      mpz_sqr(z, z);
      mpz_mod(z, z, p);

      ASSERT(mpz_cmp(z, t) == 0);

      i += 1;
    }
  }

  mpz_clear(x);
  mpz_clear(z);
  mpz_clear(t);
  mpz_clear(p);
}

static void
test_mpz_sqrtpq(mp_rng_f *rng, void *arg) {
  mpz_t z, x, p, q, n;
  int i = 0;

  printf("  - MPZ sqrtpq.\n");

  mpz_init(z);
  mpz_init(x);
  mpz_init(p);
  mpz_init(q);
  mpz_init(n);

  mpz_randprime(p, 96, rng, arg);
  mpz_randprime(q, 96, rng, arg);
  mpz_mul(n, p, q);

  for (i = 0; i < 5; i++) {
    do {
      mpz_random_nz(x, 192, rng, arg);
      mpz_mod(x, x, n);
    } while (!mpz_sqrtpq(z, x, p, q));

    mpz_sqr(z, z);
    mpz_mod(z, z, n);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(z);
  mpz_clear(x);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(n);
}

static void
test_mpz_remove(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t;
  int i, j;

  printf("  - MPZ removal of factors.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
retry:
    mpz_random_nz(x, 128, rng, arg);
    mpz_random_nz(y, 32, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_set_ui(t, 2);
    mpz_mul_2exp(x, x, 11);
    mpz_mul_2exp(y, y, 11);

    ASSERT(mpz_remove(NULL, x, t) == mpz_ctz(x));
    ASSERT(mpz_remove(NULL, y, t) == mpz_ctz(y));

    mpz_quo_2exp(x, x, 11);
    mpz_quo_2exp(y, y, 11);

    if (mpz_remove(NULL, x, y) != 0)
      goto retry;

    mpz_set(t, x);

    for (j = 0; j < 4; j++)
      mpz_mul(t, t, y);

    ASSERT(mpz_remove(z, t, y) == 4);
    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
test_mpz_factorial(void) {
  mpz_t z;

  printf("  - MPZ factorial.\n");

  mpz_init(z);

  mpz_fac_ui(z, 10);

  ASSERT(mpz_cmp_ui(z, 3628800) == 0);

  mpz_2fac_ui(z, 10);

  ASSERT(mpz_cmp_ui(z, 945) == 0);

  mpz_mfac_uiui(z, 10, 3);

  ASSERT(mpz_cmp_ui(z, 280) == 0);

  mpz_clear(z);
}

static void
test_mpz_primorial(void) {
  mpz_t z, x;

  printf("  - MPZ primorial.\n");

  mpz_init(z);
  mpz_init(x);

  /* p(16)# */
  ASSERT(mpz_set_str(x, "32589158477190044730", 10));

  mpz_primorial_ui(z, 53);

  ASSERT(mpz_cmp(z, x) == 0);

  mpz_clear(z);
  mpz_clear(x);
}

static void
test_mpz_binomial(void) {
  mpz_t n, z;

  printf("  - MPZ binomial coefficient.\n");

  mpz_init(n);
  mpz_init(z);

  mpz_set_ui(n, 9);
  mpz_bin_ui(z, n, 3);

  ASSERT(mpz_cmp_ui(z, 84) == 0);

  mpz_set_si(n, -9);
  mpz_bin_ui(z, n, 3);

  ASSERT(mpz_cmp_si(z, -165) == 0);

  mpz_bin_uiui(z, 9, 3);

  ASSERT(mpz_cmp_ui(z, 84) == 0);

  mpz_clear(n);
  mpz_clear(z);
}

static void
test_mpz_fibonacci(void) {
  mpz_t z, r;

  printf("  - MPZ fibonacci sequences.\n");

  mpz_init(z);
  mpz_init(r);

  mpz_fib_ui(z, 12);

  ASSERT(mpz_cmp_ui(z, 144) == 0);

  mpz_fib2_ui(z, r, 12);

  ASSERT(mpz_cmp_ui(z, 144) == 0);
  ASSERT(mpz_cmp_ui(r, 89) == 0);

  mpz_lucnum_ui(z, 10);

  ASSERT(mpz_cmp_ui(z, 123) == 0);

  mpz_lucnum2_ui(z, r, 10);

  ASSERT(mpz_cmp_ui(z, 123) == 0);
  ASSERT(mpz_cmp_ui(r, 76) == 0);

  mpz_clear(z);
  mpz_clear(r);
}

static void
test_mpz_mod_primorial(mp_rng_f *rng, void *arg) {
  static const mp_limb_t d1 = MP_LIMB_C(4127218095);
  static const mp_limb_t d2 = MP_LIMB_C(3948078067);
  mp_limb_t ra, rb;
  mpz_t x;
  int i;

  printf("  - MPZ mod primorial.\n");

  mpz_init(x);

  for (i = 0; i < 100; i++) {
    mpz_urandomb(x, 256, rng, arg);

    mpz_mod_primorial(&ra, &rb, x);

    ASSERT(ra == mpz_rem_ui(x, d1));
    ASSERT(rb == mpz_rem_ui(x, d2));
  }

  mpz_clear(x);
}

static void
test_mpz_primes(mp_rng_f *rng, void *arg) {
  static const mp_limb_t mr_pseudos[] = {
    /* https://oeis.org/A001262 */
    2047,
    3277,
    4033,
    4681,
    8321,
    15841,
    29341,
    42799,
    49141,
    52633,
    65281,
    74665,
    80581,
    85489,
    88357,
    90751
  };

  static const mp_limb_t lucas_pseudos[] = {
    /* https://oeis.org/A217719 */
    989,
    3239,
    5777,
    10877,
    27971,
    29681,
    30739,
    31631,
    39059,
    72389,
    73919,
    75077
  };

  arc4_t arc4;
  size_t i;
  mpz_t p;

  mpz_init(p);

  {
    printf("  - MPZ primes.\n");

    for (i = 0; i < ARRAY_SIZE(prime_vectors); i++) {
      ASSERT(mpz_set_str(p, prime_vectors[i], 10));
      ASSERT(mpz_mr_prime_p(p, 16 + 1, 1, rng, arg));
      ASSERT(mpz_mr_prime_p(p, 1, 1, rng, arg));
      ASSERT(mpz_mr_prime_p(p, 1, 0, rng, arg));
      ASSERT(mpz_mr_prime_p(p, 0, 1, rng, arg));
      ASSERT(mpz_lucas_prime_p(p, 50));
      ASSERT(mpz_probab_prime_p(p, 20, rng, arg));
    }
  }

  {
    printf("  - MPZ composites.\n");

    for (i = 0; i < ARRAY_SIZE(composite_vectors); i++) {
      ASSERT(mpz_set_str(p, composite_vectors[i], 10));

      /* Initialize deterministic RNG. */
      arc4 = arc4_initial;

      /* Miller-Rabin with a deterministic RNG. */
      ASSERT(!mpz_mr_prime_p(p, 16 + 1, 1, arc4_rng, &arc4));
      ASSERT(!mpz_mr_prime_p(p, 4, 1, arc4_rng, &arc4));
      ASSERT(!mpz_mr_prime_p(p, 4, 0, arc4_rng, &arc4));

      if (i >= 8 && i <= 42) {
        /* Lucas pseudoprime. */
        ASSERT(mpz_lucas_prime_p(p, 50));
      } else {
        ASSERT(!mpz_lucas_prime_p(p, 50));
      }

      /* No composite should ever pass Baillie-PSW, random or otherwise. */
      ASSERT(!mpz_probab_prime_p(p, 20, arc4_rng, &arc4));
      ASSERT(!mpz_probab_prime_p(p, 20, rng, arg));
    }
  }

  {
    const mp_limb_t *want = mr_pseudos;
    size_t len = ARRAY_SIZE(mr_pseudos);

    printf("  - MPZ miller-rabin pseudo-primes.\n");

    for (i = 3; i < 100000; i += 2) {
      int pseudo;

      mpz_set_ui(p, i);

      pseudo = mpz_mr_prime_p(p, 1, 1, rng, arg)
            && !mpz_lucas_prime_p(p, 50);

      if (pseudo && (len == 0 || i != want[0]))
        ASSERT(0 && "miller-rabin: want false");
      else if (!pseudo && len >= 1 && i == want[0])
        ASSERT(0 && "miller-rabin: want true");

      if (len > 0 && i == want[0]) {
        want++;
        len--;
      }
    }

    ASSERT(len == 0);
  }

  {
    const mp_limb_t *want = lucas_pseudos;
    size_t len = ARRAY_SIZE(lucas_pseudos);

    printf("  - MPZ lucas pseudo-primes.\n");

    for (i = 3; i < 100000; i += 2) {
      int pseudo;

      mpz_set_ui(p, i);

      pseudo = mpz_lucas_prime_p(p, 50)
           && !mpz_mr_prime_p(p, 1, 1, rng, arg);

      if (pseudo && (len == 0 || i != want[0]))
        ASSERT(0 && "lucas: want false");
      else if (!pseudo && len >= 1 && i == want[0])
        ASSERT(0 && "lucas: want true");

      if (len > 0 && i == want[0]) {
        want++;
        len--;
      }
    }

    ASSERT(len == 0);
  }

  mpz_clear(p);
}

static void
test_mpz_randprime(mp_rng_f *rng, void *arg) {
  mpz_t x;

  printf("  - MPZ random prime.\n");

  mpz_init(x);

  mpz_randprime(x, 128, rng, arg);

  ASSERT(mpz_probab_prime_p(x, 20, rng, arg));

  mpz_clear(x);
}

static void
test_mpz_nextprime(mp_rng_f *rng, void *arg) {
  arc4_t arc4 = arc4_initial;
  size_t i;
  mpz_t x;

  printf("  - MPZ next prime.\n");

  mpz_init(x);

  for (i = 0; i < ARRAY_SIZE(mpz_test_primes); i++) {
    mpz_nextprime(x, x, rng, arg);

    ASSERT(mpz_cmp_ui(x, mpz_test_primes[i]) == 0);
  }

  mpz_urandomb(x, 128, arc4_rng, &arc4);

  ASSERT(!mpz_probab_prime_p(x, 20, arc4_rng, &arc4));

  mpz_nextprime(x, x, arc4_rng, &arc4);

  ASSERT(mpz_probab_prime_p(x, 20, arc4_rng, &arc4));

  mpz_clear(x);
}

static void
test_mpz_findprime(void) {
  arc4_t arc4 = arc4_initial;
  mpz_t x;

  printf("  - MPZ find prime.\n");

  mpz_init(x);

  mpz_urandomb(x, 128, arc4_rng, &arc4);

  ASSERT(!mpz_probab_prime_p(x, 20, arc4_rng, &arc4));
  ASSERT(mpz_findprime(x, x, 512, arc4_rng, &arc4));
  ASSERT(mpz_probab_prime_p(x, 20, arc4_rng, &arc4));

  mpz_clear(x);
}

static void
test_mpz_helpers(void) {
  static const mp_limb_t trailp[4] = {4, 3, 2, 0};
  static const mp_limb_t oddp[4] = {3, 3, 2, 1};
  static const mp_limb_t evenp[4] = {4, 3, 2, 1};
  static const mp_limb_t tzp[4] = {0, 3, 2, 1};
  static const mpz_t tz = MPZ_ROINIT_N(tzp, 4);
  mpz_t trail, odd, even;
  mp_limb_t *xp;
  mpz_t x, y;

  printf("  - MPZ helpers.\n");

  mpz_init(x);
  mpz_init(y);

  mpz_roinit_n(trail, trailp, 4);
  mpz_roinit_n(odd, oddp, 4);
  mpz_roinit_n(even, evenp, 4);

  mpz_set_ui(x, 1);
  mpz_set_ui(y, 2);

  ASSERT(mpz_odd_p(odd));
  ASSERT(!mpz_odd_p(even));
  ASSERT(mpz_even_p(even));
  ASSERT(!mpz_even_p(odd));
  ASSERT(mpz_ctz(tz) == MP_LIMB_BITS);
  ASSERT(mpz_bitlen(trail) == 2 * MP_LIMB_BITS + 2);
  ASSERT(mpz_bytelen(trail) == ((size_t)mpz_bitlen(trail) + 7) / 8);
  ASSERT(mpz_sizeinbase(trail, 256) == (size_t)mpz_bytelen(trail));
  ASSERT(mpz_getlimbn(trail, 0) == 4);
  ASSERT(mpz_size(trail) == 3);
  ASSERT(mpz_limbs_read(trail) == trailp);
  ASSERT(!mpz_fits_ui_p(trail));
  ASSERT(!mpz_fits_si_p(trail));
  ASSERT(mpz_fits_ui_p(x));
  ASSERT(mpz_fits_si_p(x));

  mpz_swap(x, y);

  ASSERT(mpz_cmp_ui(x, 2) == 0);
  ASSERT(mpz_cmp_ui(y, 1) == 0);

  {
    mpz_limbs_write(x, 1)[0] = 3;
    mpz_limbs_finish(x, 1);

    ASSERT(mpz_cmp_ui(x, 3) == 0);

    xp = mpz_limbs_write(x, 2);
    xp[0] = 4;
    xp[1] = 0;

    mpz_limbs_finish(x, 2);

    ASSERT(mpz_size(x) == 1);
    ASSERT(mpz_cmp_ui(x, 4) == 0);
  }

  {
    mpz_limbs_modify(y, 1)[0] = 2;
    mpz_limbs_finish(y, 1);
    ASSERT(mpz_cmp_ui(y, 2) == 0);
  }

  {
    mpz_mul_2exp(x, x, MP_LIMB_BITS);

    ASSERT(x->alloc == 2);
    ASSERT(x->size == 2);

    _mpz_realloc(x, 10);

    ASSERT(x->alloc == 10);
    ASSERT(x->size == 2);

    _mpz_realloc(x, 0);

    ASSERT(x->alloc == 1);
    ASSERT(x->size == 0);
  }

  {
    mpz_mul_2exp(y, y, MP_LIMB_BITS);

    ASSERT(y->alloc == 2);
    ASSERT(y->size == 2);

    mpz_realloc2(y, MP_LIMB_BITS * 10);

    ASSERT(y->alloc == 10);
    ASSERT(y->size == 2);

    mpz_realloc2(y, 0);

    ASSERT(y->alloc == 1);
    ASSERT(y->size == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_io(mp_rng_f *rng, void *arg) {
  unsigned char raw[32];
  mpz_t x, y;
  int i;

  printf("  - MPZ I/O.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 256, rng, arg);

    mpz_set_ui(y, 0);
    mpz_export(raw, x, sizeof(raw), 1);
    mpz_import(y, raw, sizeof(raw), 1);

    ASSERT(mpz_cmp(x, y) == 0);

    mpz_set_ui(y, 0);
    mpz_export(raw, x, sizeof(raw), -1);
    mpz_import(y, raw, sizeof(raw), -1);

    ASSERT(mpz_cmp(x, y) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_io_str(mp_rng_f *rng, void *arg) {
  mpz_t x, y;
  char *str;
  int i;

  printf("  - MPZ string I/O.\n");

  mpz_init(x);
  mpz_init(y);

  ASSERT(mpz_set_str(x, "0b10101010101010101", 0));
  ASSERT(mpz_cmp_ui(x, 87381) == 0);

  ASSERT(mpz_set_str(x, "0o12345671234", 0));
  ASSERT(mpz_cmp_ui(x, 1402434204) == 0);

  ASSERT(mpz_set_str(x, "012345671234", 0));
  ASSERT(mpz_cmp_ui(x, 1402434204) == 0);

  ASSERT(mpz_set_str(x, "1234567890", 0));
  ASSERT(mpz_cmp_ui(x, 1234567890) == 0);

  ASSERT(mpz_set_str(x, "+1234567890", 0));
  ASSERT(mpz_cmp_ui(x, 1234567890) == 0);

  ASSERT(mpz_set_str(x, "0xdeadbeef", 0));
  ASSERT(mpz_cmp_ui(x, 0xdeadbeef) == 0);

  for (i = 0; i < 100; i++) {
    mpz_random_nz(x, 256, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    str = mpz_get_str(x, 10);

    ASSERT(mpz_set_str(y, str, 10));
    ASSERT(mpz_cmp(x, y) == 0);

    free(str);

    str = mpz_get_str(x, 16);

    ASSERT(mpz_set_str(y, str, 16));
    ASSERT(mpz_cmp(x, y) == 0);

    free(str);

    mpz_print(x, 10, fake_puts);
    mpz_print(x, 16, fake_puts);
  }

  for (i = 2; i <= 62; i++) {
    mpz_random_nz(x, 256, rng, arg);

    str = mpz_get_str(x, i);

    ASSERT(mpz_set_str(y, str, i));
    ASSERT(mpz_cmp(x, y) == 0);

    free(str);

    mpz_print(x, i, fake_puts);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_io_str_vectors(void) {
  char *str;
  size_t i;
  mpz_t x;

  mpz_init(x);

  printf("  - MPZ string I/O (vectors).\n");

  for (i = 0; i < ARRAY_SIZE(prime_vectors); i++) {
    ASSERT(mpz_set_str(x, prime_vectors[i], 10));

    str = mpz_get_str(x, 10);

    ASSERT(strcmp(str, prime_vectors[i]) == 0);

    free(str);
  }

  for (i = 0; i < ARRAY_SIZE(composite_vectors); i++) {
    ASSERT(mpz_set_str(x, composite_vectors[i], 10));

    str = mpz_get_str(x, 10);

    ASSERT(strcmp(str, composite_vectors[i]) == 0);

    free(str);
  }

  mpz_clear(x);
}

static void
test_mpz_random(mp_rng_f *rng, void *arg) {
  mpz_t x, y;
  int i;

  printf("  - MPZ RNG.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_urandomb(x, 256, rng, arg);

    ASSERT(mpz_bitlen(x) <= 256);
    ASSERT(mpz_sgn(x) != 0);

    mpz_urandomm(y, x, rng, arg);

    ASSERT(mpz_sgn(y) == 1);
    ASSERT(mpz_cmp(y, x) < 0);

    mpz_neg(x, x);

    mpz_urandomm(y, x, rng, arg);

    ASSERT(mpz_sgn(y) == -1);
    ASSERT(mpz_cmp(y, x) > 0);
    ASSERT(mpz_cmpabs(y, x) < 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

/*
 * Vectors
 */

typedef void mpz_op_f(mpz_t z, const mpz_t, const mpz_t);
typedef void mpz_op_ui_f(mpz_t z, const mpz_t, mp_limb_t);
typedef void mpz_op_si_f(mpz_t z, const mpz_t, mp_long_t);
typedef void mpz_ui_op_f(mpz_t z, mp_limb_t, const mpz_t);
typedef void mpz_si_op_f(mpz_t z, mp_long_t, const mpz_t);
typedef mp_limb_t mpz_op_div_ui_f(mpz_t z, const mpz_t, mp_limb_t);
typedef mp_long_t mpz_op_div_si_f(mpz_t z, const mpz_t, mp_long_t);
typedef mp_limb_t mpz_op_mod_ui_f(const mpz_t, mp_limb_t);
typedef mp_long_t mpz_op_mod_si_f(const mpz_t, mp_long_t);

static void
mpz_test(const char *name, mpz_op_f *mpz_op,
         const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, x, y);

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_test_ui(const char *name, mpz_op_ui_f *mpz_op,
            const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, x, mpz_get_ui(y));

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_test_si(const char *name, mpz_op_si_f *mpz_op,
            const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, x, mpz_get_si(y));

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_ui_test(const char *name, mpz_ui_op_f *mpz_op,
            const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, mpz_get_ui(x), y);

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_si_test(const char *name, mpz_si_op_f *mpz_op,
            const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, mpz_get_si(x), y);

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_test_div_ui(const char *name, mpz_op_div_ui_f *mpz_op,
                const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, x, mpz_get_ui(y));

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_test_div_si(const char *name, mpz_op_div_si_f *mpz_op,
                const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));

    mpz_op(t, x, mpz_get_si(y));

    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_test_mod_ui(const char *name, mpz_op_mod_ui_f *mpz_op,
                const char *vectors[][3], size_t len) {
  mpz_t x, y, z;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));
    ASSERT(mpz_cmp_ui(z, mpz_op(x, mpz_get_ui(y))) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
mpz_test_mod_si(const char *name, mpz_op_mod_si_f *mpz_op,
                const char *vectors[][3], size_t len) {
  mpz_t x, y, z;
  size_t i;

  printf("    - %s vectors.\n", name);

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));
    ASSERT(mpz_cmp_si(z, mpz_op(x, mpz_get_si(y))) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
mpz_test_gcdext(const char *vectors[][5], size_t len) {
  mpz_t x, y, g, s, t, eg, es, et;
  size_t i;

  printf("    - gcdext vectors.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(g);
  mpz_init(s);
  mpz_init(t);
  mpz_init(eg);
  mpz_init(es);
  mpz_init(et);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(eg, vectors[i][2], 10));
    ASSERT(mpz_set_str(es, vectors[i][3], 10));
    ASSERT(mpz_set_str(et, vectors[i][4], 10));

    mpz_gcdext(g, s, t, x, y);

    ASSERT(mpz_cmp(g, eg) == 0);
    ASSERT(mpz_cmp(s, es) == 0);
    ASSERT(mpz_cmp(t, et) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(g);
  mpz_clear(s);
  mpz_clear(t);
  mpz_clear(eg);
  mpz_clear(es);
  mpz_clear(et);
}

static void
mpz_test_invert(const char *vectors[][3], size_t len) {
  mpz_t x, y, z, t;
  size_t i;

  printf("    - invert vectors.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(z, vectors[i][2], 10));
    ASSERT(mpz_invert(t, x, y));
    ASSERT(mpz_cmp(t, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
}

static void
mpz_test_powm(const char *vectors[][4], size_t len) {
  mpz_t x, y, m, z, t;
  size_t i;

  printf("    - powm vectors.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(m);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < len; i++) {
    ASSERT(mpz_set_str(x, vectors[i][0], 10));
    ASSERT(mpz_set_str(y, vectors[i][1], 10));
    ASSERT(mpz_set_str(m, vectors[i][2], 10));
    ASSERT(mpz_set_str(z, vectors[i][3], 10));

    mpz_powm(t, x, y, m);

    ASSERT(mpz_cmp(t, z) == 0);

    if (mpz_sgn(y) > 0 && mpz_odd_p(m)) {
      mpz_powm_sec(t, x, y, m);

      ASSERT(mpz_cmp(t, z) == 0);
    }
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(m);
  mpz_clear(z);
  mpz_clear(t);
}

static void
test_mpz_vectors(void) {
  printf("  - MPZ vectors.\n");

  mpz_test("add", mpz_add, mpz_add_vectors, ARRAY_SIZE(mpz_add_vectors));
  mpz_test_ui("add_ui", mpz_add_ui, mpz_add_ui_vectors,
              ARRAY_SIZE(mpz_add_ui_vectors));
  mpz_test_si("add_si", mpz_add_si, mpz_add_si_vectors,
              ARRAY_SIZE(mpz_add_si_vectors));

  mpz_test("sub", mpz_sub, mpz_sub_vectors, ARRAY_SIZE(mpz_sub_vectors));
  mpz_test_ui("sub_ui", mpz_sub_ui, mpz_sub_ui_vectors,
              ARRAY_SIZE(mpz_sub_ui_vectors));
  mpz_test_si("sub_si", mpz_sub_si, mpz_sub_si_vectors,
              ARRAY_SIZE(mpz_sub_si_vectors));

  mpz_ui_test("ui_sub", mpz_ui_sub, mpz_ui_sub_vectors,
              ARRAY_SIZE(mpz_ui_sub_vectors));
  mpz_si_test("si_sub", mpz_si_sub, mpz_si_sub_vectors,
              ARRAY_SIZE(mpz_si_sub_vectors));

  mpz_test("mul", mpz_mul, mpz_mul_vectors, ARRAY_SIZE(mpz_mul_vectors));
  mpz_test_ui("mul_ui", mpz_mul_ui, mpz_mul_ui_vectors,
              ARRAY_SIZE(mpz_mul_ui_vectors));
  mpz_test_si("mul_si", mpz_mul_si, mpz_mul_si_vectors,
              ARRAY_SIZE(mpz_mul_si_vectors));

  mpz_test("quo", mpz_quo, mpz_quo_vectors, ARRAY_SIZE(mpz_quo_vectors));
  mpz_test_div_ui("quo_ui", mpz_quo_ui, mpz_quo_ui_vectors,
                  ARRAY_SIZE(mpz_quo_ui_vectors));
  mpz_test_div_si("quo_si", mpz_quo_si, mpz_quo_si_vectors,
                  ARRAY_SIZE(mpz_quo_si_vectors));

  mpz_test("rem", mpz_rem, mpz_rem_vectors, ARRAY_SIZE(mpz_rem_vectors));
  mpz_test_mod_ui("rem_ui", mpz_rem_ui, mpz_rem_ui_vectors,
                  ARRAY_SIZE(mpz_rem_ui_vectors));
  mpz_test_mod_si("rem_si", mpz_rem_si, mpz_rem_si_vectors,
                  ARRAY_SIZE(mpz_rem_si_vectors));

  mpz_test("div", mpz_div, mpz_div_vectors, ARRAY_SIZE(mpz_div_vectors));
  mpz_test_div_ui("div_ui", mpz_div_ui, mpz_div_ui_vectors,
                  ARRAY_SIZE(mpz_div_ui_vectors));
  mpz_test_div_si("div_si", mpz_div_si, mpz_div_si_vectors,
                  ARRAY_SIZE(mpz_div_si_vectors));

  mpz_test("mod", mpz_mod, mpz_mod_vectors, ARRAY_SIZE(mpz_mod_vectors));
  mpz_test_mod_ui("mod_ui", mpz_mod_ui, mpz_mod_ui_vectors,
                  ARRAY_SIZE(mpz_mod_ui_vectors));
  mpz_test_mod_si("mod_si", mpz_mod_si, mpz_mod_si_vectors,
                  ARRAY_SIZE(mpz_mod_si_vectors));

  mpz_test("divround", mpz_divround, mpz_divround_vectors,
                       ARRAY_SIZE(mpz_divround_vectors));

  mpz_test("and", mpz_and, mpz_and_vectors, ARRAY_SIZE(mpz_and_vectors));
  mpz_test_mod_ui("and_ui", mpz_and_ui, mpz_and_ui_vectors,
                  ARRAY_SIZE(mpz_and_ui_vectors));
  mpz_test_si("and_si", mpz_and_si, mpz_and_si_vectors,
              ARRAY_SIZE(mpz_and_si_vectors));

  mpz_test("ior", mpz_ior, mpz_ior_vectors, ARRAY_SIZE(mpz_ior_vectors));
  mpz_test_ui("ior_ui", mpz_ior_ui, mpz_ior_ui_vectors,
              ARRAY_SIZE(mpz_ior_ui_vectors));
  mpz_test_si("ior_si", mpz_ior_si, mpz_ior_si_vectors,
              ARRAY_SIZE(mpz_ior_si_vectors));

  mpz_test("xor", mpz_xor, mpz_xor_vectors, ARRAY_SIZE(mpz_xor_vectors));
  mpz_test_ui("xor_ui", mpz_xor_ui, mpz_xor_ui_vectors,
              ARRAY_SIZE(mpz_xor_ui_vectors));
  mpz_test_si("xor_si", mpz_xor_si, mpz_xor_si_vectors,
              ARRAY_SIZE(mpz_xor_si_vectors));

  mpz_test("gcd", mpz_gcd, mpz_gcd_vectors, ARRAY_SIZE(mpz_gcd_vectors));
  mpz_test("lcm", mpz_lcm, mpz_lcm_vectors, ARRAY_SIZE(mpz_lcm_vectors));

  mpz_test_gcdext(mpz_gcdext_vectors, ARRAY_SIZE(mpz_gcdext_vectors));
  mpz_test_invert(mpz_invert_vectors, ARRAY_SIZE(mpz_invert_vectors));
  mpz_test_invert(mpz_invert_vectors2, ARRAY_SIZE(mpz_invert_vectors2));
  mpz_test_powm(mpz_powm_vectors, ARRAY_SIZE(mpz_powm_vectors));
}

/*
 * Benchmarks
 */

static void
bench_mpn_invert(mp_start_f *start, mp_end_f *end, mp_rng_f *rng, void *arg) {
  mp_limb_t xp[MP_K256_LIMBS * 2];
  mp_limb_t zp[MP_K256_LIMBS];
  mp_limb_t mp[MP_K256_LIMBS];
  mp_limb_t scratch[MPN_INVERT_ITCH(MP_K256_LIMBS)];
  mp_size_t mn = MP_K256_LIMBS;
  uint64_t tv;
  int i;

  mpn_k256_mod(mp);

  mpn_zero(xp, mn);

  do {
    mpn_random(xp, mn * 2, rng, arg);
    mpn_mod(xp, xp, mn * 2, mp, mn);
  } while (mpn_zero_p(xp, mn));

  start(&tv, "mpn_invert");

  for (i = 0; i < 100000; i++)
    ASSERT(mpn_invert_n(zp, xp, mp, mn, scratch) == 1);

  end(&tv, i);
}

/*
 * Test
 */

void
test_mpi_internal(mp_rng_f *rng, void *arg) {
  printf("Testing internal MPI functions...\n");

  /* MP */
  test_mp_helpers(rng, arg);
  test_mp_div(rng, arg);
  test_mp_inv_2by1(rng, arg);
  test_mp_div_2by1(rng, arg);
  test_mp_inv_3by2(rng, arg);
  test_mp_div_3by2(rng, arg);
  test_mp_inv_mod(rng, arg);
  test_mp_eratosthenes();

  /* MPN */
  test_mpn_init(rng, arg);
  test_mpn_assign(rng, arg);
  test_mpn_cmp();
  test_mpn_addsub(rng, arg);
  test_mpn_addsub_1(rng, arg);
  test_mpn_add_1();
  test_mpn_muldiv(rng, arg);
  test_mpn_muldiv_1(rng, arg);
  test_mpn_muldivexact(rng, arg);
  test_mpn_muldivexact_1(rng, arg);
  test_mpn_addmul_1(rng, arg);
  test_mpn_submul_1(rng, arg);
  test_mpn_sqr(rng, arg);
  test_mpn_mod(rng, arg);
  test_mpn_mod_1(rng, arg);
  test_mpn_divround(rng, arg);
  test_mpn_divround_1(rng, arg);
  test_mpn_roots(rng, arg);
  test_mpn_and(rng, arg);
  test_mpn_ior(rng, arg);
  test_mpn_xor(rng, arg);
  test_mpn_com(rng, arg);
  test_mpn_lshift(rng, arg);
  test_mpn_rshift(rng, arg);
  test_mpn_shift(rng, arg);
  test_mpn_bits(rng, arg);
  test_mpn_scan(rng, arg);
  test_mpn_popcount(rng, arg);
  test_mpn_hamdist(rng, arg);
  test_mpn_mask(rng, arg);
  test_mpn_negate(rng, arg);
  test_mpn_mulshift(rng, arg);
  test_mpn_reduce_weak(rng, arg);
  test_mpn_barrett(rng, arg);
  test_mpn_mont(rng, arg);
  test_mpn_gcd(rng, arg);
  test_mpn_gcdext(rng, arg);
  test_mpn_invert(rng, arg);
  test_mpn_jacobi();
  test_mpn_powm(rng, arg);
  test_mpn_helpers();
  test_mpn_select(rng, arg);
  test_mpn_sec_cmp();
  test_mpn_sec_cmp_rand(rng, arg);
  test_mpn_io(rng, arg);
  test_mpn_io_str(rng, arg);
  test_mpn_random(rng, arg);

  /* MPZ */
  test_mpz_init(rng, arg);
  test_mpz_assign(rng, arg);
  test_mpz_convert();
  test_mpz_cmp();
  test_mpz_cmp_signed();
  test_mpz_addsub(rng, arg);
  test_mpz_addsub_ui(rng, arg);
  test_mpz_addsub_si(rng, arg);
  test_mpz_ui_sub(rng, arg);
  test_mpz_si_sub(rng, arg);
  test_mpz_mulquorem(rng, arg);
  test_mpz_mulquo(rng, arg);
  test_mpz_mulquo_ui(rng, arg);
  test_mpz_mulquo_si(rng, arg);
  test_mpz_muldivmod(rng, arg);
  test_mpz_muldiv(rng, arg);
  test_mpz_muldiv_ui(rng, arg);
  test_mpz_muldiv_si(rng, arg);
  test_mpz_muldivexact(rng, arg);
  test_mpz_muldivexact_ui(rng, arg);
  test_mpz_muldivexact_si(rng, arg);
  test_mpz_divisibility(rng, arg);
  test_mpz_congruence(rng, arg);
  test_mpz_sqr(rng, arg);
  test_mpz_addmul(rng, arg);
  test_mpz_submul(rng, arg);
  test_mpz_rem(rng, arg);
  test_mpz_rem_ui(rng, arg);
  test_mpz_rem_si(rng, arg);
  test_mpz_mod(rng, arg);
  test_mpz_mod_ui(rng, arg);
  test_mpz_mod_si(rng, arg);
  test_mpz_divround(rng, arg);
  test_mpz_divround_ui(rng, arg);
  test_mpz_divround_si(rng, arg);
  test_mpz_pow(rng, arg);
  test_mpz_roots(rng, arg);
  test_mpz_and(rng, arg);
  test_mpz_ior(rng, arg);
  test_mpz_xor(rng, arg);
  test_mpz_com(rng, arg);
  test_mpz_lshift(rng, arg);
  test_mpz_rshift(rng, arg);
  test_mpz_shift(rng, arg);
  test_mpz_bits(rng, arg);
  test_mpz_scan();
  test_mpz_popcount(rng, arg);
  test_mpz_hamdist(rng, arg);
  test_mpz_mask(rng, arg);
  test_mpz_negate(rng, arg);
  test_mpz_mulshift(rng, arg);
  test_mpz_gcd(rng, arg);
  test_mpz_lcm(rng, arg);
  test_mpz_gcdext(rng, arg);
  test_mpz_invert(rng, arg);
  test_mpz_jacobi();
  test_mpz_kronecker();
  test_mpz_powm(rng, arg);
  test_mpz_sqrtm(rng, arg);
  test_mpz_sqrtpq(rng, arg);
  test_mpz_remove(rng, arg);
  test_mpz_factorial();
  test_mpz_primorial();
  test_mpz_binomial();
  test_mpz_fibonacci();
  test_mpz_mod_primorial(rng, arg);
  test_mpz_primes(rng, arg);
  test_mpz_randprime(rng, arg);
  test_mpz_nextprime(rng, arg);
  test_mpz_findprime();
  test_mpz_helpers();
  test_mpz_io(rng, arg);
  test_mpz_io_str(rng, arg);
  test_mpz_io_str_vectors();
  test_mpz_random(rng, arg);
  test_mpz_vectors();
}

/*
 * Benchmarks
 */

void
bench_mpi_internal(mp_start_f *start, mp_end_f *end, mp_rng_f *rng, void *arg) {
  bench_mpn_invert(start, end, rng, arg);
}
