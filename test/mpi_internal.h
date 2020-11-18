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
#include "data/prime_vectors.h"

#undef ASSERT
#define ASSERT(expr) ASSERT_ALWAYS(expr)

/*
 * RNG
 */

static void
mpn_random_nz(mp_limb_t *zp, int zn, mp_rng_f *rng, void *arg) {
  CHECK(zn != 0);

  do {
    mpn_random(zp, zn, rng, arg);
  } while (mpn_zero_p(zp, zn));
}

static mp_limb_t
mp_random_limb(mp_rng_f *rng, void *arg) {
  mp_limb_t z;
  mpn_random(&z, 1, rng, arg);
  return z;
}

static mp_long_t
mp_random_long(mp_rng_f *rng, void *arg) {
  mp_limb_t w = mp_random_limb(rng, arg);
  mp_long_t z;

  if (w == MP_LIMB_HI)
    return MP_LONG_MIN;

  z = w & (MP_LIMB_HI - 1);

  return (w & MP_LIMB_HI) ? -z : z;
}

/*
 * P192
 */

#define MP_P192_LIMBS ((192 + MP_LIMB_BITS - 1) / MP_LIMB_BITS)

static void
mpn_p192_mod(mp_limb_t *zp) {
  /* z = 2^192 - 2^64 - 1 */
  mp_limb_t ap[MP_P192_LIMBS + 1];
  mp_limb_t bp[MP_P192_LIMBS + 1];

  mpn_zero(ap, MP_P192_LIMBS + 1);
  mpn_zero(bp, MP_P192_LIMBS + 1);

  ap[192 / MP_LIMB_BITS] |= MP_LIMB_C(1) << (192 % MP_LIMB_BITS);
  bp[64 / MP_LIMB_BITS] |= MP_LIMB_C(1) << (64 % MP_LIMB_BITS);

  mpn_sub_n(ap, ap, bp, MP_P192_LIMBS + 1);
  mpn_sub_1(zp, ap, MP_P192_LIMBS, 1);
}

static void
mpn_p192_exp(mp_limb_t *zp) {
  mp_limb_t mp[MP_P192_LIMBS + 1];
  mp_limb_t ep[MP_P192_LIMBS + 1];

  mpn_zero(mp, MP_P192_LIMBS + 1);
  mpn_zero(ep, MP_P192_LIMBS + 1);

  mpn_p192_mod(mp);

  /* e = (p + 1) / 4 */
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
 * MP
 */

static void
test_mp_helpers(void) {
  printf("  - MP helpers.\n");

  {
    mp_limb_t *ptr = mp_alloc_limbs(2);

    ptr[0] = 0;
    ptr[1] = 0;

    mp_free_limbs(ptr);
  }

  {
    mp_limb_t *ptr = mp_alloc_vla(2);

    ptr[0] = 0;
    ptr[1] = 0;

    mp_free_vla(ptr, 2);
  }

  ASSERT(mp_clz(0) == MP_LIMB_BITS);
  ASSERT(mp_ctz(0) == MP_LIMB_BITS);
  ASSERT(mp_bitlen(0) == 0);

  ASSERT(mp_clz(1) == MP_LIMB_BITS - 1);
  ASSERT(mp_ctz(1) == 0);
  ASSERT(mp_ctz(2) == 1);
  ASSERT(mp_bitlen(1) == 1);
  ASSERT(mp_bitlen(2) == 2);

  ASSERT(!mp_gt_2(0, 0, 0, 1));
  ASSERT(!mp_gt_2(1, 0, 1, 1));
  ASSERT(mp_gt_2(1, 1, 1, 0));
  ASSERT(mp_gt_2(1, 0, 0, 0));
}

/*
 * MPN
 */

static void
test_mpn_init(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];

  printf("  - MPN init.\n");

  mpn_random_nz(ap, 4, rng, arg);

  ASSERT(!mpn_zero_p(ap, 4));

  mpn_zero(ap, 4);

  ASSERT(mpn_zero_p(ap, 4));

  mpn_cleanse(ap, 4);
}

static void
test_mpn_assign(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];

  printf("  - MPN assignment.\n");

  mpn_random_nz(ap, 4, rng, arg);
  mpn_zero(bp, 4);

  ASSERT(mpn_cmp(ap, bp, 4) != 0);

  mpn_copyi(bp, ap, 4);

  ASSERT(mpn_cmp(ap, bp, 4) == 0);

  mpn_zero(bp, 4);

  ASSERT(mpn_cmp(ap, bp, 4) != 0);

  mpn_copyd(bp, ap, 4);

  ASSERT(mpn_cmp(ap, bp, 4) == 0);

  mpn_set_1(bp, 4, 1);

  ASSERT(bp[0] == 1);
  ASSERT(mpn_zero_p(bp + 1, 3));
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
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t rp[6];
  int i, an, bn, rn;

  printf("  - MPN add/sub.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(ap, 4, rng, arg);
    mpn_random_nz(bp, 4, rng, arg);

    an = mpn_strip(ap, 4);
    bn = mpn_strip(bp, 4);
    rn = MP_MAX(an, bn);

    if (an >= bn)
      rp[rn] = mpn_add(rp, ap, an, bp, bn);
    else
      rp[rn] = mpn_add(rp, bp, bn, ap, an);

    rn += (rp[rn] != 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp(rp, rn, bp, bn) > 0);

    rp[rn] = mpn_add(rp, rp, rn, bp, bn);
    rn += (rp[rn] != 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp(rp, rn, bp, bn) > 0);

    ASSERT(mpn_sub(rp, rp, rn, bp, bn) == 0);
    rn = mpn_strip(rp, rn);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp(rp, rn, bp, bn) > 0);

    ASSERT(mpn_sub(rp, rp, rn, bp, bn) == 0);
    rn = mpn_strip(rp, rn);

    ASSERT(mpv_cmp(rp, rn, ap, an) == 0);

    ASSERT(mpn_sub(rp, rp, rn, ap, an) == 0);
    rn = mpn_strip(rp, rn);

    ASSERT(rn == 0);
  }
}

static void
test_mpn_addsub_1(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t b;
  mp_limb_t rp[6];
  int i, an, rn;

  printf("  - MPN add/sub (1 limb).\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(ap, 4, rng, arg);
    mpn_random_nz(&b, 1, rng, arg);

    an = mpn_strip(ap, 4);
    rn = an;

    rp[rn] = mpn_add_1(rp, ap, an, b);
    rn += (rp[rn] != 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp_1(rp, rn, b) > 0);

    rp[rn] = mpn_add_1(rp, rp, rn, b);
    rn += (rp[rn] != 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp_1(rp, rn, b) > 0);

    ASSERT(mpn_sub_1(rp, rp, rn, b) == 0);
    rn = mpn_strip(rp, rn);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp_1(rp, rn, b) > 0);

    ASSERT(mpn_sub_1(rp, rp, rn, b) == 0);
    rn = mpn_strip(rp, rn);

    ASSERT(mpv_cmp(rp, rn, ap, an) == 0);

    ASSERT(mpn_sub(rp, rp, rn, ap, an) == 0);
    rn = mpn_strip(rp, rn);

    ASSERT(rn == 0);
  }
}

static void
test_mpn_muldiv(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t rp[8 + 1];
  mp_limb_t sp[12];
  int i, an, bn, rn, sn;

  printf("  - MPN mul/div.\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(ap, 4, rng, arg);
    mpn_random_nz(bp, 4, rng, arg);

    an = mpn_strip(ap, 4);
    bn = mpn_strip(bp, 4);

    mpn_mul(rp, ap, an, bp, bn);
    rn = an + bn - (rp[an + bn - 1] == 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp(rp, rn, bp, bn) > 0);

    mpn_mul(sp, rp, rn, bp, bn);
    sn = rn + bn - (sp[rn + bn - 1] == 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp(rp, rn, bp, bn) > 0);

    mpn_div(rp, sp, sn, bp, bn);
    rn = mpn_strip(rp, sn - bn + 1);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp(rp, rn, bp, bn) > 0);

    mpn_div(rp, rp, rn, bp, bn);
    rn = mpn_strip(rp, rn - bn + 1);

    ASSERT(mpv_cmp(rp, rn, ap, an) == 0);
  }
}

static void
test_mpn_muldiv_1(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t b;
  mp_limb_t rp[6];
  int i, an, rn;

  printf("  - MPN mul/div (1 limb).\n");

  for (i = 0; i < 100; i++) {
    mpn_random_nz(ap, 4, rng, arg);
    mpn_random_nz(&b, 1, rng, arg);

    an = mpn_strip(ap, 4);
    rn = an;

    rp[rn] = mpn_mul_1(rp, ap, an, b);
    rn += (rp[rn] != 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp_1(rp, rn, b) > 0);

    rp[rn] = mpn_mul_1(rp, rp, rn, b);
    rn += (rp[rn] != 0);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp_1(rp, rn, b) > 0);

    mpn_div_1(rp, rp, rn, b);
    rn = mpn_strip(rp, rn);

    ASSERT(mpv_cmp(rp, rn, ap, an) > 0);
    ASSERT(mpv_cmp_1(rp, rn, b) > 0);

    mpn_div_1(rp, rp, rn, b);
    rn = mpn_strip(rp, rn);

    ASSERT(mpv_cmp(rp, rn, ap, an) == 0);
  }
}

static void
test_mpn_addmul_1(mp_rng_f *rng, void *arg) {
  mp_limb_t zp[5];
  mp_limb_t xp[5];
  mp_limb_t ep[5];
  mp_limb_t tp[5];
  mp_limb_t y;
  int i;

  printf("  - MPN add+mul.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(zp, 4, rng, arg);
    mpn_random(xp, 4, rng, arg);
    mpn_random(&y, 1, rng, arg);

    zp[4] = 0;
    xp[4] = 0;

    mpn_copyi(ep, zp, 5);

    ASSERT(mpn_mul_1(tp, xp, 5, y) == 0);
    ASSERT(mpn_add_n(ep, ep, tp, 5) == 0);

    ASSERT(mpn_addmul_1(zp, xp, 5, y) == 0);

    ASSERT(mpn_cmp(zp, ep, 5) == 0);
  }
}

static void
test_mpn_submul_1(mp_rng_f *rng, void *arg) {
  mp_limb_t zp[6];
  mp_limb_t xp[6];
  mp_limb_t ep[6];
  mp_limb_t tp[6];
  mp_limb_t y;
  int i;

  printf("  - MPN sub+mul.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(zp, 6, rng, arg);
    mpn_random(xp, 4, rng, arg);
    mpn_random(&y, 1, rng, arg);

    zp[5] |= 1;
    xp[4] = 0;
    xp[5] = 0;

    mpn_copyi(ep, zp, 6);

    ASSERT(mpn_mul_1(tp, xp, 6, y) == 0);
    ASSERT(mpn_sub_n(ep, ep, tp, 6) == 0);

    ASSERT(mpn_submul_1(zp, xp, 6, y) == 0);

    ASSERT(mpn_cmp(zp, ep, 6) == 0);
  }
}

static void
test_mpn_sqr(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t rp[8];
  mp_limb_t ep[8];
  int i;

  printf("  - MPN sqr.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    mpn_sqr(rp, ap, 4, ep);
    mpn_mul_n(ep, ap, ap, 4);

    ASSERT(mpn_cmp(rp, ep, 8) == 0);
  }
}

static void
test_mpn_mod(mp_rng_f *rng, void *arg) {
  mp_limb_t np[8];
  mp_limb_t dp[4];
  mp_limb_t qp[5];
  mp_limb_t rp[4];
  mp_limb_t tp[9];
  int i, nn, dn, qn, rn, tn;

  printf("  - MPN mod.\n");

  for (i = 0; i < 100; i++) {
    do {
      mpn_random_nz(np, 8, rng, arg);
      mpn_random_nz(dp, 4, rng, arg);
    } while (mpn_strip(np, 8) < 5);

    nn = mpn_strip(np, 8);
    dn = mpn_strip(dp, 4);

    mpn_divmod(qp, rp, np, nn, dp, dn);

    qn = mpn_strip(qp, nn - dn + 1);
    rn = mpn_strip(rp, dn);

    mpn_mul(tp, qp, qn, dp, dn);
    tn = qn + dn - (tp[qn + dn - 1] == 0);

    ASSERT(mpn_sub(tp, np, nn, tp, tn) == 0);
    tn = mpn_strip(tp, tn);

    ASSERT(mpv_cmp(rp, rn, tp, tn) == 0);
  }
}

static void
test_mpn_mod_1(mp_rng_f *rng, void *arg) {
  mp_limb_t np[8];
  mp_limb_t d;
  mp_limb_t qp[8];
  mp_limb_t r;
  mp_limb_t tp[9];
  int i, nn, qn, tn;

  printf("  - MPN mod (1 limb).\n");

  for (i = 0; i < 100; i++) {
    do {
      mpn_random_nz(np, 8, rng, arg);
      mpn_random_nz(&d, 1, rng, arg);
    } while (mpn_strip(np, 8) < 2);

    nn = mpn_strip(np, 8);

    r = mpn_divmod_1(qp, np, nn, d);

    qn = mpn_strip(qp, nn);

    tn = qn;
    tp[tn] = mpn_mul_1(tp, qp, qn, d);
    tn += (tp[tn] != 0);

    ASSERT(mpn_sub(tp, np, nn, tp, tn) == 0);
    tn = mpn_strip(tp, tn);

    ASSERT(mpv_cmp_1(tp, tn, r) == 0);
  }
}

static void
test_mpn_and(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t tp[4];
  int i;

  printf("  - MPN AND.\n");

  for (i = 0; i < 4; i++)
    bp[i] = MP_LIMB_MAX;

  bp[3] = 0;

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    mpn_and_n(tp, ap, bp, 4);

    ASSERT(mpn_cmp(tp, ap, 3) == 0);
    ASSERT(tp[3] == 0);

    mpn_random(tp, 4, rng, arg);

    mpn_and(tp, ap, 4, bp, 3);

    ASSERT(mpn_cmp(tp, ap, 3) == 0);
    ASSERT(tp[3] == 0);
  }
}

static void
test_mpn_ior(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t tp[4];
  mp_limb_t ep[4];
  int i, j;

  printf("  - MPN OR.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_random(bp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ap[j] | bp[j];

    mpn_ior_n(tp, ap, bp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);

    mpn_random(tp, 4, rng, arg);

    ep[3] = ap[3];

    mpn_ior(tp, ap, 4, bp, 3);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_xor(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t tp[4];
  mp_limb_t ep[4];
  int i, j;

  printf("  - MPN XOR.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_random(bp, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ap[j] ^ bp[j];

    mpn_xor_n(tp, ap, bp, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);

    mpn_random(tp, 4, rng, arg);

    ep[3] = ap[3];

    mpn_xor(tp, ap, 4, bp, 3);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_com(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t tp[4];
  mp_limb_t ep[4];
  int i, j;

  printf("  - MPN NOT.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_random(tp, 4, rng, arg);

    for (j = 0; j < 4; j++)
      ep[j] = ~ap[j];

    mpn_com(tp, ap, 4);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_lshift(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[5];
  mp_limb_t bp[5];
  int i;

  printf("  - MPN lshift.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    bp[4] = mpn_mul_1(bp, ap, 4, 8);
    ap[4] = mpn_lshift(ap, ap, 4, 3);

    ASSERT(mpn_cmp(ap, bp, 5) == 0);
  }
}

static void
test_mpn_rshift(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  int i;

  printf("  - MPN rshift.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    mpn_div_1(bp, ap, 4, 8);
    mpn_rshift(ap, ap, 4, 3);

    ASSERT(mpn_cmp(ap, bp, 4) == 0);
  }
}

static void
test_mpn_shift(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[5];
  mp_limb_t bp[5];
  int i;

  printf("  - MPN shift.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    ap[4] = 0;
    bp[4] = mpn_lshift(bp, ap, 4, 15);

    ASSERT(mpn_cmp(bp, ap, 5) > 0);

    mpn_rshift(bp, bp, 5, 15);

    ASSERT(mpn_cmp(bp, ap, 5) == 0);
  }
}

static void
test_mpn_bits(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  int bits = 4 * MP_LIMB_BITS;
  int steps = (bits + 5 - 1) / 5;
  int i, j, b1, b2;

  printf("  - MPN bits.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_copyi(bp, ap, 4);

    for (j = 0; j < bits; j++) {
      ASSERT(mpn_get_bit(ap, 4, j) == (bp[0] & 1));

      mpn_rshift(bp, bp, 4, 1);
    }

    ASSERT(mpn_get_bit(ap, 4, j) == 0);

    for (j = steps - 1; j >= 0; j--) {
      b1 = mpn_get_bits(ap, 4, j * 5, 5);

      mpv_rshift(bp, ap, 4, j * 5);

      b2 = bp[0] & MP_MASK(5);

      ASSERT(b1 == b2);
    }

    ASSERT(mpn_get_bits(ap, 4, steps * 5, 5) == 0);
  }

  mpn_zero(ap, 4);
  mpn_zero(bp, 4);

  ap[0] = 1;

  mpv_lshift(ap, ap, 1, 87);

  mpn_set_bit(bp, 87);

  ASSERT(mpn_cmp(ap, bp, 4) == 0);

  mpn_clr_bit(bp, 87);

  ASSERT(mpn_zero_p(bp, 4));
}

static void
test_mpn_mask(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t mp[4];
  mp_limb_t tp[4];
  mp_limb_t ep[4];
  int i;

  printf("  - MPN mask.\n");

  mpn_zero(mp, 4);
  mpn_set_bit(mp, 100);
  mpn_sub_1(mp, mp, 4, 1);

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    mpn_and_n(ep, ap, mp, 4);
    mpn_mask(tp, ap, 4, 100);

    ASSERT(mpn_cmp(tp, ep, 4) == 0);
  }
}

static void
test_mpn_negate(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t cp[4];
  mp_limb_t c1, c2;
  int i;

  printf("  - MPN negation.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_copyi(bp, ap, 4);
    mpn_zero(cp, 4);

    c1 = mpn_sub_n(ap, cp, ap, 4);
    c2 = mpn_neg(bp, bp, 4);

    ASSERT(c1 == c2);
    ASSERT(mpn_cmp(ap, bp, 4) == 0);
  }
}

static void
test_mpn_mulshift(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[MP_P192_LIMBS * 2];
  mp_limb_t bp[MP_P192_LIMBS * 2];
  mp_limb_t cp[MP_P192_LIMBS];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_limb_t qp[MP_P192_LIMBS * 2];
  mp_limb_t dp[MP_P192_LIMBS * 2];
  mp_limb_t rp[MP_P192_LIMBS * 2];
  mp_limb_t scratch[MPN_MULSHIFT_ITCH(MP_P192_LIMBS)];
  int mn = MP_P192_LIMBS;
  int sn = mn * 2;
  int cn, dn, qn, rn;
  int i;

  printf("  - MPN mulshift.\n");

  mpn_p192_mod(mp);

  mpn_zero(dp, mn * 2);

  dp[196 / MP_LIMB_BITS] |= MP_LIMB_C(1) << (196 % MP_LIMB_BITS);
  dn = mpn_strip(dp, mn * 2);

  for (i = 0; i < 100; i++) {
    mpn_random(ap, mn * 2, rng, arg);
    mpn_random(bp, mn * 2, rng, arg);

    mpn_mod(ap, ap, mn * 2, mp, mn);
    mpn_mod(bp, bp, mn * 2, mp, mn);

    ASSERT(mpn_mulshift(cp, ap, bp, mn, 196, scratch) == 0);

    mpn_zero(qp, mn * 2);
    mpn_zero(rp, mn * 2);
    mpn_mul_n(sp, ap, bp, mn);
    mpn_divmod(qp, rp, sp, sn, dp, dn);

    qn = mpn_strip(qp, sn - dn + 1);
    rn = mpn_strip(rp, dn);

    if (rn != 0) {
      int cmp;

      mpn_rshift(dp, dp, dn, 1);

      cmp = mpv_cmp(rp, rn, dp, dn);

      mpn_lshift(dp, dp, dn, 1);

      qp[qn] = mpn_add_1(qp, qp, qn, cmp >= 0);
      qn += (qp[qn] != 0);
    }

    cn = mpn_strip(cp, mn);

    ASSERT(mpv_cmp(cp, cn, qp, qn) == 0);
  }
}

static void
test_mpn_reduce_weak(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t tp[4];
  mp_limb_t scratch[MPN_REDUCE_WEAK_ITCH(4)];
  int i;

  printf("  - MPN reduction (weak).\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_random(bp, 4, rng, arg);

    ap[0] |= 1;
    bp[0] |= 1;
    ap[3] |= MP_LIMB_C(1) << (MP_LIMB_BITS - 1);
    bp[3] |= MP_LIMB_C(1) << (MP_LIMB_BITS - 1);

    if (mpn_cmp(ap, bp, 4) < 0) {
      mpn_copyi(tp, ap, 4);
      mpn_copyi(ap, bp, 4);
      mpn_copyi(bp, tp, 4);
    }

    ASSERT(mpn_reduce_weak(tp, ap, bp, 4, 0, scratch));
    ASSERT(mpn_cmp(tp, bp, 4) < 0);

    ASSERT(!mpn_reduce_weak(tp, tp, bp, 4, 0, scratch));

    ASSERT(!mpn_reduce_weak(tp, bp, ap, 4, 0, scratch));
    ASSERT(mpn_reduce_weak(tp, bp, ap, 4, 1, scratch));

    if (mpn_cmp(tp, ap, 4) >= 0)
      ASSERT(mpn_reduce_weak(tp, tp, ap, 4, 0, scratch));

    ASSERT(mpn_cmp(tp, ap, 4) < 0);
  }
}

static void
test_mpn_barrett(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[MP_P192_LIMBS * 2];
  mp_limb_t bp[MP_P192_LIMBS];
  mp_limb_t rp[MP_P192_LIMBS];
  mp_limb_t np[MP_P192_LIMBS];
  mp_limb_t mp[MPN_BARRETT_RESULT(MP_P192_LIMBS, MP_P192_LIMBS * 2)];
  mp_limb_t scratch1[MPN_BARRETT_ITCH(MP_P192_LIMBS * 2)];
  mp_limb_t scratch2[MPN_REDUCE_ITCH(MP_P192_LIMBS, MP_P192_LIMBS * 2)];
  int shift = MP_P192_LIMBS * 2;
  int n = MP_P192_LIMBS;
  int i;

  printf("  - MPN reduction (barrett).\n");

  mpn_p192_mod(np);

  mpn_barrett(mp, np, n, shift, scratch1);

  for (i = 0; i < 100; i++) {
    mpn_random(ap, n * 2, rng, arg);
    mpn_mod(bp, ap, n * 2, np, n);

    mpn_reduce(rp, ap, mp, np, n, shift, scratch2);

    ASSERT(mpn_cmp(rp, bp, n) == 0);
  }
}

static void
test_mpn_mont(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[MP_P192_LIMBS * 2];
  mp_limb_t bp[MP_P192_LIMBS * 2];
  mp_limb_t xp[MP_P192_LIMBS];
  mp_limb_t yp[MP_P192_LIMBS];
  mp_limb_t rp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_limb_t ep[MP_P192_LIMBS];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_limb_t scratch1[MPN_MONT_ITCH(MP_P192_LIMBS)];
  mp_limb_t scratch2[MPN_MONTMUL_ITCH(MP_P192_LIMBS)];
  mp_limb_t k;
  int mn = MP_P192_LIMBS;
  int i;

  printf("  - MPN reduction (montgomery).\n");

  mpn_p192_mod(mp);

  mpn_mont(&k, rp, mp, mn, scratch1);

  for (i = 0; i < 100; i++) {
    mpn_random(ap, mn * 2, rng, arg);
    mpn_random(bp, mn * 2, rng, arg);

    mpn_mod(ap, ap, mn * 2, mp, mn);
    mpn_mod(bp, bp, mn * 2, mp, mn);

    mpn_mul_n(sp, ap, bp, mn);
    mpn_mod(ep, sp, mn * 2, mp, mn);

    mpn_montmul(xp, ap, rp, mp, mn, k, scratch2);
    mpn_montmul(yp, bp, rp, mp, mn, k, scratch2);
    mpn_montmul(sp, xp, yp, mp, mn, k, scratch2);
    mpn_set_1(yp, mn, 1);
    mpn_montmul(sp, sp, yp, mp, mn, k, scratch2);

    ASSERT(mpn_cmp(sp, ep, mn) == 0);

    mpn_montmul_var(xp, ap, rp, mp, mn, k, scratch2);
    mpn_montmul_var(yp, bp, rp, mp, mn, k, scratch2);
    mpn_montmul_var(sp, xp, yp, mp, mn, k, scratch2);
    mpn_set_1(yp, mn, 1);
    mpn_montmul_var(sp, sp, yp, mp, mn, k, scratch2);
    mpn_mod(sp, sp, mn, mp, mn);

    ASSERT(mpn_cmp(sp, ep, mn) == 0);
  }
}

static void
test_mpn_invert(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[MP_P192_LIMBS * 2];
  mp_limb_t bp[MP_P192_LIMBS];
  mp_limb_t cp[MP_P192_LIMBS * 2];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_limb_t scratch[MPN_INVERT_ITCH(MP_P192_LIMBS)];
  int mn = MP_P192_LIMBS;
  int i = 0;

  printf("  - MPN invert.\n");

  mpn_p192_mod(mp);

  mpn_zero(ap, mn);
  ASSERT(mpn_invert_n(bp, ap, mp, mn, scratch) == 0);

  while (i < 100) {
    mpn_random(ap, mn * 2, rng, arg);
    mpn_mod(ap, ap, mn * 2, mp, mn);

    if (mpn_zero_p(ap, mn)) {
      ASSERT(mpn_invert_n(bp, ap, mp, mn, scratch) == 0);
      continue;
    }

    ASSERT(mpn_invert_n(bp, ap, mp, mn, scratch) == 1);

    mpn_mul_n(cp, ap, bp, mn);
    mpn_mod(cp, cp, mn * 2, mp, mn);

    ASSERT(cp[0] == 1);
    ASSERT(mpn_zero_p(cp + 1, mn - 1));

    i += 1;
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

    if (v[0] < 0 && v[1] < 0)
      j = -j;

    ASSERT(j == v[2]);
  }
}

static void
test_mpn_powm(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[MP_P192_LIMBS * 2];
  mp_limb_t ep[MP_P192_LIMBS];
  mp_limb_t mp[MP_P192_LIMBS];
  mp_limb_t rp[MP_P192_LIMBS];
  mp_limb_t sp[MP_P192_LIMBS * 2];
  mp_limb_t scratch1[MPN_JACOBI_ITCH(MP_P192_LIMBS)];
  mp_limb_t scratch2[MPN_POWM_ITCH(MP_P192_LIMBS, MP_P192_LIMBS)];
  mp_limb_t scratch3[MPN_SEC_POWM_ITCH(MP_P192_LIMBS)];
  int mn = MP_P192_LIMBS;
  int i = 0;
  int k = 0;
  int j, an;

  printf("  - MPN powm.\n");

  mpn_p192_mod(mp);
  mpn_p192_exp(ep);

  mpn_zero(ap, mn);
  ASSERT(mpn_jacobi_n(ap, mp, mn, scratch1) == 0);

  mpn_copyi(ap, mp, mn);
  ASSERT(mpn_jacobi_n(ap, mp, mn, scratch1) == 0);

  while (i < 100 || k < 100) {
    mpn_random(ap, mn * 2, rng, arg);
    mpn_mod(ap, ap, mn * 2, mp, mn);

    an = mpn_strip(ap, mn);

    j = mpn_jacobi_n(ap, mp, mn, scratch1);

    if (j == 0) {
      ASSERT(mpn_zero_p(ap, mn));
      continue;
    }

    ASSERT(j == 1 || j == -1);

    mpn_powm(rp, ap, an, ep, mn, mp, mn, scratch2);
    mpn_sqr(sp, rp, mn, scratch1);
    mpn_mod(rp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(rp, ap, mn) == 0) == (j == 1));

    mpn_sec_powm(rp, ap, an, ep, mn, mp, mn, scratch3);
    mpn_sqr(sp, rp, mn, scratch1);
    mpn_mod(rp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(rp, ap, mn) == 0) == (j == 1));

    mpn_mont_powm(rp, ap, an, ep, mn, mp, mn, scratch2);
    mpn_sqr(sp, rp, mn, scratch1);
    mpn_mod(rp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(rp, ap, mn) == 0) == (j == 1));

    mpn_div_powm(rp, ap, an, ep, mn, mp, mn, scratch2);
    mpn_sqr(sp, rp, mn, scratch1);
    mpn_mod(rp, sp, mn * 2, mp, mn);

    ASSERT((mpn_cmp(rp, ap, mn) == 0) == (j == 1));

    i += (j == 1);
    k += (j == -1);
  }
}

static void
test_mpn_helpers(void) {
  static const mp_limb_t trail[4] = {4, 3, 2, 0};
  static const mp_limb_t odd[4] = {3, 3, 2, 1};
  static const mp_limb_t even[4] = {4, 3, 2, 1};
  static const mp_limb_t tz[4] = {0, 3, 2, 1};
  mp_limb_t ap[1], bp[2];
  mp_limb_t *xp = ap;
  mp_limb_t *yp = bp;
  int xn = 1;
  int yn = 2;

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

  mpn_swap(&xp, &xn, &yp, &yn);

  ASSERT(xp == bp);
  ASSERT(xn == 2);
  ASSERT(yp == ap);
  ASSERT(yn == 1);
}

static void
test_mpn_select(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  mp_limb_t rp[4];

  printf("  - MPN select.\n");

  mpn_random_nz(ap, 4, rng, arg);
  mpn_random_nz(bp, 4, rng, arg);

  mpn_select(rp, ap, bp, 4, 0);

  ASSERT(mpn_cmp(rp, ap, 4) == 0);
  ASSERT(mpn_cmp(rp, bp, 4) != 0);

  mpn_select(rp, ap, bp, 4, 1);

  ASSERT(mpn_cmp(rp, bp, 4) == 0);
  ASSERT(mpn_cmp(rp, ap, 4) != 0);

  mpn_select_zero(rp, ap, 4, 0);

  ASSERT(mpn_cmp(rp, ap, 4) == 0);
  ASSERT(!mpn_zero_p(rp, 4));

  mpn_select_zero(rp, ap, 4, 1);

  ASSERT(mpn_cmp(rp, ap, 4) != 0);
  ASSERT(mpn_zero_p(rp, 4));
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

  ASSERT(!mpn_sec_equal(mod, minus1, 4));
  ASSERT(mpn_sec_equal(mod, mod, 4));

  ASSERT(mpn_sec_lt(minus1, mod, 4));
  ASSERT(!mpn_sec_lt(mod, mod, 4));
  ASSERT(!mpn_sec_lt(plus1, mod, 4));
  ASSERT(mpn_sec_lt(mod, full, 4));
  ASSERT(!mpn_sec_lt(full, mod, 4));

  ASSERT(mpn_sec_lte(minus1, mod, 4));
  ASSERT(mpn_sec_lte(mod, mod, 4));
  ASSERT(!mpn_sec_lte(plus1, mod, 4));
  ASSERT(mpn_sec_lte(mod, full, 4));
  ASSERT(!mpn_sec_lte(full, mod, 4));

  ASSERT(!mpn_sec_gt(minus1, mod, 4));
  ASSERT(!mpn_sec_gt(mod, mod, 4));
  ASSERT(mpn_sec_gt(plus1, mod, 4));
  ASSERT(!mpn_sec_gt(mod, full, 4));
  ASSERT(mpn_sec_gt(full, mod, 4));

  ASSERT(!mpn_sec_gte(minus1, mod, 4));
  ASSERT(mpn_sec_gte(mod, mod, 4));
  ASSERT(mpn_sec_gte(plus1, mod, 4));
  ASSERT(!mpn_sec_gte(mod, full, 4));
  ASSERT(mpn_sec_gte(full, mod, 4));

  ASSERT(mpn_sec_cmp(minus1, mod, 4) == -1);
  ASSERT(mpn_sec_cmp(mod, mod, 4) == 0);
  ASSERT(mpn_sec_cmp(plus1, mod, 4) == 1);
  ASSERT(mpn_sec_cmp(mod, full, 4) == -1);
  ASSERT(mpn_sec_cmp(full, mod, 4) == 1);
}

static void
test_mpn_sec_cmp_rand(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  int i;

  printf("  - MPN secure comparison (random).\n");

  for (i = 0; i < 1000; i++) {
    mpn_random(ap, 4, rng, arg);
    mpn_random(bp, 4, rng, arg);

    ASSERT(mpn_sec_zero_p(ap, 4) == mpn_zero_p(ap, 4));
    ASSERT(mpn_sec_zero_p(bp, 4) == mpn_zero_p(bp, 4));

    ASSERT(mpn_sec_equal(ap, ap, 4));
    ASSERT(mpn_sec_equal(ap, bp, 4) == (mpn_cmp(ap, bp, 4) == 0));

    ASSERT(mpn_sec_lt(ap, ap, 4) == 0);
    ASSERT(mpn_sec_lt(bp, bp, 4) == 0);
    ASSERT(mpn_sec_lt(ap, bp, 4) == (mpn_cmp(ap, bp, 4) < 0));
    ASSERT(mpn_sec_lt(bp, ap, 4) == (mpn_cmp(bp, ap, 4) < 0));

    ASSERT(mpn_sec_lte(ap, ap, 4) == 1);
    ASSERT(mpn_sec_lte(bp, bp, 4) == 1);
    ASSERT(mpn_sec_lte(ap, bp, 4) == (mpn_cmp(ap, bp, 4) <= 0));
    ASSERT(mpn_sec_lte(bp, ap, 4) == (mpn_cmp(bp, ap, 4) <= 0));

    ASSERT(mpn_sec_gt(ap, ap, 4) == 0);
    ASSERT(mpn_sec_gt(bp, bp, 4) == 0);
    ASSERT(mpn_sec_gt(ap, bp, 4) == (mpn_cmp(ap, bp, 4) > 0));
    ASSERT(mpn_sec_gt(bp, ap, 4) == (mpn_cmp(bp, ap, 4) > 0));

    ASSERT(mpn_sec_gte(ap, ap, 4) == 1);
    ASSERT(mpn_sec_gte(bp, bp, 4) == 1);
    ASSERT(mpn_sec_gte(ap, bp, 4) == (mpn_cmp(ap, bp, 4) >= 0));
    ASSERT(mpn_sec_gte(bp, ap, 4) == (mpn_cmp(bp, ap, 4) >= 0));

    ASSERT(mpn_sec_cmp(ap, ap, 4) == 0);
    ASSERT(mpn_sec_cmp(bp, bp, 4) == 0);
    ASSERT(mpn_sec_cmp(ap, bp, 4) == mpn_cmp(ap, bp, 4));
    ASSERT(mpn_sec_cmp(bp, ap, 4) == mpn_cmp(bp, ap, 4));
  }
}

static void
test_mpn_io(mp_rng_f *rng, void *arg) {
  unsigned char raw[4 * MP_LIMB_BYTES];
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  int i;

  printf("  - MPN I/O.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    mpn_zero(bp, 4);
    mpn_export(raw, sizeof(raw), ap, 4, 1);
    mpn_import(bp, 4, raw, sizeof(raw), 1);

    ASSERT(mpn_cmp(ap, bp, 4) == 0);

    mpn_zero(bp, 4);
    mpn_export(raw, sizeof(raw), ap, 4, -1);
    mpn_import(bp, 4, raw, sizeof(raw), -1);

    ASSERT(mpn_cmp(ap, bp, 4) == 0);
  }
}

static void
test_mpn_io_str(mp_rng_f *rng, void *arg) {
  /* Base-10 size = 78 + 1 */
  /* Base-16 size = 64 + 1 */
  mp_limb_t ap[4];
  mp_limb_t bp[4];
  char str[80];
  int i;

  printf("  - MPN string I/O.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 4, rng, arg);

    mpn_zero(bp, 4);
    mpn_get_str(str, ap, 4, 10);

    ASSERT(mpn_set_str(bp, 4, str, 10));
    ASSERT(mpn_cmp(ap, bp, 4) == 0);

    mpn_zero(bp, 4);
    mpn_get_str(str, ap, 4, 16);

    ASSERT(mpn_set_str(bp, 4, str, 16));
    ASSERT(mpn_cmp(ap, bp, 4) == 0);
  }
}

static void
test_mpn_random(mp_rng_f *rng, void *arg) {
  mp_limb_t ap[8];
  int i;

  printf("  - MPN RNG.\n");

  for (i = 0; i < 100; i++) {
    mpn_random(ap, 8, rng, arg);

    ASSERT(!mpn_zero_p(ap, 8));
  }
}

/*
 * MPZ
 */

static void
test_mpz_init(mp_rng_f *rng, void *arg) {
  mpz_t x;

  printf("  - MPZ init.\n");

  mpz_init(x);
  mpz_random_bits(x, 256, rng, arg);

  ASSERT(mpz_sgn(x) != 0);

  mpz_set_ui(x, 0);

  ASSERT(mpz_sgn(x) == 0);

  mpz_clear(x);

  mpz_init(x);
  mpz_cleanse(x);
}

static void
test_mpz_assign(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;

  printf("  - MPZ assignment.\n");

  mpz_init(x);
  mpz_init(y);

  mpz_random_bits(x, 256, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long(rng, arg);

    if (y == 0)
      y = 1;

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
test_mpz_mulquorem(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z, t;
  int i;

  printf("  - MPZ mul/quorem.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);
  mpz_init(t);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long(rng, arg);

    if (y == 0)
      y = 1;

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long(rng, arg);

    if (y == 0)
      y = 1;

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);
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
    mpz_random_bits(x, 250, rng, arg);
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
  mpz_t x, z;
  mp_long_t y;
  int i;

  printf("  - MPZ mul/divexact (1 limb, signed).\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    y = mp_random_long(rng, arg);

    if (y == 0)
      y = 1;

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
test_mpz_sqr(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ sqr.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul(y, x, x);
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
    mpz_random_bits(n, 250, rng, arg);
    mpz_random_bits(d, 125, rng, arg);

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
    mpz_random_bits(n, 250, rng, arg);
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
    mpz_random_bits(n, 250, rng, arg);

    if (i & 1)
      mpz_neg(n, n);

    d = mp_random_limb(rng, arg);

    if (d == 0)
      d = 1;

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
    mpz_random_bits(n, 250, rng, arg);
    mpz_random_bits(d, 125, rng, arg);

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
    mpz_random_bits(n, 250, rng, arg);
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
    mpz_random_bits(n, 250, rng, arg);

    if (i & 1)
      mpz_neg(n, n);

    d = mp_random_long(rng, arg);

    if (d == 0)
      d = 1;

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

  const int *v;
  size_t i;
  mpz_t x, z;

  printf("  - MPZ roots.\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < ARRAY_SIZE(roots); i++) {
    v = roots[i];

    mpz_set_ui(x, v[0]);
    mpz_sqrt(x, x);

    ASSERT(mpz_cmp_ui(x, v[1]) == 0);

    mpz_set_ui(x, v[0]);

    ASSERT(mpz_perfect_square_p(x) == v[2]);
  }

  for (i = 0; i < 100; i++) {
    mpz_random_bits(z, 250, rng, arg);
    mpz_sqr(x, z);

    ASSERT(mpz_perfect_square_p(x));

    mpz_sqrt(x, x);

    ASSERT(mpz_cmp(x, z) == 0);
  }

  mpz_clear(x);
  mpz_clear(z);
}

static void
test_mpz_and(mp_rng_f *rng, void *arg) {
  mpz_t x, m, z, e;
  int i, j;

  printf("  - MPZ AND.\n");

  mpz_init(x);
  mpz_init(m);
  mpz_init(z);
  mpz_init(e);

  mpz_set_bit(m, 100);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);

    mpz_set_ui(e, 0);

    for (j = 0; j < 100; j++) {
      if (mpz_get_bit(x, j))
        mpz_set_bit(e, j);
    }

    mpz_and(z, x, m);

    ASSERT(mpz_cmp(z, e) == 0);

    ASSERT(mpz_and_ui(x, 13) == (mpz_getlimbn(x, 0) & 13));
    ASSERT(mpz_and_si(x, 13) == (mp_long_t)(mpz_getlimbn(x, 0) & 13));
  }

  mpz_clear(x);
  mpz_clear(m);
  mpz_clear(z);
  mpz_clear(e);
}

static void
test_mpz_ior(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ OR.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 31, rng, arg);

    mpz_lshift(z, x, 32);
    mpz_ior(z, z, y);

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_rshift(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_lshift(z, x, 32);
    mpz_ior_ui(z, z, mpz_get_ui(y));

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_rshift(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_lshift(z, x, 32);
    mpz_ior_si(z, z, mpz_get_si(y));

    ASSERT(mpz_and_si(z, INT32_MAX) == mpz_get_si(y));

    mpz_rshift(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_xor(mp_rng_f *rng, void *arg) {
  mpz_t x, y, z;
  int i;

  printf("  - MPZ XOR.\n");

  mpz_init(x);
  mpz_init(y);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 31, rng, arg);

    mpz_lshift(z, x, 32);
    mpz_xor(z, z, y);

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_rshift(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_lshift(z, x, 32);
    mpz_xor_ui(z, z, mpz_get_ui(y));

    ASSERT(mpz_and_ui(z, UINT32_MAX) == mpz_get_ui(y));

    mpz_rshift(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);

    mpz_lshift(z, x, 32);
    mpz_xor_si(z, z, mpz_get_si(y));

    ASSERT(mpz_and_si(z, INT32_MAX) == mpz_get_si(y));

    mpz_rshift(z, z, 32);

    ASSERT(mpz_cmp(z, x) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
}

static void
test_mpz_com(mp_rng_f *rng, void *arg) {
  mpz_t x, z;
  int i;

  printf("  - MPZ NOT.\n");

  mpz_init(x);
  mpz_init(z);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);

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
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_mul_ui(y, x, 8);
    mpz_lshift(x, x, 3);

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
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_quo_ui(y, x, 8);
    mpz_rshift(x, x, 3);

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
    mpz_random_bits(x, 250, rng, arg);

    if (i & 1)
      mpz_neg(x, x);

    mpz_lshift(y, x, 173);

    ASSERT(mpz_cmpabs(y, x) > 0);

    mpz_rshift(y, y, 173);

    ASSERT(mpz_cmp(x, y) == 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_bits(mp_rng_f *rng, void *arg) {
  int bits = 4 * MP_LIMB_BITS;
  int steps = (bits + 5 - 1) / 5;
  int i, j, b1, b2;
  mpz_t x, y;

  printf("  - MPZ bits.\n");

  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, bits, rng, arg);
    mpz_set(y, x);

    for (j = 0; j < bits; j++) {
      ASSERT(mpz_get_bit(x, j) == (mpz_getlimbn(y, 0) & 1));

      mpz_rshift(y, y, 1);
    }

    ASSERT(mpz_get_bit(x, j) == 0);

    for (j = steps - 1; j >= 0; j--) {
      b1 = mpz_get_bits(x, j * 5, 5);

      mpz_rshift(y, x, j * 5);

      b2 = mpz_getlimbn(y, 0) & MP_MASK(5);

      ASSERT(b1 == b2);
    }

    ASSERT(mpz_get_bits(x, steps * 5, 5) == 0);
  }

  mpz_set_ui(x, 1);
  mpz_lshift(x, x, 256);

  mpz_set_ui(y, 0);
  mpz_set_bit(y, 256);

  ASSERT(mpz_cmp(x, y) == 0);

  mpz_clr_bit(y, 256);

  ASSERT(mpz_sgn(y) == 0);

  mpz_clear(x);
  mpz_clear(y);
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

  mpz_set_bit(m, 100);
  mpz_sub_ui(m, m, 1);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);

    mpz_and(e, x, m);
    mpz_mask(z, x, 100);

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
    mpz_random_bits(x, 250, rng, arg);
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
test_mpz_gcd(mp_rng_f *rng, void *arg) {
  mpz_t g1, g2, x, y;
  int i;

  printf("  - MPZ gcd.\n");

  mpz_init(g1);
  mpz_init(g2);
  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_gcd(g1, x, y);
    mpz_gcdext(g2, NULL, NULL, x, y);

    ASSERT(mpz_cmp(g1, g2) == 0);
  }

  mpz_clear(g1);
  mpz_clear(g2);
  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_lcm(mp_rng_f *rng, void *arg) {
  mpz_t l1, l2, g, x, y;
  int i;

  printf("  - MPZ lcm.\n");

  mpz_init(l1);
  mpz_init(l2);
  mpz_init(g);
  mpz_init(x);
  mpz_init(y);

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 250, rng, arg);
    mpz_random_bits(y, 125, rng, arg);

    if (i & 1)
      mpz_swap(x, y);

    mpz_lcm(l1, x, y);

    mpz_gcd(g, x, y);
    mpz_divexact(l2, x, g);
    mpz_mul(l2, l2, y);

    ASSERT(mpz_cmp(l1, l2) == 0);
  }

  mpz_clear(l1);
  mpz_clear(l2);
  mpz_clear(g);
  mpz_clear(x);
  mpz_clear(y);
}

static void
test_mpz_gcdext(mp_rng_f *rng, void *arg) {
  mpz_t g, x, z, m;
  int i = 0;

  printf("  - MPZ gcdext.\n");

  mpz_init(g);
  mpz_init(x);
  mpz_init(z);
  mpz_init(m);

  mpz_p192_mod(m);

  while (i < 100) {
    mpz_random_bits(x, 384, rng, arg);
    mpz_mod(x, x, m);

    if (mpz_sgn(x) == 0)
      continue;

    mpz_gcdext(g, z, NULL, x, m);

    ASSERT(mpz_cmp_ui(g, 1) == 0);

    mpz_mod(z, z, m);
    mpz_mul(z, z, x);
    mpz_mod(z, z, m);

    ASSERT(mpz_cmp_ui(z, 1) == 0);

    i += 1;
  }

  mpz_clear(g);
  mpz_clear(x);
  mpz_clear(z);
  mpz_clear(m);
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

  while (i < 100) {
    mpz_random_bits(x, 384, rng, arg);
    mpz_mod(x, x, m);

    if (i & 1)
      mpz_neg(x, x);

    if (mpz_sgn(x) == 0) {
      ASSERT(mpz_invert(z, x, m) == 0);
      continue;
    }

    ASSERT(mpz_invert(z, x, m) == 1);

    mpz_mul(z, z, x);
    mpz_mod(z, z, m);

    ASSERT(mpz_cmp_ui(z, 1) == 0);

    i += 1;
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
    mpz_random_bits(x, 384, rng, arg);
    mpz_mod(x, x, m);

    if (i & 1) {
      mpz_neg(x, x);
      mpz_mod(t, x, m);
    } else {
      mpz_set(t, x);
    }

    j = mpz_jacobi(x, m);

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

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(z);
  mpz_clear(t);
  mpz_clear(m);
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

  size_t i;
  mpz_t p;

  mpz_init(p);

  {
    printf("  - MPZ primes.\n");

    for (i = 0; i < ARRAY_SIZE(prime_vectors); i++) {
      ASSERT(mpz_set_str(p, prime_vectors[i], 10));
      ASSERT(mpz_is_prime_mr(p, 16 + 1, 1, rng, arg));
      ASSERT(mpz_is_prime_mr(p, 1, 1, rng, arg));
      ASSERT(mpz_is_prime_mr(p, 1, 0, rng, arg));
      ASSERT(mpz_is_prime_mr(p, 0, 1, rng, arg));
      ASSERT(mpz_is_prime_lucas(p, 50));
      ASSERT(mpz_is_prime(p, 20, rng, arg));
    }
  }

  {
    printf("  - MPZ composites.\n");

    for (i = 0; i < ARRAY_SIZE(composite_vectors); i++) {
      ASSERT(mpz_set_str(p, composite_vectors[i], 10));

      /* Miller-Rabin. */
      ASSERT(!mpz_is_prime_mr(p, 16 + 1, 1, rng, arg));

      if (i >= 8 && i <= 42) {
        /* Lucas pseudoprime. */
        ASSERT(mpz_is_prime_lucas(p, 50));
      } else {
        ASSERT(!mpz_is_prime_lucas(p, 50));
      }

      /* No composite should ever pass Baillie-PSW. */
      ASSERT(!mpz_is_prime(p, 20, rng, arg));
    }
  }

  {
    const mp_limb_t *want = mr_pseudos;
    size_t len = ARRAY_SIZE(mr_pseudos);

    printf("  - MPZ miller-rabin pseudo-primes...\n");

    for (i = 3; i < 100000; i += 2) {
      int pseudo;

      mpz_set_ui(p, i);

      pseudo = mpz_is_prime_mr(p, 1, 1, rng, arg)
            && !mpz_is_prime_lucas(p, 50);

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

    printf("  - MPZ lucas pseudo-primes...\n");

    for (i = 3; i < 100000; i += 2) {
      int pseudo;

      mpz_set_ui(p, i);

      pseudo = mpz_is_prime_lucas(p, 50)
           && !mpz_is_prime_mr(p, 1, 1, rng, arg);

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
test_mpz_random_prime(mp_rng_f *rng, void *arg) {
  mpz_t x;

  printf("  - MPZ random prime.\n");

  mpz_init(x);

  mpz_random_prime(x, 128, rng, arg);

  ASSERT(mpz_is_prime(x, 20, rng, arg));

  mpz_clear(x);
}

static void
test_mpz_helpers(void) {
  mpz_t trail, odd, even, tz;
  mpz_t x, y;

  printf("  - MPZ helpers.\n");

  mpz_init(trail);
  mpz_init(odd);
  mpz_init(even);
  mpz_init(tz);
  mpz_init(x);
  mpz_init(y);

  mpz_grow(trail, 4);
  mpz_grow(odd, 4);
  mpz_grow(even, 4);
  mpz_grow(tz, 4);

  trail[0].limbs[0] = 4;
  trail[0].limbs[1] = 3;
  trail[0].limbs[2] = 2;
  trail[0].limbs[3] = 0;
  trail[0].size = 3;

  odd[0].limbs[0] = 3;
  odd[0].limbs[1] = 3;
  odd[0].limbs[2] = 2;
  odd[0].limbs[3] = 1;
  odd[0].size = 4;

  even[0].limbs[0] = 4;
  even[0].limbs[1] = 3;
  even[0].limbs[2] = 2;
  even[0].limbs[3] = 1;
  even[0].size = 4;

  tz[0].limbs[0] = 0;
  tz[0].limbs[1] = 3;
  tz[0].limbs[2] = 2;
  tz[0].limbs[3] = 1;
  tz[0].size = 4;

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

  mpz_swap(x, y);

  ASSERT(mpz_cmp_ui(x, 2) == 0);
  ASSERT(mpz_cmp_ui(y, 1) == 0);

  mpz_clear(trail);
  mpz_clear(odd);
  mpz_clear(even);
  mpz_clear(tz);
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
    mpz_random_bits(x, 256, rng, arg);

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

  for (i = 0; i < 100; i++) {
    mpz_random_bits(x, 256, rng, arg);

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
    mpz_random_bits(x, 256, rng, arg);

    ASSERT(mpz_sgn(x) != 0);

    mpz_random_int(y, x, rng, arg);

    ASSERT(mpz_sgn(y) == 1);
    ASSERT(mpz_cmp(y, x) < 0);

    mpz_neg(x, x);

    mpz_random_int(y, x, rng, arg);

    ASSERT(mpz_sgn(y) == -1);
    ASSERT(mpz_cmp(y, x) > 0);
    ASSERT(mpz_cmpabs(y, x) < 0);
  }

  mpz_clear(x);
  mpz_clear(y);
}

/*
 * Test
 */

void
test_mpi_internal(mp_rng_f *rng, void *arg) {
  printf("Testing internal MPI functions...\n");

  /* MP */
  test_mp_helpers();

  /* MPN */
  test_mpn_init(rng, arg);
  test_mpn_assign(rng, arg);
  test_mpn_cmp();
  test_mpn_addsub(rng, arg);
  test_mpn_addsub_1(rng, arg);
  test_mpn_muldiv(rng, arg);
  test_mpn_muldiv_1(rng, arg);
  test_mpn_addmul_1(rng, arg);
  test_mpn_submul_1(rng, arg);
  test_mpn_sqr(rng, arg);
  test_mpn_mod(rng, arg);
  test_mpn_mod_1(rng, arg);
  test_mpn_and(rng, arg);
  test_mpn_ior(rng, arg);
  test_mpn_xor(rng, arg);
  test_mpn_com(rng, arg);
  test_mpn_lshift(rng, arg);
  test_mpn_rshift(rng, arg);
  test_mpn_shift(rng, arg);
  test_mpn_bits(rng, arg);
  test_mpn_mask(rng, arg);
  test_mpn_negate(rng, arg);
  test_mpn_mulshift(rng, arg);
  test_mpn_reduce_weak(rng, arg);
  test_mpn_barrett(rng, arg);
  test_mpn_mont(rng, arg);
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
  test_mpz_sqr(rng, arg);
  test_mpz_rem(rng, arg);
  test_mpz_rem_ui(rng, arg);
  test_mpz_rem_si(rng, arg);
  test_mpz_mod(rng, arg);
  test_mpz_mod_ui(rng, arg);
  test_mpz_mod_si(rng, arg);
  test_mpz_roots(rng, arg);
  test_mpz_and(rng, arg);
  test_mpz_ior(rng, arg);
  test_mpz_xor(rng, arg);
  test_mpz_com(rng, arg);
  test_mpz_lshift(rng, arg);
  test_mpz_rshift(rng, arg);
  test_mpz_shift(rng, arg);
  test_mpz_bits(rng, arg);
  test_mpz_mask(rng, arg);
  test_mpz_negate(rng, arg);
  test_mpz_gcd(rng, arg);
  test_mpz_lcm(rng, arg);
  test_mpz_gcdext(rng, arg);
  test_mpz_invert(rng, arg);
  test_mpz_jacobi();
  test_mpz_powm(rng, arg);
  test_mpz_primes(rng, arg);
  test_mpz_random_prime(rng, arg);
  test_mpz_helpers();
  test_mpz_io(rng, arg);
  test_mpz_io_str(rng, arg);
  test_mpz_io_str_vectors();
  test_mpz_random(rng, arg);
}
