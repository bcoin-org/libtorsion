/*!
 * ecc-internal.h - ecc internal tests for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#undef ASSERT
#define ASSERT(expr) ASSERT_ALWAYS(expr)

/*
 * Debug Helpers
 */

static size_t
mpn_out_str(FILE *stream, int base, mp_srcptr xp, mp_size_t xn) {
  mp_size_t bytes = 0;
  mp_limb_t ch;
  mp_size_t i;

  ASSERT(base == 16);

  if (xn < 0) {
    fputc('-', stream);
    xn = -xn;
  }

  xn = mpn_normalized_size(xp, xn);

  if (xn == 0) {
    fputc('0', stream);
    return 1;
  }

  while (xn--) {
    i = MP_LIMB_BITS / 4;

    while (i--) {
      ch = (xp[xn] >> (i * 4)) & 0x0f;

      if (bytes == 0 && ch == 0)
        continue;

      if (ch < 0x0a)
        ch += '0';
      else
        ch += 'a' - 0x0a;

      fputc(ch, stream);

      bytes += 1;
    }
  }

  return bytes;
}

TORSION_UNUSED static size_t
mpz_out_str(FILE *stream, int base, const mpz_t x) {
  return mpn_out_str(stream, base, x->_mp_d, x->_mp_size);
}

TORSION_UNUSED static void
sc_print(const scalar_field_t *sc, const sc_t a) {
  mpn_out_str(stdout, 16, a, sc->limbs);
  printf("\n");
}

static void
fe_out_str(const prime_field_t *fe, const fe_t a) {
  mp_limb_t xp[MAX_FIELD_LIMBS];

  fe_get_limbs(fe, xp, a);

  mpn_out_str(stdout, 16, xp, fe->limbs);
}

TORSION_UNUSED static void
fe_print(const prime_field_t *fe, const fe_t a) {
  fe_out_str(fe, a);
  printf("\n");
}

TORSION_UNUSED static void
wge_print(const wei_t *ec, const wge_t *p) {
  const prime_field_t *fe = &ec->fe;

  if (wge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    printf("(");
    fe_out_str(fe, p->x);
    printf(", ");
    fe_out_str(fe, p->y);
    printf(")\n");
  }
}

TORSION_UNUSED static void
jge_print(const wei_t *ec, const jge_t *p) {
  const prime_field_t *fe = &ec->fe;

  if (jge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    printf("(");
    fe_out_str(fe, p->x);
    printf(", ");
    fe_out_str(fe, p->y);
    printf(", ");
    fe_out_str(fe, p->z);
    printf(")\n");
  }
}

TORSION_UNUSED static void
mge_print(const mont_t *ec, const mge_t *p) {
  const prime_field_t *fe = &ec->fe;

  if (mge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    printf("(");
    fe_out_str(fe, p->x);
    printf(", ");
    fe_out_str(fe, p->y);
    printf(")\n");
  }
}

TORSION_UNUSED static void
pge_print(const mont_t *ec, const pge_t *p) {
  const prime_field_t *fe = &ec->fe;

  if (pge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    printf("(");
    fe_out_str(fe, p->x);
    printf(", ");
    fe_out_str(fe, p->z);
    printf(")\n");
  }
}

TORSION_UNUSED static void
xge_print(const edwards_t *ec, const xge_t *p) {
  const prime_field_t *fe = &ec->fe;

  if (xge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    printf("(");
    fe_out_str(fe, p->x);
    printf(", ");
    fe_out_str(fe, p->y);
    printf(", ");
    fe_out_str(fe, p->z);
    printf(")\n");
  }
}

/*
 * Helpers
 */

static int
stupid_bytes_zero(const unsigned char *a, size_t size) {
  while (size--) {
    if (a[size] != 0)
      return 0;
  }

  return 1;
}

static int
revcmp(const unsigned char *a, const unsigned char *b, size_t size) {
  while (size--) {
    if (a[size] < b[size])
      return -1;

    if (a[size] > b[size])
      return 1;
  }

  return 0;
}

/*
 * Memcmp
 */

static void
test_memcmp(void) {
  static const unsigned char a[4] = {0, 1, 2, 3};
  static const unsigned char b[4] = {0, 1, 2, 3};
  static const unsigned char c[4] = {3, 2, 1, 0};
  static const unsigned char d[4] = {3, 2, 1, 0};

  ASSERT(torsion_memcmp(a, b, 4) == 0);
  ASSERT(torsion_memcmp(c, d, 4) == 0);
  ASSERT(torsion_memcmp(a, b, 4) >= 0);
  ASSERT(torsion_memcmp(c, d, 4) >= 0);
  ASSERT(torsion_memcmp(a, b, 4) <= 0);
  ASSERT(torsion_memcmp(c, d, 4) <= 0);
  ASSERT(torsion_memcmp(a, c, 4) != 0);
  ASSERT(torsion_memcmp(c, a, 4) != 0);
  ASSERT(torsion_memcmp(a, c, 4) < 0);
  ASSERT(torsion_memcmp(c, a, 4) > 0);
}

/*
 * Scalar
 */

static void
test_scalar(drbg_t *rng) {
  printf("  - Scalar sanity check.\n");

  {
    const unsigned char expect1[32] = {
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    const unsigned char expect2[32] = {
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01
    };

    wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
    scalar_field_t *sc = &ec->sc;
    mp_limb_t r[MAX_SCALAR_LIMBS];
    mp_limb_t t[MAX_REDUCE_LIMBS];
    unsigned char raw[MAX_SCALAR_SIZE];

    sc_export(sc, raw, sc->n);

    mpn_zero(r, ARRAY_SIZE(r));
    mpn_zero(t, ARRAY_SIZE(t));

    mpn_copyi(t, sc->n, sc->limbs);

    sc_reduce(sc, r, t);

    ASSERT(sc_is_zero(sc, r));

    raw[sc->size - 1] -= 1;

    ASSERT(sc_import(sc, r, raw));

    raw[sc->size - 1] += 1;

    ASSERT(!sc_import(sc, r, raw));

    raw[sc->size - 1] += 1;

    ASSERT(!sc_import(sc, r, raw));

    sc_export(sc, raw, sc->n);

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(torsion_memcmp(raw, expect1, 32) == 0);

    sc_export(sc, raw, sc->n);

    raw[sc->size - 1] += 1;

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(torsion_memcmp(raw, expect2, 32) == 0);

    wei_curve_destroy(ec);
  }

  {
    const unsigned char expect1[32] = {
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    const unsigned char expect2[32] = {
      0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
    scalar_field_t *sc = &ec->sc;
    mp_limb_t r[MAX_SCALAR_LIMBS];
    mp_limb_t t[MAX_REDUCE_LIMBS];
    unsigned char raw[MAX_SCALAR_SIZE];

    sc_export(sc, raw, sc->n);

    mpn_zero(r, ARRAY_SIZE(r));
    mpn_zero(t, ARRAY_SIZE(t));

    mpn_copyi(t, sc->n, sc->limbs);

    sc_reduce(sc, r, t);

    ASSERT(sc_is_zero(sc, r));

    sc_neg(sc, r, r);

    ASSERT(sc_is_zero(sc, r));

    raw[0] -= 1;

    ASSERT(sc_import(sc, r, raw));

    raw[0] += 1;

    ASSERT(!sc_import(sc, r, raw));

    raw[0] += 1;

    ASSERT(!sc_import(sc, r, raw));

    sc_export(sc, raw, sc->n);

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(torsion_memcmp(raw, expect1, 32) == 0);

    sc_export(sc, raw, sc->n);

    raw[0] += 1;

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(torsion_memcmp(raw, expect2, 32) == 0);

    edwards_curve_destroy(ec);
  }

  {
    const unsigned char expect[32] = {
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
      0x45, 0x51, 0x23, 0x19, 0x50, 0xb7, 0x5f, 0xc4,
      0x40, 0x2d, 0xa1, 0x73, 0x2f, 0xc9, 0xbe, 0xbe
    };

    wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
    scalar_field_t *sc = &ec->sc;
    mp_limb_t r[MAX_SCALAR_LIMBS];
    unsigned char max[32];

    memset(max, 0xff, 32);

    ASSERT(!sc_import_reduce(sc, r, max));

    sc_export(sc, max, r);

    ASSERT(torsion_memcmp(max, expect, 32) == 0);

    wei_curve_destroy(ec);
  }

  {
    static const sc_t one = {1, 0};
    edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED1174);
    scalar_field_t *sc = &ec->sc;
    sc_t x, y;

    do {
      sc_random(sc, x, rng);
    } while (sc_equal(sc, x, one));

    ASSERT(sc_invert(sc, y, x));

    ASSERT(!sc_is_zero(sc, x));
    ASSERT(!sc_is_zero(sc, y));

    ASSERT(!sc_equal(sc, x, one));
    ASSERT(!sc_equal(sc, y, one));

    sc_mul(sc, x, x, y);

    ASSERT(sc_equal(sc, x, one));

    edwards_curve_destroy(ec);
  }
}

/*
 * Field Element
 */

static void
test_field_element(void) {
  wei_t *ec = wei_curve_create(WEI_CURVE_P256);
  unsigned char raw[MAX_FIELD_SIZE];
  prime_field_t *fe = &ec->fe;
  fe_t t;

  printf("  - Field element sanity check.\n");

  mpn_export(raw, fe->size, fe->p, fe->limbs, fe->endian);

  raw[fe->size - 1] -= 1;

  ASSERT(fe_import(fe, t, raw));

  raw[fe->size - 1] += 1;

  ASSERT(!fe_import(fe, t, raw));

  raw[7] += 1;

  ASSERT(!fe_import(fe, t, raw));

  wei_curve_destroy(ec);
}

/*
 * Utils
 */

static void
test_zero(drbg_t *rng) {
  printf("  - Zero bytes sanity check.\n");

  {
    const unsigned char zero[4] = {0, 0, 0, 0};
    const unsigned char one[4] = {0, 0, 0, 1};
    const unsigned char full[4] = {0xff, 0xff, 0xff, 0xff};

    ASSERT(bytes_zero(zero, 4));
    ASSERT(!bytes_zero(one, 4));
    ASSERT(!bytes_zero(full, 4));
  }

  {
    size_t i;

    for (i = 0; i < 1000; i++) {
      unsigned char a[32];

      drbg_generate(rng, a, sizeof(a));

      if (i % 10 == 0)
        memset(a, 0, sizeof(a));

      ASSERT(bytes_zero(a, 32) == stupid_bytes_zero(a, 32));
    }
  }
}

static void
test_lt(drbg_t *rng) {
  printf("  - LT sanity check.\n");

  {
    const unsigned char mod[4] = {4, 3, 2, 1};
    const unsigned char minus1[4] = {3, 3, 2, 1};
    const unsigned char plus1[4] = {5, 3, 2, 1};
    const unsigned char full[4] = {0xff, 0xff, 0xff, 0xff};

    ASSERT(bytes_lt(minus1, mod, 4));
    ASSERT(!bytes_lt(mod, mod, 4));
    ASSERT(!bytes_lt(plus1, mod, 4));
    ASSERT(bytes_lt(mod, full, 4));
    ASSERT(!bytes_lt(full, mod, 4));
  }

  {
    size_t i;

    for (i = 0; i < 1000; i++) {
      unsigned char a[32];
      unsigned char b[32];

      drbg_generate(rng, a, sizeof(a));
      drbg_generate(rng, b, sizeof(b));

      ASSERT(bytes_lt(a, a, 32) == 0);
      ASSERT(bytes_lt(b, b, 32) == 0);
      ASSERT(bytes_lt(a, b, 32) == (revcmp(a, b, 32) < 0));
      ASSERT(bytes_lt(b, a, 32) == (revcmp(b, a, 32) < 0));
    }
  }
}

static void
test_mpn_cmp(drbg_t *rng) {
  printf("  - MPN comparison sanity check.\n");

  {
    const mp_limb_t mod[4] = {4, 3, 2, 1};
    const mp_limb_t minus1[4] = {3, 3, 2, 1};
    const mp_limb_t plus1[4] = {5, 3, 2, 1};
    const mp_limb_t full[4] = {0xff, 0xff, 0xff, 0xff};

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

  {
    size_t i;

    for (i = 0; i < 1000; i++) {
      mp_limb_t a[4];
      mp_limb_t b[4];
      int cab, cba;

      drbg_generate(rng, a, sizeof(a));
      drbg_generate(rng, b, sizeof(b));

      ASSERT(mpn_sec_lt(a, a, 4) == 0);
      ASSERT(mpn_sec_lt(b, b, 4) == 0);
      ASSERT(mpn_sec_lt(a, b, 4) == (mpn_cmp(a, b, 4) < 0));
      ASSERT(mpn_sec_lt(b, a, 4) == (mpn_cmp(b, a, 4) < 0));

      ASSERT(mpn_sec_lte(a, a, 4) == 1);
      ASSERT(mpn_sec_lte(b, b, 4) == 1);
      ASSERT(mpn_sec_lte(a, b, 4) == (mpn_cmp(a, b, 4) <= 0));
      ASSERT(mpn_sec_lte(b, a, 4) == (mpn_cmp(b, a, 4) <= 0));

      ASSERT(mpn_sec_gt(a, a, 4) == 0);
      ASSERT(mpn_sec_gt(b, b, 4) == 0);
      ASSERT(mpn_sec_gt(a, b, 4) == (mpn_cmp(a, b, 4) > 0));
      ASSERT(mpn_sec_gt(b, a, 4) == (mpn_cmp(b, a, 4) > 0));

      ASSERT(mpn_sec_gte(a, a, 4) == 1);
      ASSERT(mpn_sec_gte(b, b, 4) == 1);
      ASSERT(mpn_sec_gte(a, b, 4) == (mpn_cmp(a, b, 4) >= 0));
      ASSERT(mpn_sec_gte(b, a, 4) == (mpn_cmp(b, a, 4) >= 0));

      cab = mpn_cmp(a, b, 4);
      cba = mpn_cmp(b, a, 4);

      if (cab < 0)
        cab = -1;
      else if (cab > 0)
        cab = 1;

      if (cba < 0)
        cba = -1;
      else if (cba > 0)
        cba = 1;

      ASSERT(mpn_sec_cmp(a, a, 4) == 0);
      ASSERT(mpn_sec_cmp(b, b, 4) == 0);
      ASSERT(mpn_sec_cmp(a, b, 4) == cab);
      ASSERT(mpn_sec_cmp(b, a, 4) == cba);
    }
  }
}

/*
 * ECC
 */

static void
test_wei_points_p256(drbg_t *rng) {
  const unsigned char g_raw[33] = {
    0x03, 0x6b, 0x17, 0xd1, 0xf2, 0xe1, 0x2c, 0x42,
    0x47, 0xf8, 0xbc, 0xe6, 0xe5, 0x63, 0xa4, 0x40,
    0xf2, 0x77, 0x03, 0x7d, 0x81, 0x2d, 0xeb, 0x33,
    0xa0, 0xf4, 0xa1, 0x39, 0x45, 0xd8, 0x98, 0xc2,
    0x96
  };

  const unsigned char g2_raw[33] = {
    0x03, 0x7c, 0xf2, 0x7b, 0x18, 0x8d, 0x03, 0x4f,
    0x7e, 0x8a, 0x52, 0x38, 0x03, 0x04, 0xb5, 0x1a,
    0xc3, 0xc0, 0x89, 0x69, 0xe2, 0x77, 0xf2, 0x1b,
    0x35, 0xa6, 0x0b, 0x48, 0xfc, 0x47, 0x66, 0x99,
    0x78
  };

  const unsigned char g3_raw[33] = {
    0x02, 0x5e, 0xcb, 0xe4, 0xd1, 0xa6, 0x33, 0x0a,
    0x44, 0xc8, 0xf7, 0xef, 0x95, 0x1d, 0x4b, 0xf1,
    0x65, 0xe6, 0xc6, 0xb7, 0x21, 0xef, 0xad, 0xa9,
    0x85, 0xfb, 0x41, 0x66, 0x1b, 0xc6, 0xe7, 0xfd,
    0x6c
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_P256);
  jge_t jg, jp, jq, jr;
  wge_t g, p, q, r;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char p_raw[33];
  size_t p_size;

  printf("  - Testing Weierstrass group law (P256).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  wge_set(ec, &g, &ec->g);

  jge_set_wge(ec, &jg, &ec->g);

  ASSERT(wge_import(ec, &p, g_raw, 33));

  jge_set_wge(ec, &jp, &p);
  jge_set_wge(ec, &jq, &ec->g);

  ASSERT(wge_validate(ec, &p));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(wge_equal(ec, &p, &ec->g));
  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(wge_import(ec, &q, g2_raw, 33));
  ASSERT(wge_import(ec, &r, g3_raw, 33));

  jge_set_wge(ec, &jq, &q);
  jge_set_wge(ec, &jr, &r);

  wge_dbl_var(ec, &p, &ec->g);

  ASSERT(wge_equal(ec, &p, &q));

  wge_add_var(ec, &p, &p, &ec->g);

  ASSERT(wge_equal(ec, &p, &r));

  jge_dbl_var(ec, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jq));

  jge_add_var(ec, &jp, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jr));

  jge_sub_var(ec, &jp, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jq));

  jge_mixed_add_var(ec, &jp, &jp, &g);

  ASSERT(jge_equal(ec, &jp, &jr));

  jge_mixed_sub_var(ec, &jp, &jp, &g);

  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(jge_validate(ec, &jg));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(jge_validate(ec, &jr));

  ASSERT(!jge_is_zero(ec, &jg));
  ASSERT(!jge_is_zero(ec, &jp));
  ASSERT(!jge_is_zero(ec, &jq));
  ASSERT(!jge_is_zero(ec, &jr));

  wge_set_jge(ec, &p, &jp);

  ASSERT(wge_equal(ec, &p, &q));

  ASSERT(wge_export(ec, p_raw, &p_size, &p, 1));
  ASSERT(p_size == 33);

  ASSERT(torsion_memcmp(p_raw, g2_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_points_p521(drbg_t *rng) {
  const unsigned char g_raw[67] = {
    0x02, 0x00, 0xc6, 0x85, 0x8e, 0x06, 0xb7, 0x04,
    0x04, 0xe9, 0xcd, 0x9e, 0x3e, 0xcb, 0x66, 0x23,
    0x95, 0xb4, 0x42, 0x9c, 0x64, 0x81, 0x39, 0x05,
    0x3f, 0xb5, 0x21, 0xf8, 0x28, 0xaf, 0x60, 0x6b,
    0x4d, 0x3d, 0xba, 0xa1, 0x4b, 0x5e, 0x77, 0xef,
    0xe7, 0x59, 0x28, 0xfe, 0x1d, 0xc1, 0x27, 0xa2,
    0xff, 0xa8, 0xde, 0x33, 0x48, 0xb3, 0xc1, 0x85,
    0x6a, 0x42, 0x9b, 0xf9, 0x7e, 0x7e, 0x31, 0xc2,
    0xe5, 0xbd, 0x66
  };

  const unsigned char g2_raw[67] = {
    0x02, 0x00, 0x43, 0x3c, 0x21, 0x90, 0x24, 0x27,
    0x7e, 0x7e, 0x68, 0x2f, 0xcb, 0x28, 0x81, 0x48,
    0xc2, 0x82, 0x74, 0x74, 0x03, 0x27, 0x9b, 0x1c,
    0xcc, 0x06, 0x35, 0x2c, 0x6e, 0x55, 0x05, 0xd7,
    0x69, 0xbe, 0x97, 0xb3, 0xb2, 0x04, 0xda, 0x6e,
    0xf5, 0x55, 0x07, 0xaa, 0x10, 0x4a, 0x3a, 0x35,
    0xc5, 0xaf, 0x41, 0xcf, 0x2f, 0xa3, 0x64, 0xd6,
    0x0f, 0xd9, 0x67, 0xf4, 0x3e, 0x39, 0x33, 0xba,
    0x6d, 0x78, 0x3d
  };

  const unsigned char g3_raw[67] = {
    0x03, 0x01, 0xa7, 0x3d, 0x35, 0x24, 0x43, 0xde,
    0x29, 0x19, 0x5d, 0xd9, 0x1d, 0x6a, 0x64, 0xb5,
    0x95, 0x94, 0x79, 0xb5, 0x2a, 0x6e, 0x5b, 0x12,
    0x3d, 0x9a, 0xb9, 0xe5, 0xad, 0x7a, 0x11, 0x2d,
    0x7a, 0x8d, 0xd1, 0xad, 0x3f, 0x16, 0x4a, 0x3a,
    0x48, 0x32, 0x05, 0x1d, 0xa6, 0xbd, 0x16, 0xb5,
    0x9f, 0xe2, 0x1b, 0xae, 0xb4, 0x90, 0x86, 0x2c,
    0x32, 0xea, 0x05, 0xa5, 0x91, 0x9d, 0x2e, 0xde,
    0x37, 0xad, 0x7d
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_P521);
  jge_t jg, jp, jq, jr;
  wge_t g, p, q, r;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char p_raw[67];
  size_t p_size;

  printf("  - Testing Weierstrass group law (P521).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  wge_set(ec, &g, &ec->g);

  jge_set_wge(ec, &jg, &ec->g);

  ASSERT(wge_import(ec, &p, g_raw, 67));

  jge_set_wge(ec, &jp, &p);
  jge_set_wge(ec, &jq, &ec->g);

  ASSERT(wge_validate(ec, &p));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(wge_equal(ec, &p, &ec->g));
  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(wge_import(ec, &q, g2_raw, 67));
  ASSERT(wge_import(ec, &r, g3_raw, 67));

  jge_set_wge(ec, &jq, &q);
  jge_set_wge(ec, &jr, &r);

  wge_dbl_var(ec, &p, &ec->g);

  ASSERT(wge_equal(ec, &p, &q));

  wge_add_var(ec, &p, &p, &ec->g);

  ASSERT(wge_equal(ec, &p, &r));

  jge_dbl_var(ec, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jq));

  jge_add_var(ec, &jp, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jr));

  jge_sub_var(ec, &jp, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jq));

  jge_mixed_add_var(ec, &jp, &jp, &g);

  ASSERT(jge_equal(ec, &jp, &jr));

  jge_mixed_sub_var(ec, &jp, &jp, &g);

  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(jge_validate(ec, &jg));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(jge_validate(ec, &jr));

  ASSERT(!jge_is_zero(ec, &jg));
  ASSERT(!jge_is_zero(ec, &jp));
  ASSERT(!jge_is_zero(ec, &jq));
  ASSERT(!jge_is_zero(ec, &jr));

  wge_set_jge(ec, &p, &jp);

  ASSERT(wge_equal(ec, &p, &q));

  ASSERT(wge_export(ec, p_raw, &p_size, &p, 1));
  ASSERT(p_size == 67);

  ASSERT(torsion_memcmp(p_raw, g2_raw, 67) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_points_secp256k1(drbg_t *rng) {
  const unsigned char g_raw[33] = {
    0x02, 0x79, 0xbe, 0x66, 0x7e, 0xf9, 0xdc, 0xbb,
    0xac, 0x55, 0xa0, 0x62, 0x95, 0xce, 0x87, 0x0b,
    0x07, 0x02, 0x9b, 0xfc, 0xdb, 0x2d, 0xce, 0x28,
    0xd9, 0x59, 0xf2, 0x81, 0x5b, 0x16, 0xf8, 0x17,
    0x98
  };

  const unsigned char g2_raw[33] = {
    0x02, 0xc6, 0x04, 0x7f, 0x94, 0x41, 0xed, 0x7d,
    0x6d, 0x30, 0x45, 0x40, 0x6e, 0x95, 0xc0, 0x7c,
    0xd8, 0x5c, 0x77, 0x8e, 0x4b, 0x8c, 0xef, 0x3c,
    0xa7, 0xab, 0xac, 0x09, 0xb9, 0x5c, 0x70, 0x9e,
    0xe5
  };

  const unsigned char g3_raw[33] = {
    0x02, 0xf9, 0x30, 0x8a, 0x01, 0x92, 0x58, 0xc3,
    0x10, 0x49, 0x34, 0x4f, 0x85, 0xf8, 0x9d, 0x52,
    0x29, 0xb5, 0x31, 0xc8, 0x45, 0x83, 0x6f, 0x99,
    0xb0, 0x86, 0x01, 0xf1, 0x13, 0xbc, 0xe0, 0x36,
    0xf9
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  jge_t jg, jp, jq, jr;
  wge_t g, p, q, r;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char p_raw[33];
  size_t p_size;

  printf("  - Testing Weierstrass group law (SECP256K1).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  wge_set(ec, &g, &ec->g);

  jge_set_wge(ec, &jg, &ec->g);

  ASSERT(wge_import(ec, &p, g_raw, 33));

  jge_set_wge(ec, &jp, &p);
  jge_set_wge(ec, &jq, &ec->g);

  ASSERT(wge_validate(ec, &p));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(wge_equal(ec, &p, &ec->g));
  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(wge_import(ec, &q, g2_raw, 33));
  ASSERT(wge_import(ec, &r, g3_raw, 33));

  jge_set_wge(ec, &jq, &q);
  jge_set_wge(ec, &jr, &r);

  wge_dbl_var(ec, &p, &ec->g);

  ASSERT(wge_equal(ec, &p, &q));

  wge_add_var(ec, &p, &p, &ec->g);

  ASSERT(wge_equal(ec, &p, &r));

  jge_dbl_var(ec, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jq));

  jge_add_var(ec, &jp, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jr));

  jge_sub_var(ec, &jp, &jp, &jg);

  ASSERT(jge_equal(ec, &jp, &jq));

  jge_mixed_add_var(ec, &jp, &jp, &g);

  ASSERT(jge_equal(ec, &jp, &jr));

  jge_mixed_sub_var(ec, &jp, &jp, &g);

  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(jge_validate(ec, &jg));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(jge_validate(ec, &jr));

  ASSERT(!jge_is_zero(ec, &jg));
  ASSERT(!jge_is_zero(ec, &jp));
  ASSERT(!jge_is_zero(ec, &jq));
  ASSERT(!jge_is_zero(ec, &jr));

  wge_set_jge(ec, &p, &jp);

  ASSERT(wge_equal(ec, &p, &q));

  ASSERT(wge_export(ec, p_raw, &p_size, &p, 1));
  ASSERT(p_size == 33);

  ASSERT(torsion_memcmp(p_raw, g2_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_mul_g_p256(drbg_t *rng) {
  const unsigned char k_raw[32] = {
    0x38, 0xf8, 0x62, 0x0b, 0xa6, 0x0b, 0xed, 0x7c,
    0xf9, 0x0c, 0x7a, 0x99, 0xac, 0x35, 0xa4, 0x4e,
    0x39, 0x27, 0x59, 0x8e, 0x3c, 0x99, 0xbb, 0xc5,
    0xf5, 0x70, 0x75, 0x13, 0xc4, 0x0e, 0x2c, 0xe3
  };

  const unsigned char expect_raw[33] = {
    0x02, 0x1a, 0xb3, 0x49, 0x34, 0xb8, 0x11, 0xb5,
    0x5e, 0x2f, 0xa4, 0xf1, 0xcd, 0x57, 0xf1, 0x68,
    0x51, 0x3d, 0x04, 0xb9, 0x45, 0xb0, 0x43, 0xec,
    0xe9, 0x6b, 0x25, 0x53, 0x96, 0x72, 0xff, 0x52,
    0x03
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_P256);
  scalar_field_t *sc = &ec->sc;
  wge_t q, expect;
  sc_t k;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul_g (P256).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(sc_import(sc, k, k_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wei_mul_g(ec, &q, k);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_mul_p256(drbg_t *rng) {
  const unsigned char p_raw[33] = {
    0x03, 0x42, 0x67, 0xab, 0xc7, 0xde, 0x72, 0x0f,
    0x14, 0x5a, 0xbc, 0x94, 0xb9, 0x5b, 0x33, 0x50,
    0x7a, 0x37, 0x55, 0x55, 0x2b, 0xef, 0xaf, 0x57,
    0x61, 0x33, 0x7a, 0xd6, 0x7a, 0x28, 0xa9, 0x08,
    0xa1
  };

  const unsigned char k_raw[32] = {
    0xfd, 0x37, 0xfe, 0xab, 0xd9, 0xdd, 0x8d, 0xe5,
    0xfd, 0x04, 0x79, 0xf4, 0xd6, 0xea, 0xd4, 0xe6,
    0x02, 0xc7, 0x06, 0x0f, 0x43, 0x6e, 0x2b, 0xf1,
    0xc0, 0x72, 0xe9, 0x91, 0x80, 0xcb, 0x09, 0x18
  };

  const unsigned char expect_raw[33] = {
    0x02, 0x93, 0xa3, 0x55, 0xe4, 0x8f, 0x3b, 0x74,
    0xcc, 0x3b, 0xcb, 0xb4, 0x6c, 0xb2, 0x84, 0x3a,
    0xd5, 0x4e, 0xe5, 0xe0, 0x45, 0xe9, 0x17, 0x0b,
    0x00, 0x45, 0xbc, 0xc2, 0x86, 0x68, 0x8c, 0x4d,
    0x56
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_P256);
  scalar_field_t *sc = &ec->sc;
  wge_t p, q, expect;
  sc_t k;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul (P256).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(wge_import(ec, &p, p_raw, 33));
  ASSERT(sc_import(sc, k, k_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &p));
  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wei_mul(ec, &q, &p, k);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_double_mul_p256(drbg_t *rng) {
  const unsigned char p_raw[33] = {
    0x02, 0x65, 0x26, 0x45, 0xad, 0x1a, 0x36, 0x8c,
    0xdc, 0xcf, 0x81, 0x90, 0x56, 0x3b, 0x2a, 0x12,
    0xba, 0x31, 0xea, 0x33, 0x78, 0xc2, 0x23, 0x66,
    0xff, 0xf8, 0x47, 0x92, 0x63, 0x8c, 0xb8, 0xc8,
    0x94
  };

  const unsigned char k1_raw[32] = {
    0x5f, 0xd3, 0x7e, 0x3c, 0x67, 0x9e, 0xc5, 0xd0,
    0x2b, 0xb6, 0x6a, 0xa8, 0x6e, 0x56, 0xd6, 0x40,
    0x65, 0xe9, 0x47, 0x74, 0x4e, 0x50, 0xee, 0xec,
    0x80, 0xcf, 0xcc, 0xce, 0x3b, 0xd2, 0xf2, 0x1a
  };

  const unsigned char k2_raw[32] = {
    0xfb, 0x15, 0x9a, 0x7d, 0x37, 0x4d, 0x24, 0xde,
    0xde, 0x0a, 0x55, 0xb2, 0x98, 0x26, 0xe3, 0x24,
    0xf6, 0xf1, 0xd7, 0x57, 0x36, 0x53, 0xd7, 0x8a,
    0x98, 0xed, 0xa2, 0x80, 0x6d, 0xbe, 0x37, 0x98
  };

  const unsigned char expect_raw[33] = {
    0x02, 0x96, 0xf1, 0xb9, 0xe3, 0xe7, 0x0b, 0xa1,
    0x2e, 0xaf, 0x40, 0x23, 0x05, 0x64, 0x5b, 0x0f,
    0x28, 0x1b, 0xec, 0x25, 0x4f, 0xf2, 0x31, 0x8f,
    0x96, 0x9c, 0x97, 0x96, 0x0c, 0x35, 0x0b, 0x2c,
    0x6d
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_P256);
  scalar_field_t *sc = &ec->sc;
  wge_t p, q, expect;
  sc_t k1, k2;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul_double_var (P256).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(wge_import(ec, &p, p_raw, 33));
  ASSERT(sc_import(sc, k1, k1_raw));
  ASSERT(sc_import(sc, k2, k2_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &p));
  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wei_mul_double_var(ec, &q, k1, &p, k2);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_multi_mul_p256(drbg_t *rng) {
  const unsigned char p1_raw[33] = {
    0x02, 0x2b, 0xbf, 0x66, 0x0c, 0x19, 0x5a, 0xea,
    0x82, 0x82, 0x2d, 0x2e, 0x69, 0xc4, 0x02, 0x6e,
    0x57, 0x9a, 0x00, 0xbe, 0xac, 0xc8, 0xa5, 0x7d,
    0x73, 0xb7, 0x3c, 0x76, 0x7c, 0x0c, 0xdc, 0x12,
    0x2d
  };

  const unsigned char p2_raw[33] = {
    0x03, 0x5e, 0x8e, 0x2b, 0xd9, 0x9b, 0x58, 0x2f,
    0x88, 0xcd, 0x11, 0xeb, 0x27, 0x85, 0x17, 0x71,
    0x2a, 0x09, 0x14, 0x7a, 0x5b, 0x3d, 0x97, 0xbb,
    0x10, 0x4b, 0x03, 0x47, 0x42, 0x78, 0x94, 0xde,
    0xff
  };

  const unsigned char k0_raw[32] = {
    0x7e, 0x7b, 0x6a, 0x4b, 0x5b, 0x57, 0xce, 0x1e,
    0xb4, 0x68, 0x69, 0x87, 0x76, 0x63, 0xb6, 0xdc,
    0xee, 0xb7, 0x0c, 0xba, 0xa5, 0x43, 0xab, 0x87,
    0x50, 0xd5, 0x0a, 0xf7, 0x3f, 0x21, 0xa1, 0x70
  };

  const unsigned char k1_raw[32] = {
    0x31, 0x73, 0xf8, 0x6d, 0xa3, 0xd4, 0x5e, 0xe1,
    0x17, 0x7a, 0xd9, 0xde, 0x61, 0x92, 0x51, 0xa4,
    0xf7, 0x7d, 0xbf, 0x39, 0x92, 0x7f, 0x56, 0xbb,
    0x67, 0x56, 0x8c, 0xd5, 0xab, 0xc9, 0x78, 0xc9
  };

  const unsigned char k2_raw[32] = {
    0x91, 0x58, 0x61, 0x18, 0x94, 0xb5, 0xb6, 0xea,
    0x63, 0x57, 0xf0, 0xd4, 0x24, 0x36, 0x58, 0x68,
    0x86, 0x09, 0x5f, 0x51, 0xe3, 0x18, 0xaf, 0x70,
    0x40, 0x6e, 0xfc, 0xa3, 0xa0, 0x63, 0x80, 0xb9
  };

  const unsigned char expect_raw[33] = {
    0x03, 0x70, 0x76, 0xa3, 0x4c, 0x39, 0xff, 0xe2,
    0x19, 0xfa, 0x48, 0xdf, 0xb1, 0xfd, 0x19, 0xa0,
    0x6a, 0x94, 0xa3, 0xe5, 0xc6, 0x12, 0xd7, 0xf9,
    0xee, 0x19, 0xba, 0x05, 0x0e, 0x9d, 0x20, 0x0d,
    0x79
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_P256);
  wei_scratch_t *scratch = wei_scratch_create(ec, 2);
  scalar_field_t *sc = &ec->sc;
  wge_t p1, p2, q, expect;
  sc_t k0, k1, k2;
  wge_t points[2];
  sc_t coeffs[2];

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul_multi_var (P256).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(wge_import(ec, &p1, p1_raw, 33));
  ASSERT(wge_import(ec, &p2, p2_raw, 33));
  ASSERT(sc_import(sc, k0, k0_raw));
  ASSERT(sc_import(sc, k1, k1_raw));
  ASSERT(sc_import(sc, k2, k2_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &p1));
  ASSERT(wge_validate(ec, &p2));
  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wge_set(ec, &points[0], &p1);
  wge_set(ec, &points[1], &p2);

  sc_set(sc, coeffs[0], k1);
  sc_set(sc, coeffs[1], k2);

  wei_mul_multi_var(ec, &q, k0, points, (const sc_t *)coeffs, 2, scratch);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_scratch_destroy(ec, scratch);
  wei_curve_destroy(ec);
}

static void
test_wei_mul_g_secp256k1(drbg_t *rng) {
  const unsigned char k_raw[32] = {
    0xf7, 0x6c, 0xd0, 0xed, 0xfb, 0x89, 0x7f, 0x39,
    0xf6, 0x3e, 0x2c, 0x16, 0xc5, 0x73, 0x81, 0xb7,
    0x23, 0x5a, 0x1c, 0x5c, 0xe7, 0x8a, 0xab, 0xf2,
    0xbd, 0x88, 0xb0, 0xd2, 0xb2, 0x5d, 0xf8, 0x52
  };

  const unsigned char expect_raw[33] = {
    0x03, 0x44, 0x81, 0x79, 0xe0, 0x0a, 0xb1, 0xbf,
    0x6f, 0x33, 0x92, 0xeb, 0x3c, 0xdb, 0xe0, 0xff,
    0xbc, 0x69, 0x40, 0xd0, 0x0d, 0x72, 0x0d, 0xcd,
    0x1a, 0xc7, 0x80, 0xab, 0x95, 0xb0, 0xe2, 0xab,
    0x82
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  scalar_field_t *sc = &ec->sc;
  wge_t q, expect;
  sc_t k;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul_g (SECP256K1).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(sc_import(sc, k, k_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wei_mul_g(ec, &q, k);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_mul_secp256k1(drbg_t *rng) {
  const unsigned char p_raw[33] = {
    0x02, 0x0a, 0xfc, 0xf3, 0x56, 0xdb, 0x98, 0x4f,
    0xa0, 0x33, 0x98, 0x35, 0xfe, 0xb4, 0xd9, 0x65,
    0x15, 0x82, 0xee, 0xdf, 0x7a, 0x90, 0x46, 0xa1,
    0x24, 0x85, 0xf5, 0x48, 0xec, 0x75, 0x6f, 0x1a,
    0xe1
  };

  const unsigned char k_raw[32] = {
    0x86, 0xdc, 0xca, 0x82, 0x3e, 0xab, 0x1e, 0x30,
    0x65, 0xbb, 0xe2, 0x58, 0xce, 0xd4, 0xb3, 0x17,
    0xe5, 0x5c, 0x8a, 0x9e, 0x9d, 0x4e, 0xb1, 0x5c,
    0x3e, 0x6e, 0x06, 0xe4, 0x19, 0xf7, 0xcb, 0x15
  };

  const unsigned char expect_raw[33] = {
    0x02, 0xe4, 0xbe, 0xe4, 0x40, 0xb7, 0xa4, 0x75,
    0xb6, 0x8c, 0x82, 0x1a, 0xbd, 0x37, 0x3d, 0x29,
    0x05, 0xd5, 0xb2, 0x52, 0x66, 0xdd, 0x68, 0x4a,
    0xe0, 0x6f, 0x31, 0xb1, 0x6e, 0x15, 0x7c, 0xaf,
    0xd2
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  scalar_field_t *sc = &ec->sc;
  wge_t p, q, expect;
  sc_t k;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul (SECP256K1).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(wge_import(ec, &p, p_raw, 33));
  ASSERT(sc_import(sc, k, k_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &p));
  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wei_mul(ec, &q, &p, k);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_double_mul_secp256k1(drbg_t *rng) {
  const unsigned char p_raw[33] = {
    0x02, 0x08, 0x0d, 0x86, 0xe7, 0x20, 0x89, 0x5f,
    0xde, 0x96, 0x2c, 0x69, 0x75, 0xff, 0x3a, 0x34,
    0x32, 0x0f, 0x3c, 0x33, 0x56, 0x25, 0x79, 0x9a,
    0xd6, 0x47, 0xf3, 0x7f, 0x2d, 0xe5, 0x2f, 0x53,
    0x36
  };

  const unsigned char k1_raw[32] = {
    0x54, 0x69, 0xaf, 0x11, 0xd5, 0x5d, 0xe5, 0x4e,
    0x35, 0xef, 0xaa, 0x84, 0xbb, 0x2c, 0xe7, 0xa8,
    0xac, 0xce, 0xb4, 0xe3, 0x74, 0x79, 0xd6, 0xdd,
    0xea, 0x7a, 0xf6, 0x13, 0x0d, 0x2d, 0x2a, 0x0b
  };

  const unsigned char k2_raw[32] = {
    0x21, 0x6d, 0x1f, 0xc1, 0xe2, 0x3a, 0x4c, 0xae,
    0x2f, 0x35, 0xda, 0xff, 0x69, 0xb0, 0x15, 0x2b,
    0x66, 0x11, 0x46, 0xa2, 0x2e, 0x82, 0xb7, 0x81,
    0xb5, 0xd0, 0x30, 0xef, 0xb1, 0xdd, 0xa1, 0x0b
  };

  const unsigned char expect_raw[33] = {
    0x03, 0x80, 0x2c, 0x75, 0x2d, 0x04, 0xea, 0x64,
    0x0b, 0x6a, 0x8e, 0x8d, 0x11, 0x82, 0x08, 0x9c,
    0xc3, 0x12, 0x14, 0xd2, 0xc5, 0x07, 0xc0, 0xac,
    0x70, 0x11, 0x21, 0x46, 0x1d, 0x42, 0xb1, 0x73,
    0x12
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  scalar_field_t *sc = &ec->sc;
  wge_t p, q, expect;
  sc_t k1, k2;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul_double_var (SECP256K1).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(wge_import(ec, &p, p_raw, 33));
  ASSERT(sc_import(sc, k1, k1_raw));
  ASSERT(sc_import(sc, k2, k2_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &p));
  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wei_mul_double_var(ec, &q, k1, &p, k2);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_curve_destroy(ec);
}

static void
test_wei_multi_mul_secp256k1(drbg_t *rng) {
  const unsigned char p1_raw[33] = {
    0x03, 0x4a, 0x9a, 0xd8, 0x9a, 0x7e, 0xef, 0x6b,
    0xbb, 0xb5, 0x10, 0xeb, 0x3b, 0x04, 0xa8, 0x6d,
    0xa7, 0x39, 0xfa, 0xf1, 0xea, 0xf0, 0xfc, 0x9f,
    0xef, 0x50, 0x36, 0x9f, 0x06, 0x3f, 0xae, 0x07,
    0xd4
  };

  const unsigned char p2_raw[33] = {
    0x03, 0xe8, 0x61, 0x6f, 0x91, 0x4f, 0x18, 0x4e,
    0xdd, 0x67, 0xe4, 0x9b, 0xeb, 0x03, 0x53, 0xb9,
    0x9c, 0x49, 0x82, 0x9c, 0xc9, 0x5b, 0x6e, 0x1e,
    0xb4, 0xa2, 0xc3, 0x7c, 0x67, 0x31, 0x02, 0x42,
    0xde
  };

  const unsigned char k0_raw[32] = {
    0x36, 0x1e, 0x23, 0xcb, 0x88, 0x78, 0x10, 0x7e,
    0xa5, 0xf7, 0xaf, 0xc6, 0xd6, 0x69, 0x47, 0xa9,
    0x0a, 0x56, 0x3b, 0x67, 0x2f, 0x13, 0xb2, 0x5c,
    0x89, 0x98, 0xec, 0x6c, 0xa1, 0x7d, 0xc6, 0xd2
  };

  const unsigned char k1_raw[32] = {
    0x66, 0x84, 0x8d, 0x52, 0xb2, 0x58, 0xfe, 0xe1,
    0xf2, 0xa5, 0xd7, 0x44, 0x0e, 0xd4, 0x40, 0x7a,
    0x04, 0x4e, 0xde, 0x94, 0xd2, 0x69, 0x2e, 0xd7,
    0x7d, 0x6b, 0x0d, 0x22, 0x1e, 0x0a, 0x98, 0x07
  };

  const unsigned char k2_raw[32] = {
    0xcc, 0x21, 0x46, 0x94, 0x47, 0x1a, 0x81, 0xae,
    0xe8, 0xa5, 0x5a, 0x3a, 0x8b, 0x2b, 0x6f, 0xa5,
    0x34, 0x65, 0x68, 0xd5, 0xee, 0xd2, 0x63, 0x2e,
    0x82, 0xcb, 0x7e, 0xf2, 0xb1, 0x96, 0x15, 0x29
  };

  const unsigned char expect_raw[33] = {
    0x03, 0xbc, 0xc1, 0x34, 0x04, 0xd2, 0x3e, 0x37,
    0x33, 0x7f, 0x15, 0xfe, 0xf0, 0x1d, 0xb4, 0xe1,
    0x31, 0xca, 0x1c, 0x0d, 0x14, 0x42, 0x3d, 0xbd,
    0x2e, 0xff, 0x72, 0xda, 0xd5, 0x60, 0xcb, 0xad,
    0xd6
  };

  wei_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  wei_scratch_t *scratch = wei_scratch_create(ec, 2);
  scalar_field_t *sc = &ec->sc;
  wge_t p1, p2, q, expect;
  sc_t k0, k1, k2;
  wge_t points[2];
  sc_t coeffs[2];

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[33];
  size_t q_size;

  printf("  - Testing wei_mul_multi_var (SECP256K1).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  ASSERT(wge_import(ec, &p1, p1_raw, 33));
  ASSERT(wge_import(ec, &p2, p2_raw, 33));
  ASSERT(sc_import(sc, k0, k0_raw));
  ASSERT(sc_import(sc, k1, k1_raw));
  ASSERT(sc_import(sc, k2, k2_raw));
  ASSERT(wge_import(ec, &expect, expect_raw, 33));

  ASSERT(wge_validate(ec, &p1));
  ASSERT(wge_validate(ec, &p2));
  ASSERT(wge_validate(ec, &expect));
  ASSERT(wge_equal(ec, &expect, &expect));
  ASSERT(!wge_equal(ec, &expect, &ec->g));

  wge_set(ec, &points[0], &p1);
  wge_set(ec, &points[1], &p2);

  sc_set(sc, coeffs[0], k1);
  sc_set(sc, coeffs[1], k2);

  wei_mul_multi_var(ec, &q, k0, points, (const sc_t *)coeffs, 2, scratch);

  ASSERT(wge_equal(ec, &q, &expect));

  ASSERT(wge_export(ec, q_raw, &q_size, &q, 1));
  ASSERT(q_size == 33);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 33) == 0);

  wei_scratch_destroy(ec, scratch);
  wei_curve_destroy(ec);
}

static void
test_mont_points_x25519(void) {
  const unsigned char g_raw[32] = {
    0x09, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  };

  const unsigned char g2_raw[32] = {
    0xfb, 0x4e, 0x68, 0xdd, 0x9c, 0x46, 0xae, 0x5c,
    0x5c, 0x0b, 0x35, 0x1e, 0xed, 0x5c, 0x3f, 0x8f,
    0x14, 0x71, 0x15, 0x7d, 0x68, 0x0c, 0x75, 0xd9,
    0xb7, 0xf1, 0x73, 0x18, 0xd5, 0x42, 0xd3, 0x20
  };

  const unsigned char g3_raw[32] = {
    0x12, 0x3c, 0x71, 0xfb, 0xaf, 0x03, 0x0a, 0xc0,
    0x59, 0x08, 0x1c, 0x62, 0x67, 0x4e, 0x82, 0xf8,
    0x64, 0xba, 0x1b, 0xc2, 0x91, 0x4d, 0x53, 0x45,
    0xe6, 0xab, 0x57, 0x6d, 0x1a, 0xbc, 0x12, 0x1c
  };

  mont_t *ec = mont_curve_create(MONT_CURVE_X25519);
  pge_t jg, jp, jq, jr;
  mge_t g, p, q, r;

  unsigned char p_raw[32];

  printf("  - Testing Montgomery group law (X25519).\n");

  mge_set(ec, &g, &ec->g);

  pge_set_mge(ec, &jg, &ec->g);

  ASSERT(mge_import(ec, &p, g_raw, 0));

  pge_set_mge(ec, &jp, &p);
  pge_set_mge(ec, &jq, &ec->g);

  ASSERT(mge_validate(ec, &p));
  ASSERT(pge_validate(ec, &jp));
  ASSERT(pge_validate(ec, &jq));
  ASSERT(mge_equal(ec, &p, &ec->g));
  ASSERT(pge_equal(ec, &jp, &jq));

  ASSERT(mge_import(ec, &q, g2_raw, 0));
  ASSERT(mge_import(ec, &r, g3_raw, 0));

  pge_set_mge(ec, &jq, &q);
  pge_set_mge(ec, &jr, &r);

  mge_dbl(ec, &p, &ec->g);

  ASSERT(mge_equal(ec, &p, &q));

  mge_add(ec, &p, &p, &ec->g);

  ASSERT(mge_equal(ec, &p, &r));

  pge_dbl(ec, &jp, &jg);

  ASSERT(pge_equal(ec, &jp, &jq));

  ASSERT(pge_validate(ec, &jg));
  ASSERT(pge_validate(ec, &jp));
  ASSERT(pge_validate(ec, &jq));
  ASSERT(pge_validate(ec, &jr));

  ASSERT(!pge_is_zero(ec, &jg));
  ASSERT(!pge_is_zero(ec, &jp));
  ASSERT(!pge_is_zero(ec, &jq));
  ASSERT(!pge_is_zero(ec, &jr));

  mge_set_pge(ec, &p, &jp, 0);

  ASSERT(mge_equal(ec, &p, &q));

  ASSERT(mge_export(ec, p_raw, &p));

  ASSERT(torsion_memcmp(p_raw, g2_raw, 32) == 0);

  mont_curve_destroy(ec);
}

static void
test_edwards_points_ed25519(drbg_t *rng) {
  const unsigned char g_raw[32] = {
    0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
  };

  const unsigned char g2_raw[32] = {
    0xc9, 0xa3, 0xf8, 0x6a, 0xae, 0x46, 0x5f, 0x0e,
    0x56, 0x51, 0x38, 0x64, 0x51, 0x0f, 0x39, 0x97,
    0x56, 0x1f, 0xa2, 0xc9, 0xe8, 0x5e, 0xa2, 0x1d,
    0xc2, 0x29, 0x23, 0x09, 0xf3, 0xcd, 0x60, 0x22
  };

  const unsigned char g3_raw[32] = {
    0xd4, 0xb4, 0xf5, 0x78, 0x48, 0x68, 0xc3, 0x02,
    0x04, 0x03, 0x24, 0x67, 0x17, 0xec, 0x16, 0x9f,
    0xf7, 0x9e, 0x26, 0x60, 0x8e, 0xa1, 0x26, 0xa1,
    0xab, 0x69, 0xee, 0x77, 0xd1, 0xb1, 0x67, 0x12
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  xge_t jg, jp, jq, jr;
  xge_t g, p, q, r;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char p_raw[32];

  printf("  - Testing Edwards group law (ED25519).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  edwards_randomize(ec, entropy);

  xge_set(ec, &g, &ec->g);
  xge_set(ec, &jg, &ec->g);

  ASSERT(xge_import(ec, &p, g_raw));

  xge_set(ec, &jp, &p);
  xge_set(ec, &jq, &ec->g);

  ASSERT(xge_validate(ec, &p));
  ASSERT(xge_validate(ec, &jp));
  ASSERT(xge_validate(ec, &jq));
  ASSERT(xge_equal(ec, &p, &ec->g));
  ASSERT(xge_equal(ec, &jp, &jq));

  ASSERT(xge_import(ec, &q, g2_raw));
  ASSERT(xge_import(ec, &r, g3_raw));

  xge_set(ec, &jq, &q);
  xge_set(ec, &jr, &r);

  xge_dbl(ec, &p, &ec->g);

  ASSERT(xge_equal(ec, &p, &q));

  xge_add(ec, &p, &p, &ec->g);

  ASSERT(xge_equal(ec, &p, &r));

  xge_dbl(ec, &jp, &jg);

  ASSERT(xge_equal(ec, &jp, &jq));

  xge_add(ec, &jp, &jp, &jg);

  ASSERT(xge_equal(ec, &jp, &jr));

  xge_sub(ec, &jp, &jp, &jg);

  ASSERT(xge_equal(ec, &jp, &jq));

  xge_add(ec, &jp, &jp, &jg);

  ASSERT(xge_equal(ec, &jp, &jr));

  xge_sub(ec, &jp, &jp, &jg);

  ASSERT(xge_equal(ec, &jp, &jq));

  ASSERT(xge_validate(ec, &jg));
  ASSERT(xge_validate(ec, &jp));
  ASSERT(xge_validate(ec, &jq));
  ASSERT(xge_validate(ec, &jr));

  ASSERT(!xge_is_zero(ec, &jg));
  ASSERT(!xge_is_zero(ec, &jp));
  ASSERT(!xge_is_zero(ec, &jq));
  ASSERT(!xge_is_zero(ec, &jr));

  xge_set(ec, &p, &jp);

  ASSERT(xge_equal(ec, &p, &q));

  xge_export(ec, p_raw, &p);

  ASSERT(torsion_memcmp(p_raw, g2_raw, 32) == 0);

  edwards_curve_destroy(ec);
}

static void
test_edwards_mul_g_ed25519(drbg_t *rng) {
  const unsigned char k_raw[32] = {
    0x69, 0x33, 0xe6, 0x67, 0x60, 0xa1, 0x4f, 0x93,
    0xa7, 0x75, 0x24, 0xb9, 0x07, 0xe4, 0x03, 0x75,
    0xf5, 0x4f, 0xd6, 0x0d, 0xec, 0x03, 0x3f, 0x1b,
    0x74, 0xf1, 0x54, 0x6d, 0xf2, 0xd7, 0xe8, 0x09
  };

  const unsigned char expect_raw[32] = {
    0x3c, 0xcf, 0xe3, 0x87, 0xa1, 0x5e, 0x19, 0x62,
    0x7e, 0xa5, 0x2f, 0x01, 0x3c, 0x42, 0x8e, 0xb3,
    0x50, 0x7c, 0xa8, 0x2f, 0x19, 0x90, 0x61, 0xd2,
    0xcb, 0xe1, 0x07, 0x8f, 0x01, 0x21, 0xaf, 0x21
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  scalar_field_t *sc = &ec->sc;
  xge_t q, expect;
  sc_t k;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[32];

  printf("  - Testing edwards_mul_g (ED25519).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  edwards_randomize(ec, entropy);

  ASSERT(sc_import(sc, k, k_raw));
  ASSERT(xge_import(ec, &expect, expect_raw));

  ASSERT(xge_validate(ec, &expect));
  ASSERT(xge_equal(ec, &expect, &expect));
  ASSERT(!xge_equal(ec, &expect, &ec->g));

  edwards_mul_g(ec, &q, k);

  ASSERT(xge_equal(ec, &q, &expect));

  xge_export(ec, q_raw, &q);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 32) == 0);

  edwards_curve_destroy(ec);
}

static void
test_edwards_mul_ed25519(drbg_t *rng) {
  const unsigned char p_raw[32] = {
    0xe1, 0xb0, 0x38, 0x25, 0x32, 0xd6, 0x5a, 0x8f,
    0x8d, 0xb3, 0x9a, 0x07, 0xff, 0x2a, 0x80, 0x9c,
    0x10, 0xbe, 0x39, 0xee, 0xf8, 0x6f, 0x85, 0x3b,
    0x18, 0xcb, 0x8d, 0xff, 0x52, 0x5d, 0x18, 0xea
  };

  const unsigned char k_raw[32] = {
    0x80, 0x41, 0xd8, 0xb0, 0xc9, 0x39, 0x74, 0x5f,
    0xd7, 0x86, 0xa0, 0x39, 0x50, 0x5b, 0x60, 0x14,
    0xf0, 0xa7, 0x2f, 0xf6, 0x54, 0x17, 0x82, 0x50,
    0x02, 0xde, 0x22, 0xae, 0x95, 0xd6, 0x6f, 0x04
  };

  const unsigned char expect_raw[32] = {
    0xca, 0x61, 0x5f, 0x44, 0xdf, 0xb4, 0xae, 0x00,
    0x53, 0xff, 0x3d, 0xb7, 0x1f, 0x7a, 0xbe, 0xfd,
    0x52, 0xf5, 0x8f, 0x36, 0x2c, 0xd5, 0x7c, 0x04,
    0x68, 0x2a, 0xee, 0xb0, 0x12, 0xec, 0x90, 0xbe
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  scalar_field_t *sc = &ec->sc;
  xge_t p, q, expect;
  sc_t k;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[32];

  printf("  - Testing edwards_mul (ED25519).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  edwards_randomize(ec, entropy);

  ASSERT(xge_import(ec, &p, p_raw));
  ASSERT(sc_import(sc, k, k_raw));
  ASSERT(xge_import(ec, &expect, expect_raw));

  ASSERT(xge_validate(ec, &p));
  ASSERT(xge_validate(ec, &expect));
  ASSERT(xge_equal(ec, &expect, &expect));
  ASSERT(!xge_equal(ec, &expect, &ec->g));

  edwards_mul(ec, &q, &p, k);

  ASSERT(xge_equal(ec, &q, &expect));

  xge_export(ec, q_raw, &q);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 32) == 0);

  edwards_curve_destroy(ec);
}

static void
test_edwards_double_mul_ed25519(drbg_t *rng) {
  const unsigned char p_raw[32] = {
    0x3e, 0x5f, 0x12, 0xb3, 0xec, 0x9e, 0x60, 0x1c,
    0x2b, 0x6c, 0x5b, 0x5b, 0xe6, 0x31, 0x97, 0xc0,
    0xc0, 0xef, 0x70, 0x61, 0x0a, 0xa6, 0x59, 0xef,
    0x1c, 0xdd, 0x6d, 0x0d, 0xac, 0x5b, 0x75, 0x5d
  };

  const unsigned char k1_raw[32] = {
    0xaa, 0x78, 0x33, 0xea, 0x1b, 0x8c, 0xd3, 0x7e,
    0x8c, 0x4b, 0xbc, 0xea, 0x75, 0x88, 0x53, 0x1d,
    0x45, 0x11, 0x5c, 0x03, 0x4f, 0xee, 0xfc, 0x98,
    0x51, 0x54, 0x37, 0x3e, 0xe5, 0x59, 0xe8, 0x03
  };

  const unsigned char k2_raw[32] = {
    0x87, 0xd3, 0x3b, 0x4a, 0x4c, 0x11, 0x08, 0xf2,
    0x2d, 0x6f, 0x53, 0xd5, 0xfd, 0x7a, 0x14, 0x8e,
    0xc3, 0x68, 0x63, 0x6d, 0xc1, 0x80, 0x0d, 0xd9,
    0x28, 0xd7, 0x05, 0x49, 0x18, 0x0c, 0x9c, 0x01
  };

  const unsigned char expect_raw[32] = {
    0x1c, 0x54, 0xae, 0xd8, 0xa6, 0x76, 0x2d, 0xde,
    0x34, 0x28, 0xe8, 0xdc, 0x59, 0xdf, 0x54, 0x29,
    0x92, 0x96, 0xb3, 0xb6, 0xa5, 0xb6, 0x7c, 0x78,
    0xf9, 0xc3, 0x7b, 0xb1, 0x26, 0xa3, 0x7f, 0x50
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  scalar_field_t *sc = &ec->sc;
  xge_t p, q, expect;
  sc_t k1, k2;

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[32];

  printf("  - Testing edwards_mul_double_var (ED25519).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  edwards_randomize(ec, entropy);

  ASSERT(xge_import(ec, &p, p_raw));
  ASSERT(sc_import(sc, k1, k1_raw));
  ASSERT(sc_import(sc, k2, k2_raw));
  ASSERT(xge_import(ec, &expect, expect_raw));

  ASSERT(xge_validate(ec, &p));
  ASSERT(xge_validate(ec, &expect));
  ASSERT(xge_equal(ec, &expect, &expect));
  ASSERT(!xge_equal(ec, &expect, &ec->g));

  edwards_mul_double_var(ec, &q, k1, &p, k2);

  ASSERT(xge_equal(ec, &q, &expect));

  xge_export(ec, q_raw, &q);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 32) == 0);

  edwards_curve_destroy(ec);
}

static void
test_edwards_multi_mul_ed25519(drbg_t *rng) {
  const unsigned char p1_raw[32] = {
    0xfa, 0xa1, 0xa7, 0xbc, 0xc6, 0x40, 0x81, 0xe9,
    0xc4, 0x9b, 0x28, 0x49, 0x8a, 0x99, 0x50, 0xbb,
    0x58, 0x11, 0x4a, 0xd0, 0xf6, 0xf9, 0x1e, 0x11,
    0xfb, 0xf1, 0xce, 0xb3, 0x28, 0xf6, 0xac, 0xce
  };

  const unsigned char p2_raw[32] = {
    0x70, 0x21, 0x0f, 0xd8, 0xd5, 0x32, 0x0f, 0xcc,
    0xd1, 0x9e, 0x97, 0x96, 0x4f, 0x58, 0xd2, 0xf8,
    0x1e, 0x3c, 0x3e, 0x14, 0xfa, 0xed, 0xaa, 0xe2,
    0x32, 0x8f, 0x5d, 0x9a, 0x4e, 0x86, 0x5d, 0xc1
  };

  const unsigned char k0_raw[32] = {
    0x1b, 0xaf, 0x38, 0x8f, 0x81, 0x39, 0x2c, 0xb6,
    0x9d, 0xfe, 0xd8, 0xca, 0x76, 0x3d, 0xfd, 0xa0,
    0xb0, 0x3a, 0x90, 0xfa, 0xf3, 0x76, 0xd2, 0x94,
    0xff, 0x2a, 0x57, 0x44, 0xf3, 0x36, 0x59, 0x0f
  };

  const unsigned char k1_raw[32] = {
    0xa9, 0x4e, 0x67, 0x1e, 0xd0, 0xee, 0xa4, 0xe0,
    0x96, 0xd4, 0xa9, 0x10, 0x84, 0x28, 0xbb, 0x03,
    0x58, 0x00, 0x1c, 0xc0, 0x4a, 0x02, 0xd9, 0x02,
    0x75, 0xf9, 0x88, 0x94, 0x56, 0x06, 0x49, 0x0f
  };

  const unsigned char k2_raw[32] = {
    0xc3, 0x0f, 0xd1, 0x50, 0x6f, 0x29, 0x89, 0x7f,
    0x6e, 0xc8, 0x52, 0x04, 0x9e, 0xa4, 0x9f, 0xb0,
    0x46, 0x7e, 0xb2, 0x17, 0xb8, 0x4c, 0xd2, 0xbe,
    0x4e, 0xac, 0xd1, 0xe3, 0xeb, 0x7f, 0xda, 0x02
  };

  const unsigned char expect_raw[32] = {
    0x4f, 0x05, 0x88, 0x00, 0xfc, 0x45, 0xcb, 0xfe,
    0x75, 0x06, 0x48, 0x21, 0x23, 0xc6, 0xf7, 0x4e,
    0xbd, 0x44, 0x4f, 0x62, 0x4b, 0x6e, 0x28, 0x6a,
    0x36, 0x1d, 0x25, 0xea, 0x2a, 0xd5, 0x48, 0x86
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  edwards_scratch_t *scratch = edwards_scratch_create(ec, 2);
  scalar_field_t *sc = &ec->sc;
  xge_t p1, p2, q, expect;
  sc_t k0, k1, k2;
  xge_t points[2];
  sc_t coeffs[2];

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char q_raw[32];

  printf("  - Testing edwards_mul_multi_var (ED25519).\n");

  drbg_generate(rng, entropy, sizeof(entropy));

  edwards_randomize(ec, entropy);

  ASSERT(xge_import(ec, &p1, p1_raw));
  ASSERT(xge_import(ec, &p2, p2_raw));
  ASSERT(sc_import(sc, k0, k0_raw));
  ASSERT(sc_import(sc, k1, k1_raw));
  ASSERT(sc_import(sc, k2, k2_raw));
  ASSERT(xge_import(ec, &expect, expect_raw));

  ASSERT(xge_validate(ec, &p1));
  ASSERT(xge_validate(ec, &p2));
  ASSERT(xge_validate(ec, &expect));
  ASSERT(xge_equal(ec, &expect, &expect));
  ASSERT(!xge_equal(ec, &expect, &ec->g));

  xge_set(ec, &points[0], &p1);
  xge_set(ec, &points[1], &p2);

  sc_set(sc, coeffs[0], k1);
  sc_set(sc, coeffs[1], k2);

  edwards_mul_multi_var(ec, &q, k0, points, (const sc_t *)coeffs, 2, scratch);

  ASSERT(xge_equal(ec, &q, &expect));

  xge_export(ec, q_raw, &q);

  ASSERT(torsion_memcmp(q_raw, expect_raw, 32) == 0);

  edwards_scratch_destroy(ec, scratch);
  edwards_curve_destroy(ec);
}

/*
 * Ristretto
 */

static void
test_ristretto_basepoint_multiples(drbg_t *unused) {
  /* https://ristretto.group/test_vectors/ristretto255.html */
  static const unsigned char multiples[][32] = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0},
    {226, 242, 174, 10, 106, 188, 78, 113, 168, 132, 169, 97, 197, 0, 81, 95,
     88, 227, 11, 106, 165, 130, 221, 141, 182, 166, 89, 69, 224, 141, 45, 118},
    {106, 73, 50, 16, 247, 73, 156, 209, 127, 236, 181, 16, 174, 12, 234, 35,
     161, 16, 232, 213, 185, 1, 248, 172, 173, 211, 9, 92, 115, 163, 185, 25},
    {148, 116, 31, 93, 93, 82, 117, 94, 206, 79, 35, 240, 68, 238, 39, 213,
     209, 234, 30, 43, 209, 150, 180, 98, 22, 107, 22, 21, 42, 157, 2, 89},
    {218, 128, 134, 39, 115, 53, 139, 70, 111, 250, 223, 224, 179, 41, 58, 179,
     217, 253, 83, 197, 234, 108, 149, 83, 88, 245, 104, 50, 45, 175, 106, 87},
    {232, 130, 177, 49, 1, 107, 82, 193, 211, 51, 112, 128, 24, 124, 247, 104,
     66, 62, 252, 203, 181, 23, 187, 73, 90, 184, 18, 196, 22, 15, 244, 78},
    {246, 71, 70, 211, 201, 43, 19, 5, 14, 216, 216, 2, 54, 167, 240, 0, 124,
     59, 63, 150, 47, 91, 167, 147, 209, 154, 96, 30, 187, 29, 244, 3},
    {68, 245, 53, 32, 146, 110, 200, 31, 189, 90, 56, 120, 69, 190, 183, 223,
     133, 169, 106, 36, 236, 225, 135, 56, 189, 207, 166, 167, 130, 42, 23,
     109},
    {144, 50, 147, 216, 242, 40, 126, 190, 16, 226, 55, 77, 193, 165, 62, 11,
     200, 135, 229, 146, 105, 159, 2, 208, 119, 213, 38, 60, 221, 85, 96, 28},
    {2, 98, 42, 206, 143, 115, 3, 163, 28, 175, 198, 63, 143, 196, 143, 220,
     22, 225, 200, 200, 210, 52, 178, 240, 214, 104, 82, 130, 169, 7, 96, 49},
    {32, 112, 111, 215, 136, 178, 114, 10, 30, 210, 165, 218, 212, 149, 43, 1,
     244, 19, 188, 240, 231, 86, 77, 232, 205, 200, 22, 104, 158, 45, 185, 95},
    {188, 232, 63, 139, 165, 221, 47, 165, 114, 134, 76, 36, 186, 24, 16, 249,
     82, 43, 198, 0, 74, 254, 149, 135, 122, 199, 50, 65, 202, 253, 171, 66},
    {228, 84, 158, 225, 107, 154, 160, 48, 153, 202, 32, 140, 103, 173, 175,
     202, 250, 76, 63, 62, 78, 83, 3, 222, 96, 38, 227, 202, 143, 248, 68, 96},
    {170, 82, 224, 0, 223, 46, 22, 245, 95, 177, 3, 47, 195, 59, 196, 39, 66,
     218, 214, 189, 90, 143, 192, 190, 1, 103, 67, 108, 89, 72, 80, 31},
    {70, 55, 107, 128, 244, 9, 178, 157, 194, 181, 246, 240, 197, 37, 145, 153,
     8, 150, 229, 113, 111, 65, 71, 124, 211, 0, 133, 171, 127, 16, 48, 30},
    {224, 196, 24, 247, 200, 217, 196, 205, 215, 57, 91, 147, 234, 18, 79, 58,
     217, 144, 33, 187, 104, 29, 252, 51, 2, 169, 217, 154, 46, 83, 230, 78}
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char tmp[32];
  rge_t p, q;
  xge_t a, b;
  size_t i;

  (void)unused;

  printf("  - Testing ristretto_basepoint_multiples (ED25519).\n");

  xge_zero(ec, &p);

  for (i = 0; i < ARRAY_SIZE(multiples); i++) {
    const unsigned char *raw = multiples[i];

    ASSERT(rge_import(ec, &q, raw));
    ASSERT(rge_equal(ec, &q, &p));

    ASSERT(xge_validate(ec, &q));

    xge_mulh(ec, &a, &p);
    xge_mulh(ec, &b, &q);

    ASSERT(xge_equal(ec, &a, &b));

    rge_export(ec, tmp, &p);

    ASSERT(torsion_memcmp(tmp, raw, 32) == 0);

    rge_export(ec, tmp, &q);

    ASSERT(torsion_memcmp(tmp, raw, 32) == 0);

    xge_add(ec, &p, &p, &ec->g);
  }

  edwards_curve_destroy(ec);
}

static void
test_ristretto_elligator(drbg_t *unused) {
  /* https://ristretto.group/test_vectors/ristretto255.html */
  static const unsigned char bytes[][32] = {
    {184, 249, 135, 49, 253, 123, 89, 113, 67, 160, 6, 239, 7, 105, 211, 41,
     192, 249, 185, 57, 9, 102, 70, 198, 15, 127, 7, 26, 160, 102, 134, 71},
    {229, 14, 241, 227, 75, 9, 118, 60, 128, 153, 226, 21, 183, 217, 91, 136,
     98, 0, 231, 156, 124, 77, 82, 139, 142, 134, 164, 169, 169, 62, 250, 52},
    {115, 109, 36, 220, 180, 223, 99, 6, 204, 169, 19, 29, 169, 68, 84, 23, 21,
     109, 189, 149, 127, 205, 91, 102, 172, 35, 112, 35, 134, 69, 186, 34},
    {16, 49, 96, 107, 171, 199, 164, 9, 129, 16, 64, 62, 241, 63, 132, 173,
     209, 160, 112, 215, 105, 50, 157, 81, 253, 105, 1, 154, 229, 25, 120, 83},
    {156, 131, 161, 162, 236, 251, 5, 187, 167, 171, 17, 178, 148, 210, 90,
     207, 86, 21, 79, 161, 167, 215, 234, 1, 136, 242, 182, 248, 38, 85, 79,
     86},
    {251, 177, 124, 54, 18, 101, 75, 235, 245, 186, 19, 46, 133, 157, 229, 64,
     10, 136, 181, 185, 78, 144, 254, 167, 137, 49, 107, 10, 61, 10, 21, 25},
    {232, 193, 20, 68, 240, 77, 186, 77, 183, 40, 44, 86, 150, 31, 198, 212,
     76, 81, 3, 217, 197, 8, 126, 128, 126, 152, 164, 208, 153, 44, 189, 77},
    {173, 229, 149, 177, 37, 230, 30, 69, 61, 56, 172, 190, 219, 115, 167, 194,
     71, 134, 59, 75, 28, 244, 118, 26, 162, 97, 64, 16, 15, 189, 30, 64},
    {106, 71, 61, 107, 250, 117, 42, 151, 91, 202, 212, 100, 52, 188, 190, 21,
     125, 218, 31, 18, 253, 241, 160, 133, 57, 242, 3, 164, 189, 68, 111, 75},
    {112, 204, 182, 90, 220, 198, 120, 73, 173, 107, 193, 17, 227, 40, 162, 36,
     150, 141, 235, 55, 172, 183, 12, 39, 194, 136, 43, 153, 244, 118, 91, 89},
    {111, 24, 203, 123, 254, 189, 11, 162, 51, 196, 163, 136, 204, 143, 10,
     222, 33, 112, 81, 205, 34, 35, 8, 66, 90, 6, 164, 58, 170, 177, 34, 25},
    {225, 183, 30, 52, 236, 82, 6, 183, 109, 25, 227, 181, 25, 82, 41, 193, 80,
     77, 161, 80, 242, 203, 79, 204, 136, 245, 131, 110, 237, 106, 3, 58},
    {207, 246, 38, 56, 30, 86, 176, 90, 27, 200, 61, 42, 221, 27, 56, 210, 79,
     178, 189, 120, 68, 193, 120, 167, 77, 185, 53, 197, 124, 128, 191, 126},
    {1, 136, 215, 80, 240, 46, 63, 147, 16, 244, 230, 207, 82, 189, 74, 50,
     106, 169, 138, 86, 30, 131, 214, 202, 166, 125, 251, 228, 98, 24, 36, 21},
    {210, 207, 228, 56, 155, 116, 207, 54, 84, 195, 251, 215, 249, 199, 116,
     75, 109, 239, 196, 251, 194, 246, 252, 228, 70, 146, 156, 35, 25, 39, 241,
     4},
    {34, 116, 123, 9, 8, 40, 93, 189, 9, 103, 57, 103, 66, 227, 3, 2, 157, 107,
     134, 219, 202, 74, 230, 154, 78, 107, 219, 195, 214, 14, 84, 80}
  };

  static const unsigned char images[][32] = {
    {176, 157, 237, 97, 66, 29, 140, 166, 168, 94, 26, 157, 212, 216, 229, 160,
     195, 246, 232, 239, 169, 112, 63, 193, 64, 32, 152, 69, 11, 190, 246, 86},
    {234, 141, 77, 203, 181, 225, 250, 74, 171, 62, 15, 118, 78, 212, 150, 19,
     131, 14, 188, 238, 194, 244, 141, 138, 166, 162, 83, 122, 228, 201, 19,
     26},
    {232, 231, 51, 92, 5, 168, 80, 36, 173, 179, 104, 68, 186, 149, 68, 40,
     140, 170, 27, 103, 99, 140, 21, 242, 43, 62, 250, 134, 208, 255, 61, 89},
    {208, 120, 140, 129, 177, 179, 237, 159, 252, 160, 28, 13, 206, 5, 211,
     241, 192, 218, 1, 97, 130, 241, 20, 169, 119, 46, 246, 29, 79, 80, 77, 84},
    {202, 11, 236, 145, 58, 12, 181, 157, 209, 6, 213, 88, 75, 147, 11, 119,
     191, 139, 47, 142, 33, 36, 153, 193, 223, 183, 178, 8, 205, 120, 248, 110},
    {26, 66, 231, 67, 203, 175, 116, 130, 32, 136, 62, 253, 215, 46, 5, 214,
     166, 248, 108, 237, 216, 71, 244, 173, 72, 133, 82, 6, 143, 240, 104, 41},
    {40, 157, 102, 96, 201, 223, 200, 197, 150, 181, 106, 83, 103, 126, 143,
     33, 145, 230, 78, 6, 171, 146, 210, 143, 112, 5, 245, 23, 183, 138, 18,
     120},
    {220, 37, 27, 203, 239, 196, 176, 131, 37, 66, 188, 243, 185, 250, 113, 23,
     167, 211, 154, 243, 168, 215, 54, 171, 159, 36, 195, 81, 13, 150, 43, 43},
    {232, 121, 176, 222, 183, 196, 159, 90, 238, 193, 105, 52, 101, 167, 244,
     170, 121, 114, 196, 6, 67, 152, 80, 185, 221, 7, 83, 105, 176, 208, 224,
     121},
    {226, 181, 183, 52, 241, 163, 61, 179, 221, 207, 220, 73, 245, 242, 25,
     236, 67, 84, 179, 222, 167, 62, 167, 182, 32, 9, 92, 30, 165, 127, 204,
     68},
    {226, 119, 16, 242, 200, 139, 240, 87, 11, 222, 92, 146, 156, 243, 46, 119,
     65, 59, 1, 248, 92, 183, 50, 175, 87, 40, 206, 53, 208, 220, 148, 13},
    {70, 240, 79, 112, 54, 157, 228, 146, 74, 122, 216, 88, 232, 62, 158, 13,
     14, 146, 115, 117, 176, 222, 90, 225, 244, 23, 94, 190, 150, 7, 136, 96},
    {22, 71, 241, 103, 45, 193, 195, 144, 183, 101, 154, 50, 39, 68, 49, 110,
     51, 44, 62, 0, 229, 113, 72, 81, 168, 29, 73, 106, 102, 40, 132, 24},
    {196, 133, 107, 11, 130, 105, 74, 33, 204, 171, 133, 221, 174, 193, 241,
     36, 38, 179, 196, 107, 219, 185, 181, 253, 228, 47, 155, 42, 231, 73, 41,
     78},
    {58, 255, 225, 197, 115, 208, 160, 143, 39, 197, 82, 69, 143, 235, 92, 170,
     74, 40, 57, 11, 171, 227, 26, 185, 217, 207, 90, 185, 197, 190, 35, 60},
    {88, 43, 92, 118, 223, 136, 105, 145, 238, 186, 115, 8, 214, 112, 153, 253,
     38, 108, 205, 230, 157, 130, 11, 66, 101, 85, 253, 110, 110, 14, 148, 112}
  };

  static const unsigned int hints[] = {
    0,
    0,
    0,
    0,
    0,
    0,
    4,
    4,
    4,
    4,
    0,
    4,
    4,
    0,
    0,
    4
  };

  static const size_t totals[] = {
    3,
    5,
    2,
    4,
    4,
    5,
    5,
    4,
    7,
    5,
    4,
    5,
    6,
    4,
    5,
    6
  };

  edwards_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  prime_field_t *fe = &ec->fe;
  unsigned char tmp[32];
  size_t i, j, total;
  fe_t r0, r1;
  rge_t p, q;

  (void)unused;

  printf("  - Testing ristretto_elligator (ED25519).\n");

  for (i = 0; i < ARRAY_SIZE(images); i++) {
    ASSERT(fe_import(fe, r0, bytes[i]));

    ristretto_elligator(ec, &p, r0);

    rge_export(ec, tmp, &p);

    ASSERT(torsion_memcmp(tmp, images[i], 32) == 0);

    ASSERT(ristretto_invert(ec, r1, &p, hints[i]));

    fe_set_odd(fe, r1, r1, fe_is_odd(fe, r0));

    ASSERT(fe_equal(fe, r1, r0));

    total = 0;

    for (j = 0; j < 8; j++) {
      if (ristretto_invert(ec, r1, &p, j)) {
        ristretto_elligator(ec, &q, r1);

        ASSERT(rge_equal(ec, &q, &p));

        total += 1;
      }
    }

    ASSERT(total == totals[i]);
  }

  edwards_curve_destroy(ec);
}

void
test_ecc_internal(drbg_t *rng) {
  printf("Testing internal ECC functions...\n");

  /* Memcmp */
  test_memcmp();

  /* Scalar */
  test_scalar(rng);

  /* Field Element */
  test_field_element();

  /* Utils */
  test_zero(rng);
  test_lt(rng);
  test_mpn_cmp(rng);

  /* ECC */
  test_wei_points_p256(rng);
  test_wei_points_p521(rng);
  test_wei_points_secp256k1(rng);
  test_wei_mul_g_p256(rng);
  test_wei_mul_p256(rng);
  test_wei_double_mul_p256(rng);
  test_wei_multi_mul_p256(rng);
  test_wei_mul_g_secp256k1(rng);
  test_wei_mul_secp256k1(rng);
  test_wei_double_mul_secp256k1(rng);
  test_wei_multi_mul_secp256k1(rng);
  test_mont_points_x25519();
  test_edwards_points_ed25519(rng);
  test_edwards_mul_g_ed25519(rng);
  test_edwards_mul_ed25519(rng);
  test_edwards_double_mul_ed25519(rng);
  test_edwards_multi_mul_ed25519(rng);
  test_ristretto_basepoint_multiples(rng);
  test_ristretto_elligator(rng);
}
