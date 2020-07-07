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
 * Scalar
 */

static void
test_scalar(void) {
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

    memcpy(raw, sc->raw, sc->size);

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

    memcpy(raw, sc->raw, sc->size);

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(memcmp(raw, expect1, 32) == 0);

    memcpy(raw, sc->raw, sc->size);

    raw[sc->size - 1] += 1;

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(memcmp(raw, expect2, 32) == 0);

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

    memcpy(raw, sc->raw, sc->size);

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

    memcpy(raw, sc->raw, sc->size);

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(memcmp(raw, expect1, 32) == 0);

    memcpy(raw, sc->raw, sc->size);

    raw[0] += 1;

    ASSERT(!sc_import_reduce(sc, r, raw));

    sc_export(sc, raw, r);

    ASSERT(memcmp(raw, expect2, 32) == 0);

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

    ASSERT(memcmp(max, expect, 32) == 0);

    wei_curve_destroy(ec);
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

  memcpy(raw, fe->raw, fe->size);

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
    const unsigned char mod[4] = {1, 2, 3, 4};
    const unsigned char minus1[4] = {1, 2, 3, 3};
    const unsigned char plus1[4] = {1, 2, 3, 5};
    const unsigned char full[4] = {0xff, 0xff, 0xff, 0xff};

    ASSERT(bytes_lt(minus1, mod, 4, 1));
    ASSERT(!bytes_lt(mod, mod, 4, 1));
    ASSERT(!bytes_lt(plus1, mod, 4, 1));
    ASSERT(bytes_lt(mod, full, 4, 1));
    ASSERT(!bytes_lt(full, mod, 4, 1));
  }

  {
    const unsigned char mod[4] = {4, 3, 2, 1};
    const unsigned char minus1[4] = {3, 3, 2, 1};
    const unsigned char plus1[4] = {5, 3, 2, 1};
    const unsigned char full[4] = {0xff, 0xff, 0xff, 0xff};

    ASSERT(bytes_lt(minus1, mod, 4, -1));
    ASSERT(!bytes_lt(mod, mod, 4, -1));
    ASSERT(!bytes_lt(plus1, mod, 4, -1));
    ASSERT(bytes_lt(mod, full, 4, -1));
    ASSERT(!bytes_lt(full, mod, 4, -1));
  }

  {
    size_t i;

    for (i = 0; i < 1000; i++) {
      unsigned char a[32];
      unsigned char b[32];

      drbg_generate(rng, a, sizeof(a));
      drbg_generate(rng, b, sizeof(b));

      ASSERT(bytes_lt(a, a, 32, 1) == 0);
      ASSERT(bytes_lt(b, b, 32, 1) == 0);
      ASSERT(bytes_lt(a, b, 32, 1) == (memcmp(a, b, 32) < 0));
      ASSERT(bytes_lt(b, a, 32, 1) == (memcmp(b, a, 32) < 0));

      ASSERT(bytes_lt(a, a, 32, -1) == 0);
      ASSERT(bytes_lt(b, b, 32, -1) == 0);
      ASSERT(bytes_lt(a, b, 32, -1) == (revcmp(a, b, 32) < 0));
      ASSERT(bytes_lt(b, a, 32, -1) == (revcmp(b, a, 32) < 0));
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
  wge_to_jge(ec, &jg, &ec->g);

  ASSERT(wge_import(ec, &p, g_raw, 33));

  wge_to_jge(ec, &jp, &p);
  wge_to_jge(ec, &jq, &ec->g);

  ASSERT(wge_validate(ec, &p));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(wge_equal(ec, &p, &ec->g));
  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(wge_import(ec, &q, g2_raw, 33));
  ASSERT(wge_import(ec, &r, g3_raw, 33));

  wge_to_jge(ec, &jq, &q);
  wge_to_jge(ec, &jr, &r);

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

  jge_to_wge(ec, &p, &jp);

  ASSERT(wge_equal(ec, &p, &q));

  ASSERT(wge_export(ec, p_raw, &p_size, &p, 1));
  ASSERT(p_size == 33);

  ASSERT(memcmp(p_raw, g2_raw, 33) == 0);

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
  wge_to_jge(ec, &jg, &ec->g);

  ASSERT(wge_import(ec, &p, g_raw, 67));

  wge_to_jge(ec, &jp, &p);
  wge_to_jge(ec, &jq, &ec->g);

  ASSERT(wge_validate(ec, &p));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(wge_equal(ec, &p, &ec->g));
  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(wge_import(ec, &q, g2_raw, 67));
  ASSERT(wge_import(ec, &r, g3_raw, 67));

  wge_to_jge(ec, &jq, &q);
  wge_to_jge(ec, &jr, &r);

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

  jge_to_wge(ec, &p, &jp);

  ASSERT(wge_equal(ec, &p, &q));

  ASSERT(wge_export(ec, p_raw, &p_size, &p, 1));
  ASSERT(p_size == 67);

  ASSERT(memcmp(p_raw, g2_raw, 67) == 0);

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
  wge_to_jge(ec, &jg, &ec->g);

  ASSERT(wge_import(ec, &p, g_raw, 33));

  wge_to_jge(ec, &jp, &p);
  wge_to_jge(ec, &jq, &ec->g);

  ASSERT(wge_validate(ec, &p));
  ASSERT(jge_validate(ec, &jp));
  ASSERT(jge_validate(ec, &jq));
  ASSERT(wge_equal(ec, &p, &ec->g));
  ASSERT(jge_equal(ec, &jp, &jq));

  ASSERT(wge_import(ec, &q, g2_raw, 33));
  ASSERT(wge_import(ec, &r, g3_raw, 33));

  wge_to_jge(ec, &jq, &q);
  wge_to_jge(ec, &jr, &r);

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

  jge_to_wge(ec, &p, &jp);

  ASSERT(wge_equal(ec, &p, &q));

  ASSERT(wge_export(ec, p_raw, &p_size, &p, 1));
  ASSERT(p_size == 33);

  ASSERT(memcmp(p_raw, g2_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 33) == 0);

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
  mge_to_pge(ec, &jg, &ec->g);

  ASSERT(mge_import(ec, &p, g_raw, 0));

  mge_to_pge(ec, &jp, &p);
  mge_to_pge(ec, &jq, &ec->g);

  ASSERT(mge_validate(ec, &p));
  ASSERT(pge_validate(ec, &jp));
  ASSERT(pge_validate(ec, &jq));
  ASSERT(mge_equal(ec, &p, &ec->g));
  ASSERT(pge_equal(ec, &jp, &jq));

  ASSERT(mge_import(ec, &q, g2_raw, 0));
  ASSERT(mge_import(ec, &r, g3_raw, 0));

  mge_to_pge(ec, &jq, &q);
  mge_to_pge(ec, &jr, &r);

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

  pge_to_mge(ec, &p, &jp, 0);

  ASSERT(mge_equal(ec, &p, &q));

  ASSERT(mge_export(ec, p_raw, &p));

  ASSERT(memcmp(p_raw, g2_raw, 32) == 0);

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

  ASSERT(memcmp(p_raw, g2_raw, 32) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 32) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 32) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 32) == 0);

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

  ASSERT(memcmp(q_raw, expect_raw, 32) == 0);

  edwards_scratch_destroy(ec, scratch);
  edwards_curve_destroy(ec);
}

void
test_ecc_internal(drbg_t *rng) {
  printf("Testing internal ECC functions...\n");

  /* Scalar */
  test_scalar();

  /* Field Element */
  test_field_element();

  /* Utils */
  test_zero(rng);
  test_lt(rng);

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
}
