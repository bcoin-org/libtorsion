#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include "field.h"

#define NAF_WIDTH 4
#define NAF_WIDTH_PRE 8

#define HASH_SHA256 0
#define HASH_SHA384 1
#define HASH_SHA512 2
#define HASH_SHAKE256 3

/*
 * Structs
 */

typedef void clamp_func(unsigned char *raw);

/* ege = edwards group element (affine) */
typedef struct ege_s {
  fe_t x;
  fe_t y;
} ege_t;

/* xge = extended group element (projective) */
typedef struct xge_s {
  fe_t x;
  fe_t y;
  fe_t z;
  fe_t t;
} xge_t;

typedef struct edwards_s {
  int hash;
  int context;
  const char *prefix;
  prime_field_t fe;
  scalar_field_t sc;
  fe_t a;
  fe_t d;
  fe_t k;
  ege_t g;
  sc_t blind;
  xge_t unblind;
  xge_t points[(1 << NAF_WIDTH_PRE) - 1];
  clamp_func *clamp;
} edwards_t;

typedef struct edwards_def_s {
  const char *id;
  int hash;
  int context;
  const char *prefix;
  const prime_def_t *fe;
  const scalar_def_t *sc;
  const unsigned char a[MAX_FIELD_SIZE];
  const unsigned char d[MAX_FIELD_SIZE];
  unsigned int h;
  int z;
  const unsigned char c[MAX_FIELD_SIZE];
  const unsigned char x[MAX_FIELD_SIZE];
  const unsigned char y[MAX_FIELD_SIZE];
  clamp_func *clamp;
} edwards_def_t;

static void
edwards_mul_a(edwards_t *ec, fe_t r, const fe_t x);

/*
 * Affine Point
 */

static void
ege_zero(edwards_t *ec, ege_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_zero(fe, r->x);
  fe_set(fe, r->y, fe->one);
}

static void
ege_cleanse(edwards_t *ec, ege_t *r) {
  prime_field_t *fe = &ec->fe;
  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
}

static int
ege_validate(edwards_t *ec, const ege_t *p) {
  /* [TWISTED] Definition 2.1, Page 3, Section 2. */
  /*           Page 11, Section 6. */
  /* a * x^2 + y^2 = 1 + d * x^2 * y^2 */
  prime_field_t *fe = &ec->fe;
  fe_t x2, y2, dxy, lhs, rhs;

  fe_sqr(fe, x2, p->x);
  fe_sqr(fe, y2, p->y);
  fe_mul(fe, dxy, ec->d, x2);
  fe_mul(fe, dxy, dxy, y2);
  fe_mul(fe, lhs, ec->a, x2);
  fe_add(fe, lhs, lhs, y2);
  fe_add(fe, rhs, fe->one, dxy);

  return fe_equal(fe, lhs, rhs);
}

static int
ege_set_y(edwards_t *ec, ege_t *r, const fe_t y, int sign) {
  /* [RFC8032] Section 5.1.3 & 5.2.3. */
  /* x^2 = (y^2 - 1) / (d * y^2 - a) */
  prime_field_t *fe = &ec->fe;
  fe_t x, y2, lhs, rhs;
  int ret;

  fe_sqr(fe, y2, y);
  fe_sub(fe, lhs, y2, fe->one);
  fe_mul(fe, rhs, ec->d, y2);
  fe_sub(fe, rhs, rhs, ec->a);

  ret = fe_isqrt(fe, x, lhs, rhs);

  if (sign != -1)
    fe_set_odd(fe, x, x, sign);

  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);

  return ret;
}

static void
ege_set_xy(edwards_t *ec, ege_t *r, const fe_t x, const fe_t y) {
  prime_field_t *fe = &ec->fe;
  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
}

static int
ege_import(edwards_t *ec, ege_t *r, const unsigned char *raw) {
  prime_field_t *fe = &ec->fe;
  int sign;

  /* Quirk: we need an extra byte (p448). */
  if ((fe->bits & 7) == 0) {
    if ((raw[fe->size] & 0x7f) != 0x00)
      return 0;

    if (!fe_import(fe, r->y, raw))
      return 0;

    sign = (raw[fe->size] & 0x80) != 0;
  } else {
    unsigned char tmp[MAX_FIELD_SIZE];

    memcpy(tmp, raw, fe->size);

    tmp[fe->size - 1] &= 0x7f;

    if (!fe_import(fe, r->y, tmp))
      return 0;

    sign = (raw[fe->size - 1] & 0x80) != 0;
  }

  return ege_set_y(ec, r, r->y, sign);
}

static int
ege_export(edwards_t *ec,
          unsigned char *raw,
          const ege_t *p) {
  /* [RFC8032] Section 5.1.2. */
  prime_field_t *fe = &ec->fe;

  fe_export(fe, raw, p->y);

  /* Quirk: we need an extra byte (p448). */
  if ((fe->bits & 7) == 0)
    raw[fe->size] = fe_is_odd(fe, p->x) << 7;
  else
    raw[fe->size - 1] |= fe_is_odd(fe, p->x) << 7;

  return 1;
}

static void
ege_swap(edwards_t *ec, ege_t *a, ege_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_swap(fe, a->x, b->x, flag != 0);
  fe_swap(fe, a->y, b->y, flag != 0);
}

static void
ege_set(edwards_t *ec, ege_t *r, const ege_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
}

static int
ege_equal(edwards_t *ec, const ege_t *a, const ege_t *b) {
  prime_field_t *fe = &ec->fe;

  /* X1 = X2, Y1 = Y2 */
  return fe_equal(fe, a->x, b->x)
       & fe_equal(fe, a->y, b->y);
}

static int
ege_is_zero(edwards_t *ec, const ege_t *a) {
  prime_field_t *fe = &ec->fe;

  return fe_is_zero(fe, a->x)
       & fe_equal(fe, a->y, fe->one);
}

static void
ege_neg(edwards_t *ec, ege_t *r, const ege_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_neg(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
}

static void
ege_add(edwards_t *ec, ege_t *r, const ege_t *a, const ege_t *b) {
  /* Affine twisted addition formula:
   *
   *   x3 = (x1 * y2 + y1 * x2) / (1 + d * x1 * x2 * y1 * y2)
   *   y3 = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)
   */
  prime_field_t *fe = &ec->fe;
  fe_t x3, y3, x1x2, y1y2, x1y2, y1x2, z, z1, z2;

  fe_mul(fe, x1x2, a->x, b->x);
  fe_mul(fe, y1y2, a->y, b->y);
  fe_mul(fe, x1y2, a->x, b->y);
  fe_mul(fe, y1x2, a->y, b->x);

  fe_add(fe, x3, x1y2, y1x2);
  edwards_mul_a(ec, y3, x1x2);
  fe_sub(fe, y3, y1y2, y3);

  fe_mul(fe, z, x1x2, y1y2);
  fe_mul(fe, z, z, ec->d);
  fe_add(fe, z1, fe->one, z);
  fe_sub(fe, z2, fe->one, z);
  fe_mul(fe, z, z1, z2);

  fe_invert_var(fe, z, z);

  fe_mul(fe, x3, x3, z);
  fe_mul(fe, y3, y3, z);

  fe_mul(fe, r->x, x3, z2);
  fe_mul(fe, r->y, y3, z1);
}

static void
ege_dbl(edwards_t *ec, ege_t *r, const ege_t *a) {
  ege_add(ec, r, a, a);
}

static void
ege_sub(edwards_t *ec, ege_t *r, const ege_t *a, const ege_t *b) {
  ege_t c;
  ege_neg(ec, &c, b);
  ege_add(ec, r, a, &c);
}

static void
ege_to_xge(edwards_t *ec, xge_t *r, const ege_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, fe->one);
  fe_mul(fe, r->t, r->x, r->y);
}

static void
ege_print(edwards_t *ec, const ege_t *p) {
  prime_field_t *fe = &ec->fe;

  if (ege_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    mp_limb_t xp[MAX_FIELD_LIMBS];
    mp_limb_t yp[MAX_FIELD_LIMBS];

    fe_get_limbs(fe, xp, p->x);
    fe_get_limbs(fe, yp, p->y);

    printf("(");
    mpn_print(xp, fe->limbs, 16);
    printf(", ");
    mpn_print(yp, fe->limbs, 16);
    printf(")\n");
  }
}

/*
 * Extended Point
 */

static void
xge_zero(edwards_t *ec, xge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_zero(fe, r->x);
  fe_set(fe, r->y, fe->one);
  fe_set(fe, r->z, fe->one);
  fe_zero(fe, r->t);
}

static void
xge_cleanse(edwards_t *ec, xge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
  fe_cleanse(fe, r->z);
  fe_cleanse(fe, r->t);
}

static void
xge_swap(edwards_t *ec, xge_t *a, xge_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_swap(fe, a->x, b->x, flag);
  fe_swap(fe, a->y, b->y, flag);
  fe_swap(fe, a->z, b->z, flag);
  fe_swap(fe, a->t, b->t, flag);
}

static void
xge_select(edwards_t *ec,
           xge_t *r,
           const xge_t *a,
           const xge_t *b,
           unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, b->x, flag);
  fe_select(fe, r->y, a->y, b->y, flag);
  fe_select(fe, r->z, a->z, b->z, flag);
  fe_select(fe, r->t, a->t, b->t, flag);
}

static void
xge_set(edwards_t *ec, xge_t *r, const xge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
  fe_set(fe, r->t, a->t);
}

static int
xge_is_zero(edwards_t *ec, const xge_t *a) {
  prime_field_t *fe = &ec->fe;

  return fe_is_zero(fe, a->x)
       & fe_equal(fe, a->y, a->z);
}

static int
xge_equal(edwards_t *ec, const xge_t *a, const xge_t *b) {
  prime_field_t *fe = &ec->fe;
  fe_t e1, e2;
  int ret = 1;

  /* X1 * Z2 == X2 * Z1 */
  fe_mul(fe, e1, a->x, b->z);
  fe_mul(fe, e2, b->x, a->z);

  ret &= fe_equal(fe, e1, e2);

  /* Y1 * Z2 == Y2 * Z1 */
  fe_mul(fe, e1, a->y, b->z);
  fe_mul(fe, e2, b->y, a->z);

  ret &= fe_equal(fe, e1, e2);

  return ret;
}

static void
xge_neg(edwards_t *ec, xge_t *r, const xge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_neg(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
  fe_neg(fe, r->t, a->t);
}

static void
xge_zero_cond(edwards_t *ec, xge_t *r, const xge_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, fe->zero, flag);
  fe_select(fe, r->y, a->y, fe->one, flag);
  fe_select(fe, r->z, a->z, fe->one, flag);
  fe_select(fe, r->t, a->t, fe->zero, flag);
}

static void
xge_neg_cond(edwards_t *ec, xge_t *r, const xge_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_neg_cond(fe, r->x, a->x, flag);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
  fe_neg_cond(fe, r->t, a->t, flag);
}

static void
xge_dbl(edwards_t *ec, xge_t *r, const xge_t *a) {
  /* https://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html#doubling-dbl-2008-hwcd
   * 4M + 4S + 6A + 1*a + 1*2
   */
  prime_field_t *fe = &ec->fe;
  fe_t A, b, c, d, e, g, f, h;

  /* A = X1^2 */
  fe_sqr(fe, A, a->x);

  /* B = Y1^2 */
  fe_sqr(fe, b, a->y);

  /* C = 2 * Z1^2 */
  fe_sqr(fe, c, a->z);
  fe_mulw(fe, c, c, 2);

  /* D = a * A */
  edwards_mul_a(ec, d, A);

  /* E = (X1 + Y1)^2 - A - B */
  fe_add(fe, e, a->x, a->y);
  fe_sqr(fe, e, e);
  fe_sub(fe, e, e, A);
  fe_sub(fe, e, e, b);

  /* G = D + B */
  fe_add(fe, g, d, b);

  /* F = G - C */
  fe_sub(fe, f, g, c);

  /* H = D - B */
  fe_sub(fe, h, d, b);

  /* X3 = E * F */
  fe_mul(fe, r->x, e, f);

  /* Y3 = G * H */
  fe_mul(fe, r->y, g, h);

  /* T3 = E * H */
  fe_mul(fe, r->t, e, h);

  /* Z3 = F * G */
  fe_mul(fe, r->z, f, g);
}

static void
xge_add_a(edwards_t *ec, xge_t *r, const xge_t *a, const xge_t *b) {
  /* https://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html#addition-add-2008-hwcd
   * 9M + 7A + 1*a + 1*d
   */
  prime_field_t *fe = &ec->fe;
  fe_t A, B, c, d, e, f, g, h;

  /* A = X1 * X2 */
  fe_mul(fe, A, a->x, b->x);

  /* B = Y1 * Y2 */
  fe_mul(fe, B, a->y, b->y);

  /* C = T1 * d * T2 */
  fe_mul(fe, c, a->t, b->t);
  fe_mul(fe, c, c, ec->d);

  /* D = Z1 * Z2 */
  fe_mul(fe, d, a->z, b->z);

  /* E = (X1 + Y1) * (X2 + Y2) - A - B */
  fe_add(fe, f, a->x, a->y);
  fe_add(fe, g, b->x, b->y);
  fe_mul(fe, e, f, g);
  fe_sub(fe, e, e, A);
  fe_sub(fe, e, e, B);

  /* F = D - C */
  fe_sub(fe, f, d, c);

  /* G = D + C */
  fe_add(fe, g, d, c);

  /* H = B - a * A */
  edwards_mul_a(ec, h, A);
  fe_sub(fe, h, B, h);

  /* X3 = E * F */
  fe_mul(fe, r->x, e, f);

  /* Y3 = G * H */
  fe_mul(fe, r->y, g, h);

  /* T3 = E * H */
  fe_mul(fe, r->t, e, h);

  /* Z3 = F * G */
  fe_mul(fe, r->z, f, g);
}

static void
xge_add_m1(edwards_t *ec, xge_t *r, const xge_t *a, const xge_t *b) {
  /* Assumes a = -1.
   *
   * https://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#addition-add-2008-hwcd-3
   * 8M + 8A + 1*k + 1*2
   */
  prime_field_t *fe = &ec->fe;
  fe_t A, B, c, d, e, f, g, h;

  /* A = (Y1 - X1) * (Y2 - X2) */
  fe_sub(fe, c, a->y, a->x);
  fe_sub(fe, d, b->y, b->x);
  fe_mul(fe, A, c, d);

  /* B = (Y1 + X1) * (Y2 + X2) */
  fe_add(fe, c, a->y, a->x);
  fe_add(fe, d, b->y, b->x);
  fe_mul(fe, B, c, d);

  /* C = T1 * k * T2 */
  fe_mul(fe, c, a->t, b->t);
  fe_mul(fe, c, c, ec->k);

  /* D = Z1 * 2 * Z2 */
  fe_mul(fe, d, a->z, b->z);
  fe_mulw(fe, d, d, 2);

  /* E = B - A */
  fe_sub(fe, e, B, A);

  /* F = D - C */
  fe_sub(fe, f, d, c);

  /* G = D + C */
  fe_add(fe, g, d, c);

  /* H = B + A */
  fe_add(fe, h, B, A);

  /* X3 = E * F */
  fe_mul(fe, r->x, e, f);

  /* Y3 = G * H */
  fe_mul(fe, r->y, g, h);

  /* T3 = E * H */
  fe_mul(fe, r->t, e, h);

  /* Z3 = F * G */
  fe_mul(fe, r->z, f, g);
}

static void
xge_add(edwards_t *ec, xge_t *r, const xge_t *a, const xge_t *b) {
  if (ec->fe.bits == 255)
    xge_add_m1(ec, r, a, b);
  else
    xge_add_a(ec, r, a, b);
}

static void
xge_sub(edwards_t *ec, xge_t *r, const xge_t *a, const xge_t *b) {
  xge_t c;
  xge_neg(ec, &c, b);
  xge_add(ec, r, a, &c);
}

static void
xge_dblp(edwards_t *ec, xge_t *r, const xge_t *a, size_t pow) {
  size_t i;

  xge_set(ec, r, a);

  for (i = 0; i < pow; i++)
    xge_dbl(ec, r, r);
}

static void
xge_to_ege(edwards_t *ec, ege_t *r, const xge_t *p) {
  /* https://hyperelliptic.org/EFD/g1p/auto-edwards-projective.html#scaling-z
   * 1I + 2M (+ 1M if extended)
   */
  prime_field_t *fe = &ec->fe;
  fe_t a;

  /* A = 1 / Z1 */
  fe_invert(fe, a, p->z);

  /* X3 = X1 * A */
  fe_mul(fe, r->x, p->x, a);

  /* Y3 = Y1 * A */
  fe_mul(fe, r->y, p->y, a);
}

static void
xge_to_ege_var(edwards_t *ec, ege_t *r, const xge_t *p) {
  /* https://hyperelliptic.org/EFD/g1p/auto-edwards-projective.html#scaling-z
   * 1I + 2M (+ 1M if extended)
   */
  prime_field_t *fe = &ec->fe;
  fe_t a;

  /* Z = 1 */
  if (fe_equal(fe, p->z, fe->one)) {
    fe_set(fe, r->x, p->x);
    fe_set(fe, r->y, p->y);
    return;
  }

  /* A = 1 / Z1 */
  fe_invert_var(fe, a, p->z);

  /* X3 = X1 * A */
  fe_mul(fe, r->x, p->x, a);

  /* Y3 = Y1 * A */
  fe_mul(fe, r->y, p->y, a);
}

static int
xge_validate(edwards_t *ec, const xge_t *p) {
  /* [TWISTED] Definition 2.1, Page 3, Section 2. */
  /*           Page 11, Section 6. */
  /* (a * x^2 + y^2) * z^2 = z^4 + d * x^2 * y^2 */
  prime_field_t *fe = &ec->fe;
  fe_t lhs, rhs, x2, y2, ax2, z2, z4;

  fe_sqr(fe, x2, p->x);
  fe_sqr(fe, y2, p->y);
  fe_sqr(fe, z2, p->z);
  fe_sqr(fe, z4, z2);

  fe_mul(fe, ax2, ec->a, x2);
  fe_add(fe, lhs, ax2, y2);
  fe_mul(fe, lhs, lhs, z2);

  fe_mul(fe, rhs, x2, y2);
  fe_mul(fe, rhs, rhs, ec->d);
  fe_add(fe, rhs, rhs, z4);

  return fe_equal(fe, lhs, rhs);
}

static void
xge_naf_points(edwards_t *ec, xge_t *points, const ege_t *p, size_t width) {
  size_t size = (1 << width) - 1;
  xge_t dbl;
  size_t i;

  ege_to_xge(ec, &points[0], p);
  xge_dbl(ec, &dbl, &points[0]);

  for (i = 1; i < size; i++)
    xge_add(ec, &points[i], &points[i - 1], &dbl);
}

static void
xge_print(edwards_t *ec, const xge_t *p) {
  prime_field_t *fe = &ec->fe;

  if (xge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    mp_limb_t xp[MAX_FIELD_LIMBS];
    mp_limb_t yp[MAX_FIELD_LIMBS];
    mp_limb_t zp[MAX_FIELD_LIMBS];

    fe_get_limbs(fe, xp, p->x);
    fe_get_limbs(fe, yp, p->y);
    fe_get_limbs(fe, zp, p->z);

    printf("(");
    mpn_print(xp, fe->limbs, 16);
    printf(", ");
    mpn_print(yp, fe->limbs, 16);
    printf(", ");
    mpn_print(zp, fe->limbs, 16);
    printf(")\n");
  }
}

/*
 * Curve
 */

static void
edwards_init(edwards_t *ec, const edwards_def_t *def) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;
  mpz_t p, n;

  ec->hash = def->hash;
  ec->context = def->context;
  ec->prefix = def->prefix;
  ec->clamp = def->clamp;

  prime_field_init(fe, def->fe, -1);
  scalar_field_init(sc, def->sc, -1);

  mpz_roinit_n(p, fe->p, fe->limbs);
  mpz_roinit_n(n, sc->n, sc->limbs);

  fe_import_be(fe, ec->a, def->a);
  fe_import_be(fe, ec->d, def->d);
  fe_mulw(fe, ec->k, ec->d, 2);

  fe_import_be(fe, ec->g.x, def->x);
  fe_import_be(fe, ec->g.y, def->y);

  sc_zero(sc, ec->blind);
  xge_zero(ec, &ec->unblind);

  xge_naf_points(ec, ec->points, &ec->g, NAF_WIDTH_PRE);
}

static void
edwards_clamp(edwards_t *ec, unsigned char *raw) {
  ec->clamp(raw);
}

static void
edwards_mul_a(edwards_t *ec, fe_t r, const fe_t x) {
  if (ec->fe.bits == 255)
    fe_neg(&ec->fe, r, x); /* a = -1 */
  else if (ec->fe.bits == 448)
    fe_set(&ec->fe, r, x); /* a = 1 */
  else
    fe_mul(&ec->fe, r, x, ec->a);
}

static void
edwards_jmul_g_var(edwards_t *ec, xge_t *r, const sc_t k) {
  /* Window NAF method for point multiplication.
   *
   * [GECC] Algorithm 3.36, Page 100, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  xge_t *points = ec->points;
  int32_t naf[MAX_SCALAR_BITS + 1];
  size_t max;
  int32_t i;
  sc_t k0;
  xge_t acc;

  /* Blind if available. */
  sc_add(sc, k0, k, ec->blind);

  /* Calculate max size. */
  max = sc_bitlen(sc, k0) + 1;

  /* Get NAF form. */
  sc_naf(sc, naf, k0, NAF_WIDTH_PRE, max);

  /* Add `this`*(N+1) for every w-NAF index. */
  xge_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    /* Count zeroes. */
    size_t k = 0;
    int32_t z;

    for (; i >= 0 && naf[i] == 0; i--)
      k += 1;

    if (i >= 0)
      k += 1;

    xge_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z = naf[i];

    assert(z != 0);

    if (z > 0)
      xge_add(ec, &acc, &acc, &points[(z - 1) >> 1]);
    else
      xge_sub(ec, &acc, &acc, &points[(-z - 1) >> 1]);
  }

  /* Unblind. */
  xge_add(ec, &acc, &acc, &ec->unblind);

  xge_set(ec, r, &acc);

  sc_cleanse(sc, k0);
}

static void
edwards_jmul_var(edwards_t *ec, xge_t *r, const ege_t *p, const sc_t k) {
  /* Window NAF method for point multiplication.
   *
   * [GECC] Algorithm 3.36, Page 100, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  xge_t points[(1 << NAF_WIDTH) - 1];
  int32_t naf[MAX_SCALAR_BITS + 1];
  size_t max = sc_bitlen(sc, k) + 1;
  int32_t i;
  xge_t acc;

  /* Precompute window. */
  xge_naf_points(ec, points, p, NAF_WIDTH);

  /* Get NAF form. */
  sc_naf(sc, naf, k, NAF_WIDTH, max);

  /* Add `this`*(N+1) for every w-NAF index. */
  xge_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    /* Count zeroes. */
    size_t k = 0;
    int32_t z;

    for (; i >= 0 && naf[i] == 0; i--)
      k += 1;

    if (i >= 0)
      k += 1;

    xge_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z = naf[i];

    assert(z != 0);

    if (z > 0)
      xge_add(ec, &acc, &acc, &points[(z - 1) >> 1]);
    else
      xge_sub(ec, &acc, &acc, &points[(-z - 1) >> 1]);
  }

  xge_set(ec, r, &acc);
}

static void
edwards_jmul_double_var(edwards_t *ec,
                        xge_t *r,
                        const sc_t k1,
                        const ege_t *p2,
                        const sc_t k2) {
  /* Multiple point multiplication, also known
   * as "Shamir's trick" (with interleaved NAFs).
   *
   * [GECC] Algorithm 3.48, Page 109, Section 3.3.3.
   *        Algorithm 3.51, Page 112, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  xge_t *wnd1 = ec->points;
  xge_t wnd2[(1 << NAF_WIDTH) - 1];
  int32_t naf1[MAX_SCALAR_BITS + 1];
  int32_t naf2[MAX_SCALAR_BITS + 1];
  size_t max1 = sc_bitlen(sc, k1) + 1;
  size_t max2 = sc_bitlen(sc, k2) + 1;
  size_t max = max1 > max2 ? max1 : max2;
  int32_t i;
  xge_t acc;

  sc_naf(sc, naf1, k1, NAF_WIDTH_PRE, max);
  sc_naf(sc, naf2, k2, NAF_WIDTH, max);

  xge_naf_points(ec, wnd2, p2, NAF_WIDTH);

  /* Multiply and add. */
  xge_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    size_t k = 0;
    int32_t z1, z2;

    while (i >= 0) {
      if (naf1[i] != 0 || naf2[i] != 0)
        break;

      k += 1;
      i -= 1;
    }

    if (i >= 0)
      k += 1;

    xge_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z1 = naf1[i];
    z2 = naf2[i];

    if (z1 > 0)
      xge_add(ec, &acc, &acc, &wnd1[(z1 - 1) >> 1]);
    else if (z1 < 0)
      xge_sub(ec, &acc, &acc, &wnd1[(-z1 - 1) >> 1]);

    if (z2 > 0)
      xge_add(ec, &acc, &acc, &wnd2[(z2 - 1) >> 1]);
    else if (z2 < 0)
      xge_sub(ec, &acc, &acc, &wnd2[(-z2 - 1) >> 1]);
  }

  xge_set(ec, r, &acc);
}

static void
edwards_jmul(edwards_t *ec, xge_t *r, const ege_t *p, const sc_t k) {
  edwards_jmul_var(ec, r, p, k);
}

static void
edwards_mul_g_var(edwards_t *ec, ege_t *r, const sc_t k) {
  xge_t j;
  edwards_jmul_g_var(ec, &j, k);
  xge_to_ege_var(ec, r, &j);
}

static void
edwards_mul_var(edwards_t *ec, ege_t *r, const ege_t *p, const sc_t k) {
  xge_t j;
  edwards_jmul_var(ec, &j, p, k);
  xge_to_ege_var(ec, r, &j);
}

static void
edwards_mul_double_var(edwards_t *ec,
                     ege_t *r,
                     const sc_t k1,
                     const ege_t *p2,
                     const sc_t k2) {
  xge_t j;
  edwards_jmul_double_var(ec, &j, k1, p2, k2);
  xge_to_ege_var(ec, r, &j);
}

static void
edwards_jmul_g(edwards_t *ec, xge_t *r, const sc_t k) {
  scalar_field_t *sc = &ec->sc;
  sc_t k0;

  /* Blind if available. */
  sc_add(sc, k0, k, ec->blind);

  /* Multiply in constant time. */
  edwards_jmul(ec, r, &ec->g, k0);

  /* Unblind. */
  xge_add(ec, r, r, &ec->unblind);

  /* Cleanse. */
  sc_cleanse(sc, k0);
}

static void
edwards_mul(edwards_t *ec, ege_t *r, const ege_t *p, const sc_t k) {
  xge_t j;
  edwards_jmul(ec, &j, p, k);
  xge_to_ege(ec, r, &j);
}

static void
edwards_mul_g(edwards_t *ec, ege_t *r, const sc_t k) {
  xge_t j;
  edwards_jmul_g(ec, &j, k);
  xge_to_ege(ec, r, &j);
}

static void
edwards_randomize(edwards_t *ec, const unsigned char *entropy) {
  scalar_field_t *sc = &ec->sc;
  sc_t blind;
  xge_t unblind;

  sc_import_lax(sc, blind, entropy);
  edwards_jmul_g(ec, &unblind, blind);
  xge_neg(ec, &unblind, &unblind);

  sc_set(sc, ec->blind, blind);
  xge_set(ec, &ec->unblind, &unblind);

  sc_cleanse(sc, blind);
  xge_cleanse(ec, &unblind);
}

/*
 * Curves
 */

static const edwards_def_t curve_ed25519 = {
  .id = "ED25519",
  .hash = HASH_SHA512,
  .context = 0,
  .prefix = "SigEd25519 no Ed25519 collisions",
  .fe = &field_p25519,
  .sc = &field_q25519,
  /* -1 mod p */
  .a = {
    0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xec
  },
  /* -121665 / 121666 mod p */
  .d = {
    0x52, 0x03, 0x6c, 0xee, 0x2b, 0x6f, 0xfe, 0x73,
    0x8c, 0xc7, 0x40, 0x79, 0x77, 0x79, 0xe8, 0x98,
    0x00, 0x70, 0x0a, 0x4d, 0x41, 0x41, 0xd8, 0xab,
    0x75, 0xeb, 0x4d, 0xca, 0x13, 0x59, 0x78, 0xa3
  },
  .h = 8,
  /* Elligator 2 */
  .z = 2,
  .x = {
    0x21, 0x69, 0x36, 0xd3, 0xcd, 0x6e, 0x53, 0xfe,
    0xc0, 0xa4, 0xe2, 0x31, 0xfd, 0xd6, 0xdc, 0x5c,
    0x69, 0x2c, 0xc7, 0x60, 0x95, 0x25, 0xa7, 0xb2,
    0xc9, 0x56, 0x2d, 0x60, 0x8f, 0x25, 0xd5, 0x1a
  },
  /* 4/5 */
  .y = {
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
    0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x58
  },
  .clamp = p25519_clamp
};

static const edwards_def_t curve_ed448 = {
  .id = "ED448",
  .hash = HASH_SHAKE256,
  .context = 1,
  .prefix = "SigEd448",
  .fe = &field_p448,
  .sc = &field_q448,
  .a = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01
  },
  /* -39081 mod p */
  .d = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x67, 0x56
  },
  .h = 4,
  /* Elligator 2 */
  .z = -1,
  .x = {
    0x4f, 0x19, 0x70, 0xc6, 0x6b, 0xed, 0x0d, 0xed,
    0x22, 0x1d, 0x15, 0xa6, 0x22, 0xbf, 0x36, 0xda,
    0x9e, 0x14, 0x65, 0x70, 0x47, 0x0f, 0x17, 0x67,
    0xea, 0x6d, 0xe3, 0x24, 0xa3, 0xd3, 0xa4, 0x64,
    0x12, 0xae, 0x1a, 0xf7, 0x2a, 0xb6, 0x65, 0x11,
    0x43, 0x3b, 0x80, 0xe1, 0x8b, 0x00, 0x93, 0x8e,
    0x26, 0x26, 0xa8, 0x2b, 0xc7, 0x0c, 0xc0, 0x5e
  },
  .y = {
    0x69, 0x3f, 0x46, 0x71, 0x6e, 0xb6, 0xbc, 0x24,
    0x88, 0x76, 0x20, 0x37, 0x56, 0xc9, 0xc7, 0x62,
    0x4b, 0xea, 0x73, 0x73, 0x6c, 0xa3, 0x98, 0x40,
    0x87, 0x78, 0x9c, 0x1e, 0x05, 0xa0, 0xc2, 0xd7,
    0x3a, 0xd3, 0xff, 0x1c, 0xe6, 0x7c, 0x39, 0xc4,
    0xfd, 0xbd, 0x13, 0x2c, 0x4e, 0xd7, 0xc8, 0xad,
    0x98, 0x08, 0x79, 0x5b, 0xf2, 0x30, 0xfa, 0x14
  },
  .clamp = p448_clamp
};
