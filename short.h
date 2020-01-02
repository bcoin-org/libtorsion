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

typedef struct ge_s {
  fe_t x;
  fe_t y;
  int inf;
} ge_t;

typedef struct gej_s {
  fe_t x;
  fe_t y;
  fe_t z;
} gej_t;

typedef struct curve_s {
  int hash;
  prime_field_t fe;
  scalar_field_t sc;
  mp_limb_t pmodn[MAX_SCALAR_LIMBS * 4];
  fe_t red_n;
  fe_t a;
  fe_t b;
  ge_t g;
  sc_t blind;
  ge_t unblind;
  ge_t points[(1 << NAF_WIDTH_PRE) - 1];
} curve_t;

typedef struct curve_def_s {
  const char *id;
  int hash;
  const prime_def_t *fe;
  const scalar_def_t *sc;
  const unsigned char a[MAX_FIELD_SIZE];
  const unsigned char b[MAX_FIELD_SIZE];
  unsigned int h;
  int z;
  const unsigned char c[MAX_FIELD_SIZE];
  const unsigned char x[MAX_FIELD_SIZE];
  const unsigned char y[MAX_FIELD_SIZE];
} curve_def_t;

/*
 * Affine Point
 */

static void
ge_zero(curve_t *ec, ge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_zero(fe, r->x);
  fe_zero(fe, r->y);
  r->inf = 1;
}

static void
ge_cleanse(curve_t *ec, ge_t *r) {
  prime_field_t *fe = &ec->fe;
  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
  r->inf = 1;
}

static int
ge_validate(curve_t *ec, const ge_t *p) {
  /* [GECC] Page 89, Section 3.2.2. */
  prime_field_t *fe = &ec->fe;
  fe_t lhs, rhs, ax;

  /* y^2 = x^3 + a * x + b */
  fe_sqr(fe, lhs, p->y);
  fe_sqr(fe, rhs, p->x);
  fe_mul(fe, rhs, rhs, p->x);
  fe_mul(fe, ax, ec->a, p->x);
  fe_add(fe, rhs, rhs, ax);
  fe_add(fe, rhs, rhs, ec->b);

  return fe_equal(fe, lhs, rhs) | p->inf;
}

static int
ge_set_x(curve_t *ec, ge_t *r, const fe_t x, int sign) {
  /* [GECC] Page 89, Section 3.2.2. */
  prime_field_t *fe = &ec->fe;
  fe_t y, ax;
  int ret;

  /* y^2 = x^3 + a * x + b */
  fe_sqr(fe, y, x);
  fe_mul(fe, y, y, x);
  fe_mul(fe, ax, ec->a, x);
  fe_add(fe, y, y, ax);
  fe_add(fe, y, y, ec->b);

  ret = fe_sqrt(fe, y, y);

  if (sign != -1)
    fe_set_odd(fe, y, y, sign);

  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
  r->inf = 0;

  return ret;
}

static void
ge_set_xy(curve_t *ec, ge_t *r, const fe_t x, const fe_t y) {
  prime_field_t *fe = &ec->fe;
  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
  r->inf = 0;
}

static int
ge_import(curve_t *ec, ge_t *r, const unsigned char *raw, size_t len) {
  /* [SEC1] Page 11, Section 2.3.4. */
  prime_field_t *fe = &ec->fe;
  int form;

  if (len == 0)
    return 0;

  form = raw[0];

  switch (form) {
    case 0x02:
    case 0x03: {
      if (len != 1 + fe->size)
        return 0;

      if (!fe_import(fe, r->x, raw + 1))
        return 0;

      if (!ge_set_x(ec, r, r->x, form & 1))
        return 0;

      return 1;
    }
    case 0x04:
    case 0x06:
    case 0x07: {
      if (len != 1 + fe->size * 2)
        return 0;

      if (!fe_import(fe, r->x, raw + 1))
        return 0;

      if (!fe_import(fe, r->y, raw + 1 + fe->size))
        return 0;

      r->inf = 0;

      if (form != 0x04 && form != (0x06 | fe_is_odd(fe, r->y)))
        return 0;

      if (!ge_validate(ec, r))
        return 0;

      return 1;
    }
    default: {
      return 0;
    }
  }
}

static int
ge_export(curve_t *ec,
          unsigned char *raw,
          size_t *len,
          const ge_t *p,
          int compact) {
  /* [SEC1] Page 10, Section 2.3.3. */
  prime_field_t *fe = &ec->fe;

  if (p->inf)
    return 0;

  if (compact) {
    raw[0] = 0x02 | fe_is_odd(fe, p->y);
    fe_export(fe, raw + 1, p->x);

    if (len != NULL)
      *len = 1 + fe->size;
  } else {
    raw[0] = 0x04;
    fe_export(fe, raw + 1, p->x);
    fe_export(fe, raw + 1 + fe->size, p->y);

    if (len != NULL)
      *len = 1 + fe->size * 2;
  }

  return 1;
}

static int
ge_import_x(curve_t *ec, ge_t *r, const unsigned char *raw) {
  /* [SCHNORR] "Specification". */
  prime_field_t *fe = &ec->fe;

  if (!fe_import(fe, r->x, raw))
    return 0;

  return ge_set_x(ec, r, r->x, -1);
}

static int
ge_export_x(curve_t *ec, unsigned char *raw, const ge_t *p) {
  /* [SCHNORR] "Specification". */
  prime_field_t *fe = &ec->fe;

  if (p->inf)
    return 0;

  fe_export(fe, raw, p->x);

  return 1;
}

static void
ge_swap(curve_t *ec, ge_t *a, ge_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;
  int cond = (flag != 0);
  int inf1 = a->inf;
  int inf2 = b->inf;

  fe_swap(fe, a->x, b->x, flag);
  fe_swap(fe, a->y, b->y, flag);

  a->inf = (inf1 & (cond ^ 1)) | (inf2 & cond);
  b->inf = (inf2 & (cond ^ 1)) | (inf1 & cond);
}

static void
ge_set(curve_t *ec, ge_t *r, const ge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  r->inf = a->inf;
}

static int
ge_equal(curve_t *ec, const ge_t *a, const ge_t *b) {
  prime_field_t *fe = &ec->fe;
  int both = a->inf & b->inf;
  int ret = 1;

  /* P = O, Q = O */
  ret &= (a->inf ^ b->inf) ^ 1;

  /* X1 = X2 */
  ret &= fe_equal(fe, a->x, b->x) | both;

  /* Y1 = Y2 */
  ret &= fe_equal(fe, a->y, b->y) | both;

  return ret;
}

static int
ge_is_zero(curve_t *ec, const ge_t *a) {
  return a->inf;
}

static void
ge_neg(curve_t *ec, ge_t *r, const ge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg(fe, r->y, a->y);
  r->inf = a->inf;
}

static void
ge_dbl(curve_t *ec, ge_t *r, const ge_t *a) {
  /* [GECC] Page 80, Section 3.1.2.
   *
   * Addition Law (doubling):
   *
   *   l = (3 * x1^2 + a) / (2 * y1)
   *   x3 = l^2 - 2 * x1
   *   y3 = l * (x1 - x3) - y1
   *
   * 1I + 2M + 2S + 3A + 2*2 + 1*3
   */
  prime_field_t *fe = &ec->fe;
  fe_t l, t, x3, y3;

  /* P = O */
  if (a->inf) {
    ge_zero(ec, r);
    return;
  }

  /* Y1 = 0 */
  if (fe_is_zero(fe, a->y)) {
    ge_zero(ec, r);
    return;
  }

  /* L = (3 * X1^2 + a) / (2 * Y1) */
  fe_sqr(fe, l, a->x);
  fe_mulw(fe, l, l, 3);
  fe_add(fe, l, l, ec->a);
  fe_mulw(fe, t, a->y, 2);
  fe_invert_var(fe, t, t);
  fe_mul(fe, l, l, t);

  /* X3 = L^2 - 2 * X1 */
  fe_sqr(fe, x3, l);
  fe_sub(fe, x3, x3, a->x);
  fe_sub(fe, x3, x3, a->x);

  /* Y3 = L * (X1 - X3) - Y1 */
  fe_sub(fe, t, a->x, x3);
  fe_mul(fe, y3, l, t);
  fe_sub(fe, y3, y3, a->y);

  fe_set(fe, r->x, x3);
  fe_set(fe, r->y, y3);
  r->inf = 0;
}

static void
ge_add(curve_t *ec, ge_t *r, const ge_t *a, const ge_t *b) {
  /* [GECC] Page 80, Section 3.1.2.
   *
   * Addition Law:
   *
   *   l = (y1 - y2) / (x1 - x2)
   *   x3 = l^2 - x1 - x2
   *   y3 = l * (x1 - x3) - y1
   *
   * 1I + 2M + 1S + 6A
   */
  prime_field_t *fe = &ec->fe;
  fe_t l, t, x3, y3;

  /* O + P = P */
  if (a->inf) {
    ge_set(ec, r, b);
    return;
  }

  /* P + O = P */
  if (b->inf) {
    ge_set(ec, r, a);
    return;
  }

  /* P + P, P + -P */
  if (fe_equal(fe, a->x, b->x)) {
    /* P + -P = O */
    if (!fe_equal(fe, a->y, b->y)) {
      ge_zero(ec, r);
      return;
    }

    /* P + P = 2P */
    ge_dbl(ec, r, a);
    return;
  }

  /* X1 != X2, Y1 = Y2 */
  if (fe_equal(fe, a->y, b->y)) {
    /* X3 = -X1 - X2 */
    fe_neg(fe, x3, a->x);
    fe_sub(fe, x3, x3, b->x);

    /* Y3 = -Y1 */
    fe_neg(fe, y3, a->y);

    /* Skip the inverse. */
    return;
  }

  /* L = (Y1 - Y2) / (X1 - X2) */
  fe_sub(fe, l, a->y, b->y);
  fe_sub(fe, t, a->x, b->x);
  fe_invert_var(fe, t, t);
  fe_mul(fe, l, l, t);

  /* X3 = L^2 - X1 - X2 */
  fe_sqr(fe, x3, l);
  fe_sub(fe, x3, x3, a->x);
  fe_sub(fe, x3, x3, b->x);

  /* Y3 = L * (X1 - X3) - Y1 */
  fe_sub(fe, t, a->x, x3);
  fe_mul(fe, y3, l, t);
  fe_sub(fe, y3, y3, a->y);

  fe_set(fe, r->x, x3);
  fe_set(fe, r->y, y3);
  r->inf = 0;
}

static void
ge_sub(curve_t *ec, ge_t *r, const ge_t *a, const ge_t *b) {
  ge_t c;
  ge_neg(ec, &c, b);
  ge_add(ec, r, a, &c);
}

static void
ge_to_gej(curve_t *ec, gej_t *r, const ge_t *a) {
  prime_field_t *fe = &ec->fe;

  if (a->inf) {
    fe_set(fe, r->x, fe->one);
    fe_set(fe, r->y, fe->one);
    fe_zero(fe, r->z);
    return;
  }

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, fe->one);
}

static void
ge_naf_points(curve_t *ec, ge_t *points, const ge_t *p, size_t width) {
  size_t size = (1 << width) - 1;
  ge_t dbl;
  size_t i;

  ge_dbl(ec, &dbl, p);
  ge_set(ec, &points[0], p);

  for (i = 1; i < size; i++)
    ge_add(ec, &points[i], &points[i - 1], &dbl);
}

static void
ge_print(curve_t *ec, const ge_t *p) {
  prime_field_t *fe = &ec->fe;

  if (ge_is_zero(ec, p)) {
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
 * Jacobian Point
 */

static void
gej_zero(curve_t *ec, gej_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, fe->one);
  fe_set(fe, r->y, fe->one);
  fe_zero(fe, r->z);
}

static void
gej_cleanse(curve_t *ec, gej_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
  fe_cleanse(fe, r->z);
}

static void
gej_swap(curve_t *ec, gej_t *a, gej_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_swap(fe, a->x, b->x, flag);
  fe_swap(fe, a->y, b->y, flag);
  fe_swap(fe, a->z, b->z, flag);
}

static void
gej_select(curve_t *ec,
           gej_t *r,
           const gej_t *a,
           const gej_t *b,
           unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, b->x, flag);
  fe_select(fe, r->y, a->y, b->y, flag);
  fe_select(fe, r->z, a->z, b->z, flag);
}

static void
gej_set(curve_t *ec, gej_t *r, const gej_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
}

static int
gej_is_zero(curve_t *ec, const gej_t *a) {
  prime_field_t *fe = &ec->fe;

  return fe_is_zero(fe, a->z);
}

static int
gej_equal(curve_t *ec, const gej_t *a, const gej_t *b) {
  prime_field_t *fe = &ec->fe;
  fe_t z1, z2, e1, e2;
  int inf1 = gej_is_zero(ec, a);
  int inf2 = gej_is_zero(ec, b);
  int both = inf1 & inf2;
  int ret = 1;

  /* P = O, Q = O */
  ret &= (inf1 ^ inf2) ^ 1;

  /* X1 * Z2^2 == X2 * Z1^2 */
  fe_sqr(fe, z1, a->z);
  fe_sqr(fe, z2, b->z);
  fe_mul(fe, e1, a->x, z2);
  fe_mul(fe, e2, b->x, z1);

  ret &= fe_equal(fe, e1, e2) | both;

  /* Y1 * Z2^3 == Y2 * Z1^3 */
  fe_mul(fe, z1, z1, a->z);
  fe_mul(fe, z2, z2, b->z);
  fe_mul(fe, e1, a->y, z2);
  fe_mul(fe, e2, b->y, z1);

  ret &= fe_equal(fe, e1, e2) | both;

  return ret;
}

static int
gej_equal_r(curve_t *ec, const gej_t *p, const sc_t x) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;
  mp_limb_t cp[MAX_FIELD_LIMBS];
  mp_size_t cn = fe->limbs;
  fe_t zz, rx, rn;

  assert(fe->limbs >= sc->limbs);

  if (gej_is_zero(ec, p))
    return 0;

  fe_sqr(fe, zz, p->z);

  fe_set_sc(fe, sc, rx, x);
  fe_mul(fe, rx, rx, zz);

  if (fe_equal(fe, p->x, rx))
    return 1;

  mpn_zero(cp, cn);
  mpn_copyi(cp, x, sc->limbs);

  fe_mul(fe, rn, ec->red_n, zz);

  for (;;) {
    mpn_add_n(cp, cp, sc->n, cn);

    if (mpn_cmp(cp, fe->p, cn) >= 0)
      return 0;

    fe_add(fe, rx, rx, rn);

    if (fe_equal(fe, p->x, rx))
      break;
  }

  return 1;
}

static void
gej_neg(curve_t *ec, gej_t *r, const gej_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
}

static void
gej_zero_cond(curve_t *ec, gej_t *r, const gej_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, fe->one, flag);
  fe_select(fe, r->y, a->y, fe->one, flag);
  fe_select(fe, r->z, a->z, fe->zero, flag);
}

static void
gej_neg_cond(curve_t *ec, gej_t *r, const gej_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg_cond(fe, r->y, a->y, flag);
  fe_set(fe, r->z, a->z);
}

static void
gej_to_ge(curve_t *ec, ge_t *r, const gej_t *p);

static void
gej_dbl(curve_t *ec, gej_t *r, const gej_t *a) {
  /* https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-1998-cmo-2
   * 3M + 6S + 4A + 1*a + 2*2 + 1*3 + 1*4 + 1*8
   * (implemented as: 3M + 6S + 5A + 1*a + 1*2 + 1*3 + 1*4 + 1*8)
   */
  prime_field_t *fe = &ec->fe;
  fe_t xx, yy, zz, s, m, t, t0, t1, x3, y3, z3;

  if (gej_is_zero(ec, a)) {
    gej_zero(ec, r);
    return;
  }

  if (fe_is_zero(fe, a->y)) {
    gej_zero(ec, r);
    return;
  }

  /* XX = X1^2 */
  fe_sqr(fe, xx, a->x);

  /* YY = Y1^2 */
  fe_sqr(fe, yy, a->y);

  /* ZZ = Z1^2 */
  fe_sqr(fe, zz, a->z);

  /* S = 4 * X1 * YY */
  fe_mul(fe, s, a->x, yy);
  fe_mulw(fe, s, s, 4);

  /* M = 3 * XX + a * ZZ^2 */
  fe_mulw(fe, m, xx, 3);
  fe_sqr(fe, t, zz);
  fe_mul(fe, t, t, ec->a);
  fe_add(fe, m, m, t);

  /* T = M^2 - 2 * S */
  fe_sqr(fe, t, m);
  fe_sub(fe, t, t, s);
  fe_sub(fe, t, t, s);

  /* X3 = T */
  fe_set(fe, x3, t);

  /* Y3 = M * (S - T) - 8 * YY^2 */
  fe_sub(fe, t0, s, t);
  fe_sqr(fe, t1, yy);
  fe_mulw(fe, t1, t1, 8);
  fe_mul(fe, y3, m, t0);
  fe_sub(fe, y3, y3, t1);

  /* Z3 = 2 * Y1 * Z1 */
  fe_mul(fe, z3, a->y, a->z);
  fe_mulw(fe, z3, z3, 2);

  fe_set(fe, r->x, x3);
  fe_set(fe, r->y, y3);
  fe_set(fe, r->z, z3);
}

static void
gej_add(curve_t *ec, gej_t *r, const gej_t *a, const gej_t *b) {
  /* No assumptions.
   * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-1998-cmo-2
   * 12M + 4S + 6A + 1*2 (implemented as: 12M + 4S + 7A)
   */
  prime_field_t *fe = &ec->fe;
  fe_t z1z1, z2z2, u1, u2, s1, s2, h, r0, hh, hhh, v, x3, y3, z3;

  if (gej_is_zero(ec, a)) {
    gej_set(ec, r, b);
    return;
  }

  if (gej_is_zero(ec, b)) {
    gej_set(ec, r, a);
    return;
  }

  /* Z1Z1 = Z1^2 */
  fe_sqr(fe, z1z1, a->z);

  /* Z2Z2 = Z2^2 */
  fe_sqr(fe, z2z2, b->z);

  /* U1 = X1 * Z2Z2 */
  fe_mul(fe, u1, a->x, z2z2);

  /* U2 = X2 * Z1Z1 */
  fe_mul(fe, u2, b->x, z1z1);

  /* S1 = Y1 * Z2 * Z2Z2 */
  fe_mul(fe, s1, a->y, b->z);
  fe_mul(fe, s1, s1, z2z2);

  /* S2 = Y2 * Z1 * Z1Z1 */
  fe_mul(fe, s2, b->y, a->z);
  fe_mul(fe, s2, s2, z1z1);

  /* H = U2 - U1 */
  fe_sub(fe, h, u2, u1);

  /* r = S2 - S1 */
  fe_sub(fe, r0, s2, s1);

  /* H = 0 */
  if (fe_is_zero(fe, h)) {
    if (!fe_is_zero(fe, r0)) {
      gej_zero(ec, r);
      return;
    }

    gej_dbl(ec, r, a);
    return;
  }

  /* HH = H^2 */
  fe_sqr(fe, hh, h);

  /* HHH = H * HH */
  fe_mul(fe, hhh, h, hh);

  /* V = U1 * HH */
  fe_mul(fe, v, u1, hh);

  /* X3 = r^2 - HHH - 2 * V */
  fe_sqr(fe, x3, r0);
  fe_sub(fe, x3, x3, hhh);
  fe_sub(fe, x3, x3, v);
  fe_sub(fe, x3, x3, v);

  /* Y3 = r * (V - X3) - S1 * HHH */
  fe_sub(fe, u1, v, x3);
  fe_mul(fe, u2, s1, hhh);
  fe_mul(fe, y3, r0, u1);
  fe_sub(fe, y3, y3, u2);

  /* Z3 = Z1 * Z2 * H */
  fe_mul(fe, z3, a->z, b->z);
  fe_mul(fe, z3, z3, h);

  fe_set(fe, r->x, x3);
  fe_set(fe, r->y, y3);
  fe_set(fe, r->z, z3);
}

static void
gej_sub(curve_t *ec, gej_t *r, const gej_t *a, const gej_t *b) {
  gej_t c;
  gej_neg(ec, &c, b);
  gej_add(ec, r, a, &c);
}

static void
gej_mixed_add(curve_t *ec, gej_t *r, const gej_t *a, const ge_t *b) {
  /* Assumes Z2 = 1.
   * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-madd
   * 8M + 3S + 6A + 5*2 (implemented as: 8M + 3S + 7A + 4*2)
   */
  prime_field_t *fe = &ec->fe;
  fe_t z1z1, u2, s2, h, r0, i, j, v, x3, y3, z3;

  if (gej_is_zero(ec, a)) {
    ge_to_gej(ec, r, b);
    return;
  }

  if (ge_is_zero(ec, b)) {
    gej_set(ec, r, a);
    return;
  }

  /* Z1Z1 = Z1^2 */
  fe_sqr(fe, z1z1, a->z);

  /* U2 = X2 * Z1Z1 */
  fe_mul(fe, u2, b->x, z1z1);

  /* S2 = Y2 * Z1 * Z1Z1 */
  fe_mul(fe, s2, b->y, a->z);
  fe_mul(fe, s2, s2, z1z1);

  /* H = U2 - X1 */
  fe_sub(fe, h, u2, a->x);

  /* r = 2 * (S2 - Y1) */
  fe_sub(fe, r0, s2, a->y);
  fe_mulw(fe, r0, r0, 2);

  /* H = 0 */
  if (fe_is_zero(fe, h)) {
    if (!fe_is_zero(fe, r0)) {
      gej_zero(ec, r);
      return;
    }

    gej_dbl(ec, r, a);
    return;
  }

  /* I = (2 * H)^2 */
  fe_mulw(fe, i, h, 2);
  fe_sqr(fe, i, i);

  /* J = H * I */
  fe_mul(fe, j, h, i);

  /* V = X1 * I */
  fe_mul(fe, v, a->x, i);

  /* X3 = r^2 - J - 2 * V */
  fe_sqr(fe, x3, r0);
  fe_sub(fe, x3, x3, j);
  fe_sub(fe, x3, x3, v);
  fe_sub(fe, x3, x3, v);

  /* Y3 = r * (V - X3) - 2 * Y1 * J */
  fe_sub(fe, u2, v, x3);
  fe_mul(fe, s2, a->y, j);
  fe_mulw(fe, s2, s2, 2);
  fe_mul(fe, y3, r0, u2);
  fe_sub(fe, y3, y3, s2);

  /* Z3 = 2 * Z1 * H */
  fe_mul(fe, z3, a->z, h);
  fe_mulw(fe, z3, z3, 2);

  fe_set(fe, r->x, x3);
  fe_set(fe, r->y, y3);
  fe_set(fe, r->z, z3);
}

static void
gej_mixed_sub(curve_t *ec, gej_t *r, const gej_t *a, const ge_t *b) {
  ge_t c;
  ge_neg(ec, &c, b);
  gej_mixed_add(ec, r, a, &c);
}

static void
gej_zaddu(curve_t *ec, gej_t *r, gej_t *p, const gej_t *a, const gej_t *b) {
  /* Co-Z addition with update (ZADDU) */
  /* Algorithm 19, Page 15, Appendix C */
  prime_field_t *fe = &ec->fe;
  fe_t t1, t2, t3, t4, t5, t6;

  /* T1 = X1 */
  fe_set(fe, t1, a->x);

  /* T2 = Y1 */
  fe_set(fe, t2, a->y);

  /* T3 = Z */
  fe_set(fe, t3, a->z);

  /* T4 = X2 */
  fe_set(fe, t4, b->x);

  /* T5 = Y2 */
  fe_set(fe, t5, b->y);

  /* T6 = T1 - T4 */
  fe_sub(fe, t6, t1, t4);

  /* T3 = T3 * T6 */
  fe_mul(fe, t3, t3, t6);

  /* T6 = T6^2 */
  fe_sqr(fe, t6, t6);

  /* T1 = T1 * T6 */
  fe_mul(fe, t1, t1, t6);

  /* T6 = T6 * T4 */
  fe_mul(fe, t6, t6, t4);

  /* T5 = T2 - T5 */
  fe_sub(fe, t5, t2, t5);

  /* T4 = T5^2 */
  fe_sqr(fe, t4, t5);

  /* T4 = T4 - T1 */
  fe_sub(fe, t4, t4, t1);

  /* T4 = T4 - T6 */
  fe_sub(fe, t4, t4, t6);

  /* T6 = T1 - T6 */
  fe_sub(fe, t6, t1, t6);

  /* T2 = T2 * T6 */
  fe_mul(fe, t2, t2, t6);

  /* T6 = T1 - T4 */
  fe_sub(fe, t6, t1, t4);

  /* T5 = T5 * T6 */
  fe_mul(fe, t5, t5, t6);

  /* T5 = T5 - T2 */
  fe_sub(fe, t5, t5, t2);

  /* R = (T4, T5, T3) */
  fe_set(fe, r->x, t4);
  fe_set(fe, r->y, t5);
  fe_set(fe, r->z, t3);

  /* P = (T1, T2, T3) */
  fe_set(fe, p->x, t1);
  fe_set(fe, p->y, t2);
  fe_set(fe, p->z, t3);
}

static void
gej_zaddc(curve_t *ec, gej_t *r, gej_t *s, const gej_t *a, const gej_t *b) {
  /* Conjugate co-Z addition (ZADDC) */
  /* Algorithm 20, Page 15, Appendix C */
  prime_field_t *fe = &ec->fe;
  fe_t t1, t2, t3, t4, t5, t6, t7;

  /* T1 = X1 */
  fe_set(fe, t1, a->x);

  /* T2 = Y1 */
  fe_set(fe, t2, a->y);

  /* T3 = Z */
  fe_set(fe, t3, a->z);

  /* T4 = X2 */
  fe_set(fe, t4, b->x);

  /* T5 = Y2 */
  fe_set(fe, t5, b->y);

  /* T6 = T1 - T4 */
  fe_sub(fe, t6, t1, t4);

  /* T3 = T3 * T6 */
  fe_mul(fe, t3, t3, t6);

  /* T6 = T6^2 */
  fe_sqr(fe, t6, t6);

  /* T7 = T1 * T6 */
  fe_mul(fe, t7, t1, t6);

  /* T6 = T6 * T4 */
  fe_mul(fe, t6, t6, t4);

  /* T1 = T2 + T5 */
  fe_add(fe, t1, t2, t5);

  /* T4 = T1^2 */
  fe_sqr(fe, t4, t1);

  /* T4 = T4 - T7 */
  fe_sub(fe, t4, t4, t7);

  /* T4 = T4 - T6 */
  fe_sub(fe, t4, t4, t6);

  /* T1 = T2 - T5 */
  fe_sub(fe, t1, t2, t5);

  /* T1 = T1^2 */
  fe_sqr(fe, t1, t1);

  /* T1 = T1 - T7 */
  fe_sub(fe, t1, t1, t7);

  /* T1 = T1 - T6 */
  fe_sub(fe, t1, t1, t6);

  /* T6 = T6 - T7 */
  fe_sub(fe, t6, t6, t7);

  /* T6 = T6 * T2 */
  fe_mul(fe, t6, t6, t2);

  /* T2 = T2 - T5 */
  fe_sub(fe, t2, t2, t5);

  /* T5 = 2 * T5 */
  fe_mulw(fe, t5, t5, 2);

  /* T5 = T2 + T5 */
  fe_add(fe, t5, t2, t5);

  /* T7 = T7 - T4 */
  fe_sub(fe, t7, t7, t4);

  /* T5 = T5 * T7 */
  fe_mul(fe, t5, t5, t7);

  /* T5 = T5 + T6 */
  fe_add(fe, t5, t5, t6);

  /* T7 = T4 + T7 */
  fe_add(fe, t7, t4, t7);

  /* T7 = T7 - T1 */
  fe_sub(fe, t7, t7, t1);

  /* T2 = T2 * T7 */
  fe_mul(fe, t2, t2, t7);

  /* T2 = T2 + T6 */
  fe_add(fe, t2, t2, t6);

  /* R = (T1, T2, T3) */
  fe_set(fe, r->x, t1);
  fe_set(fe, r->y, t2);
  fe_set(fe, r->z, t3);

  /* S = (T4, T5, T3) */
  fe_set(fe, s->x, t4);
  fe_set(fe, s->y, t5);
  fe_set(fe, s->z, t3);
}

static void
gej_zdblu(curve_t *ec, gej_t *r, gej_t *p, const gej_t *a) {
  /* Co-Z doubling with update (DBLU) */
  /* Algorithm 21, Page 15, Appendix C */
  prime_field_t *fe = &ec->fe;
  fe_t t0, t1, t2, t3, t4, t5;

  /* T0 = a */
  fe_set(fe, t0, ec->a);

  /* T1 = X1 */
  fe_set(fe, t1, a->x);

  /* T2 = Y1 */
  fe_set(fe, t2, a->y);

  /* T3 = 2 * T2 */
  fe_mulw(fe, t3, t2, 2);

  /* T2 = T2^2 */
  fe_sqr(fe, t2, t2);

  /* T4 = T1 + T2 */
  fe_add(fe, t4, t1, t2);

  /* T4 = T4^2 */
  fe_sqr(fe, t4, t4);

  /* T5 = T1^2 */
  fe_sqr(fe, t5, t1);

  /* T4 = T4 - T5 */
  fe_sub(fe, t4, t4, t5);

  /* T2 = T2^2 */
  fe_sqr(fe, t2, t2);

  /* T4 = T4 - T2 */
  fe_sub(fe, t4, t4, t2);

  /* T1 = 2 * T4 */
  fe_mulw(fe, t1, t4, 2);

  /* T0 = T0 + T5 */
  fe_add(fe, t0, t0, t5);

  /* T5 = 2 * T5 */
  fe_mulw(fe, t5, t5, 2);

  /* T0 = T0 + T5 */
  fe_add(fe, t0, t0, t5);

  /* T4 = T0^2 */
  fe_sqr(fe, t4, t0);

  /* T5 = 2 * T1 */
  fe_mulw(fe, t5, t1, 2);

  /* T4 = T4 - T5 */
  fe_sub(fe, t4, t4, t5);

  /* T2 = 8 * T2 */
  fe_mulw(fe, t2, t2, 8);

  /* T5 = T1 - T4 */
  fe_sub(fe, t5, t1, t4);

  /* T5 = T5 * T0 */
  fe_mul(fe, t5, t5, t0);

  /* T5 = T5 - T2 */
  fe_sub(fe, t5, t5, t2);

  /* R = (T4, T5, T3) */
  fe_set(fe, r->x, t4);
  fe_set(fe, r->y, t5);
  fe_set(fe, r->z, t3);

  /* P = (T1, T2, T3) */
  fe_set(fe, p->x, t1);
  fe_set(fe, p->y, t2);
  fe_set(fe, p->z, t3);
}

static void
gej_dblp(curve_t *ec, gej_t *r, const gej_t *a, size_t pow) {
  size_t i;

  gej_set(ec, r, a);

  for (i = 0; i < pow; i++)
    gej_dbl(ec, r, r);
}

static void
gej_to_ge(curve_t *ec, ge_t *r, const gej_t *p) {
  /* https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#scaling-z
   * 1I + 3M + 1S
   */
  prime_field_t *fe = &ec->fe;
  fe_t a, aa;

  /* A = 1 / Z1 */
  fe_invert(fe, a, p->z);

  /* AA = A^2 */
  fe_sqr(fe, aa, a);

  /* X3 = X1 * AA */
  fe_mul(fe, r->x, p->x, aa);

  /* Y3 = Y1 * AA * A */
  fe_mul(fe, r->y, p->y, aa);
  fe_mul(fe, r->y, r->y, a);
  r->inf = fe_is_zero(fe, p->z);
}

static void
gej_to_ge_var(curve_t *ec, ge_t *r, const gej_t *p) {
  /* https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#scaling-z
   * 1I + 3M + 1S
   */
  prime_field_t *fe = &ec->fe;
  fe_t a, aa;

  /* P = O */
  if (gej_is_zero(ec, p)) {
    r->inf = 1;
    return;
  }

  /* Z = 1 */
  if (fe_equal(fe, p->z, fe->one)) {
    fe_set(fe, r->x, p->x);
    fe_set(fe, r->y, p->y);
    r->inf = 0;
    return;
  }

  /* A = 1 / Z1 */
  fe_invert_var(fe, a, p->z);

  /* AA = A^2 */
  fe_sqr(fe, aa, a);

  /* X3 = X1 * AA */
  fe_mul(fe, r->x, p->x, aa);

  /* Y3 = Y1 * AA * A */
  fe_mul(fe, r->y, p->y, aa);
  fe_mul(fe, r->y, r->y, a);
  r->inf = 0;
}

static int
gej_validate(curve_t *ec, const gej_t *p) {
  /* [GECC] Example 3.20, Page 88, Section 3. */
  prime_field_t *fe = &ec->fe;
  fe_t lhs, x3, z2, z4, z6, rhs;

  /* y^2 = x^3 + a * x * z^4 + b * z^6 */
  fe_sqr(fe, lhs, p->y);
  fe_sqr(fe, x3, p->x);
  fe_mul(fe, x3, x3, p->x);
  fe_sqr(fe, z2, p->z);
  fe_sqr(fe, z4, z2);
  fe_mul(fe, z6, z4, z2);
  fe_mul(fe, rhs, ec->b, z6);
  fe_add(fe, rhs, rhs, x3);
  fe_mul(fe, x3, ec->a, z4);
  fe_mul(fe, x3, x3, p->x);
  fe_add(fe, rhs, rhs, x3);

  return fe_equal(fe, lhs, rhs)
       | gej_is_zero(ec, p);
}

static void
gej_naf_points(curve_t *ec, gej_t *points, const ge_t *p, size_t width) {
  size_t size = (1 << width) - 1;
  gej_t dbl;
  size_t i;

  ge_to_gej(ec, &points[0], p);
  gej_dbl(ec, &dbl, &points[0]);

  for (i = 1; i < size; i++)
    gej_add(ec, &points[i], &points[i - 1], &dbl);
}

static void
gej_print(curve_t *ec, const gej_t *p) {
  prime_field_t *fe = &ec->fe;

  if (gej_is_zero(ec, p)) {
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
curve_init(curve_t *ec, const curve_def_t *def) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;

  ec->hash = def->hash;

  prime_field_init(fe, def->fe, 1);
  scalar_field_init(sc, def->sc, 1);

  sc_reduce(sc, ec->pmodn, fe->p);

  fe_set_limbs(fe, ec->red_n, sc->n, sc->limbs);
  fe_import(fe, ec->a, def->a);
  fe_import(fe, ec->b, def->b);

  fe_import(fe, ec->g.x, def->x);
  fe_import(fe, ec->g.y, def->y);
  ec->g.inf = 0;

  sc_zero(sc, ec->blind);
  ge_zero(ec, &ec->unblind);

  ge_naf_points(ec, ec->points, &ec->g, NAF_WIDTH_PRE);
}

static void
curve_jmul_g_var(curve_t *ec, gej_t *r, const sc_t k) {
  /* Window NAF method for point multiplication.
   *
   * [GECC] Algorithm 3.36, Page 100, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  ge_t *points = ec->points;
  int32_t naf[MAX_SCALAR_BITS + 1];
  size_t max;
  int32_t i;
  sc_t k0;
  gej_t acc;

  /* Blind if available. */
  sc_add(sc, k0, k, ec->blind);

  /* Calculate max size. */
  max = sc_bitlen(sc, k0) + 1;

  /* Get NAF form. */
  sc_naf(sc, naf, k0, NAF_WIDTH_PRE, max);

  /* Add `this`*(N+1) for every w-NAF index. */
  gej_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    /* Count zeroes. */
    size_t k = 0;
    int32_t z;

    for (; i >= 0 && naf[i] == 0; i--)
      k += 1;

    if (i >= 0)
      k += 1;

    gej_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z = naf[i];

    assert(z != 0);

    if (z > 0)
      gej_mixed_add(ec, &acc, &acc, &points[(z - 1) >> 1]);
    else
      gej_mixed_sub(ec, &acc, &acc, &points[(-z - 1) >> 1]);
  }

  /* Unblind. */
  gej_mixed_add(ec, &acc, &acc, &ec->unblind);

  gej_set(ec, r, &acc);

  sc_cleanse(sc, k0);
}

static void
curve_jmul_var(curve_t *ec, gej_t *r, const ge_t *p, const sc_t k) {
  /* Window NAF method for point multiplication.
   *
   * [GECC] Algorithm 3.36, Page 100, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  ge_t points[(1 << NAF_WIDTH) - 1];
  int32_t naf[MAX_SCALAR_BITS + 1];
  size_t max = sc_bitlen(sc, k) + 1;
  int32_t i;
  gej_t acc;

  /* Precompute window. */
  ge_naf_points(ec, points, p, NAF_WIDTH);

  /* Get NAF form. */
  sc_naf(sc, naf, k, NAF_WIDTH, max);

  /* Add `this`*(N+1) for every w-NAF index. */
  gej_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    /* Count zeroes. */
    size_t k = 0;
    int32_t z;

    for (; i >= 0 && naf[i] == 0; i--)
      k += 1;

    if (i >= 0)
      k += 1;

    gej_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z = naf[i];

    assert(z != 0);

    if (z > 0)
      gej_mixed_add(ec, &acc, &acc, &points[(z - 1) >> 1]);
    else
      gej_mixed_sub(ec, &acc, &acc, &points[(-z - 1) >> 1]);
  }

  gej_set(ec, r, &acc);
}

static void
curve_jmul_double_var(curve_t *ec,
                      gej_t *r,
                      const sc_t k1,
                      const ge_t *p2,
                      const sc_t k2) {
  /* Multiple point multiplication, also known
   * as "Shamir's trick" (with interleaved NAFs).
   *
   * [GECC] Algorithm 3.48, Page 109, Section 3.3.3.
   *        Algorithm 3.51, Page 112, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  ge_t *wnd1 = ec->points;
  ge_t wnd2[(1 << NAF_WIDTH) - 1];
  int32_t naf1[MAX_SCALAR_BITS + 1];
  int32_t naf2[MAX_SCALAR_BITS + 1];
  size_t max1 = sc_bitlen(sc, k1) + 1;
  size_t max2 = sc_bitlen(sc, k2) + 1;
  size_t max = max1 > max2 ? max1 : max2;
  int32_t i;
  gej_t acc;

  sc_naf(sc, naf1, k1, NAF_WIDTH_PRE, max);
  sc_naf(sc, naf2, k2, NAF_WIDTH, max);

  ge_naf_points(ec, wnd2, p2, NAF_WIDTH);

  /* Multiply and add. */
  gej_zero(ec, &acc);

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

    gej_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z1 = naf1[i];
    z2 = naf2[i];

    if (z1 > 0)
      gej_mixed_add(ec, &acc, &acc, &wnd1[(z1 - 1) >> 1]);
    else if (z1 < 0)
      gej_mixed_sub(ec, &acc, &acc, &wnd1[(-z1 - 1) >> 1]);

    if (z2 > 0)
      gej_mixed_add(ec, &acc, &acc, &wnd2[(z2 - 1) >> 1]);
    else if (z2 < 0)
      gej_mixed_sub(ec, &acc, &acc, &wnd2[(-z2 - 1) >> 1]);
  }

  gej_set(ec, r, &acc);
}

static void
curve_jmul(curve_t *ec, gej_t *r, const ge_t *p, const sc_t k) {
  /* Co-Z Montgomery Ladder.
   *
   * [COZ] Algorithm 9, Page 6, Section 4.
   */
  scalar_field_t *sc = &ec->sc;
  gej_t a, b, c;
  mp_limb_t swap = 0;
  sc_t u, v;
  uint32_t ub, vb, negated, zero, minus1;
  int i, bits;

  /* Negate scalar. */
  sc_set(sc, u, k);
  sc_neg(sc, v, k);

  /* Get bit lengths. */
  ub = sc_bitlen(sc, u);
  vb = sc_bitlen(sc, v);

  /* Negate if ceil(log2(k)) < ceil(log2(-k)). */
  negated = (ub - vb) >> 31;

  /* Possibly negate. */
  sc_swap(sc, u, v, negated);

  /* Calculate the new scalar's length. */
  bits = sc_bitlen(sc, u);

  /* Edge case (k = 0). */
  zero = sc_is_zero(sc, u);

  /* Edge case (k = -1). */
  sc_set_word(sc, v, 1);
  sc_neg(sc, v, v);
  minus1 = sc_equal(sc, u, v);

  /* Multiply with Co-Z arithmetic. */
  ge_to_gej(ec, &c, p);
  gej_zdblu(ec, &a, &b, &c);

  for (i = bits - 2; i >= 0; i--) {
    mp_limb_t bit = (u[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS)) & 1;

    gej_swap(ec, &a, &b, swap ^ bit);
    gej_zaddc(ec, &a, &b, &b, &a);
    gej_zaddu(ec, &b, &a, &a, &b);

    swap = bit;
  }

  /* Finalize loop. */
  gej_swap(ec, &a, &b, swap);

  /* Handle edge case (k = 0). */
  gej_zero_cond(ec, &b, &b, zero);

  /* Handle edge case (k = -1). */
  gej_neg(ec, &c, &c);
  gej_swap(ec, &b, &c, minus1);

  /* Adjust sign. */
  gej_neg_cond(ec, &b, &b, negated);

  /* Result. */
  gej_set(ec, r, &b);

  /* Zero scalars. */
  sc_cleanse(sc, u);
  sc_cleanse(sc, v);
}

static void
curve_mul_g_var(curve_t *ec, ge_t *r, const sc_t k) {
  gej_t j;
  curve_jmul_g_var(ec, &j, k);
  gej_to_ge_var(ec, r, &j);
}

static void
curve_mul_var(curve_t *ec, ge_t *r, const ge_t *p, const sc_t k) {
  gej_t j;
  curve_jmul_var(ec, &j, p, k);
  gej_to_ge_var(ec, r, &j);
}

static void
curve_mul_double_var(curve_t *ec,
                     ge_t *r,
                     const sc_t k1,
                     const ge_t *p2,
                     const sc_t k2) {
  gej_t j;
  curve_jmul_double_var(ec, &j, k1, p2, k2);
  gej_to_ge_var(ec, r, &j);
}

static void
curve_jmul_g(curve_t *ec, gej_t *r, const sc_t k) {
  scalar_field_t *sc = &ec->sc;
  sc_t k0;

  /* Blind if available. */
  sc_add(sc, k0, k, ec->blind);

  /* Multiply in constant time. */
  curve_jmul(ec, r, &ec->g, k0);

  /* Unblind. */
  gej_mixed_add(ec, r, r, &ec->unblind);

  /* Cleanse. */
  sc_cleanse(sc, k0);
}

static void
curve_mul(curve_t *ec, ge_t *r, const ge_t *p, const sc_t k) {
  gej_t j;
  curve_jmul(ec, &j, p, k);
  gej_to_ge(ec, r, &j);
}

static void
curve_mul_g(curve_t *ec, ge_t *r, const sc_t k) {
  gej_t j;
  curve_jmul_g(ec, &j, k);
  gej_to_ge(ec, r, &j);
}

static void
curve_randomize(curve_t *ec, const unsigned char *entropy) {
  scalar_field_t *sc = &ec->sc;
  sc_t blind;
  ge_t unblind;

  sc_import_lax(sc, blind, entropy);
  curve_mul_g(ec, &unblind, blind);
  ge_neg(ec, &unblind, &unblind);

  sc_set(sc, ec->blind, blind);
  ge_set(ec, &ec->unblind, &unblind);

  sc_cleanse(sc, blind);
  ge_cleanse(ec, &unblind);
}

/*
 * Curves
 */

static const curve_def_t curve_p192 = {
  .id = "P192",
  .hash = HASH_SHA256,
  .fe = &field_p192,
  .sc = &field_q192,
  /* -3 mod p */
  .a = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfc
  },
  .b = {
    0x64, 0x21, 0x05, 0x19, 0xe5, 0x9c, 0x80, 0xe7,
    0x0f, 0xa7, 0xe9, 0xab, 0x72, 0x24, 0x30, 0x49,
    0xfe, 0xb8, 0xde, 0xec, 0xc1, 0x46, 0xb9, 0xb1
  },
  .h = 1,
  /* Icart */
  .z = -5,
  .x = {
    0x18, 0x8d, 0xa8, 0x0e, 0xb0, 0x30, 0x90, 0xf6,
    0x7c, 0xbf, 0x20, 0xeb, 0x43, 0xa1, 0x88, 0x00,
    0xf4, 0xff, 0x0a, 0xfd, 0x82, 0xff, 0x10, 0x12
  },
  .y = {
    0x07, 0x19, 0x2b, 0x95, 0xff, 0xc8, 0xda, 0x78,
    0x63, 0x10, 0x11, 0xed, 0x6b, 0x24, 0xcd, 0xd5,
    0x73, 0xf9, 0x77, 0xa1, 0x1e, 0x79, 0x48, 0x11
  }
};

static const curve_def_t curve_p224 = {
  .id = "P224",
  .hash = HASH_SHA256,
  .fe = &field_p224,
  .sc = &field_q224,
  /* -3 mod p */
  .a = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xfe
  },
  .b = {
    0xb4, 0x05, 0x0a, 0x85, 0x0c, 0x04, 0xb3, 0xab,
    0xf5, 0x41, 0x32, 0x56, 0x50, 0x44, 0xb0, 0xb7,
    0xd7, 0xbf, 0xd8, 0xba, 0x27, 0x0b, 0x39, 0x43,
    0x23, 0x55, 0xff, 0xb4
  },
  .h = 1,
  /* SSWU */
  .z = 31,
  .x = {
    0xb7, 0x0e, 0x0c, 0xbd, 0x6b, 0xb4, 0xbf, 0x7f,
    0x32, 0x13, 0x90, 0xb9, 0x4a, 0x03, 0xc1, 0xd3,
    0x56, 0xc2, 0x11, 0x22, 0x34, 0x32, 0x80, 0xd6,
    0x11, 0x5c, 0x1d, 0x21
  },
  .y = {
    0xbd, 0x37, 0x63, 0x88, 0xb5, 0xf7, 0x23, 0xfb,
    0x4c, 0x22, 0xdf, 0xe6, 0xcd, 0x43, 0x75, 0xa0,
    0x5a, 0x07, 0x47, 0x64, 0x44, 0xd5, 0x81, 0x99,
    0x85, 0x00, 0x7e, 0x34
  }
};

static const curve_def_t curve_p256 = {
  .id = "P256",
  .hash = HASH_SHA256,
  .fe = &field_p256,
  .sc = &field_q256,
  /* -3 mod p */
  .a = {
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfc
  },
  .b = {
    0x5a, 0xc6, 0x35, 0xd8, 0xaa, 0x3a, 0x93, 0xe7,
    0xb3, 0xeb, 0xbd, 0x55, 0x76, 0x98, 0x86, 0xbc,
    0x65, 0x1d, 0x06, 0xb0, 0xcc, 0x53, 0xb0, 0xf6,
    0x3b, 0xce, 0x3c, 0x3e, 0x27, 0xd2, 0x60, 0x4b
  },
  .h = 1,
  /* SSWU */
  .z = -10,
  .x = {
    0x6b, 0x17, 0xd1, 0xf2, 0xe1, 0x2c, 0x42, 0x47,
    0xf8, 0xbc, 0xe6, 0xe5, 0x63, 0xa4, 0x40, 0xf2,
    0x77, 0x03, 0x7d, 0x81, 0x2d, 0xeb, 0x33, 0xa0,
    0xf4, 0xa1, 0x39, 0x45, 0xd8, 0x98, 0xc2, 0x96
  },
  .y = {
    0x4f, 0xe3, 0x42, 0xe2, 0xfe, 0x1a, 0x7f, 0x9b,
    0x8e, 0xe7, 0xeb, 0x4a, 0x7c, 0x0f, 0x9e, 0x16,
    0x2b, 0xce, 0x33, 0x57, 0x6b, 0x31, 0x5e, 0xce,
    0xcb, 0xb6, 0x40, 0x68, 0x37, 0xbf, 0x51, 0xf5
  }
};

static const curve_def_t curve_p384 = {
  .id = "P384",
  .hash = HASH_SHA384,
  .fe = &field_p384,
  .sc = &field_q384,
  /* -3 mod p */
  .a = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xfc
  },
  .b = {
    0xb3, 0x31, 0x2f, 0xa7, 0xe2, 0x3e, 0xe7, 0xe4,
    0x98, 0x8e, 0x05, 0x6b, 0xe3, 0xf8, 0x2d, 0x19,
    0x18, 0x1d, 0x9c, 0x6e, 0xfe, 0x81, 0x41, 0x12,
    0x03, 0x14, 0x08, 0x8f, 0x50, 0x13, 0x87, 0x5a,
    0xc6, 0x56, 0x39, 0x8d, 0x8a, 0x2e, 0xd1, 0x9d,
    0x2a, 0x85, 0xc8, 0xed, 0xd3, 0xec, 0x2a, 0xef
  },
  .h = 1,
  /* Icart */
  .z = -12,
  .x = {
    0xaa, 0x87, 0xca, 0x22, 0xbe, 0x8b, 0x05, 0x37,
    0x8e, 0xb1, 0xc7, 0x1e, 0xf3, 0x20, 0xad, 0x74,
    0x6e, 0x1d, 0x3b, 0x62, 0x8b, 0xa7, 0x9b, 0x98,
    0x59, 0xf7, 0x41, 0xe0, 0x82, 0x54, 0x2a, 0x38,
    0x55, 0x02, 0xf2, 0x5d, 0xbf, 0x55, 0x29, 0x6c,
    0x3a, 0x54, 0x5e, 0x38, 0x72, 0x76, 0x0a, 0xb7
  },
  .y = {
    0x36, 0x17, 0xde, 0x4a, 0x96, 0x26, 0x2c, 0x6f,
    0x5d, 0x9e, 0x98, 0xbf, 0x92, 0x92, 0xdc, 0x29,
    0xf8, 0xf4, 0x1d, 0xbd, 0x28, 0x9a, 0x14, 0x7c,
    0xe9, 0xda, 0x31, 0x13, 0xb5, 0xf0, 0xb8, 0xc0,
    0x0a, 0x60, 0xb1, 0xce, 0x1d, 0x7e, 0x81, 0x9d,
    0x7a, 0x43, 0x1d, 0x7c, 0x90, 0xea, 0x0e, 0x5f
  }
};

static const curve_def_t curve_p521 = {
  .id = "P521",
  .hash = HASH_SHA512,
  .fe = &field_p521,
  .sc = &field_q521,
  /* -3 mod p */
  .a = {
    0x01, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xfc
  },
  .b = {
    0x00, 0x51, 0x95, 0x3e, 0xb9, 0x61, 0x8e, 0x1c,
    0x9a, 0x1f, 0x92, 0x9a, 0x21, 0xa0, 0xb6, 0x85,
    0x40, 0xee, 0xa2, 0xda, 0x72, 0x5b, 0x99, 0xb3,
    0x15, 0xf3, 0xb8, 0xb4, 0x89, 0x91, 0x8e, 0xf1,
    0x09, 0xe1, 0x56, 0x19, 0x39, 0x51, 0xec, 0x7e,
    0x93, 0x7b, 0x16, 0x52, 0xc0, 0xbd, 0x3b, 0xb1,
    0xbf, 0x07, 0x35, 0x73, 0xdf, 0x88, 0x3d, 0x2c,
    0x34, 0xf1, 0xef, 0x45, 0x1f, 0xd4, 0x6b, 0x50,
    0x3f, 0x00
  },
  .h = 1,
  /* SSWU */
  .z = -4,
  .x = {
    0x00, 0xc6, 0x85, 0x8e, 0x06, 0xb7, 0x04, 0x04,
    0xe9, 0xcd, 0x9e, 0x3e, 0xcb, 0x66, 0x23, 0x95,
    0xb4, 0x42, 0x9c, 0x64, 0x81, 0x39, 0x05, 0x3f,
    0xb5, 0x21, 0xf8, 0x28, 0xaf, 0x60, 0x6b, 0x4d,
    0x3d, 0xba, 0xa1, 0x4b, 0x5e, 0x77, 0xef, 0xe7,
    0x59, 0x28, 0xfe, 0x1d, 0xc1, 0x27, 0xa2, 0xff,
    0xa8, 0xde, 0x33, 0x48, 0xb3, 0xc1, 0x85, 0x6a,
    0x42, 0x9b, 0xf9, 0x7e, 0x7e, 0x31, 0xc2, 0xe5,
    0xbd, 0x66
  },
  .y = {
    0x01, 0x18, 0x39, 0x29, 0x6a, 0x78, 0x9a, 0x3b,
    0xc0, 0x04, 0x5c, 0x8a, 0x5f, 0xb4, 0x2c, 0x7d,
    0x1b, 0xd9, 0x98, 0xf5, 0x44, 0x49, 0x57, 0x9b,
    0x44, 0x68, 0x17, 0xaf, 0xbd, 0x17, 0x27, 0x3e,
    0x66, 0x2c, 0x97, 0xee, 0x72, 0x99, 0x5e, 0xf4,
    0x26, 0x40, 0xc5, 0x50, 0xb9, 0x01, 0x3f, 0xad,
    0x07, 0x61, 0x35, 0x3c, 0x70, 0x86, 0xa2, 0x72,
    0xc2, 0x40, 0x88, 0xbe, 0x94, 0x76, 0x9f, 0xd1,
    0x66, 0x50
  }
};

static const curve_def_t curve_secp256k1 = {
  .id = "SECP256K1",
  .hash = HASH_SHA256,
  .fe = &field_p256k1,
  .sc = &field_q256k1,
  .a = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  .b = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x07
  },
  .h = 1,
  /* SVDW */
  .z = 1,
  /* sqrt(-3) */
  .c = {
    0x0a, 0x2d, 0x2b, 0xa9, 0x35, 0x07, 0xf1, 0xdf,
    0x23, 0x37, 0x70, 0xc2, 0xa7, 0x97, 0x96, 0x2c,
    0xc6, 0x1f, 0x6d, 0x15, 0xda, 0x14, 0xec, 0xd4,
    0x7d, 0x8d, 0x27, 0xae, 0x1c, 0xd5, 0xf8, 0x52
  },
  .x = {
    0x79, 0xbe, 0x66, 0x7e, 0xf9, 0xdc, 0xbb, 0xac,
    0x55, 0xa0, 0x62, 0x95, 0xce, 0x87, 0x0b, 0x07,
    0x02, 0x9b, 0xfc, 0xdb, 0x2d, 0xce, 0x28, 0xd9,
    0x59, 0xf2, 0x81, 0x5b, 0x16, 0xf8, 0x17, 0x98
  },
  .y = {
    0x48, 0x3a, 0xda, 0x77, 0x26, 0xa3, 0xc4, 0x65,
    0x5d, 0xa4, 0xfb, 0xfc, 0x0e, 0x11, 0x08, 0xa8,
    0xfd, 0x17, 0xb4, 0x48, 0xa6, 0x85, 0x54, 0x19,
    0x9c, 0x47, 0xd0, 0x8f, 0xfb, 0x10, 0xd4, 0xb8
  }
};
