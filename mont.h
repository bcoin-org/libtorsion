#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include "field.h"

/*
 * Structs
 */

/* mge = montgomery group element (affine) */
typedef struct mge_s {
  fe_t x;
  fe_t y;
  int inf;
} mge_t;

/* pge = projective group element (x/z) */
typedef struct pge_s {
  fe_t x;
  fe_t z;
} pge_t;

typedef struct mont_s {
  const char *prefix;
  prime_field_t fe;
  scalar_field_t sc;
  fe_t a;
  fe_t b;
  fe_t bi;
  mge_t g;
} mont_t;

typedef struct mont_def_s {
  const char *id;
  const char *prefix;
  const prime_def_t *fe;
  const scalar_def_t *sc;
  const unsigned char a[MAX_FIELD_SIZE];
  const unsigned char b[MAX_FIELD_SIZE];
  unsigned int h;
  int z;
  const unsigned char c[MAX_FIELD_SIZE];
  const unsigned char x[MAX_FIELD_SIZE];
  const unsigned char y[MAX_FIELD_SIZE];
} mont_def_t;

/*
 * Affine Point
 */

static void
mge_zero(mont_t *ec, mge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_zero(fe, r->x);
  fe_zero(fe, r->y);
  r->inf = 1;
}

static void
mge_cleanse(mont_t *ec, mge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
  r->inf = 1;
}

static int
mge_validate(mont_t *ec, const mge_t *p) {
  /* [MONT3] Page 3, Section 2. */
  /* B * y^2 = x^3 + A * x^2 + x */
  prime_field_t *fe = &ec->fe;
  fe_t lhs, rhs, x2, x3;

  fe_sqr(fe, lhs, p->y);
  fe_mul(fe, lhs, ec->b, lhs);
  fe_sqr(fe, x2, p->x);
  fe_mul(fe, x3, x2, p->x);
  fe_mul(fe, x2, ec->a, x2);
  fe_add(fe, rhs, x3, x2);
  fe_add(fe, rhs, rhs, p->x);

  return fe_equal(fe, lhs, rhs) | p->inf;
}

static int
mge_set_x(mont_t *ec, mge_t *r, const fe_t x, int sign) {
  /* [MONT3] Page 3, Section 2. */
  /* B * y^2 = x^3 + A * x^2 + x */
  prime_field_t *fe = &ec->fe;
  fe_t y, x2, x3;
  int ret;

  fe_sqr(fe, x2, x);
  fe_mul(fe, x3, x2, x);
  fe_mul(fe, x2, ec->a, x2);
  fe_add(fe, y, x3, x2);
  fe_add(fe, y, y, x);
  fe_mul(fe, y, y, ec->bi);

  ret = fe_sqrt(fe, y, y);

  if (sign != -1)
    fe_set_odd(fe, y, y, sign);

  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
  r->inf = 0;

  return ret;
}

static void
mge_set_xy(mont_t *ec, mge_t *r, const fe_t x, const fe_t y) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
  r->inf = 0;
}

static int
mge_import(mont_t *ec, mge_t *r, const unsigned char *raw) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;
  pge_t p;

  pge_import(ec, &p, raw);

  return pge_to_mge(ec, r, &p);
}

static int
mge_export(mont_t *ec,
          unsigned char *raw,
          const mge_t *p) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;

  if (p->inf)
    return 0;

  fe_export(fe, raw, p->x);

  return 1;
}

static void
mge_swap(mont_t *ec, mge_t *a, mge_t *b, unsigned int flag) {
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
mge_set(mont_t *ec, mge_t *r, const mge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  r->inf = a->inf;
}

static int
mge_equal(mont_t *ec, const mge_t *a, const mge_t *b) {
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
mge_is_zero(mont_t *ec, const mge_t *a) {
  return a->inf;
}

static void
mge_neg(mont_t *ec, mge_t *r, const mge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg(fe, r->y, a->y);
  r->inf = a->inf;
}

static void
mge_dbl(mont_t *ec, mge_t *r, const mge_t *a) {
  /* [MONT1] Page 8, Section 4.3.2.
   *
   * Addition Law (doubling):
   *
   *   l = (3 * x1^2 + 2 * a * x1 + 1) / (2 * b * y1)
   *   x3 = b * l^2 - a - 2 * x1
   *   y3 = l * (x1 - x3) - y1
   *
   * 1I + 3M + 2S + 7A + 1*a + 1*b + 1*b + 2*2 + 1*3
   */
  prime_field_t *fe = &ec->fe;
  fe_t l, t, x3, y3;

  /* P = O */
  if (a->inf) {
    mge_zero(ec, r);
    return;
  }

  /* Y1 = 0 */
  if (fe_is_zero(fe, a->y)) {
    mge_zero(ec, r);
    return;
  }

  /* L = (3 * X1^2 + 2 * a * X1 + 1) / (2 * b * Y1) */
  fe_mulw(fe, x3, ec->a, 2);
  fe_mul(fe, x3, x3, a->x);
  fe_add(fe, x3, x3, fe->one);
  fe_sqr(fe, l, a->x);
  fe_mulw(fe, l, l, 3);
  fe_add(fe, l, l, x3);
  fe_mulw(fe, t, a->y, 2);
  fe_mul(fe, t, t, ec->b);
  fe_invert_var(fe, t, t);
  fe_mul(fe, l, l, t);

  /* X3 = b * L^2 - a - 2 * X1 */
  fe_sqr(fe, x3, l);
  fe_mul(fe, x3, x3, ec->b);
  fe_sub(fe, x3, x3, ec->a);
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
mge_add(mont_t *ec, mge_t *r, const mge_t *a, const mge_t *b) {
  /* [MONT1] Page 8, Section 4.3.2.
   *
   * Addition Law:
   *
   *   l = (y2 - y1) / (x2 - x1)
   *   x3 = b * l^2 - a - x1 - x2
   *   y3 = l * (x1 - x3) - y1
   *
   * 1I + 2M + 1S + 7A + 1*b
   */
  prime_field_t *fe = &ec->fe;
  fe_t l, t, x3, y3;

  /* O + P = P */
  if (a->inf) {
    mge_set(ec, r, b);
    return;
  }

  /* P + O = P */
  if (b->inf) {
    mge_set(ec, r, a);
    return;
  }

  /* P + P, P + -P */
  if (fe_equal(fe, a->x, b->x)) {
    /* P + -P = O */
    if (!fe_equal(fe, a->y, b->y)) {
      mge_zero(ec, r);
      return;
    }

    /* P + P = 2P */
    mge_dbl(ec, r, a);
    return;
  }

  /* L = (Y2 - Y1) / (X2 - X1) */
  fe_sub(fe, l, b->y, a->y);
  fe_sub(fe, t, b->x, a->x);
  fe_invert_var(fe, t, t);
  fe_mul(fe, l, l, t);

  /* X3 = b * L^2 - a - X1 - X2 */
  fe_sqr(fe, x3, l);
  fe_mul(fe, x3, x3, ec->b);
  fe_sub(fe, x3, x3, ec->a);
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
mge_sub(mont_t *ec, mge_t *r, const mge_t *a, const mge_t *b) {
  mge_t c;
  mge_neg(ec, &c, b);
  mge_add(ec, r, a, &c);
}

static void
mge_to_pge(mont_t *ec, pge_t *r, const mge_t *a) {
  prime_field_t *fe = &ec->fe;

  if (a->inf) {
    fe_set(fe, r->x, fe->one);
    fe_zero(fe, r->z);
    return;
  }

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->z, fe->one);
}

static void
mge_print(mont_t *ec, const mge_t *p) {
  prime_field_t *fe = &ec->fe;

  if (mge_is_zero(ec, p)) {
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
 * Projective Point
 */

static void
pge_zero(mont_t *ec, pge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, fe->one);
  fe_zero(fe, r->z);
}

static void
pge_cleanse(mont_t *ec, pge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->z);
}

static int
pge_validate(mont_t *ec, const pge_t *p) {
  prime_field_t *fe = &ec->fe;
  fe_t x2, x3, z2, ax2, xz2, y2;

  // B * y^2 * z = x^3 + A * x^2 * z + x * z^2
  fe_sqr(fe, x2, p->x);
  fe_mul(fe, x3, x2, p->x);
  fe_sqr(fe, z2, p->z);
  fe_mul(fe, ax2, ec->a, x2);
  fe_mul(fe, ax2, x2, p->z);
  fe_mul(fe, xz2, p->x, z2);
  fe_add(fe, y2, x3, ax2);
  fe_add(fe, y2, y2, xz2);
  fe_mul(fe, y2, y2, ec->bi);
  fe_mul(fe, y2, y2, p->z);

  // sqrt(y^2 * z^4) = y * z^2
  return fe_is_square(fe, y2);
}

static void
pge_import(mont_t *ec, pge_t *r, const unsigned char *raw) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;

  if ((fe->bits & 7) != 0) {
    /* Ignore the hi bit for curve25519. */
    unsigned char tmp[MAX_FIELD_SIZE];
    unsigned int ignore = fe->size * 8 - fe->bits;
    unsigned int mask = (1 << (8 - ignore)) - 1;

    memcpy(tmp, raw, fe->size);

    tmp[fe->size - 1] &= mask;

    fe_import(fe, r->x, tmp);
  } else {
    fe_import(fe, r->x, raw);
  }
}

static int
pge_export(mont_t *ec,
          unsigned char *raw,
          const pge_t *p) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;
  int ret;
  fe_t x;

  ret = fe_invert(fe, x, p->z);

  fe_mul(fe, x, x, p->x);
  fe_export(fe, raw, x);

  return ret;
}

static void
pge_swap(mont_t *ec, pge_t *a, pge_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_swap(fe, a->x, b->x, flag);
  fe_swap(fe, a->z, b->z, flag);
}

static void
pge_set(mont_t *ec, pge_t *r, const pge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->z, a->z);
}

static int
pge_is_zero(mont_t *ec, const pge_t *a) {
  prime_field_t *fe = &ec->fe;

  return fe_is_zero(fe, a->z);
}

static int
pge_equal(mont_t *ec, const pge_t *a, const pge_t *b) {
  prime_field_t *fe = &ec->fe;
  fe_t e1, e2;

  /* X1 * Z2 == X2 * Z1 */
  fe_mul(fe, e1, a->x, b->z);
  fe_mul(fe, e2, b->x, a->z);

  return fe_equal(fe, e1, e2);
}

static void
pge_dbl(mont_t *ec, pge_t *p3, const pge_t *p1) {
  /* https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-3
   * 2M + 2S + 4A + 1*a24
   */
  prime_field_t *fe = &ec->fe;
  fe_t a, aa, b, bb, c, t;

  /* A = X1 + Z1 */
  fe_add(fe, a, p1->x, p1->z);

  /* AA = A^2 */
  fe_sqr(fe, aa, a);

  /* B = X1 - Z1 */
  fe_sub(fe, b, p1->x, p1->z);

  /* BB = B^2 */
  fe_sqr(fe, bb, b);

  /* C = AA - BB */
  fe_sub(fe, c, aa, bb);

  /* X3 = AA * BB */
  fe_mul(fe, p3->x, aa, bb);

  /* Z3 = C * (BB + a24 * C) */
  fe_mul(fe, t, c, ec->a24);
  fe_add(fe, t, bb, t);
  fe_mul(fe, p3->z, c, t);
}

static void
pge_dad(mont_t *ec, pge_t *p5, pge_t *p4, const pge_t *p1, const pge_t *p3, const pge_t *p2) {
  /* https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#ladder-ladd-1987-m-3
   * 6M + 4S + 8A + 1*a24
   */
  prime_field_t *fe = &ec->fe;
  fe_t a, aa, b, bb, e, c, d, da, cb, t;

  assert(p5 != p1);

  /* A = X2 + Z2 */
  fe_add(fe, a, p2->x, p2->z);

  /* AA = A^2 */
  fe_sqr(fe, aa, a);

  /* B = X2 - Z2 */
  fe_sub(fe, b, p2->x, p2->z);

  /* BB = B^2 */
  fe_sqr(fe, bb, b);

  /* E = AA - BB */
  fe_sub(fe, aa, bb);

  /* C = X3 + Z3 */
  fe_add(fe, c, p3->x, p3->z);

  /* D = X3 - Z3 */
  fe_sub(fe, d, p3->x, p3->z);

  /* DA = D * A */
  fe_mul(fe, da, d, a);

  /* CB = C * B */
  fe_mul(fe, cb, c, b);

  /* X5 = Z1 * (DA + CB)^2 */
  fe_add(fe, t, da, cb);
  fe_sqr(fe, t, t);
  fe_mul(fe, p5->x, p1->z, t);

  /* Z5 = X1 * (DA - CB)^2 */
  fe_sub(fe, t, da, cb);
  fe_sqr(fe, t, t);
  fe_mul(fe, p5->z, p1->x, t);

  /* X4 = AA * BB */
  fe_mul(fe, p4->x, aa, bb);

  /* Z4 = E * (BB + a24 * E) */
  fe_mul(fe, t, ec->a24, e);
  fe_add(fe, t, bb, t);
  fe_mul(fe, p4->z, e, t);
}

static int
pge_to_mge(mont_t *ec, mge_t *r, const pge_t *p) {
  /* https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#scaling-scale
   * 1I + 1M
   */
  prime_field_t *fe = &ec->fe;
  fe_t a;

  /* A = 1 / Z1 */
  fe_invert(fe, a, p->z);

  /* X3 = X1 * A */
  fe_mul(fe, r->x, p->x, a);

  /* Take the principle square root. */
  return mge_set_x(ec, r, r->x, -1);
}

static void
pge_print(mont_t *ec, const pge_t *p) {
  prime_field_t *fe = &ec->fe;

  if (pge_is_zero(ec, p)) {
    printf("(infinity)\n");
  } else {
    mp_limb_t xp[MAX_FIELD_LIMBS];
    mp_limb_t zp[MAX_FIELD_LIMBS];

    fe_get_limbs(fe, xp, p->x);
    fe_get_limbs(fe, zp, p->z);

    printf("(");
    mpn_print(xp, fe->limbs, 16);
    printf(", ");
    mpn_print(zp, fe->limbs, 16);
    printf(")\n");
  }
}

/*
 * Curve
 */

static void
mont_init(mont_t *ec, const mont_def_t *def) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;

  prime_field_init(fe, def->fe, -1);
  scalar_field_init(sc, def->sc, -1);

  fe_import_be(fe, ec->a, def->a);
  fe_import_be(fe, ec->b, def->b);
  fe_invert_var(fe, ec->bi, ec->b);

  fe_import_be(fe, ec->g.x, def->x);
  fe_import_be(fe, ec->g.y, def->y);
  ec->g.inf = 0;
}

static void
mont_mul(mont_t *ec, pge_t *r, const pge_t *p, const sc_t k) {
  /* Multiply with the Montgomery Ladder.
   *
   * [MONT3] Algorithm 7, Page 16, Section 5.3.
   *         Algorithm 8, Page 16, Section 5.3.
   *
   * [RFC7748] Page 7, Section 5.
   *
   * Note that any clamping is meant to
   * be done _outside_ of this function.
   */
  pge_t a, b;
  int swap = 0;
  mp_size_t i;

  /* Clone points (for safe swapping). */
  pge_set(ec, a, p);
  pge_zero(ec, b, b);

  /* Climb the ladder. */
  for (i = fe->bits - 1; i >= 0; i--) {
    mp_limb_t bit = (k[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS)) & 1;

    /* Maybe swap. */
    pge_swap(ec, a, b, swap ^ bit);

    /* Single coordinate add+double. */
    pge_dad(ec, a, b, p, a, b);

    swap = bit;
  }

  /* Finalize loop. */
  pge_swap(ec, a, b, swap);

  pge_set(ec, r, b);
}

static void
mont_mul_g(mont_t *ec, pge_t *r, const sc_t k) {
  pge_t g;
  mge_to_pge(ec, &g, &ec->g);
  mont_ladder(ec, r, &g, k);
}

/*
 * Curves
 */
