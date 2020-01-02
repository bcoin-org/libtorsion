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

typedef struct ge_s {
  fe_t x;
  fe_t y;
  int inf;
} ge_t;

typedef struct gej_s {
  fe_t x;
  fe_t z;
} gej_t;

typedef struct curve_s {
  const char *prefix;
  prime_field_t fe;
  scalar_field_t sc;
  fe_t a;
  fe_t b;
  fe_t bi;
  ge_t g;
} curve_t;

typedef struct curve_def_s {
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
ge_set_x(curve_t *ec, ge_t *r, const fe_t x, int sign) {
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
ge_set_xy(curve_t *ec, ge_t *r, const fe_t x, const fe_t y) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
  r->inf = 0;
}

static int
ge_import(curve_t *ec, ge_t *r, const unsigned char *raw) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;

  if ((fe->bits & 7) != 0) {
    /* Ignore the hi bit for curve25519. */
    unsigned char tmp[MAX_FIELD_SIZE];
    uint32_t ignore = fe->size * 8 - fe->bits;
    uint32_t mask = (1 << (8 - ignore)) - 1;

    memcpy(tmp, raw, fe->size);

    tmp[fe->size - 1] &= mask;

    if (!fe_import(fe, r->x, tmp))
      return 0;
  } else {
    if (!fe_import(fe, r->x, raw))
      return 0;
  }

  /* Take the principle square root. */
  return ge_set_x(ec, r, r->x, -1);
}

static int
ge_export(curve_t *ec,
          unsigned char *raw,
          const ge_t *p) {
  /* [RFC7748] Section 5. */
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
    ge_zero(ec, r);
    return;
  }

  /* Y1 = 0 */
  if (fe_is_zero(fe, a->y)) {
    ge_zero(ec, r);
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
ge_add(curve_t *ec, ge_t *r, const ge_t *a, const ge_t *b) {
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
    fe_zero(fe, r->z);
    return;
  }

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->z, fe->one);
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
 * Projective Point
 */

static void
gej_zero(curve_t *ec, gej_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, fe->one);
  fe_zero(fe, r->z);
}

static void
gej_cleanse(curve_t *ec, gej_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->z);
}

static int
gej_validate(curve_t *ec, const gej_t *p) {
  prime_field_t *fe = &ec->fe;

  // B * y^2 * z = x^3 + A * x^2 * z + x * z^2
  const {x, z} = this;
  const x2 = x.redSqr();
  const x3 = x2.redMul(x);
  const z2 = z.redSqr();
  const ax2 = this.curve.a.redMul(x2).redMul(z);
  const by2 = x3.redIAdd(ax2).redIAdd(x.redMul(z2));
  const y2 = by2.redMul(this.curve.bi);

  // sqrt(y^2 * z^4) = y * z^2
  return y2.redMul(z).redJacobi() !== -1;
}

static int
gej_import(curve_t *ec, gej_t *r, const unsigned char *raw) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;

  if ((fe->bits & 7) != 0) {
    /* Ignore the hi bit for curve25519. */
    unsigned char tmp[MAX_FIELD_SIZE];
    uint32_t ignore = fe->size * 8 - fe->bits;
    uint32_t mask = (1 << (8 - ignore)) - 1;

    memcpy(tmp, raw, fe->size);

    tmp[fe->size - 1] &= mask;

    return fe_import(fe, r->x, tmp);
  }

  return fe_import(fe, r->x, raw);
}

static int
gej_export(curve_t *ec,
          unsigned char *raw,
          const gej_t *p) {
  /* [RFC7748] Section 5. */
  prime_field_t *fe = &ec->fe;
  fe_t x;

  fe_invert(fe, x, p->z);
  fe_mul(fe, x, x, p->x);
  fe_export(fe, raw, x);

  return fe_is_zero(fe, p->z) ^ 1;
}

static void
gej_swap(curve_t *ec, gej_t *a, gej_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_swap(fe, a->x, b->x, flag);
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
  fe_select(fe, r->z, a->z, b->z, flag);
}

static void
gej_set(curve_t *ec, gej_t *r, const gej_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
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
  fe_t e1, e2;

  /* X1 * Z2 == X2 * Z1 */
  fe_mul(fe, e1, a->x, b->z);
  fe_mul(fe, e2, b->x, a->z);

  return fe_equal(fe, e1, e2);
}

static void
gej_zero_cond(curve_t *ec, gej_t *r, const gej_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, fe->one, flag);
  fe_select(fe, r->z, a->z, fe->zero, flag);
}

static void
gej_dbl(curve_t *ec, gej_t *r, const gej_t *a) {
  prime_field_t *fe = &ec->fe;

  // https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-3
  // 2M + 2S + 4A + 1*a24

  // A = X1 + Z1
  const a = this.x.redAdd(this.z);

  // AA = A^2
  const aa = a.redSqr();

  // B = X1 - Z1
  const b = this.x.redSub(this.z);

  // BB = B^2
  const bb = b.redSqr();

  // C = AA - BB
  const c = aa.redSub(bb);

  // X3 = AA * BB
  const nx = aa.redMul(bb);

  // Z3 = C * (BB + a24 * C)
  const nz = c.redMul(bb.redIAdd(this.curve.a24.redMul(c)));

  return this.curve.xpoint(nx, nz);
}

diffAddDbl(p, q) {
  // https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#ladder-ladd-1987-m-3
  // Note that we swap P2 and P3 here (for consistency).
  // 6M + 4S + 8A + 1*a24
  assert(p instanceof XPoint);
  assert(q instanceof XPoint);

  // A = X2 + Z2
  const a = q.x.redAdd(q.z);

  // AA = A^2
  const aa = a.redSqr();

  // B = X2 - Z2
  const b = q.x.redSub(q.z);

  // BB = B^2
  const bb = b.redSqr();

  // E = AA - BB
  const e = aa.redSub(bb);

  // C = X3 + Z3
  const c = p.x.redAdd(p.z);

  // D = X3 - Z3
  const d = p.x.redSub(p.z);

  // DA = D * A
  const da = d.redMul(a);

  // CB = C * B
  const cb = c.redMul(b);

  // X5 = Z1 * (DA + CB)^2
  const nx = this.z.redMul(da.redAdd(cb).redSqr());

  // Z5 = X1 * (DA - CB)^2
  const nz = this.x.redMul(da.redISub(cb).redSqr());

  // X4 = AA * BB
  const dx = aa.redMul(bb);

  // Z4 = E * (BB + a24 * E)
  const dz = e.redMul(bb.redIAdd(this.curve.a24.redMul(e)));

  return [
    this.curve.xpoint(nx, nz),
    this.curve.xpoint(dx, dz)
  ];
}

ladderConst(k, rng) {
  // Multiply with the Montgomery Ladder.
  //
  // [MONT3] Algorithm 7, Page 16, Section 5.3.
  //         Algorithm 8, Page 16, Section 5.3.
  //
  // [RFC7748] Page 7, Section 5.
  //
  // Note that any clamping is meant to
  // be done _outside_ of this function.
  assert(k instanceof BN);
  assert(!k.red);

  const bits = Math.max(k.bitLength(), this.curve.fieldBits);
  const bytes = (bits + 7) >>> 3;

  // Recode scalar to base256.
  const exp = k.toArray('le', bytes);

  // Randomize if available.
  const point = rng ? this.randomize(rng) : this;

  // Clone points (for safe swapping).
  let a = point.clone();
  let b = this.curve.xpoint().clone();
  let swap = 0;

  // Climb the ladder.
  for (let i = bits - 1; i >= 0; i--) {
    const bit = (exp[i >> 3] >> (i & 7)) & 1;

    // Maybe swap.
    a.swap(b, swap ^ bit);

    // Single coordinate add+double.
    [a, b] = point.diffAddDbl(a, b);

    swap = bit;
  }

  // Finalize loop.
  a.swap(b, swap);

  return [a, b];
}



static void
gej_to_ge(curve_t *ec, ge_t *r, const gej_t *p) {
  // https://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#scaling-scale
  // 1I + 1M
  prime_field_t *fe = &ec->fe;
  fe_t a;

  /* A = 1 / Z1 */
  fe_invert(fe, a, p->z);

  /* X3 = X1 * A */
  fe_mul(fe, r->x, p->x, a);

  ge_set_x(ec, r, r->x, -1);
}

static void
gej_print(curve_t *ec, const gej_t *p) {
  prime_field_t *fe = &ec->fe;

  if (gej_is_zero(ec, p)) {
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
curve_init(curve_t *ec, const curve_def_t *def) {
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

/*
 * Curves
 */
