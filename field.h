#define BCRYPTO_HAS_GMP
#define BCRYPTO_EC_64BIT

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>

#include "util.h"
#include "mpn.h"
#include "mpz.h"

#include "fields/p224.h"
#include "fields/p256.h"
#include "fields/p384.h"
#include "fields/p521.h"
#include "fields/secp256k1.h"
#include "fields/p25519.h"
#include "fields/p448.h"

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

#define MAX_FIELD_BITS 521
#define MAX_SCALAR_BITS 521
#define MAX_FIELD_SIZE 66
#define MAX_SCALAR_SIZE 66

#define MAX_SCALAR_LIMBS \
  ((MAX_SCALAR_BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

#define MAX_FIELD_LIMBS \
  ((MAX_FIELD_BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

#ifdef BCRYPTO_EC_64BIT
typedef uint64_t fe_word_t;
#define FIELD_WORD_SIZE 64
#define MAX_FIELD_WORDS 9
#else
typedef uint32_t fe_word_t;
#define FIELD_WORD_SIZE 32
#define MAX_FIELD_WORDS 17
#endif

/*
 * Structs
 */

typedef void fe_add_func(fe_word_t *out1, const fe_word_t *arg1, const fe_word_t *arg2);
typedef void fe_sub_func(fe_word_t *out1, const fe_word_t *arg1, const fe_word_t *arg2);
typedef void fe_opp_func(fe_word_t *out1, const fe_word_t *arg1);
typedef void fe_mul_func(fe_word_t *out1, const fe_word_t *arg1, const fe_word_t *arg2);
typedef void fe_sqr_func(fe_word_t *out1, const fe_word_t *arg1);
typedef void fe_from_montgomery_func(fe_word_t *out1, const fe_word_t *arg1);
typedef void fe_nonzero_func(fe_word_t *out1, const fe_word_t *arg1);
typedef void fe_selectznz_func(fe_word_t *out1, unsigned char arg1, const fe_word_t *arg2, const fe_word_t *arg3);
typedef void fe_to_bytes_func(uint8_t *out1, const fe_word_t *arg1);
typedef void fe_from_bytes_func(fe_word_t *out1, const uint8_t *arg1);
typedef void fe_carry_func(fe_word_t *out1, const fe_word_t *arg1);
typedef void fe_invert_func(fe_word_t *out, const fe_word_t *in);
typedef int fe_sqrt_func(fe_word_t *out, const fe_word_t *in);
typedef int fe_isqrt_func(fe_word_t *out, const fe_word_t *u, const fe_word_t *v);
typedef void fe_scmul_121666(fe_word_t *out1, const fe_word_t *arg1);

typedef fe_word_t fe_t[MAX_FIELD_WORDS];
typedef fe_word_t *fe_ptr;

typedef mp_limb_t sc_t[MAX_SCALAR_LIMBS];
typedef mp_limb_t *sc_ptr;

typedef struct scalar_field_s {
  int endian;
  size_t size;
  size_t bits;
  size_t shift;
  mp_limb_t n[MAX_SCALAR_LIMBS * 4];
  mp_limb_t nh[MAX_SCALAR_LIMBS * 4];
  mp_limb_t m[MAX_SCALAR_LIMBS * 4];
  mp_size_t limbs;
  unsigned char raw[MAX_SCALAR_SIZE];
} scalar_field_t;

typedef struct scalar_def_s {
  size_t bits;
  const unsigned char n[MAX_FIELD_SIZE];
} scalar_def_t;

typedef struct prime_field_s {
  int endian;
  size_t size;
  size_t bits;
  size_t shift;
  size_t words;
  mp_limb_t p[MAX_SCALAR_LIMBS * 4];
  mp_size_t limbs;
  unsigned char raw[MAX_FIELD_SIZE];
  scalar_field_t sc;
  fe_add_func *add;
  fe_sub_func *sub;
  fe_opp_func *opp;
  fe_mul_func *mul;
  fe_sqr_func *square;
  fe_from_montgomery_func *from_montgomery;
  fe_nonzero_func *nonzero;
  fe_selectznz_func *selectznz;
  fe_to_bytes_func *to_bytes;
  fe_from_bytes_func *from_bytes;
  fe_carry_func *carry;
  fe_invert_func *invert;
  fe_sqrt_func *sqrt;
  fe_isqrt_func *isqrt;
  fe_scmul_121666 *scmul_121666;
  fe_t zero;
  fe_t one;
  fe_t two;
  fe_t three;
} prime_field_t;

typedef struct prime_def_s {
  size_t bits;
  size_t words;
  const unsigned char p[MAX_FIELD_SIZE];
  fe_add_func *add;
  fe_sub_func *sub;
  fe_opp_func *opp;
  fe_mul_func *mul;
  fe_sqr_func *square;
  fe_from_montgomery_func *from_montgomery;
  fe_nonzero_func *nonzero;
  fe_selectznz_func *selectznz;
  fe_to_bytes_func *to_bytes;
  fe_from_bytes_func *from_bytes;
  fe_carry_func *carry;
  fe_invert_func *invert;
  fe_sqrt_func *sqrt;
  fe_isqrt_func *isqrt;
  fe_scmul_121666 *scmul_121666;
} prime_def_t;

/*
 * Scalar
 */

static void
sc_zero(scalar_field_t *sc, sc_t r) {
  mpn_zero(r, sc->limbs);
}

static void
sc_cleanse(scalar_field_t *sc, sc_t r) {
  mpn_cleanse(r, sc->limbs);
}

static int
sc_import(scalar_field_t *sc, sc_t r, const unsigned char *raw) {
  int zero = is_zero(raw, sc->size);
  int lt = less_than(raw, sc->raw, sc->size, sc->endian);

  mpn_import(r, sc->limbs, raw, sc->size, sc->endian);

  return (zero ^ 1) & lt;
}

static void
sc_reduce(scalar_field_t *sc, sc_t r, const sc_t ap);

static int
sc_import_lax(scalar_field_t *sc, sc_t r, const unsigned char *raw) {
  int zero = is_zero(raw, sc->size);
  int lt = less_than(raw, sc->raw, sc->size, sc->endian);
  mp_limb_t tmp[MAX_SCALAR_LIMBS * 4];

  mpn_import(tmp, sc->limbs * 4, raw, sc->size, sc->endian);

  sc_reduce(sc, r, tmp);

  cleanse(tmp, sizeof(tmp));

  return (zero ^ 1) & lt;
}

static void
sc_export(scalar_field_t *sc, unsigned char *raw, const sc_t a) {
  mpn_export(raw, sc->size, a, sc->limbs, sc->endian);
}

static void
sc_set(scalar_field_t *sc, sc_t r, const sc_t a) {
  mpn_copyi(r, a, sc->limbs);
}

static void
sc_swap(scalar_field_t *sc, sc_t a, sc_t b, unsigned int flag) {
  cnd_swap(flag != 0, a, b, sc->limbs);
}

static void
sc_select(scalar_field_t *sc, sc_t r,
          const sc_t a, const sc_t b,
          unsigned int flag) {
  cnd_select(flag != 0, r, a, b, sc->limbs);
}

static void
sc_print(scalar_field_t *sc, const sc_t a) {
  mpn_print(a, sc->limbs, 16);
  printf("\n");
}

static void
sc_set_word(scalar_field_t *sc, sc_t r, mp_limb_t word) {
  mpn_zero(r, sc->limbs);
  r[0] = word;
}

static int
sc_equal(scalar_field_t *sc, const sc_t a, const sc_t b) {
  mp_limb_t z = 0;
  mp_size_t i;

  for (i = 0; i < sc->limbs; i++)
    z |= a[i] ^ b[i];

  return z == 0;
}

static int
sc_cmp_var(scalar_field_t *sc, const sc_t a, const sc_t b) {
  return mpn_cmp(a, b, sc->limbs);
}

static int
sc_is_zero(scalar_field_t *sc, const sc_t a) {
  mp_limb_t z = 0;
  mp_size_t i;

  for (i = 0; i < sc->limbs; i++)
    z |= a[i];

  return z == 0;
}

static int
sc_is_high_var(scalar_field_t *sc, const sc_t a) {
  return mpn_cmp(a, sc->nh, sc->limbs) > 0;
}

static void
sc_neg(scalar_field_t *sc, sc_t r, const sc_t a) {
  const mp_limb_t *np = sc->n;
  mp_size_t nn = sc->limbs;
  mp_limb_t rp[MAX_SCALAR_LIMBS];
  mp_limb_t zero = sc_is_zero(sc, a);
  mp_size_t i;

  mpn_zero(rp, nn);
  mpn_sub_n(rp, np, a, nn);

  for (i = 0; i < nn; i++)
    r[i] = rp[i] & -(zero ^ 1);
}

static void
sc_add(scalar_field_t *sc, sc_t r, const sc_t ap, const sc_t bp) {
  const mp_limb_t *np = sc->n;
  mp_size_t nn = sc->limbs + 1;
  mp_limb_t up[MAX_SCALAR_LIMBS + 1];
  mp_limb_t vp[MAX_SCALAR_LIMBS + 1];
  mp_limb_t c;

  assert(np[nn - 1] == 0);

  mpn_copyi(up, ap, nn - 1);
  mpn_copyi(vp, bp, nn - 1);

  up[nn - 1] = 0;
  vp[nn - 1] = 0;

  /* r = a + b */
  mpn_add_n(up, up, vp, nn);

  /* r = r - n if u >= n */
  c = mpn_sub_n(vp, up, np, nn);
  cnd_swap(c == 0, up, vp, nn);

  mpn_copyi(r, up, nn - 1);
}

static void
sc_sub(scalar_field_t *sc, sc_t r, const sc_t a, const sc_t b) {
  sc_t mb;
  sc_neg(sc, mb, b);
  sc_add(sc, r, a, mb);
  sc_cleanse(sc, mb);
}

static void
sc_reduce(scalar_field_t *sc, sc_t r, const mp_limb_t *ap) {
  /* Barrett reduction. */
  const mp_limb_t *np = sc->n;
  const mp_limb_t *mp = sc->m;
  mp_size_t shift = sc->shift;
  mp_size_t nn = sc->limbs;
  mp_limb_t qp[MAX_SCALAR_LIMBS * 4];
  mp_limb_t up[MAX_SCALAR_LIMBS * 4];
  mp_limb_t vp[MAX_SCALAR_LIMBS * 4];
  mp_limb_t *hp = qp;
  mp_limb_t c;

  mpn_zero(qp, nn * 4);
  mpn_zero(up, nn * 4);
  mpn_zero(vp, nn * 4);

  /* q = a * m */
  mpn_mul_n(qp, ap, mp, nn * 2);

  /* h = q >> k */
  hp += shift / GMP_NUMB_BITS;
  shift %= GMP_NUMB_BITS;

  /* Realign (necessary for curves like p521). */
  if (shift != 0)
    mpn_rshift(hp, hp, nn * 2, shift);

  /* u = r - h * n */
  mpn_mul_n(vp, hp, np, nn * 2);
  c = mpn_sub_n(up, ap, vp, nn * 4);
  assert(c == 0);

  /* u = u - n if u >= n */
  c = mpn_sub_n(vp, up, np, nn * 2);
  cnd_swap(c == 0, up, vp, nn * 2);

  mpn_copyi(r, up, nn);
}

static void
sc_mul(scalar_field_t *sc, sc_t r, const sc_t a, const sc_t b) {
  mp_limb_t ap[MAX_SCALAR_LIMBS * 4];

  mpn_zero(ap, sc->limbs * 4);
  mpn_mul_n(ap, a, b, sc->limbs);

  sc_reduce(sc, r, ap);
}

static void
sc_sqr(scalar_field_t *sc, sc_t r, const sc_t a) {
  mp_limb_t ap[MAX_SCALAR_LIMBS * 4];

  mpn_zero(ap, sc->limbs * 4);
  mpn_sqr(ap, a, sc->limbs);

  sc_reduce(sc, r, ap);
}

static int
sc_invert_var(scalar_field_t *sc, sc_t r, const sc_t a) {
  mpz_t rn, an, nn;

  if (sc_is_zero(sc, a)) {
    sc_zero(sc, r);
    return 0;
  }

  mpz_init(rn); /* Warning: allocation. */
  mpz_roinit_n(an, a, sc->limbs);
  mpz_roinit_n(nn, sc->n, sc->limbs);

  assert(mpz_invert(rn, an, nn));

  mpn_set_mpz(r, rn, sc->limbs);

  mpz_clear(rn);

  return 1;
}

static int
sc_invert(scalar_field_t *sc, sc_t r, const sc_t a) {
  mp_limb_t zero = sc_is_zero(sc, a);
  mp_size_t i;
  sc_t e, b;

  mpn_copyi(e, sc->n, sc->limbs);
  mpn_sub_1(e, e, sc->limbs, 2);

  sc_set_word(sc, b, 1);

  for (i = sc->bits - 1; i >= 0; i--) {
    mp_limb_t bit = e[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS);

    sc_sqr(sc, b, b);

    if (bit & 1)
      sc_mul(sc, b, b, a);
  }

  for (i = 0; i < sc->limbs; i++)
    r[i] = b[i] & -(zero ^ 1);

  return zero ^ 1;
}

static size_t
sc_bitlen(scalar_field_t *sc, const sc_t a) {
  return mpn_bitlen(a, sc->limbs);
}

static void
sc_naf(scalar_field_t *sc, int32_t *naf,
       const sc_t x, size_t width, size_t max) {
  /* Computing the NAF of a positive integer.
   *
   * [GECC] Algorithm 3.30, Page 98, Section 3.3.
   *        Algorithm 3.35, Page 100, Section 3.3.
   */
  mp_size_t nn = sc->limbs;
  mp_limb_t k[MAX_SCALAR_LIMBS];
  int32_t size = 1 << (width + 1);
  size_t i = 0;

  assert(width + 1 < GMP_LIMB_BITS);

  mpn_zero(k, MAX_SCALAR_LIMBS);
  mpn_copyi(k, x, nn);

  while (k[0] != 0) {
    int32_t z = 0;
    size_t shift = 1;
    size_t j;

    if (k[0] & 1) {
      int32_t mod = k[0] & (size - 1);

      if (mod > (size >> 1) - 1)
        z = (size >> 1) - mod;
      else
        z = mod;

      if (z < 0)
        mpn_add_1(k, k, nn, -z);
      else
        mpn_sub_1(k, k, nn, z);
    }

    naf[i++] = z;

    /* Optimization: shift by word if possible. */
    if (k[0] != 0 && (k[0] & (size - 1)) == 0)
      shift = width + 1;

    for (j = 1; j < shift; j++)
      naf[i++] = 0;

    mpn_rshift(k, k, nn, shift);
  }

  assert(i <= max);

  for (; i < max; i++)
    naf[i] = 0;

  mpn_zero(k, nn);
}

/*
 * Field Element
 */

static void
fe_zero(prime_field_t *fe, fe_t r) {
  memset(r, 0, fe->words * sizeof(fe_word_t));
}

static void
fe_cleanse(prime_field_t *fe, fe_t r) {
  cleanse(r, fe->words * sizeof(fe_word_t));
}

static int
fe_import(prime_field_t *fe, fe_t r, const unsigned char *raw) {
  unsigned char tmp[MAX_FIELD_SIZE];
  size_t i;

  if (fe->from_montgomery) {
    /* Hack: Use a constant time barret reduction
     *       to montgomerize the field element.
     */
    mp_limb_t xp[MAX_FIELD_LIMBS * 4];

    /* We must be aligned to the limb size. */
    /* Every montgomerized field satisfies this requirement. */
    assert((fe->shift % GMP_NUMB_BITS) == 0);

    /* x = (x << shift) mod p */
    mpn_zero(xp, fe->limbs * 4);
    mpn_import(xp + fe->limbs, fe->limbs, raw, fe->size, fe->endian);
    sc_reduce(&fe->sc, xp, xp);

    /* Export as little endian. */
    mpn_export_le(tmp, fe->size, xp, fe->limbs);

    fe->from_bytes(r, tmp);
  } else {
    if (fe->endian == 1) {
      /* Swap endianness. */
      for (i = 0; i < fe->size; i++)
        tmp[i] = raw[fe->size - 1 - i];

      fe->from_bytes(r, tmp);
    } else {
      fe->from_bytes(r, raw);
    }
  }

  return less_than(raw, fe->raw, fe->size, fe->endian);
}

static int
fe_import_be(prime_field_t *fe, fe_t r, const unsigned char *raw) {
  if (fe->endian == -1) {
    unsigned char tmp[MAX_FIELD_SIZE];
    size_t i;

    for (i = 0; i < fe->size; i++)
      tmp[i] = raw[fe->size - 1 - i];

    return fe_import(fe, r, tmp);
  }

  return fe_import(fe, r, raw);
}

static void
fe_export(prime_field_t *fe, unsigned char *raw, const fe_t a) {
  int i = 0;
  int j = fe->size - 1;

  if (fe->from_montgomery) {
    fe_t tmp;
    fe->from_montgomery(tmp, a);
    fe->to_bytes(raw, tmp);
  } else {
    fe->to_bytes(raw, a);
  }

  if (fe->endian == 1) {
    while (i < j) {
      unsigned char t = raw[j];

      raw[j] = raw[i];
      raw[i] = t;

      i += 1;
      j -= 1;
    }
  }
}

static void
fe_swap(prime_field_t *fe, fe_t a, fe_t b, unsigned int flag) {
  fe_word_t cond = (flag != 0);
  fe_word_t mask = -cond;
  size_t i;

  for (i = 0; i < fe->words; i++) {
    fe_word_t word = (a[i] ^ b[i]) & mask;

    a[i] ^= word;
    b[i] ^= word;
  }
}

static void
fe_select(prime_field_t *fe, fe_t r,
          const fe_t a, const fe_t b,
          unsigned int flag) {
  fe->selectznz(r, flag != 0, a, b);
}

static void
fe_set(prime_field_t *fe, fe_t r, const fe_t a) {
  memcpy(r, a, sizeof(fe_t));
}

static int
fe_set_limbs(prime_field_t *fe, fe_t r, const mp_limb_t *p, mp_size_t n) {
  unsigned char tmp[MAX_FIELD_SIZE];

  assert(n <= fe->limbs);

  mpn_export(tmp, fe->size, p, n, fe->endian);

  return fe_import(fe, r, tmp);
}

static void
fe_get_limbs(prime_field_t *fe, mp_limb_t *r, const fe_t a) {
  unsigned char tmp[MAX_FIELD_SIZE];

  fe_export(fe, tmp, a);

  mpn_import(r, fe->limbs, tmp, fe->size, fe->endian);
}

static int
fe_set_num(prime_field_t *fe, fe_t r, const mpz_t a) {
  return fe_set_limbs(fe, r, mpz_limbs_read(a), mpz_size(a));
}

static void
fe_print(prime_field_t *fe, const fe_t a) {
  mp_limb_t xp[MAX_FIELD_LIMBS];

  fe_get_limbs(fe, xp, a);

  mpn_print(xp, fe->limbs, 16);
  printf("\n");
}

static int
fe_set_sc(prime_field_t *fe, scalar_field_t *sc, fe_t r, const sc_t a) {
  unsigned char tmp[MAX_SCALAR_SIZE];
  int ret;

  sc_export(sc, tmp, a);

  ret = fe_import(fe, r, tmp);

  cleanse(tmp, sizeof(tmp));

  return ret;
}

static int
fe_get_sc(prime_field_t *fe, scalar_field_t *sc, sc_t r, const fe_t a) {
  unsigned char tmp[MAX_FIELD_SIZE];
  fe_export(fe, tmp, a);
  return sc_import_lax(sc, r, tmp);
}

static void
fe_set_word(prime_field_t *fe, fe_t r, uint32_t word) {
  if (fe->from_montgomery) {
    unsigned char tmp[MAX_FIELD_SIZE];

    memset(tmp, 0, sizeof(tmp));

    if (fe->endian == 1) {
      tmp[fe->size - 4] = (word >> 24) & 0xff;
      tmp[fe->size - 3] = (word >> 16) & 0xff;
      tmp[fe->size - 2] = (word >> 8) & 0xff;
      tmp[fe->size - 1] = word & 0xff;
    } else {
      tmp[0] = word & 0xff;
      tmp[1] = (word >> 8) & 0xff;
      tmp[2] = (word >> 16) & 0xff;
      tmp[3] = (word >> 24) & 0xff;
    }

    fe_import(fe, r, tmp);
  } else {
    /* Maybe need to deserialize? */
    fe_zero(fe, r);
    r[0] = word;
  }
}

static int
fe_is_zero(prime_field_t *fe, const fe_t a) {
  fe_word_t z = 0;

  if (fe->nonzero) {
    fe->nonzero(&z, a);
  } else {
    unsigned char tmp[MAX_FIELD_SIZE];
    size_t i;

    fe->to_bytes(tmp, a);

    for (i = 0; i < fe->size; i++)
      z |= (fe_word_t)tmp[i];
  }

  return z == 0;
}

static int
fe_equal(prime_field_t *fe, const fe_t a, const fe_t b) {
  fe_t c;

  fe->sub(c, a, b);

  if (fe->carry)
    fe->carry(c, c);

  return fe_is_zero(fe, c);
}

static int
fe_is_odd(prime_field_t *fe, const fe_t a) {
  int sign;

  if (fe->from_montgomery) {
    fe_t tmp;
    fe->from_montgomery(tmp, a);
    sign = tmp[0] & 1;
  } else {
    /* Maybe need to serialize? */
    sign = a[0] & 1;
  }

  return sign;
}

static void
fe_neg(prime_field_t *fe, fe_t r, const fe_t a) {
  fe->opp(r, a);

  if (fe->carry)
    fe->carry(r, r);
}

static void
fe_neg_noreduce(prime_field_t *fe, fe_t r, const fe_t a) {
  fe->opp(r, a);
}

static void
fe_neg_cond(prime_field_t *fe, fe_t r, const fe_t a, unsigned int flag) {
  fe_t b;
  fe_neg(fe, b, a);
  fe_select(fe, r, a, b, flag);
}

static void
fe_set_odd(prime_field_t *fe, fe_t r, const fe_t a, unsigned int odd) {
  fe_neg_cond(fe, r, a, fe_is_odd(fe, a) ^ (odd != 0));
}

static void
fe_add(prime_field_t *fe, fe_t r, const fe_t a, const fe_t b) {
  fe->add(r, a, b);

  if (fe->carry)
    fe->carry(r, r);
}

static void
fe_add_noreduce(prime_field_t *fe, fe_t r, const fe_t a, const fe_t b) {
  fe->add(r, a, b);
}

static void
fe_sub(prime_field_t *fe, fe_t r, const fe_t a, const fe_t b) {
  fe->sub(r, a, b);

  if (fe->carry)
    fe->carry(r, r);
}

static void
fe_sub_noreduce(prime_field_t *fe, fe_t r, const fe_t a, const fe_t b) {
  fe->sub(r, a, b);
}

static void
fe_reduce(prime_field_t *fe, fe_t r, const fe_t a) {
  if (fe->carry)
    fe->carry(r, a);
}

static void
fe_mulw(prime_field_t *fe, fe_t r, const fe_t a, long b) {
  fe_t c;
  long bits = 0;
  int neg = (b < 0);
  long tmp, i;

  if (neg)
    b = -b;

  tmp = b;

  while (tmp != 0) {
    bits += 1;
    tmp >>= 1;
  }

  fe_set(fe, c, a);
  fe_zero(fe, r);

  for (i = bits - 1; i >= 0; i--) {
    fe_add(fe, r, r, r);

    if ((b >> i) & 1)
      fe_add(fe, r, r, c);
  }

  if (neg)
    fe_neg(fe, r, r);
}

static void
fe_mul(prime_field_t *fe, fe_t r, const fe_t a, const fe_t b) {
  fe->mul(r, a, b);
}

static void
fe_sqr(prime_field_t *fe, fe_t r, const fe_t a) {
  fe->square(r, a);
}

static void
fe_mul121666(prime_field_t *fe, fe_t r, const fe_t a) {
  assert(fe->scmul_121666 != NULL);
  fe->scmul_121666(r, a);
}

static int
fe_invert_var(prime_field_t *fe, fe_t r, const fe_t a) {
  int ret = !fe_is_zero(fe, a);

  if (ret) {
#ifdef BCRYPTO_HAS_GMP
    mp_limb_t gp[MAX_FIELD_LIMBS + 1];
    mp_limb_t sp[MAX_FIELD_LIMBS + 1];
    mp_limb_t up[MAX_FIELD_LIMBS + 1];
    mp_limb_t vp[MAX_FIELD_LIMBS + 1];
    mp_size_t sn = fe->limbs + 1;
    mp_size_t gn;

    fe_get_limbs(fe, up, a);
    mpn_copyi(vp, fe->p, fe->limbs);

    gn = mpn_gcdext(gp, sp, &sn, up, fe->limbs, vp, fe->limbs);

    assert(gn == 1);
    assert(gp[0] == 1);

    if (sn < 0) {
      mpn_sub(sp, fe->p, fe->limbs, sp, -sn);
      sn = fe->limbs;
    }

    fe_set_limbs(fe, r, sp, sn);
#else
    mp_limb_t ap[MAX_FIELD_LIMBS];
    mpz_t rn, an, pn;

    fe_get_limbs(fe, ap, a);

    mpz_init(rn); /* Warning: allocation. */
    mpz_roinit_n(an, ap, fe->limbs);
    mpz_roinit_n(pn, fe->p, fe->limbs);

    assert(mpz_invert(rn, an, pn));

    fe_set_num(fe, r, rn);

    mpz_clear(rn);
#endif
  } else {
    fe_zero(fe, r);
  }

  return ret;
}

static int
fe_invert_slow(prime_field_t *fe, fe_t r, const fe_t a) {
  mp_limb_t e[MAX_FIELD_LIMBS];
  fe_word_t zero = fe_is_zero(fe, a);
  mp_size_t i;
  size_t j;
  fe_t b;

  mpn_copyi(e, fe->p, fe->limbs);
  mpn_sub_1(e, e, fe->limbs, 2);

  fe_set(fe, b, fe->one);

  for (i = fe->bits - 1; i >= 0; i--) {
    mp_limb_t bit = e[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS);

    fe_sqr(fe, b, b);

    if (bit & 1)
      fe_mul(fe, b, b, a);
  }

  for (j = 0; j < fe->words; j++)
    r[j] = b[j] & -(zero ^ 1);

  return zero ^ 1;
}

static int
fe_invert(prime_field_t *fe, fe_t r, const fe_t a) {
  if (fe->invert) {
    /* Faster inversion chain. */
    int zero = fe_is_zero(fe, a);
    fe->invert(r, a);
    return zero ^ 1;
  }

  return fe_invert_slow(fe, r, a);
}

static int
fe_sqrt(prime_field_t *fe, fe_t r, const fe_t a) {
  int ret;

  if (fe->sqrt) {
    /* Faster square root chain. */
    ret = fe->sqrt(r, a);
  } else {
    mp_limb_t ap[MAX_FIELD_LIMBS];
    mpz_t rn, an, pn;

    fe_get_limbs(fe, ap, a);

    mpz_init(rn); /* Warning: allocation. */
    mpz_roinit_n(an, ap, fe->limbs);
    mpz_roinit_n(pn, fe->p, fe->limbs);

    ret = mpz_sqrtm(rn, an, pn);

    if (ret)
      fe_set_num(fe, r, rn);

    mpz_clear(rn);
  }

  return ret;
}

static int
fe_is_square_var(prime_field_t *fe, const fe_t a) {
  mp_limb_t ap[MAX_FIELD_LIMBS];
  mpz_t an, pn;

  fe_get_limbs(fe, ap, a);

  mpz_roinit_n(an, ap, fe->limbs);
  mpz_roinit_n(pn, fe->p, fe->limbs);

  return mpz_jacobi(an, pn) >= 0;
}

static int
fe_is_square(prime_field_t *fe, const fe_t a) {
  int ret;

  if (fe->sqrt) {
    fe_t tmp;
    ret = fe->sqrt(tmp, a);
  } else {
    mp_limb_t ap[MAX_FIELD_LIMBS];
    mpz_t an, pn;

    fe_get_limbs(fe, ap, a);

    mpz_roinit_n(an, ap, fe->limbs);
    mpz_roinit_n(pn, fe->p, fe->limbs);

    ret = mpz_legendre(an, pn) >= 0;
  }

  return ret;
}

static int
fe_isqrt(prime_field_t *fe, fe_t r, const fe_t u, const fe_t v) {
  int ret = 1;

  if (fe->isqrt) {
    /* Fast inverse square root chain. */
    ret = fe->isqrt(r, u, v);
  } else {
    fe_t z;

    ret &= fe_invert(fe, z, v);

    fe_mul(fe, z, z, u);

    ret &= fe_sqrt(fe, r, z);
  }

  return ret;
}

/*
 * Scalar Field
 */

static void
scalar_field_set(scalar_field_t *sc,
                 const unsigned char *modulus,
                 size_t bits,
                 int endian) {
  sc->endian = endian;
  sc->limbs = (bits + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
  sc->size = (bits + 7) / 8;
  sc->bits = bits;
  sc->shift = bits * 2;

  mpn_import_be(sc->n, ARRAY_SIZE(sc->n), modulus, sc->size);

  mpn_rshift(sc->nh, sc->n, ARRAY_SIZE(sc->n), 1);

  /* Compute the barrett reduction constant `m`:
   *
   *   m = (1 << (bits * 2)) / n
   */
#ifndef BCRYPTO_HAS_GMP
  {
    mpz_t m, n;
    mpz_init(m);
    mpz_roinit_n(n, sc->n, sc->limbs);
    mpz_set_ui(m, 0);
    mpz_setbit(m, sc->shift);
    mpz_tdiv_q(m, m, n);
    mpn_set_mpz(sc->m, m, ARRAY_SIZE(sc->m));
    mpz_clear(m);
  }
#else
  {
    /* Maintain the philosophy of zero allocations. */
    mp_limb_t x[MAX_SCALAR_LIMBS * 2 + 1];
    mp_size_t index = sc->shift / GMP_NUMB_BITS;
    mp_limb_t bit = sc->shift % GMP_NUMB_BITS;

    mpn_zero(sc->m, ARRAY_SIZE(sc->m));
    mpn_zero(x, ARRAY_SIZE(x));

    x[index] = (mp_limb_t)1 << bit;

    mpn_tdiv_qr(sc->m, x, 0, x, index + 1, sc->n, sc->limbs);
  }
#endif

  mpn_export(sc->raw, sc->size, sc->n, ARRAY_SIZE(sc->n), sc->endian);
}

static void
scalar_field_init(scalar_field_t *sc, const scalar_def_t *def, int endian) {
  scalar_field_set(sc, def->n, def->bits, endian);
}

/*
 * Prime Field
 */

static void
prime_field_init(prime_field_t *fe, const prime_def_t *def, int endian) {
  fe->endian = endian;
  fe->limbs = (def->bits + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
  fe->size = (def->bits + 7) / 8;
  fe->bits = def->bits;
  fe->shift = def->bits;
  fe->words = def->words;

  if ((fe->shift % FIELD_WORD_SIZE) != 0)
    fe->shift += FIELD_WORD_SIZE - (fe->shift % FIELD_WORD_SIZE);

  mpn_import_be(fe->p, ARRAY_SIZE(fe->p), def->p, fe->size);

  mpn_export(fe->raw, fe->size, fe->p, ARRAY_SIZE(fe->p), fe->endian);

  scalar_field_set(&fe->sc, def->p, def->bits, endian);

  fe->add = def->add;
  fe->sub = def->sub;
  fe->opp = def->opp;
  fe->mul = def->mul;
  fe->square = def->square;
  fe->from_montgomery = def->from_montgomery;
  fe->nonzero = def->nonzero;
  fe->selectznz = def->selectznz;
  fe->to_bytes = def->to_bytes;
  fe->from_bytes = def->from_bytes;
  fe->carry = def->carry;
  fe->invert = def->invert;
  fe->sqrt = def->sqrt;
  fe->isqrt = def->isqrt;
  fe->scmul_121666 = def->scmul_121666;

  fe_set_word(fe, fe->zero, 0);
  fe_set_word(fe, fe->one, 1);
  fe_set_word(fe, fe->two, 2);
  fe_set_word(fe, fe->three, 3);
}

/*
 * Fields
 */

/*
 * P192
 */

static const prime_def_t field_p192 = {
  .bits = 192,
  .words = 192 / FIELD_WORD_SIZE,
  /* 2^192 - 2^64 - 1 (= 3 mod 4) */
  .p = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
  },
  .add = NULL,
  .sub = NULL,
  .opp = NULL,
  .mul = NULL,
  .square = NULL,
  .from_montgomery = NULL,
  .nonzero = NULL,
  .selectznz = NULL,
  .to_bytes = NULL,
  .from_bytes = NULL,
  .carry = NULL,
  .invert = NULL,
  .sqrt = NULL,
  .isqrt = NULL,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q192 = {
  .bits = 192,
  .n = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0x99, 0xde, 0xf8, 0x36,
    0x14, 0x6b, 0xc9, 0xb1, 0xb4, 0xd2, 0x28, 0x31
  }
};

/*
 * P224
 */

static const prime_def_t field_p224 = {
  .bits = 224,
  .words = P224_FIELD_WORDS,
  /* 2^224 - 2^96 + 1 (no congruence) */
  .p = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01
  },
  .add = fiat_p224_add,
  .sub = fiat_p224_sub,
  .opp = fiat_p224_opp,
  .mul = fiat_p224_mul,
  .square = fiat_p224_square,
  .from_montgomery = fiat_p224_from_montgomery,
  .nonzero = fiat_p224_nonzero,
  .selectznz = fiat_p224_selectznz,
  .to_bytes = fiat_p224_to_bytes,
  .from_bytes = fiat_p224_from_bytes,
  .carry = NULL,
  .invert = NULL,
  .sqrt = NULL,
  .isqrt = NULL,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q224 = {
  .bits = 224,
  .n = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x16, 0xa2,
    0xe0, 0xb8, 0xf0, 0x3e, 0x13, 0xdd, 0x29, 0x45,
    0x5c, 0x5c, 0x2a, 0x3d
  }
};

/*
 * P256
 */

static const prime_def_t field_p256 = {
  .bits = 256,
  .words = P256_FIELD_WORDS,
  /* 2^256 - 2^224 + 2^192 + 2^96 - 1 (= 3 mod 4) */
  .p = {
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
  },
  .add = fiat_p256_add,
  .sub = fiat_p256_sub,
  .opp = fiat_p256_opp,
  .mul = fiat_p256_mul,
  .square = fiat_p256_square,
  .from_montgomery = fiat_p256_from_montgomery,
  .nonzero = fiat_p256_nonzero,
  .selectznz = fiat_p256_selectznz,
  .to_bytes = fiat_p256_to_bytes,
  .from_bytes = fiat_p256_from_bytes,
  .carry = NULL,
  .invert = p256_fe_invert,
  .sqrt = p256_fe_sqrt,
  .isqrt = NULL,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q256 = {
  .bits = 256,
  .n = {
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xbc, 0xe6, 0xfa, 0xad, 0xa7, 0x17, 0x9e, 0x84,
    0xf3, 0xb9, 0xca, 0xc2, 0xfc, 0x63, 0x25, 0x51
  }
};

/*
 * P384
 */

static const prime_def_t field_p384 = {
  .bits = 384,
  .words = P384_FIELD_WORDS,
  /* 2^384 - 2^128 - 2^96 + 2^32 - 1 (= 3 mod 4) */
  .p = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xff, 0xff, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff
  },
  .add = fiat_p384_add,
  .sub = fiat_p384_sub,
  .opp = fiat_p384_opp,
  .mul = fiat_p384_mul,
  .square = fiat_p384_square,
  .from_montgomery = fiat_p384_from_montgomery,
  .nonzero = fiat_p384_nonzero,
  .selectznz = fiat_p384_selectznz,
  .to_bytes = fiat_p384_to_bytes,
  .from_bytes = fiat_p384_from_bytes,
  .carry = NULL,
  .invert = p384_fe_invert,
  .sqrt = NULL,
  .isqrt = NULL,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q384 = {
  .bits = 384,
  .n = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xc7, 0x63, 0x4d, 0x81, 0xf4, 0x37, 0x2d, 0xdf,
    0x58, 0x1a, 0x0d, 0xb2, 0x48, 0xb0, 0xa7, 0x7a,
    0xec, 0xec, 0x19, 0x6a, 0xcc, 0xc5, 0x29, 0x73
  }
};

/*
 * P521
 */

static const prime_def_t field_p521 = {
  .bits = 521,
  .words = P521_FIELD_WORDS,
  /* 2^521 - 1 (= 3 mod 4) */
  .p = {
    0x01, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff
  },
  .add = fiat_p521_add,
  .sub = fiat_p521_sub,
  .opp = fiat_p521_opp,
  .mul = fiat_p521_carry_mul,
  .square = fiat_p521_carry_square,
  .from_montgomery = NULL,
  .nonzero = NULL,
  .selectznz = fiat_p521_selectznz,
  .to_bytes = fiat_p521_to_bytes,
  .from_bytes = fiat_p521_from_bytes,
  .carry = fiat_p521_carry,
  .invert = p521_fe_invert,
  .sqrt = NULL,
  .isqrt = NULL,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q521 = {
  .bits = 521,
  .n = {
    0x01, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xfa, 0x51, 0x86, 0x87, 0x83, 0xbf, 0x2f,
    0x96, 0x6b, 0x7f, 0xcc, 0x01, 0x48, 0xf7, 0x09,
    0xa5, 0xd0, 0x3b, 0xb5, 0xc9, 0xb8, 0x89, 0x9c,
    0x47, 0xae, 0xbb, 0x6f, 0xb7, 0x1e, 0x91, 0x38,
    0x64, 0x09
  }
};

/*
 * P256K1
 */

static const prime_def_t field_p256k1 = {
  .bits = 256,
  .words = SECP256K1_FIELD_WORDS,
  /* 2^256 - 2^32 - 977 (= 3 mod 4) */
  .p = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xfc, 0x2f
  },
  .add = fiat_secp256k1_add,
  .sub = fiat_secp256k1_sub,
  .opp = fiat_secp256k1_opp,
  .mul = fiat_secp256k1_mul,
  .square = fiat_secp256k1_square,
  .from_montgomery = fiat_secp256k1_from_montgomery,
  .nonzero = fiat_secp256k1_nonzero,
  .selectznz = fiat_secp256k1_selectznz,
  .to_bytes = fiat_secp256k1_to_bytes,
  .from_bytes = fiat_secp256k1_from_bytes,
  .carry = NULL,
  .invert = secp256k1_fe_invert,
  .sqrt = secp256k1_fe_sqrt,
  .isqrt = secp256k1_fe_isqrt,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q256k1 = {
  .bits = 256,
  .n = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe,
    0xba, 0xae, 0xdc, 0xe6, 0xaf, 0x48, 0xa0, 0x3b,
    0xbf, 0xd2, 0x5e, 0x8c, 0xd0, 0x36, 0x41, 0x41
  }
};

/*
 * P25519
 */

static const prime_def_t field_p25519 = {
  .bits = 255,
  .words = P25519_FIELD_WORDS,
  /* 2^255 - 19 (= 5 mod 8) */
  .p = {
    0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xed
  },
  .add = fiat_25519_add,
  .sub = fiat_25519_sub,
  .opp = fiat_25519_opp,
  .mul = fiat_25519_carry_mul,
  .square = fiat_25519_carry_square,
  .from_montgomery = NULL,
  .nonzero = NULL,
  .selectznz = fiat_25519_selectznz,
  .to_bytes = fiat_25519_to_bytes,
  .from_bytes = fiat_25519_from_bytes,
  .carry = fiat_25519_carry,
  .invert = p25519_fe_invert,
  .sqrt = p25519_fe_sqrt,
  .isqrt = p25519_fe_isqrt,
  .scmul_121666 = fiat_25519_carry_scmul_121666
};

static const scalar_def_t field_q25519 = {
  .bits = 253,
  .n = {
    0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x14, 0xde, 0xf9, 0xde, 0xa2, 0xf7, 0x9c, 0xd6,
    0x58, 0x12, 0x63, 0x1a, 0x5c, 0xf5, 0xd3, 0xed
  }
};

/*
 * P448
 */

static const prime_def_t field_p448 = {
  .bits = 448,
  .words = P448_FIELD_WORDS,
  /* 2^448 - 2^224 - 1 (= 3 mod 4) */
  .p = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
  },
  .add = fiat_p448_add,
  .sub = fiat_p448_sub,
  .opp = fiat_p448_opp,
  .mul = fiat_p448_carry_mul,
  .square = fiat_p448_carry_square,
  .from_montgomery = NULL,
  .nonzero = NULL,
  .selectznz = fiat_p448_selectznz,
  .to_bytes = fiat_p448_to_bytes,
  .from_bytes = fiat_p448_from_bytes,
  .carry = fiat_p448_carry,
  .invert = p448_fe_invert,
  .sqrt = p448_fe_sqrt,
  .isqrt = p448_fe_isqrt,
  .scmul_121666 = NULL
};

static const scalar_def_t field_q448 = {
  .bits = 446,
  .n = {
    0x3f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0x7c, 0xca, 0x23, 0xe9,
    0xc4, 0x4e, 0xdb, 0x49, 0xae, 0xd6, 0x36, 0x90,
    0x21, 0x6c, 0xc2, 0x72, 0x8d, 0xc5, 0x8f, 0x55,
    0x23, 0x78, 0xc2, 0x92, 0xab, 0x58, 0x44, 0xf3
  }
};
