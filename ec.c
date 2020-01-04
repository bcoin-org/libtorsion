#define BCRYPTO_HAS_GMP
#define BCRYPTO_EC_64BIT

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef _WIN32
/* For SecureZeroMemory (actually defined in winbase.h). */
#include <windows.h>
#endif

#include "fields/p224.h"
#include "fields/p256.h"
#include "fields/p384.h"
#include "fields/p521.h"
#include "fields/secp256k1.h"
#include "fields/p25519.h"
#include "fields/p448.h"
#include "hash.h"

#ifdef BCRYPTO_HAS_GMP
#include <gmp.h>
#else
#include "mini-gmp.h"
#endif

#ifndef BCRYPTO_HAS_GMP
#define GMP_LIMB_BITS (sizeof(mp_limb_t) * CHAR_BIT)
#define GMP_NAIL_BITS 0
#define GMP_NUMB_BITS GMP_LIMB_BITS
#define GMP_NUMB_MASK (~((mp_limb_t)0))
#define GMP_NUMB_MAX GMP_NUMB_MASK
#define GMP_NAIL_MASK 0
#endif

/* Nails probably break our code. */
#if GMP_NAIL_BITS != 0 || GMP_LIMB_BITS != GMP_NUMB_BITS
#error "please use a build of gmp without nails"
#endif

#if (GMP_NUMB_BITS & 31) != 0
#error "invalid gmp bit alignment"
#endif

#if CHAR_BIT != 8
#error "sane char widths please"
#endif

#ifdef BCRYPTO_EC_64BIT
typedef uint64_t fe_word_t;
#define FIELD_WORD_SIZE 64
#define MAX_FIELD_WORDS 9
#else
typedef uint32_t fe_word_t;
#define FIELD_WORD_SIZE 32
#define MAX_FIELD_WORDS 17
#endif

#define MAX_FIELD_BITS 521
#define MAX_SCALAR_BITS 521
#define MAX_FIELD_SIZE 66
#define MAX_SCALAR_SIZE 66

#define MAX_SCALAR_LIMBS \
  ((MAX_SCALAR_BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

#define MAX_FIELD_LIMBS \
  ((MAX_FIELD_BITS + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)

#define NAF_WIDTH 4
#define NAF_WIDTH_PRE 8

/*
 * Scalar Field
 */

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

/*
 * Prime Field
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
  fe_t four;
  fe_t mone;
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
 * Short Weierstrass
 */

/* wge = weierstrass group element (affine) */
typedef struct wge_s {
  fe_t x;
  fe_t y;
  int inf;
} wge_t;

/* jge = jacobian group element */
typedef struct jge_s {
  fe_t x;
  fe_t y;
  fe_t z;
} jge_t;

typedef struct wei_s {
  int hash;
  prime_field_t fe;
  scalar_field_t sc;
  mp_limb_t pmodn[MAX_SCALAR_LIMBS * 4];
  fe_t red_n;
  fe_t a;
  fe_t b;
  wge_t g;
  sc_t blind;
  wge_t unblind;
  wge_t points[(1 << NAF_WIDTH_PRE) - 1];
} wei_t;

typedef struct wei_def_s {
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
} wei_def_t;

/*
 * Montgomery
 */

typedef void clamp_func(unsigned char *raw);

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
  fe_t i4;
  fe_t a24;
  fe_t a0;
  fe_t b0;
  mge_t g;
  clamp_func *clamp;
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
  clamp_func *clamp;
} mont_def_t;

/*
 * Edwards
 */

/* ege = edwards group element (affine) */
typedef struct ege_s {
  fe_t x;
  fe_t y;
} ege_t;

/* xge = extended group element */
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

/*
 * Helpers
 */

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))
#endif

static void
cleanse(void *ptr, size_t len) {
#if defined(_WIN32)
  /* https://github.com/jedisct1/libsodium/blob/3b26a5c/src/libsodium/sodium/utils.c#L112 */
  SecureZeroMemory(ptr, len);
#elif defined(__GNUC__)
  /* https://github.com/torvalds/linux/blob/37d4e84/include/linux/string.h#L233 */
  /* https://github.com/torvalds/linux/blob/37d4e84/include/linux/compiler-gcc.h#L21 */
  memset(ptr, 0, len);
  __asm__ __volatile__("": :"r"(ptr) :"memory");
#else
  /* http://www.daemonology.net/blog/2014-09-04-how-to-zero-a-buffer.html */
  static void *(*const volatile memset_ptr)(void *, int, size_t) = memset;
  (memset_ptr)(ptr, 0, len);
#endif
}

static uint32_t
less_than(const unsigned char *a,
          const unsigned char *b,
          int size,
          int endian) {
  /* Compare a < b in constant time. */
  int le = (endian == -1);
  int32_t eq = -1;
  int32_t gt = 0;
  int i = le ? size - 1 : 0;

  for (; le ? i >= 0 : i < size; le ? i-- : i++) {
    int32_t x = a[i];
    int32_t y = b[i];

    gt = (~eq & gt) | (eq & ((x - y) >> 31));
    eq = eq & (((x ^ y) - 1) >> 31);
  }

  return (uint32_t)(~eq & 1 & gt);
}

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

#if __has_builtin(__builtin_clz)
#define count_bits(x) (sizeof(unsigned int) * CHAR_BIT - __builtin_clz(x))
#else
static int
count_bits(unsigned int x) {
  int bits = 0;

  while (x != 0) {
    bits += 1;
    x >>= 1;
  }

  return bits;
}
#endif

/*
 * GMP Extras
 */

static mp_size_t
mpn_bitlen(const mp_limb_t *xp, mp_size_t n) {
  mp_size_t i, b;
  mp_limb_t w;

  for (i = n - 1; i >= 0; i--) {
    if (xp[i] != 0)
      break;
  }

  if (i < 0)
    return 0;

  w = xp[i];
  b = 0;

  while (w != 0) {
    w >>= 1;
    b += 1;
  }

  return i * GMP_NUMB_BITS + b;
}

static void
mpn_cleanse(mp_limb_t *p, mp_size_t n) {
  cleanse(p, n * sizeof(mp_limb_t));
}

static void
mpn_set_mpz(mp_limb_t *xp, mpz_srcptr x, mp_size_t n) {
  mp_size_t xn = mpz_size(x);

  assert(xn <= n);

  mpn_copyi(xp, mpz_limbs_read(x), xn);

  if (xn < n)
    mpn_zero(xp + xn, n - xn);
}

static void
mpz_set_mpn(mpz_t r, const mp_limb_t *xp, mp_size_t xn) {
  mpn_copyi(mpz_limbs_write(r, xn), xp, xn);
  mpz_limbs_finish(r, xn);
}

static void
cnd_swap(mp_limb_t cnd, mp_limb_t *ap, mp_limb_t *bp, mp_size_t n) {
  mp_limb_t mask = -(mp_limb_t)(cnd != 0);
  mp_size_t i;

  for (i = 0; i < n; i++) {
    mp_limb_t a = ap[i];
    mp_limb_t b = bp[i];
    mp_limb_t w = (a ^ b) & mask;

    ap[i] = a ^ w;
    bp[i] = b ^ w;
  }
}

static void
cnd_select(mp_limb_t cnd,
           mp_limb_t *rp,
           const mp_limb_t *ap,
           const mp_limb_t *bp,
           mp_size_t n) {
  mp_limb_t cond = (cnd != 0);
  mp_limb_t mask0 = cond - 1;
  mp_limb_t mask1 = ~mask0;
  mp_size_t i;

  for (i = 0; i < n; i++)
    rp[i] = (ap[i] & mask0) | (bp[i] & mask1);
}

static void
mpn_import_be(mp_limb_t *rp, mp_size_t rn,
              const unsigned char *xp, size_t xn) {
  size_t xi;
  mp_limb_t out;
  unsigned bits;

  for (xi = xn, out = bits = 0; xi > 0 && rn > 0;) {
    mp_limb_t in = xp[--xi];

    out |= (in << bits) & GMP_NUMB_MASK;
    bits += 8;

    if (bits >= GMP_NUMB_BITS) {
      *rp++ = out;
      rn--;

      bits -= GMP_NUMB_BITS;
      out = in >> (8 - bits);
    }
  }

  if (rn > 0) {
    *rp++ = out;
    if (--rn > 0)
      mpn_zero(rp, rn);
  }
}

static void
mpn_import_le(mp_limb_t *rp, mp_size_t rn,
              const unsigned char *xp, size_t xn) {
  size_t xi;
  mp_limb_t out;
  unsigned bits;

  for (xi = 0, out = bits = 0; xi < xn && rn > 0; ) {
    mp_limb_t in = xp[xi++];

    out |= (in << bits) & GMP_NUMB_MASK;
    bits += 8;

    if (bits >= GMP_NUMB_BITS) {
      *rp++ = out;
      rn--;

      bits -= GMP_NUMB_BITS;
      out = in >> (8 - bits);
    }
  }

  if (rn > 0) {
    *rp++ = out;
    if (--rn > 0)
      mpn_zero(rp, rn);
  }
}

static void
mpn_import(mp_limb_t *rp, mp_size_t rn,
           const unsigned char *xp, size_t xn, int endian) {
  if (endian == 1)
    mpn_import_be(rp, rn, xp, xn);
  else
    mpn_import_le(rp, rn, xp, xn);
}

static void
mpn_export_be(unsigned char *rp, size_t rn,
              const mp_limb_t *xp, mp_size_t xn) {
  unsigned bits;
  mp_limb_t in;

  for (bits = in = 0; xn > 0 && rn > 0;) {
    if (bits >= 8) {
      rp[--rn] = in;
      in >>= 8;
      bits -= 8;
    } else {
      unsigned char old = in;
      in = *xp++;
      xn--;
      rp[--rn] = old | (in << bits);
      in >>= (8 - bits);
      bits += GMP_NUMB_BITS - 8;
    }
  }

  while (rn > 0) {
    rp[--rn] = in;
    in >>= 8;
  }
}

static void
mpn_export_le(unsigned char *rp, size_t rn,
              const mp_limb_t *xp, mp_size_t xn) {
  unsigned bits;
  mp_limb_t in;

  for (bits = in = 0; xn > 0 && rn > 0;) {
    if (bits >= 8) {
      *rp++ = in;
      rn--;
      in >>= 8;
      bits -= 8;
    } else {
      unsigned char old = in;
      in = *xp++;
      xn--;
      *rp++ = old | (in << bits);
      rn--;
      in >>= (8 - bits);
      bits += GMP_NUMB_BITS - 8;
    }
  }

  while (rn > 0) {
    *rp++ = in;
    rn--;
    in >>= 8;
  }
}

static void
mpn_export(unsigned char *rp, size_t rn,
           const mp_limb_t *xp, mp_size_t xn, int endian) {
  if (endian == 1)
    mpn_export_be(rp, rn, xp, xn);
  else
    mpn_export_le(rp, rn, xp, xn);
}

static int
mpn_invert_n(mp_limb_t *rp,
             const mp_limb_t *xp,
             const mp_limb_t *yp,
             mp_size_t n) {
#ifdef BCRYPTO_HAS_GMP
#define MAX_EGCD_LIMBS ((521 + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS)
  mp_limb_t gp[MAX_EGCD_LIMBS + 1];
  mp_limb_t sp[MAX_EGCD_LIMBS + 1];
  mp_limb_t up[MAX_EGCD_LIMBS + 1];
  mp_limb_t vp[MAX_EGCD_LIMBS + 1];
  mp_size_t sn = n + 1;
  mp_size_t gn;

  assert(n <= MAX_EGCD_LIMBS);

  if (mpn_zero_p(xp, n)) {
    mpn_zero(rp, n);
    return 0;
  }

  mpn_copyi(up, xp, n);
  mpn_copyi(vp, yp, n);

  gn = mpn_gcdext(gp, sp, &sn, up, n, vp, n);

  assert(gn == 1);
  assert(gp[0] == 1);

  if (sn < 0) {
    mpn_sub(sp, yp, n, sp, -sn);
    sn = n;
  }

  assert(sn <= n);

  mpn_zero(rp + sn, n - sn);
  mpn_copyi(rp, sp, sn);

  return 1;
#undef MAX_EGCD_LIMBS
#else
  mpz_t rn, un, vn;

  if (mpn_zero_p(xp, n)) {
    mpn_zero(rp, n);
    return 0;
  }

  mpz_init(rn);
  mpz_roinit_n(un, xp, n);
  mpz_roinit_n(vn, yp, n);

  assert(mpz_invert(rn, un, vn));

  mpn_set_mpz(rp, rn, n);

  mpz_clear(rn);

  return 1;
#endif
}

static void
mpn_print(const mp_limb_t *p, mp_size_t n, int base) {
  mpz_t x;
  mpz_roinit_n(x, p, n);
  mpz_out_str(stdout, base, x);
}

#ifndef BCRYPTO_HAS_GMP
/* `mpz_jacobi` is not implemented in mini-gmp. */
/* https://github.com/golang/go/blob/aadaec5/src/math/big/int.go#L754 */
static int
mpz_jacobi(const mpz_t x, const mpz_t y) {
  mpz_t a, b, c;
  unsigned long s, bmod8;
  int j = 1;

  assert(mpz_sgn(x) >= 0);
  assert(mpz_sgn(y) > 0);
  assert(mpz_odd_p(y));

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  /* a = x */
  mpz_set(a, x);

  /* b = y */
  mpz_set(b, y);

  for (;;) {
    /* if b == 1 */
    if (mpz_cmp_ui(b, 1) == 0)
      break;

    /* if a == 0 */
    if (mpz_sgn(a) == 0) {
      j = 0;
      break;
    }

    /* a = a mod b */
    mpz_mod(a, a, b);

    /* if a == 0 */
    if (mpz_sgn(a) == 0) {
      j = 0;
      break;
    }

    /* s = a factors of 2 */
    s = mpz_scan1(a, 0);

    if (s & 1) {
      /* bmod8 = b mod 8 */
      bmod8 = mpz_getlimbn(b, 0) & 7;

      if (bmod8 == 3 || bmod8 == 5)
        j = -j;
    }

    /* c = a >> s */
    mpz_tdiv_q_2exp(c, a, s);

    /* if b mod 4 == 3 and c mod 4 == 3 */
    if ((mpz_getlimbn(b, 0) & 3) == 3 && (mpz_getlimbn(c, 0) & 3) == 3)
      j = -j;

    /* a = b */
    mpz_set(a, b);

    /* b = c */
    mpz_set(b, c);
  }

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);

  return j;
}
#endif

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
  mpn_import(r, sc->limbs, raw, sc->size, sc->endian);
  return less_than(raw, sc->raw, sc->size, sc->endian);
}

static void
sc_reduce(scalar_field_t *sc, sc_t r, const sc_t ap);

static int
sc_import_reduce(scalar_field_t *sc, sc_t r, const unsigned char *raw) {
  mp_limb_t tmp[MAX_SCALAR_LIMBS * 4];

  mpn_import(tmp, sc->limbs * 4, raw, sc->size, sc->endian);

  sc_reduce(sc, r, tmp);

  cleanse(tmp, sizeof(tmp));

  return less_than(raw, sc->raw, sc->size, sc->endian);
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
sc_set_word(scalar_field_t *sc, sc_t r, uint32_t word) {
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
  return mpn_invert_n(r, a, sc->n, sc->limbs);
}

static void
sc_pow(scalar_field_t *sc, sc_t r, const sc_t a, const sc_t e) {
  /* Used for inversion if not available otherwise. */
  mp_size_t i;
  sc_t b;

  sc_set_word(sc, b, 1);

  for (i = sc->bits - 1; i >= 0; i--) {
    mp_limb_t bit = e[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS);

    sc_sqr(sc, b, b);

    if (bit & 1)
      sc_mul(sc, b, b, a);
  }

  sc_set(sc, r, b);
}

static int
sc_invert(scalar_field_t *sc, sc_t r, const sc_t a) {
  mp_limb_t zero = sc_is_zero(sc, a);
  mp_size_t i;
  sc_t e;

  /* e = n - 2 */
  mpn_copyi(e, sc->n, sc->limbs);
  mpn_sub_1(e, e, sc->limbs, 2);

  sc_pow(sc, r, a, e);

  for (i = 0; i < sc->limbs; i++)
    r[i] &= -(zero ^ 1);

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
    /* Use a constant time barret reduction
     * to montgomerize the field element.
     *
     * We must be aligned to the limb size.
     * Normally we could just shift left by
     * the remainder after the import, but we
     * can only handle 2*max limbs in sc_reduce.
     *
     * This should only be an issue when
     *
     *   GMP_NUMB_BITS != FIELD_WORD_SIZE
     *
     * but this should never occur with a
     * proper build of the library.
     *
     * In particular, an assertion error
     * will be triggered with P224 when:
     *
     *   GMP_NUM_BITS == 64 && FIELD_WORD_SIZE == 32
     *
     * As the P224 `shift` is already aligned
     * to the field words, it does not align
     * the to the gmp limbs. Note:
     *
     *   224 mod 64 == 32
     *   224 mod 32 == 0
     */
    mp_limb_t xp[MAX_FIELD_LIMBS * 4];
    mp_size_t shift = fe->shift / GMP_NUMB_BITS;

    /* Check alignment. */
    assert((fe->shift % GMP_NUMB_BITS) == 0);

    /* We can only handle 2*max limbs. */
    assert(shift <= fe->limbs);

    /* x = (x << shift) mod p */
    mpn_zero(xp, fe->limbs * 4);
    mpn_import(xp + shift, fe->limbs, raw, fe->size, fe->endian);
    sc_reduce(&fe->sc, xp, xp);

#if GMP_NUMB_BITS == FIELD_WORD_SIZE
    /* Import directly. */
    assert(sizeof(mp_limb_t) == sizeof(fe_word_t));
    assert((size_t)fe->limbs == fe->words);
    memcpy(r, xp, fe->limbs * sizeof(mp_limb_t));
#else
    /* Export as little endian. */
    mpn_export_le(tmp, fe->size, xp, fe->limbs);
    fe->from_bytes(r, tmp);
#endif
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
  return sc_import_reduce(sc, r, tmp);
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
fe_sub(prime_field_t *fe, fe_t r, const fe_t a, const fe_t b) {
  fe->sub(r, a, b);

  if (fe->carry)
    fe->carry(r, r);
}

static void
fe_mulw(prime_field_t *fe, fe_t r, const fe_t a, unsigned int b) {
  int bits = count_bits(b);
  int i;

  if (b > 1 && (b & (b - 1)) == 0) {
    fe_add(fe, r, a, a);

    for (i = 1; i < bits - 1; i++)
      fe_add(fe, r, r, r);
  } else {
    fe_t c;

    fe_set(fe, c, a);
    fe_zero(fe, r);

    for (i = bits - 1; i >= 0; i--) {
      fe_add(fe, r, r, r);

      if ((b >> i) & 1)
        fe_add(fe, r, r, c);
    }
  }
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

static void
fe_mulm3(prime_field_t *fe, fe_t r, const fe_t a) {
  fe_t c;
  fe_add(fe, c, a, a);
  fe_add(fe, c, c, a);
  fe_neg(fe, r, c);
}

static void
fe_pow(prime_field_t *fe, fe_t r, const fe_t a, const mp_limb_t *e) {
  /* Used for inversion and legendre if not available otherwise. */
  mp_size_t i;
  fe_t b;

  fe_set(fe, b, fe->one);

  for (i = fe->bits - 1; i >= 0; i--) {
    mp_limb_t bit = e[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS);

    fe_sqr(fe, b, b);

    if (bit & 1)
      fe_mul(fe, b, b, a);
  }

  fe_set(fe, r, b);
}

static int
fe_invert_var(prime_field_t *fe, fe_t r, const fe_t a) {
  mp_limb_t rp[MAX_FIELD_LIMBS];
  int ret;

  fe_get_limbs(fe, rp, a);

  ret = mpn_invert_n(rp, rp, fe->p, fe->limbs);

  fe_set_limbs(fe, r, rp, fe->limbs);

  return ret;
}

static int
fe_invert(prime_field_t *fe, fe_t r, const fe_t a) {
  int zero = fe_is_zero(fe, a);
  int ret = zero ^ 1;

  if (fe->invert) {
    /* Fast inversion chain. */
    fe->invert(r, a);
  } else {
    /* Fermat's little theorem. */
    mp_limb_t e[MAX_FIELD_LIMBS];
    fe_word_t zero = fe_is_zero(fe, a);
    size_t i;

    /* e = p - 2 */
    mpn_copyi(e, fe->p, fe->limbs);
    mpn_sub_1(e, e, fe->limbs, 2);

    fe_pow(fe, r, a, e);

    for (i = 0; i < fe->words; i++)
      r[i] &= -(zero ^ 1);
  }

  return ret;
}

static int
fe_sqrt(prime_field_t *fe, fe_t r, const fe_t a) {
  int ret;

  if (fe->sqrt) {
    /* Fast square root chain. */
    ret = fe->sqrt(r, a);
  } else {
    /* Handle p = 3 mod 4 and p = 5 mod 8. */
    mp_limb_t e[MAX_FIELD_LIMBS + 1];
    fe_t b, b2;

    if ((fe->p[0] & 3) == 3) {
      /* b = a^((p + 1) / 4) mod p */
      mpn_copyi(e, fe->p, fe->limbs + 1);
      mpn_add_1(e, e, fe->limbs + 1, 1);
      mpn_rshift(e, e, fe->limbs + 1, 2);
      fe_pow(fe, b, a, e);
    } else if ((fe->p[0] & 7) == 5) {
      fe_t a2, c;

      /* a2 = a * 2 mod p */
      fe_mulw(fe, a2, a, 2);

      /* c = a2^((p - 5) / 8) mod p */
      mpn_copyi(e, fe->p, fe->limbs);
      mpn_sub_1(e, e, fe->limbs, 5);
      mpn_rshift(e, e, fe->limbs, 3);
      fe_pow(fe, c, a2, e);

      /* b = (c^2 * a2 - 1) * a * c mod p */
      fe_sqr(fe, b, c);
      fe_mul(fe, b, b, a2);
      fe_sub(fe, b, b, fe->one);
      fe_mul(fe, b, b, a);
      fe_mul(fe, b, b, c);
    } else {
      assert(0 && "no sqrt implementation");
    }

    /* b2 = b^2 mod p */
    fe_sqr(fe, b2, b);

    ret = fe_equal(fe, b2, a);

    fe_set(fe, r, b);
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

  if (fe->sqrt && fe->bits != 224) {
    /* Fast square root chain. */
    fe_t tmp;
    ret = fe->sqrt(tmp, a);
  } else {
    /* Euler's criterion. */
    mp_limb_t e[MAX_FIELD_LIMBS];
    int x, y, z;
    fe_t b;

    /* e = (p - 1) / 2 */
    mpn_copyi(e, fe->p, fe->limbs);
    mpn_sub_1(e, e, fe->limbs, 1);
    mpn_rshift(e, e, fe->limbs, 1);

    fe_pow(fe, b, a, e);

    x = fe_is_zero(fe, a);
    y = fe_equal(fe, b, fe->one);
    z = fe_equal(fe, b, fe->mone);

    assert(x + y + z == 1);

    ret = x | y;
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
#ifdef BCRYPTO_HAS_GMP
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
#else
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
  fe_set_word(fe, fe->four, 4);
  fe_neg(fe, fe->mone, fe->one);
}

/*
 * Short Weierstrass
 */

/*
 * Short Weierstrass Affine Point
 */

static void
wge_zero(wei_t *ec, wge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_zero(fe, r->x);
  fe_zero(fe, r->y);
  r->inf = 1;
}

static void
wge_cleanse(wei_t *ec, wge_t *r) {
  prime_field_t *fe = &ec->fe;
  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
  r->inf = 1;
}

static int
wge_validate(wei_t *ec, const wge_t *p) {
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
wge_set_x(wei_t *ec, wge_t *r, const fe_t x, int sign) {
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
wge_set_xy(wei_t *ec, wge_t *r, const fe_t x, const fe_t y) {
  prime_field_t *fe = &ec->fe;
  fe_set(fe, r->x, x);
  fe_set(fe, r->y, y);
  r->inf = 0;
}

static int
wge_import(wei_t *ec, wge_t *r, const unsigned char *raw, size_t len) {
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

      if (!wge_set_x(ec, r, r->x, form & 1))
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

      if (!wge_validate(ec, r))
        return 0;

      return 1;
    }
    default: {
      return 0;
    }
  }
}

static int
wge_export(wei_t *ec,
          unsigned char *raw,
          size_t *len,
          const wge_t *p,
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
wge_import_x(wei_t *ec, wge_t *r, const unsigned char *raw) {
  /* [SCHNORR] "Specification". */
  prime_field_t *fe = &ec->fe;

  if (!fe_import(fe, r->x, raw))
    return 0;

  return wge_set_x(ec, r, r->x, -1);
}

static int
wge_export_x(wei_t *ec, unsigned char *raw, const wge_t *p) {
  /* [SCHNORR] "Specification". */
  prime_field_t *fe = &ec->fe;

  if (p->inf)
    return 0;

  fe_export(fe, raw, p->x);

  return 1;
}

static void
wge_swap(wei_t *ec, wge_t *a, wge_t *b, unsigned int flag) {
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
wge_set(wei_t *ec, wge_t *r, const wge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  r->inf = a->inf;
}

static int
wge_equal(wei_t *ec, const wge_t *a, const wge_t *b) {
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
wge_is_zero(wei_t *ec, const wge_t *a) {
  return a->inf;
}

static void
wge_neg(wei_t *ec, wge_t *r, const wge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg(fe, r->y, a->y);
  r->inf = a->inf;
}

static void
wge_dbl(wei_t *ec, wge_t *r, const wge_t *a) {
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
    wge_zero(ec, r);
    return;
  }

  /* Y1 = 0 */
  if (fe_is_zero(fe, a->y)) {
    wge_zero(ec, r);
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
wge_add(wei_t *ec, wge_t *r, const wge_t *a, const wge_t *b) {
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
    wge_set(ec, r, b);
    return;
  }

  /* P + O = P */
  if (b->inf) {
    wge_set(ec, r, a);
    return;
  }

  /* P + P, P + -P */
  if (fe_equal(fe, a->x, b->x)) {
    /* P + -P = O */
    if (!fe_equal(fe, a->y, b->y)) {
      wge_zero(ec, r);
      return;
    }

    /* P + P = 2P */
    wge_dbl(ec, r, a);
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
wge_sub(wei_t *ec, wge_t *r, const wge_t *a, const wge_t *b) {
  wge_t c;
  wge_neg(ec, &c, b);
  wge_add(ec, r, a, &c);
}

static void
wge_to_jge(wei_t *ec, jge_t *r, const wge_t *a) {
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
wge_naf_points(wei_t *ec, wge_t *points, const wge_t *p, size_t width) {
  size_t size = (1 << width) - 1;
  wge_t dbl;
  size_t i;

  wge_dbl(ec, &dbl, p);
  wge_set(ec, &points[0], p);

  for (i = 1; i < size; i++)
    wge_add(ec, &points[i], &points[i - 1], &dbl);
}

static void
wge_print(wei_t *ec, const wge_t *p) {
  prime_field_t *fe = &ec->fe;

  if (wge_is_zero(ec, p)) {
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
 * Short Weierstrass Jacobian Point
 */

static void
jge_zero(wei_t *ec, jge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, fe->one);
  fe_set(fe, r->y, fe->one);
  fe_zero(fe, r->z);
}

static void
jge_cleanse(wei_t *ec, jge_t *r) {
  prime_field_t *fe = &ec->fe;

  fe_cleanse(fe, r->x);
  fe_cleanse(fe, r->y);
  fe_cleanse(fe, r->z);
}

static void
jge_swap(wei_t *ec, jge_t *a, jge_t *b, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_swap(fe, a->x, b->x, flag);
  fe_swap(fe, a->y, b->y, flag);
  fe_swap(fe, a->z, b->z, flag);
}

static void
jge_select(wei_t *ec,
           jge_t *r,
           const jge_t *a,
           const jge_t *b,
           unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, b->x, flag);
  fe_select(fe, r->y, a->y, b->y, flag);
  fe_select(fe, r->z, a->z, b->z, flag);
}

static void
jge_set(wei_t *ec, jge_t *r, const jge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_set(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
}

static int
jge_is_zero(wei_t *ec, const jge_t *a) {
  prime_field_t *fe = &ec->fe;

  return fe_is_zero(fe, a->z);
}

static int
jge_equal(wei_t *ec, const jge_t *a, const jge_t *b) {
  prime_field_t *fe = &ec->fe;
  fe_t z1, z2, e1, e2;
  int inf1 = jge_is_zero(ec, a);
  int inf2 = jge_is_zero(ec, b);
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
jge_equal_r(wei_t *ec, const jge_t *p, const sc_t x) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;
  mp_limb_t cp[MAX_FIELD_LIMBS + 1];
  mp_size_t cn = fe->limbs + 1;
  fe_t zz, rx, rn;

  assert(fe->limbs >= sc->limbs);

  if (jge_is_zero(ec, p))
    return 0;

  fe_sqr(fe, zz, p->z);

  fe_set_sc(fe, sc, rx, x);
  fe_mul(fe, rx, rx, zz);

  if (fe_equal(fe, p->x, rx))
    return 1;

  mpn_zero(cp + sc->limbs, cn - sc->limbs);
  mpn_copyi(cp, x, sc->limbs);

  fe_mul(fe, rn, ec->red_n, zz);

  assert(sc->n[cn - 1] == 0);
  assert(fe->p[cn - 1] == 0);

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
jge_neg(wei_t *ec, jge_t *r, const jge_t *a) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg(fe, r->y, a->y);
  fe_set(fe, r->z, a->z);
}

static void
jge_zero_cond(wei_t *ec, jge_t *r, const jge_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_select(fe, r->x, a->x, fe->one, flag);
  fe_select(fe, r->y, a->y, fe->one, flag);
  fe_select(fe, r->z, a->z, fe->zero, flag);
}

static void
jge_neg_cond(wei_t *ec, jge_t *r, const jge_t *a, unsigned int flag) {
  prime_field_t *fe = &ec->fe;

  fe_set(fe, r->x, a->x);
  fe_neg_cond(fe, r->y, a->y, flag);
  fe_set(fe, r->z, a->z);
}

static void
jge_to_wge(wei_t *ec, wge_t *r, const jge_t *p);

static void
jge_dbl(wei_t *ec, jge_t *r, const jge_t *a) {
  /* https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-1998-cmo-2
   * 3M + 6S + 4A + 1*a + 2*2 + 1*3 + 1*4 + 1*8
   * (implemented as: 3M + 6S + 5A + 1*a + 1*2 + 1*3 + 1*4 + 1*8)
   */
  prime_field_t *fe = &ec->fe;
  fe_t xx, yy, zz, s, m, t, t0, t1, x3, y3, z3;

  if (jge_is_zero(ec, a)) {
    jge_zero(ec, r);
    return;
  }

  if (fe_is_zero(fe, a->y)) {
    jge_zero(ec, r);
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
jge_add(wei_t *ec, jge_t *r, const jge_t *a, const jge_t *b) {
  /* No assumptions.
   * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-1998-cmo-2
   * 12M + 4S + 6A + 1*2 (implemented as: 12M + 4S + 7A)
   */
  prime_field_t *fe = &ec->fe;
  fe_t z1z1, z2z2, u1, u2, s1, s2, h, r0, hh, hhh, v, x3, y3, z3;

  if (jge_is_zero(ec, a)) {
    jge_set(ec, r, b);
    return;
  }

  if (jge_is_zero(ec, b)) {
    jge_set(ec, r, a);
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
      jge_zero(ec, r);
      return;
    }

    jge_dbl(ec, r, a);
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
jge_sub(wei_t *ec, jge_t *r, const jge_t *a, const jge_t *b) {
  jge_t c;
  jge_neg(ec, &c, b);
  jge_add(ec, r, a, &c);
}

static void
jge_mixed_add(wei_t *ec, jge_t *r, const jge_t *a, const wge_t *b) {
  /* Assumes Z2 = 1.
   * https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-madd
   * 8M + 3S + 6A + 5*2 (implemented as: 8M + 3S + 7A + 4*2)
   */
  prime_field_t *fe = &ec->fe;
  fe_t z1z1, u2, s2, h, r0, i, j, v, x3, y3, z3;

  if (jge_is_zero(ec, a)) {
    wge_to_jge(ec, r, b);
    return;
  }

  if (wge_is_zero(ec, b)) {
    jge_set(ec, r, a);
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
      jge_zero(ec, r);
      return;
    }

    jge_dbl(ec, r, a);
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
jge_mixed_sub(wei_t *ec, jge_t *r, const jge_t *a, const wge_t *b) {
  wge_t c;
  wge_neg(ec, &c, b);
  jge_mixed_add(ec, r, a, &c);
}

static void
jge_zaddu(wei_t *ec, jge_t *r, jge_t *p, const jge_t *a, const jge_t *b) {
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
jge_zaddc(wei_t *ec, jge_t *r, jge_t *s, const jge_t *a, const jge_t *b) {
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
jge_zdblu(wei_t *ec, jge_t *r, jge_t *p, const jge_t *a) {
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
jge_dblp(wei_t *ec, jge_t *r, const jge_t *a, size_t pow) {
  size_t i;

  jge_set(ec, r, a);

  for (i = 0; i < pow; i++)
    jge_dbl(ec, r, r);
}

static void
jge_to_wge(wei_t *ec, wge_t *r, const jge_t *p) {
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
jge_to_wge_var(wei_t *ec, wge_t *r, const jge_t *p) {
  /* https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#scaling-z
   * 1I + 3M + 1S
   */
  prime_field_t *fe = &ec->fe;
  fe_t a, aa;

  /* P = O */
  if (jge_is_zero(ec, p)) {
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
jge_validate(wei_t *ec, const jge_t *p) {
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
       | jge_is_zero(ec, p);
}

static void
jge_naf_points(wei_t *ec, jge_t *points, const wge_t *p, size_t width) {
  size_t size = (1 << width) - 1;
  jge_t dbl;
  size_t i;

  wge_to_jge(ec, &points[0], p);
  jge_dbl(ec, &dbl, &points[0]);

  for (i = 1; i < size; i++)
    jge_add(ec, &points[i], &points[i - 1], &dbl);
}

static void
jge_print(wei_t *ec, const jge_t *p) {
  prime_field_t *fe = &ec->fe;

  if (jge_is_zero(ec, p)) {
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
 * Short Weierstrass Curve
 */

static void
wei_init(wei_t *ec, const wei_def_t *def) {
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
  wge_zero(ec, &ec->unblind);

  wge_naf_points(ec, ec->points, &ec->g, NAF_WIDTH_PRE);
}

static void
wei_jmul_g_var(wei_t *ec, jge_t *r, const sc_t k) {
  /* Window NAF method for point multiplication.
   *
   * [GECC] Algorithm 3.36, Page 100, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  wge_t *points = ec->points;
  int32_t naf[MAX_SCALAR_BITS + 1];
  size_t max;
  int32_t i;
  sc_t k0;
  jge_t acc;

  /* Blind if available. */
  sc_add(sc, k0, k, ec->blind);

  /* Calculate max size. */
  max = sc_bitlen(sc, k0) + 1;

  /* Get NAF form. */
  sc_naf(sc, naf, k0, NAF_WIDTH_PRE, max);

  /* Add `this`*(N+1) for every w-NAF index. */
  jge_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    /* Count zeroes. */
    size_t k = 0;
    int32_t z;

    for (; i >= 0 && naf[i] == 0; i--)
      k += 1;

    if (i >= 0)
      k += 1;

    jge_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z = naf[i];

    assert(z != 0);

    if (z > 0)
      jge_mixed_add(ec, &acc, &acc, &points[(z - 1) >> 1]);
    else
      jge_mixed_sub(ec, &acc, &acc, &points[(-z - 1) >> 1]);
  }

  /* Unblind. */
  jge_mixed_add(ec, &acc, &acc, &ec->unblind);

  jge_set(ec, r, &acc);

  sc_cleanse(sc, k0);
}

static void
wei_jmul_var(wei_t *ec, jge_t *r, const wge_t *p, const sc_t k) {
  /* Window NAF method for point multiplication.
   *
   * [GECC] Algorithm 3.36, Page 100, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  wge_t points[(1 << NAF_WIDTH) - 1];
  int32_t naf[MAX_SCALAR_BITS + 1];
  size_t max = sc_bitlen(sc, k) + 1;
  int32_t i;
  jge_t acc;

  /* Precompute window. */
  wge_naf_points(ec, points, p, NAF_WIDTH);

  /* Get NAF form. */
  sc_naf(sc, naf, k, NAF_WIDTH, max);

  /* Add `this`*(N+1) for every w-NAF index. */
  jge_zero(ec, &acc);

  for (i = max - 1; i >= 0; i--) {
    /* Count zeroes. */
    size_t k = 0;
    int32_t z;

    for (; i >= 0 && naf[i] == 0; i--)
      k += 1;

    if (i >= 0)
      k += 1;

    jge_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z = naf[i];

    assert(z != 0);

    if (z > 0)
      jge_mixed_add(ec, &acc, &acc, &points[(z - 1) >> 1]);
    else
      jge_mixed_sub(ec, &acc, &acc, &points[(-z - 1) >> 1]);
  }

  jge_set(ec, r, &acc);
}

static void
wei_jmul_double_var(wei_t *ec,
                      jge_t *r,
                      const sc_t k1,
                      const wge_t *p2,
                      const sc_t k2) {
  /* Multiple point multiplication, also known
   * as "Shamir's trick" (with interleaved NAFs).
   *
   * [GECC] Algorithm 3.48, Page 109, Section 3.3.3.
   *        Algorithm 3.51, Page 112, Section 3.3.
   */
  scalar_field_t *sc = &ec->sc;
  wge_t *wnd1 = ec->points;
  wge_t wnd2[(1 << NAF_WIDTH) - 1];
  int32_t naf1[MAX_SCALAR_BITS + 1];
  int32_t naf2[MAX_SCALAR_BITS + 1];
  size_t max1 = sc_bitlen(sc, k1) + 1;
  size_t max2 = sc_bitlen(sc, k2) + 1;
  size_t max = max1 > max2 ? max1 : max2;
  int32_t i;
  jge_t acc;

  sc_naf(sc, naf1, k1, NAF_WIDTH_PRE, max);
  sc_naf(sc, naf2, k2, NAF_WIDTH, max);

  wge_naf_points(ec, wnd2, p2, NAF_WIDTH);

  /* Multiply and add. */
  jge_zero(ec, &acc);

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

    jge_dblp(ec, &acc, &acc, k);

    if (i < 0)
      break;

    z1 = naf1[i];
    z2 = naf2[i];

    if (z1 > 0)
      jge_mixed_add(ec, &acc, &acc, &wnd1[(z1 - 1) >> 1]);
    else if (z1 < 0)
      jge_mixed_sub(ec, &acc, &acc, &wnd1[(-z1 - 1) >> 1]);

    if (z2 > 0)
      jge_mixed_add(ec, &acc, &acc, &wnd2[(z2 - 1) >> 1]);
    else if (z2 < 0)
      jge_mixed_sub(ec, &acc, &acc, &wnd2[(-z2 - 1) >> 1]);
  }

  jge_set(ec, r, &acc);
}

static void
wei_jmul(wei_t *ec, jge_t *r, const wge_t *p, const sc_t k) {
  /* Co-Z Montgomery Ladder.
   *
   * [COZ] Algorithm 9, Page 6, Section 4.
   */
  scalar_field_t *sc = &ec->sc;
  jge_t a, b, c;
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
  wge_to_jge(ec, &c, p);
  jge_zdblu(ec, &a, &b, &c);

  for (i = bits - 2; i >= 0; i--) {
    mp_limb_t bit = (u[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS)) & 1;

    jge_swap(ec, &a, &b, swap ^ bit);
    jge_zaddc(ec, &a, &b, &b, &a);
    jge_zaddu(ec, &b, &a, &a, &b);

    swap = bit;
  }

  /* Finalize loop. */
  jge_swap(ec, &a, &b, swap);

  /* Handle edge case (k = 0). */
  jge_zero_cond(ec, &b, &b, zero);

  /* Handle edge case (k = -1). */
  jge_neg(ec, &c, &c);
  jge_swap(ec, &b, &c, minus1);

  /* Adjust sign. */
  jge_neg_cond(ec, &b, &b, negated);

  /* Result. */
  jge_set(ec, r, &b);

  /* Zero scalars. */
  sc_cleanse(sc, u);
  sc_cleanse(sc, v);
}

static void
wei_mul_g_var(wei_t *ec, wge_t *r, const sc_t k) {
  jge_t j;
  wei_jmul_g_var(ec, &j, k);
  jge_to_wge_var(ec, r, &j);
}

static void
wei_mul_var(wei_t *ec, wge_t *r, const wge_t *p, const sc_t k) {
  jge_t j;
  wei_jmul_var(ec, &j, p, k);
  jge_to_wge_var(ec, r, &j);
}

static void
wei_mul_double_var(wei_t *ec,
                     wge_t *r,
                     const sc_t k1,
                     const wge_t *p2,
                     const sc_t k2) {
  jge_t j;
  wei_jmul_double_var(ec, &j, k1, p2, k2);
  jge_to_wge_var(ec, r, &j);
}

static void
wei_jmul_g(wei_t *ec, jge_t *r, const sc_t k) {
  scalar_field_t *sc = &ec->sc;
  sc_t k0;

  /* Blind if available. */
  sc_add(sc, k0, k, ec->blind);

  /* Multiply in constant time. */
  wei_jmul(ec, r, &ec->g, k0);

  /* Unblind. */
  jge_mixed_add(ec, r, r, &ec->unblind);

  /* Cleanse. */
  sc_cleanse(sc, k0);
}

static void
wei_mul(wei_t *ec, wge_t *r, const wge_t *p, const sc_t k) {
  jge_t j;
  wei_jmul(ec, &j, p, k);
  jge_to_wge(ec, r, &j);
}

static void
wei_mul_g(wei_t *ec, wge_t *r, const sc_t k) {
  jge_t j;
  wei_jmul_g(ec, &j, k);
  jge_to_wge(ec, r, &j);
}

static void
wei_randomize(wei_t *ec, const unsigned char *entropy) {
  scalar_field_t *sc = &ec->sc;
  sc_t blind;
  wge_t unblind;

  sc_import_reduce(sc, blind, entropy);
  wei_mul_g(ec, &unblind, blind);
  wge_neg(ec, &unblind, &unblind);

  sc_set(sc, ec->blind, blind);
  wge_set(ec, &ec->unblind, &unblind);

  sc_cleanse(sc, blind);
  wge_cleanse(ec, &unblind);
}

/*
 * Montgomery
 */

/*
 * Montgomery Affine Point
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

static void
pge_import(mont_t *ec, pge_t *r, const unsigned char *raw);

static int
pge_to_mge(mont_t *ec, mge_t *r, const pge_t *p);

static int
mge_import(mont_t *ec, mge_t *r, const unsigned char *raw) {
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
 * Montgomery Projective Point
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

  /* B * y^2 * z = x^3 + A * x^2 * z + x * z^2 */
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

  /* sqrt(y^2 * z^4) = y * z^2 */
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
  fe_sub(fe, e, aa, bb);

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
 * Montgomery Curve
 */

static void
mont_init(mont_t *ec, const mont_def_t *def) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;

  ec->clamp = def->clamp;

  prime_field_init(fe, def->fe, -1);
  scalar_field_init(sc, def->sc, -1);

  fe_import_be(fe, ec->a, def->a);
  fe_import_be(fe, ec->b, def->b);
  fe_invert_var(fe, ec->bi, ec->b);
  fe_invert_var(fe, ec->i4, fe->four);

  /* a24 = (a + 2) / 4 */
  fe_add(fe, ec->a24, ec->a, fe->two);
  fe_mul(fe, ec->a24, ec->a24, ec->i4);

  /* a0 = a / b */
  fe_mul(fe, ec->a0, ec->a, ec->bi);

  /* b0 = 1 / b^2 */
  fe_sqr(fe, ec->b0, ec->bi);

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
  prime_field_t *fe = &ec->fe;
  pge_t a, b;
  int swap = 0;
  mp_size_t i;

  /* Clone points (for safe swapping). */
  pge_set(ec, &a, p);
  pge_zero(ec, &b);

  /* Climb the ladder. */
  for (i = fe->bits - 1; i >= 0; i--) {
    mp_limb_t bit = (k[i / GMP_NUMB_BITS] >> (i % GMP_NUMB_BITS)) & 1;

    /* Maybe swap. */
    pge_swap(ec, &a, &b, swap ^ bit);

    /* Single coordinate add+double. */
    pge_dad(ec, &a, &b, p, &a, &b);

    swap = bit;
  }

  /* Finalize loop. */
  pge_swap(ec, &a, &b, swap);

  pge_set(ec, r, &b);
}

static void
mont_mul_g(mont_t *ec, pge_t *r, const sc_t k) {
  pge_t g;
  mge_to_pge(ec, &g, &ec->g);
  mont_mul(ec, r, &g, k);
}

/*
 * Edwards
 */

static void
edwards_mul_a(edwards_t *ec, fe_t r, const fe_t x);

/*
 * Edwards Affine Point
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
 * Edwards Extended Point
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
 * Edwards Curve
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

  sc_import_reduce(sc, blind, entropy);
  edwards_jmul_g(ec, &unblind, blind);
  xge_neg(ec, &unblind, &unblind);

  sc_set(sc, ec->blind, blind);
  xge_set(ec, &ec->unblind, &unblind);

  sc_cleanse(sc, blind);
  xge_cleanse(ec, &unblind);
}

/*
 * ECDSA
 */

static int
ecdsa_reduce(wei_t *ec, sc_t r, const unsigned char *msg, size_t msg_len) {
  scalar_field_t *sc = &ec->sc;
  mp_limb_t *np = sc->n;
  mp_limb_t mp[MAX_SCALAR_LIMBS * 4];
  mp_size_t nn = sc->limbs;
  mp_size_t mn = nn * 4;
  long shift;
  int zero, cmp;

  if (msg_len > sc->size)
    msg_len = sc->size;

  mpn_import(mp, mn, msg, msg_len, sc->endian);

  shift = (long)msg_len * 8 - (long)sc->bits;

  if (shift > 0)
    mpn_rshift(mp, mp, mn, shift);

  zero = mpn_zero_p(mp, nn);
  cmp = mpn_cmp(mp, np, nn);

  if (cmp >= 0)
    sc_reduce(sc, r, mp);
  else
    sc_set(sc, r, mp);

  return !zero && cmp < 0;
}

static int
ecdsa_pubkey_create(wei_t *ec,
                    unsigned char *pub,
                    size_t *pub_len,
                    const unsigned char *priv,
                    int compact) {
  scalar_field_t *sc = &ec->sc;
  sc_t a;
  wge_t A;
  int ret = 0;

  if (!sc_import(sc, a, priv))
    goto fail;

  if (sc_is_zero(sc, a))
    goto fail;

  wei_mul_g(ec, &A, a);

  if (!wge_export(ec, pub, pub_len, &A, compact))
    goto fail;

  ret = 1;
fail:
  sc_cleanse(sc, a);
  wge_cleanse(ec, &A);
  return ret;
}

static int
ecdsa_sign(wei_t *ec,
           unsigned char *sig,
           unsigned int *param,
           const unsigned char *msg,
           size_t msg_len,
           const unsigned char *priv) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;
  drbg_t rng;
  sc_t a, m, k, r, s;
  wge_t R;
  unsigned char bytes[MAX_SCALAR_SIZE * 2];
  unsigned int hint;
  int ret = 0;

  if (param != NULL)
    *param = 0;

  if (!sc_import(sc, a, priv))
    goto fail;

  if (sc_is_zero(sc, a))
    goto fail;

  ecdsa_reduce(ec, m, msg, msg_len);

  memcpy(bytes, priv, sc->size);
  sc_export(sc, bytes + sc->size, m);

  drbg_init(&rng, ec->hash, bytes, sc->size * 2);

  for (;;) {
    drbg_generate(&rng, bytes, sc->size);

    if (!ecdsa_reduce(ec, k, bytes, sc->size))
      continue;

    wei_mul_g(ec, &R, k);

    if (wge_is_zero(ec, &R))
      continue;

    hint = fe_is_odd(fe, R.y);

    if (!fe_get_sc(fe, sc, r, R.x))
      hint |= 2;

    if (sc_is_zero(sc, r))
      continue;

    sc_invert(sc, k, k);
    sc_mul(sc, s, r, a);
    sc_add(sc, s, s, m);
    sc_mul(sc, s, s, k);

    if (sc_is_high_var(sc, s)) {
      sc_neg(sc, s, s);
      hint ^= 1;
    }

    sc_export(sc, sig, r);
    sc_export(sc, sig + sc->size, s);

    if (param != NULL)
      *param = hint;

    break;
  }

  ret = 1;
fail:
  cleanse(&rng, sizeof(drbg_t));
  sc_cleanse(sc, a);
  sc_cleanse(sc, m);
  sc_cleanse(sc, k);
  sc_cleanse(sc, r);
  sc_cleanse(sc, s);
  wge_cleanse(ec, &R);
  cleanse(bytes, sizeof(bytes));
  return ret;
}

static int
ecdsa_verify(wei_t *ec,
             const unsigned char *msg,
             size_t msg_len,
             const unsigned char *sig,
             const unsigned char *pub,
             size_t pub_len) {
  scalar_field_t *sc = &ec->sc;
  sc_t m, r, s, u1, u2;
  wge_t A;
#ifdef WITH_TRICK
  jge_t R;
#else
  prime_field_t *fe = &ec->fe;
  wge_t R;
  sc_t re;
#endif

  ecdsa_reduce(ec, m, msg, msg_len);

  if (!wge_import(ec, &A, pub, pub_len))
    return 0;

  if (!sc_import(sc, r, sig))
    return 0;

  if (!sc_import(sc, s, sig + sc->size))
    return 0;

  if (sc_is_zero(sc, r) || sc_is_zero(sc, s))
    return 0;

  sc_invert_var(sc, s, s);
  sc_mul(sc, u1, m, s);
  sc_mul(sc, u2, r, s);
#ifdef WITH_TRICK
  wei_jmul_double_var(ec, &R, u1, &A, u2);

  return jge_equal_r(ec, &R, r);
#else
  wei_mul_double_var(ec, &R, u1, &A, u2);

  fe_get_sc(fe, sc, re, R.x);

  return sc_equal(sc, r, re);
#endif
}

static int
ecdsa_recover(wei_t *ec,
              unsigned char *pub,
              size_t *pub_len,
              const unsigned char *msg,
              size_t msg_len,
              const unsigned char *sig,
              unsigned int param,
              int compact) {
  prime_field_t *fe = &ec->fe;
  scalar_field_t *sc = &ec->sc;
  unsigned int sign = param & 1;
  unsigned int high = param >> 1;
  sc_t m, r, s, s1, s2;
  fe_t x;
  wge_t R, A;

  ecdsa_reduce(ec, m, msg, msg_len);

  if (!sc_import(sc, r, sig))
    return 0;

  if (!sc_import(sc, s, sig + sc->size))
    return 0;

  if (sc_is_zero(sc, r) || sc_is_zero(sc, s))
    return 0;

  /* Assumes n < p. */
  fe_set_sc(fe, sc, x, r);

  if (high) {
    if (sc_cmp_var(sc, r, ec->pmodn) >= 0)
      return 0;

    fe_add(fe, x, x, ec->red_n);
  }

  if (!wge_set_x(ec, &R, x, sign))
    return 0;

  sc_invert_var(sc, r, r);
  sc_mul(sc, s1, m, r);
  sc_mul(sc, s2, s, r);
  sc_neg(sc, s1, s1);
  wei_mul_double_var(ec, &A, s1, &R, s2);

  return wge_export(ec, pub, pub_len, &A, compact);
}

static int
ecdsa_derive(wei_t *ec,
             unsigned char *secret,
             size_t *secret_len,
             const unsigned char *pub,
             const size_t pub_len,
             const unsigned char *priv,
             int compact) {
  scalar_field_t *sc = &ec->sc;
  sc_t a;
  wge_t A, P;
  int ret = 0;

  if (!sc_import(sc, a, priv))
    goto fail;

  if (sc_is_zero(sc, a))
    goto fail;

  if (!wge_import(ec, &A, pub, pub_len))
    goto fail;

  wei_mul(ec, &P, &A, a);

  if (!wge_export(ec, secret, secret_len, &P, compact))
    goto fail;

  ret = 1;
fail:
  sc_cleanse(sc, a);
  wge_cleanse(ec, &A);
  wge_cleanse(ec, &P);
  return ret;
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
  .invert = p224_fe_invert,
  .sqrt = p224_fe_sqrt_var,
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
  .sqrt = p384_fe_sqrt,
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
  .sqrt = p521_fe_sqrt,
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

/*
 * Short Weierstrass Curves
 */

static const wei_def_t curve_p192 = {
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

static const wei_def_t curve_p224 = {
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

static const wei_def_t curve_p256 = {
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

static const wei_def_t curve_p384 = {
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

static const wei_def_t curve_p521 = {
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

static const wei_def_t curve_secp256k1 = {
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

/*
 * Edwards Curves
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
