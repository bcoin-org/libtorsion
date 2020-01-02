#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>

static size_t
mpn_bitlen(const mp_limb_t *xp, mp_size_t n) {
  mpz_t x;

  mpz_roinit_n(x, xp, n);

  if (mpz_sgn(x) == 0)
    return 0;

  return mpz_sizeinbase(x, 2);
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
  mp_limb_t mask1 = ~mask1;
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

static void
mpn_invert_n(mp_limb_t *rp,
             const mp_limb_t *xp,
             const mp_limb_t *yp,
             mp_size_t n) {
#ifdef BCRYPTO_HAS_GMP
#ifdef BCRYPTO_EC_64BIT
#define MAX_EGCD_LIMBS 9
#else
#define MAX_EGCD_LIMBS 17
#endif
    mp_limb_t gp[MAX_EGCD_LIMBS + 1];
    mp_limb_t sp[MAX_EGCD_LIMBS + 1];
    mp_limb_t up[MAX_EGCD_LIMBS + 1];
    mp_limb_t vp[MAX_EGCD_LIMBS + 1];
    mp_size_t sn = n + 1;
    mp_size_t gn;

    assert(n <= MAX_EGCD_LIMBS);

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
#undef MAX_EGCD_LIMBS
#else
    mpz_t rn, un, vn;

    mpz_init(rn);
    mpz_roinit_n(un, xp, n);
    mpz_roinit_n(vn, yp, n);

    assert(mpz_invert(rn, un, vn));

    mpn_set_mpz(rp, rn, n);

    mpz_clear(rn);
#endif
}

static void
mpn_print(const mp_limb_t *p, mp_size_t n, int base) {
  mpz_t x;
  mpz_roinit_n(x, p, n);
  mpz_out_str(stdout, base, x);
}
