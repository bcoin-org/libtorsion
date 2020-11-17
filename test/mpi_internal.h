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

#undef ASSERT
#define ASSERT(expr) ASSERT_ALWAYS(expr)

/*
 * Debug Helpers
 */

TORSION_UNUSED static size_t
mpn_out_str(FILE *stream, int base, const mp_limb_t *xp, int xn) {
  mp_limb_t ch;
  int bytes = 0;
  int i;

  ASSERT(base == 16);

  if (xn < 0) {
    fputc('-', stream);
    xn = -xn;
  }

  xn = mpn_strip(xp, xn);

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

/*
 * MPN
 */

static void
test_mpn_sec_cmp(mp_rng_f *rng, void *arg) {
  printf("  - MPN comparison sanity check.\n");

  {
    static const mp_limb_t mod[4] = {4, 3, 2, 1};
    static const mp_limb_t minus1[4] = {3, 3, 2, 1};
    static const mp_limb_t plus1[4] = {5, 3, 2, 1};
    static const mp_limb_t full[4] = {0xff, 0xff, 0xff, 0xff};

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
    int i;

    for (i = 0; i < 1000; i++) {
      mp_limb_t a[4];
      mp_limb_t b[4];
      int cab, cba;

      rng(a, sizeof(a), arg);
      rng(b, sizeof(b), arg);

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
 * MPZ
 */

static void
test_mpz_primes(mp_rng_f *rng, void *arg) {
  /* TODO */
  (void)rng;
  (void)arg;
}

/*
 * Test
 */

void
test_mpi_internal(mp_rng_f *rng, void *arg) {
  printf("Testing internal MPI functions...\n");

  /* MPN */
  test_mpn_sec_cmp(rng, arg);

  /* MPZ */
  test_mpz_primes(rng, arg);
}
