#ifdef BCRYPTO_EC_64BIT
typedef uint64_t p224_fe_word_t;
#define P224_FIELD_WORDS 4
#include "p224_64.h"
#else
typedef uint32_t p224_fe_word_t;
#define P224_FIELD_WORDS 7
#include "p224_32.h"
#endif

typedef p224_fe_word_t p224_fe_t[P224_FIELD_WORDS];

#define p224_fe_add fiat_p224_add
#define p224_fe_sub fiat_p224_sub
#define p224_fe_neg fiat_p224_opp
#define p224_fe_mul fiat_p224_mul
#define p224_fe_sqr fiat_p224_square
#define p224_fe_nonzero fiat_p224_nonzero

static void
p224_fe_set(p224_fe_t out, const p224_fe_t in) {
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
  out[3] = in[3];
#if P224_FIELD_WORDS == 7
  out[4] = in[4];
  out[5] = in[5];
  out[6] = in[6];
#endif
}

static void
p224_fe_sqrn(p224_fe_t out, const p224_fe_t in, int rounds) {
  int i;

  p224_fe_set(out, in);

  for (i = 0; i < rounds; i++)
    p224_fe_sqr(out, out);
}

/* https://github.com/openssl/openssl/blob/master/crypto/ec/ecp_nistp224.c#L701 */
/* TODO: optimize */
static void
p224_fe_invert(p224_fe_t out, const p224_fe_t in) {
  p224_fe_t ftmp, ftmp2, ftmp3, ftmp4;
  p224_fe_t tmp;
  unsigned int i;

  p224_fe_sqr(tmp, in);
  p224_fe_set(ftmp, tmp); /* 2 */
  p224_fe_mul(tmp, in, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^2 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^3 - 2 */
  p224_fe_mul(tmp, in, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^3 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp2, tmp); /* 2^4 - 2 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp2, tmp); /* 2^5 - 4 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp2, tmp); /* 2^6 - 8 */
  p224_fe_mul(tmp, ftmp2, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^6 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp2, tmp); /* 2^7 - 2 */

  for (i = 0; i < 5; ++i) { /* 2^12 - 2^6 */
    p224_fe_sqr(tmp, ftmp2);
    p224_fe_set(ftmp2, tmp);
  }

  p224_fe_mul(tmp, ftmp2, ftmp);
  p224_fe_set(ftmp2, tmp); /* 2^12 - 1 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp3, tmp); /* 2^13 - 2 */

  for (i = 0; i < 11; ++i) { /* 2^24 - 2^12 */
    p224_fe_sqr(tmp, ftmp3);
    p224_fe_set(ftmp3, tmp);
  }

  p224_fe_mul(tmp, ftmp3, ftmp2);
  p224_fe_set(ftmp2, tmp); /* 2^24 - 1 */
  p224_fe_sqr(tmp, ftmp2);
  p224_fe_set(ftmp3, tmp); /* 2^25 - 2 */

  for (i = 0; i < 23; ++i) { /* 2^48 - 2^24 */
    p224_fe_sqr(tmp, ftmp3);
    p224_fe_set(ftmp3, tmp);
  }

  p224_fe_mul(tmp, ftmp3, ftmp2);
  p224_fe_set(ftmp3, tmp); /* 2^48 - 1 */
  p224_fe_sqr(tmp, ftmp3);
  p224_fe_set(ftmp4, tmp); /* 2^49 - 2 */

  for (i = 0; i < 47; ++i) { /* 2^96 - 2^48 */
    p224_fe_sqr(tmp, ftmp4);
    p224_fe_set(ftmp4, tmp);
  }

  p224_fe_mul(tmp, ftmp3, ftmp4);
  p224_fe_set(ftmp3, tmp); /* 2^96 - 1 */
  p224_fe_sqr(tmp, ftmp3);
  p224_fe_set(ftmp4, tmp); /* 2^97 - 2 */

  for (i = 0; i < 23; ++i) { /* 2^120 - 2^24 */
    p224_fe_sqr(tmp, ftmp4);
    p224_fe_set(ftmp4, tmp);
  }

  p224_fe_mul(tmp, ftmp2, ftmp4);
  p224_fe_set(ftmp2, tmp); /* 2^120 - 1 */

  for (i = 0; i < 6; ++i) { /* 2^126 - 2^6 */
    p224_fe_sqr(tmp, ftmp2);
    p224_fe_set(ftmp2, tmp);
  }

  p224_fe_mul(tmp, ftmp2, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^126 - 1 */
  p224_fe_sqr(tmp, ftmp);
  p224_fe_set(ftmp, tmp); /* 2^127 - 2 */
  p224_fe_mul(tmp, ftmp, in);
  p224_fe_set(ftmp, tmp); /* 2^127 - 1 */

  for (i = 0; i < 97; ++i) { /* 2^224 - 2^97 */
    p224_fe_sqr(tmp, ftmp);
    p224_fe_set(ftmp, tmp);
  }

  p224_fe_mul(tmp, ftmp, ftmp3);
  p224_fe_set(out, tmp); /* 2^224 - 2^96 - 1 */
}
