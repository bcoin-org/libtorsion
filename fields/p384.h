#ifdef BCRYPTO_EC_64BIT
typedef uint64_t p384_fe_word_t;
#define P384_FIELD_WORDS 6
#include "p384_64.h"
#else
typedef uint32_t p384_fe_word_t;
#define P384_FIELD_WORDS 12
#include "p384_32.h"
#endif

typedef p384_fe_word_t p384_fe_t[P384_FIELD_WORDS];

#define p384_fe_add fiat_p384_add
#define p384_fe_sub fiat_p384_sub
#define p384_fe_neg fiat_p384_opp
#define p384_fe_mul fiat_p384_mul
#define p384_fe_sqr fiat_p384_square
#define p384_fe_nonzero fiat_p384_nonzero

static void
p384_fe_set(p384_fe_t out, const p384_fe_t in) {
  out[0] = in[0];
  out[1] = in[1];
  out[2] = in[2];
  out[3] = in[3];
  out[4] = in[4];
  out[5] = in[5];
#if P384_FIELD_WORDS == 12
  out[4] = in[4];
  out[5] = in[5];
  out[6] = in[6];
  out[7] = in[7];
  out[8] = in[8];
  out[9] = in[9];
  out[10] = in[10];
  out[11] = in[11];
#endif
}

static void
p384_fe_sqrn(p384_fe_t out, const p384_fe_t in, int rounds) {
  int i;

  p384_fe_set(out, in);

  for (i = 0; i < rounds; i++)
    p384_fe_sqr(out, out);
}

/* https://briansmith.org/ecc-inversion-addition-chains-01 */
/* https://github.com/randombit/botan/blob/master/src/lib/pubkey/ec_group/curve_gfp.cpp */
static void
p384_fe_invert(p384_fe_t out, const p384_fe_t in) {
  p384_fe_t r, x2, x3, x15, x30, rl;

  p384_fe_set(r, in);
  p384_fe_sqr(r, r);
  p384_fe_mul(r, r, in);
  p384_fe_set(x2, r);

  p384_fe_sqr(r, r);
  p384_fe_mul(r, r, in);

  p384_fe_set(x3, r);

  p384_fe_sqrn(r, r, 3);
  p384_fe_mul(r, r, x3);

  p384_fe_set(rl, r);
  p384_fe_sqrn(r, r, 6);
  p384_fe_mul(r, r, rl);

  p384_fe_sqrn(r, r, 3);
  p384_fe_mul(r, r, x3);

  p384_fe_set(x15, r);
  p384_fe_sqrn(r, r, 15);
  p384_fe_mul(r, r, x15);

  p384_fe_set(x30, r);
  p384_fe_sqrn(r, r, 30);
  p384_fe_mul(r, r, x30);

  p384_fe_set(rl, r);
  p384_fe_sqrn(r, r, 60);
  p384_fe_mul(r, r, rl);

  p384_fe_set(rl, r);
  p384_fe_sqrn(r, r, 120);
  p384_fe_mul(r, r, rl);

  p384_fe_sqrn(r, r, 15);
  p384_fe_mul(r, r, x15);

  p384_fe_sqrn(r, r, 31);
  p384_fe_mul(r, r, x30);

  p384_fe_sqrn(r, r, 2);
  p384_fe_mul(r, r, x2);

  p384_fe_sqrn(r, r, 94);
  p384_fe_mul(r, r, x30);

  p384_fe_sqrn(r, r, 2);

  p384_fe_mul(r, r, in);
  p384_fe_set(out, r);
}
