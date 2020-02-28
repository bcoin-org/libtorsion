/*!
 * scalar.h - scalar optimizations for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 *
 * Parts of this software are based on bitcoin-core/secp256k1:
 *   Copyright (c) 2013 Pieter Wuille
 *   https://github.com/bitcoin-core/secp256k1
 */

static void
sc_sqrn(const scalar_field_t *sc, sc_t r, const sc_t x, int rounds) {
  int i;

  sc_set(sc, r, x);

  for (i = 0; i < rounds; i++)
    sc_sqr(sc, r, r);
}

static void
secp256k1_sc_invert(const scalar_field_t *sc, sc_t r, const sc_t x) {
  /* https://github.com/bitcoin-core/secp256k1/blob/master/src/scalar_impl.h */
  sc_t x2, x3, x6, x8, x14, x28, x56, x112, x126;
  sc_t u2, u5, u9, u11, u13;
  sc_t t;

  sc_sqr(sc, u2, x);
  sc_mul(sc, x2, u2, x);
  sc_mul(sc, u5, u2,x2);
  sc_mul(sc, x3, u5, u2);
  sc_mul(sc, u9, x3, u2);
  sc_mul(sc, u11, u9, u2);
  sc_mul(sc, u13, u11, u2);

  sc_sqr(sc, x6, u13);
  sc_sqr(sc, x6, x6);
  sc_mul(sc, x6, x6, u11);

  sc_sqr(sc, x8, x6);
  sc_sqr(sc, x8, x8);
  sc_mul(sc, x8, x8,  x2);

  sc_sqr(sc, x14, x8);
  sc_sqrn(sc, x14, x14, 5);
  sc_mul(sc, x14, x14, x6);

  sc_sqr(sc, x28, x14);
  sc_sqrn(sc, x28, x28, 13);
  sc_mul(sc, x28, x28, x14);

  sc_sqr(sc, x56, x28);
  sc_sqrn(sc, x56, x56, 27);
  sc_mul(sc, x56, x56, x28);

  sc_sqr(sc, x112, x56);
  sc_sqrn(sc, x112, x112, 55);
  sc_mul(sc, x112, x112, x56);

  sc_sqr(sc, x126, x112);
  sc_sqrn(sc, x126, x126, 13);
  sc_mul(sc, x126, x126, x14);

  sc_set(sc, t, x126);
  sc_sqrn(sc, t, t, 3);
  sc_mul(sc, t, t, u5); /* 101 */
  sc_sqrn(sc, t, t, 4); /* 0 */
  sc_mul(sc, t, t, x3); /* 111 */
  sc_sqrn(sc, t, t, 4); /* 0 */
  sc_mul(sc, t, t, u5); /* 101 */
  sc_sqrn(sc, t, t, 5); /* 0 */
  sc_mul(sc, t, t, u11); /* 1011 */
  sc_sqrn(sc, t, t, 4);
  sc_mul(sc, t, t, u11); /* 1011 */
  sc_sqrn(sc, t, t, 4); /* 0 */
  sc_mul(sc, t, t, x3); /* 111 */
  sc_sqrn(sc, t, t, 5); /* 00 */
  sc_mul(sc, t, t, x3); /* 111 */
  sc_sqrn(sc, t, t, 6); /* 00 */
  sc_mul(sc, t, t, u13); /* 1101 */
  sc_sqrn(sc, t, t, 4); /* 0 */
  sc_mul(sc, t, t, u5); /* 101 */
  sc_sqrn(sc, t, t, 3);
  sc_mul(sc, t, t, x3); /* 111 */
  sc_sqrn(sc, t, t, 5); /* 0 */
  sc_mul(sc, t, t, u9); /* 1001 */
  sc_sqrn(sc, t, t, 6); /* 000 */
  sc_mul(sc, t, t, u5); /* 101 */
  sc_sqrn(sc, t, t, 10); /* 0000000 */
  sc_mul(sc, t, t, x3); /* 111 */
  sc_sqrn(sc, t, t, 4); /* 0 */
  sc_mul(sc, t, t, x3); /* 111 */
  sc_sqrn(sc, t, t, 9); /* 0 */
  sc_mul(sc, t, t, x8); /* 11111111 */
  sc_sqrn(sc, t, t, 5); /* 0 */
  sc_mul(sc, t, t, u9); /* 1001 */
  sc_sqrn(sc, t, t, 6); /* 00 */
  sc_mul(sc, t, t, u11); /* 1011 */
  sc_sqrn(sc, t, t, 4);
  sc_mul(sc, t, t, u13); /* 1101 */
  sc_sqrn(sc, t, t, 5);
  sc_mul(sc, t, t, x2); /* 11 */
  sc_sqrn(sc, t, t, 6); /* 00 */
  sc_mul(sc, t, t, u13); /* 1101 */
  sc_sqrn(sc, t, t, 10); /* 000000 */
  sc_mul(sc, t, t, u13); /* 1101 */
  sc_sqrn(sc, t, t, 4);
  sc_mul(sc, t, t, u9); /* 1001 */
  sc_sqrn(sc, t, t, 6); /* 00000 */
  sc_mul(sc, t, t, x); /* 1 */
  sc_sqrn(sc, t, t, 8); /* 00 */
  sc_mul(sc, r, t, x6); /* 111111 */
}
