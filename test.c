#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "ec.c"

int RAND_status(void);
int RAND_poll();
int RAND_bytes(unsigned char *buf, int num);

static int
random_bytes(unsigned char *dst, size_t len) {
  memset(dst, 0x00, len);

  if (len > (size_t)INT_MAX)
    return 0;

  for (;;) {
    int status = RAND_status();

    assert(status >= 0);

    if (status != 0)
      break;

    if (RAND_poll() == 0)
      break;
  }

  return RAND_bytes(dst, (int)len) == 1;
}

int
main(void) {
  wei_t *ec = malloc(sizeof(wei_t));
  sc_t k, k1, k2;
  wge_t g, p, q, r, expect;
  jge_t jg, jp, jq, jr;
  unsigned char entropy[32];
  unsigned char p_raw[33], q_raw[33];
  size_t p_size, q_size;
  prime_field_t *fe;
  scalar_field_t *sc;

  assert(ec != NULL);

  wei_init(ec, &curve_p256);
  fe = &ec->fe;
  sc = &ec->sc;

  random_bytes(entropy, sizeof(entropy));

  wei_randomize(ec, entropy);

  wge_set(ec, &g, &ec->g);
  wge_to_jge(ec, &jg, &ec->g);

#if 0
  fe_print(fe, ec->a);
  fe_print(fe, ec->b);
  fe_print(fe, ec->two);

  {
    sc_t x, y, z;
    sc_set_word(sc, x, 2);
    sc_set_word(sc, y, 3);
    sc_print(sc, x);
    sc_print(sc, y);
    sc_mul(sc, x, x, y);
    sc_print(sc, x);
    sc_add(sc, x, x, y);
    sc_print(sc, x);
    sc_set_word(sc, x, 1);
    sc_set_word(sc, y, 1);
    sc_neg(sc, x, x);
    sc_neg(sc, y, y);
    sc_print(sc, x);
    sc_mul(sc, z, x, y);
    sc_print(sc, z);
    sc_set_word(sc, x, 4);
    sc_set_word(sc, y, 2);
    sc_invert_var(sc, y, y);
    sc_print(sc, y);
    sc_mul(sc, z, x, y);
    sc_print(sc, z);
    assert(!sc_is_zero(sc, z));

    sc_set_word(sc, x, 4);
    sc_set_word(sc, y, 2);
    sc_invert(sc, y, y);
    sc_print(sc, y);
    sc_mul(sc, z, x, y);
    sc_print(sc, z);
  }
#endif

  {
    mp_limb_t r[MAX_SCALAR_LIMBS];
    mp_limb_t t[MAX_SCALAR_LIMBS * 4];

    mpn_zero(r, ARRAY_SIZE(r));
    mpn_zero(t, ARRAY_SIZE(t));

    mpn_copyi(t, sc->n, sc->limbs);

    sc_reduce(sc, r, t);

    assert(sc_is_zero(sc, r));

    assert(!sc_import(sc, r, sc->raw));
  }

  {
    fe_t t;

    assert(!fe_import(fe, t, fe->raw));
  }


  {
    const unsigned char g_raw[33] = {
      0x03, 0x6b, 0x17, 0xd1, 0xf2, 0xe1, 0x2c, 0x42,
      0x47, 0xf8, 0xbc, 0xe6, 0xe5, 0x63, 0xa4, 0x40,
      0xf2, 0x77, 0x03, 0x7d, 0x81, 0x2d, 0xeb, 0x33,
      0xa0, 0xf4, 0xa1, 0x39, 0x45, 0xd8, 0x98, 0xc2,
      0x96
    };

    const unsigned char g2_raw[33] = {
      0x03, 0x7c, 0xf2, 0x7b, 0x18, 0x8d, 0x03, 0x4f,
      0x7e, 0x8a, 0x52, 0x38, 0x03, 0x04, 0xb5, 0x1a,
      0xc3, 0xc0, 0x89, 0x69, 0xe2, 0x77, 0xf2, 0x1b,
      0x35, 0xa6, 0x0b, 0x48, 0xfc, 0x47, 0x66, 0x99,
      0x78
    };

    const unsigned char g3_raw[33] = {
      0x02, 0x5e, 0xcb, 0xe4, 0xd1, 0xa6, 0x33, 0x0a,
      0x44, 0xc8, 0xf7, 0xef, 0x95, 0x1d, 0x4b, 0xf1,
      0x65, 0xe6, 0xc6, 0xb7, 0x21, 0xef, 0xad, 0xa9,
      0x85, 0xfb, 0x41, 0x66, 0x1b, 0xc6, 0xe7, 0xfd,
      0x6c
    };

    assert(wge_import(ec, &p, g_raw, 33));

    wge_to_jge(ec, &jp, &p);
    wge_to_jge(ec, &jq, &ec->g);

    assert(wge_validate(ec, &p));
    assert(jge_validate(ec, &jp));
    assert(jge_validate(ec, &jq));
    assert(wge_equal(ec, &p, &ec->g));
    assert(jge_equal(ec, &jp, &jq));

    assert(wge_import(ec, &q, g2_raw, 33));
    assert(wge_import(ec, &r, g3_raw, 33));

    wge_to_jge(ec, &jq, &q);
    wge_to_jge(ec, &jr, &r);

    wge_dbl(ec, &p, &ec->g);

    assert(wge_equal(ec, &p, &q));

    wge_add(ec, &p, &p, &ec->g);

    assert(wge_equal(ec, &p, &r));

    jge_dbl(ec, &jp, &jg);

    assert(jge_equal(ec, &jp, &jq));

    jge_add(ec, &jp, &jp, &jg);

    assert(jge_equal(ec, &jp, &jr));

    jge_sub(ec, &jp, &jp, &jg);

    assert(jge_equal(ec, &jp, &jq));

    jge_mixed_add(ec, &jp, &jp, &g);

    assert(jge_equal(ec, &jp, &jr));

    jge_mixed_sub(ec, &jp, &jp, &g);

    assert(jge_equal(ec, &jp, &jq));

    assert(jge_validate(ec, &jg));
    assert(jge_validate(ec, &jp));
    assert(jge_validate(ec, &jq));
    assert(jge_validate(ec, &jr));

    assert(!jge_is_zero(ec, &jg));
    assert(!jge_is_zero(ec, &jp));
    assert(!jge_is_zero(ec, &jq));
    assert(!jge_is_zero(ec, &jr));

    jge_to_wge(ec, &p, &jp);

    assert(wge_equal(ec, &p, &q));

    assert(wge_export(ec, p_raw, &p_size, &p, 1));
    assert(p_size == 33);

    assert(memcmp(p_raw, g2_raw, 33) == 0);
  }

  {
    const unsigned char k_raw[32] = {
      0x38, 0xf8, 0x62, 0x0b, 0xa6, 0x0b, 0xed, 0x7c,
      0xf9, 0x0c, 0x7a, 0x99, 0xac, 0x35, 0xa4, 0x4e,
      0x39, 0x27, 0x59, 0x8e, 0x3c, 0x99, 0xbb, 0xc5,
      0xf5, 0x70, 0x75, 0x13, 0xc4, 0x0e, 0x2c, 0xe3
    };

    const unsigned char expect_raw[33] = {
      0x02, 0x1a, 0xb3, 0x49, 0x34, 0xb8, 0x11, 0xb5,
      0x5e, 0x2f, 0xa4, 0xf1, 0xcd, 0x57, 0xf1, 0x68,
      0x51, 0x3d, 0x04, 0xb9, 0x45, 0xb0, 0x43, 0xec,
      0xe9, 0x6b, 0x25, 0x53, 0x96, 0x72, 0xff, 0x52,
      0x03
    };

    assert(sc_import(sc, k, k_raw));
    assert(wge_import(ec, &expect, expect_raw, 33));

    assert(wge_validate(ec, &expect));
    assert(wge_equal(ec, &expect, &expect));
    assert(!wge_equal(ec, &expect, &ec->g));

    wei_mul_g(ec, &q, k);

    assert(wge_equal(ec, &q, &expect));

    assert(wge_export(ec, q_raw, &q_size, &q, 1));
    assert(q_size == 33);

    assert(memcmp(q_raw, expect_raw, 33) == 0);

    wei_mul_g_var(ec, &q, k);

    assert(wge_equal(ec, &q, &expect));

    assert(wge_export(ec, q_raw, &q_size, &q, 1));
    assert(q_size == 33);

    assert(memcmp(q_raw, expect_raw, 33) == 0);
  }

  {
    const unsigned char p_raw[33] = {
      0x03, 0x42, 0x67, 0xab, 0xc7, 0xde, 0x72, 0x0f,
      0x14, 0x5a, 0xbc, 0x94, 0xb9, 0x5b, 0x33, 0x50,
      0x7a, 0x37, 0x55, 0x55, 0x2b, 0xef, 0xaf, 0x57,
      0x61, 0x33, 0x7a, 0xd6, 0x7a, 0x28, 0xa9, 0x08,
      0xa1
    };

    const unsigned char k_raw[32] = {
      0xfd, 0x37, 0xfe, 0xab, 0xd9, 0xdd, 0x8d, 0xe5,
      0xfd, 0x04, 0x79, 0xf4, 0xd6, 0xea, 0xd4, 0xe6,
      0x02, 0xc7, 0x06, 0x0f, 0x43, 0x6e, 0x2b, 0xf1,
      0xc0, 0x72, 0xe9, 0x91, 0x80, 0xcb, 0x09, 0x18
    };

    const unsigned char expect_raw[33] = {
      0x02, 0x93, 0xa3, 0x55, 0xe4, 0x8f, 0x3b, 0x74,
      0xcc, 0x3b, 0xcb, 0xb4, 0x6c, 0xb2, 0x84, 0x3a,
      0xd5, 0x4e, 0xe5, 0xe0, 0x45, 0xe9, 0x17, 0x0b,
      0x00, 0x45, 0xbc, 0xc2, 0x86, 0x68, 0x8c, 0x4d,
      0x56
    };

    assert(wge_import(ec, &p, p_raw, 33));
    assert(sc_import(sc, k, k_raw));
    assert(wge_import(ec, &expect, expect_raw, 33));

    assert(wge_validate(ec, &p));
    assert(wge_validate(ec, &expect));
    assert(wge_equal(ec, &expect, &expect));
    assert(!wge_equal(ec, &expect, &ec->g));

    wei_mul(ec, &q, &p, k);

    assert(wge_equal(ec, &q, &expect));

    assert(wge_export(ec, q_raw, &q_size, &q, 1));
    assert(q_size == 33);

    assert(memcmp(q_raw, expect_raw, 33) == 0);

    wei_mul_var(ec, &q, &p, k);

    assert(wge_equal(ec, &q, &expect));

    assert(wge_export(ec, q_raw, &q_size, &q, 1));
    assert(q_size == 33);

    assert(memcmp(q_raw, expect_raw, 33) == 0);
  }

  {
    const unsigned char p_raw[33] = {
      0x02, 0x65, 0x26, 0x45, 0xad, 0x1a, 0x36, 0x8c,
      0xdc, 0xcf, 0x81, 0x90, 0x56, 0x3b, 0x2a, 0x12,
      0xba, 0x31, 0xea, 0x33, 0x78, 0xc2, 0x23, 0x66,
      0xff, 0xf8, 0x47, 0x92, 0x63, 0x8c, 0xb8, 0xc8,
      0x94
    };

    const unsigned char k1_raw[32] = {
      0x5f, 0xd3, 0x7e, 0x3c, 0x67, 0x9e, 0xc5, 0xd0,
      0x2b, 0xb6, 0x6a, 0xa8, 0x6e, 0x56, 0xd6, 0x40,
      0x65, 0xe9, 0x47, 0x74, 0x4e, 0x50, 0xee, 0xec,
      0x80, 0xcf, 0xcc, 0xce, 0x3b, 0xd2, 0xf2, 0x1a
    };

    const unsigned char k2_raw[32] = {
      0xfb, 0x15, 0x9a, 0x7d, 0x37, 0x4d, 0x24, 0xde,
      0xde, 0x0a, 0x55, 0xb2, 0x98, 0x26, 0xe3, 0x24,
      0xf6, 0xf1, 0xd7, 0x57, 0x36, 0x53, 0xd7, 0x8a,
      0x98, 0xed, 0xa2, 0x80, 0x6d, 0xbe, 0x37, 0x98
    };

    const unsigned char expect_raw[33] = {
      0x02, 0x96, 0xf1, 0xb9, 0xe3, 0xe7, 0x0b, 0xa1,
      0x2e, 0xaf, 0x40, 0x23, 0x05, 0x64, 0x5b, 0x0f,
      0x28, 0x1b, 0xec, 0x25, 0x4f, 0xf2, 0x31, 0x8f,
      0x96, 0x9c, 0x97, 0x96, 0x0c, 0x35, 0x0b, 0x2c,
      0x6d
    };

    assert(wge_import(ec, &p, p_raw, 33));
    assert(sc_import(sc, k1, k1_raw));
    assert(sc_import(sc, k2, k2_raw));
    assert(wge_import(ec, &expect, expect_raw, 33));

    assert(wge_validate(ec, &p));
    assert(wge_validate(ec, &expect));
    assert(wge_equal(ec, &expect, &expect));
    assert(!wge_equal(ec, &expect, &ec->g));

    wei_mul_double_var(ec, &q, k1, &p, k2);

    assert(wge_equal(ec, &q, &expect));

    assert(wge_export(ec, q_raw, &q_size, &q, 1));
    assert(q_size == 33);

    assert(memcmp(q_raw, expect_raw, 33) == 0);
  }

  {
    const unsigned char priv[32] = {
      0x73, 0x23, 0x53, 0xb1, 0x1e, 0x72, 0x8d, 0xa3,
      0x44, 0x96, 0x8b, 0xa4, 0xeb, 0xb9, 0xad, 0xf5,
      0x69, 0x1f, 0x32, 0x5c, 0xc6, 0x76, 0x0d, 0x20,
      0x64, 0xb3, 0xbd, 0x08, 0x20, 0xf7, 0xf5, 0x5a
    };

    const unsigned char pub[33] = {
      0x02, 0x20, 0xb8, 0x61, 0x12, 0xd4, 0x30, 0x37,
      0xf0, 0xef, 0xdf, 0x0d, 0xd8, 0xc6, 0x4e, 0xbf,
      0x33, 0x4b, 0xed, 0x32, 0x32, 0x27, 0xd1, 0x33,
      0xf0, 0x01, 0xac, 0xa8, 0xd6, 0x83, 0x40, 0xdb,
      0x11
    };

    const unsigned char msg[32] = {
      0xc1, 0x4a, 0xd0, 0x9e, 0x00, 0x4a, 0x27, 0x78,
      0x98, 0x0b, 0x0b, 0x11, 0x60, 0x76, 0xec, 0x75,
      0x30, 0xf5, 0x61, 0x27, 0x39, 0x3e, 0x5d, 0xf1,
      0x14, 0x50, 0x61, 0xf7, 0xeb, 0x1a, 0x2c, 0xf7
    };

    const unsigned char sig[64] = {
      0x08, 0xbc, 0x4c, 0x1b, 0xe2, 0x2c, 0xc1, 0xc7,
      0xd9, 0xed, 0xf6, 0xd4, 0xc0, 0x40, 0x20, 0x5b,
      0x98, 0x95, 0xb9, 0x32, 0xac, 0xc1, 0x8c, 0xba,
      0x9c, 0x26, 0x3c, 0x24, 0x61, 0xa2, 0xaa, 0xe2,
      0x1c, 0xff, 0x3d, 0x0e, 0xe4, 0xdc, 0xed, 0x73,
      0xbb, 0x43, 0x82, 0x56, 0x80, 0x69, 0xe4, 0x6f,
      0x81, 0x05, 0x70, 0x04, 0xf9, 0xfd, 0x8e, 0x2c,
      0x89, 0x0f, 0x2e, 0x34, 0x42, 0x7a, 0x72, 0x7d
    };

    unsigned int param = 0;
    unsigned char pub2[33];
    size_t pub2_len;

    assert(ecdsa_pubkey_create(ec, pub2, &pub2_len, priv, 1));
    assert(pub2_len == 33);
    assert(memcmp(pub, pub2, 33) == 0);
    assert(ecdsa_verify(ec, msg, 32, sig, pub, 33));
    assert(ecdsa_recover(ec, pub2, &pub2_len, msg, 32, sig, param, 1));
    assert(pub2_len == 33);
    assert(memcmp(pub, pub2, 33) == 0);
  }

  {
    const unsigned char priv[32] = {
      0x5f, 0xd3, 0x7e, 0x3c, 0x67, 0x9e, 0xc5, 0xd0,
      0x2b, 0xb6, 0x6a, 0xa8, 0x6e, 0x56, 0xd6, 0x40,
      0x65, 0xe9, 0x47, 0x74, 0x4e, 0x50, 0xee, 0xec,
      0x80, 0xcf, 0xcc, 0xce, 0x3b, 0xd2, 0xf2, 0x1a
    };

    const unsigned char msg[32] = {
      0xfb, 0x15, 0x9a, 0x7d, 0x37, 0x4d, 0x24, 0xde,
      0xde, 0x0a, 0x55, 0xb2, 0x98, 0x26, 0xe3, 0x24,
      0xf6, 0xf1, 0xd7, 0x57, 0x36, 0x53, 0xd7, 0x8a,
      0x98, 0xed, 0xa2, 0x80, 0x6d, 0xbe, 0x37, 0x98
    };

    unsigned char sig[64];
    unsigned char pub[33];
    size_t pub_len;
    unsigned char pub2[33];
    size_t pub2_len;
    unsigned int param;

    assert(ecdsa_sign(ec, sig, &param, msg, 32, priv));
    assert(ecdsa_pubkey_create(ec, pub, &pub_len, priv, 1));
    assert(pub_len == 33);
    assert(ecdsa_verify(ec, msg, 32, sig, pub, 33));
    assert(ecdsa_recover(ec, pub2, &pub2_len, msg, 32, sig, param, 1));
    assert(pub2_len == 33);
    assert(memcmp(pub, pub2, 33) == 0);
  }

  {
    wei_t *ec = malloc(sizeof(wei_t));
    unsigned char entropy[32];
    unsigned char priv[32];
    unsigned char msg[32];
    unsigned char sig[64];
    unsigned char pub[33];
    size_t pub_len;
    unsigned char pub2[33];
    size_t pub2_len;
    unsigned int param;

    wei_init(ec, &curve_secp256k1);

    random_bytes(entropy, sizeof(entropy));
    random_bytes(priv, sizeof(priv));
    random_bytes(msg, sizeof(msg));

    priv[0] &= 0x7f;

    wei_randomize(ec, entropy);

    assert(ecdsa_sign(ec, sig, &param, msg, 32, priv));
    assert(ecdsa_pubkey_create(ec, pub, &pub_len, priv, 1));
    assert(pub_len == 33);
    assert(ecdsa_verify(ec, msg, 32, sig, pub, 33));
    assert(ecdsa_recover(ec, pub2, &pub2_len, msg, 32, sig, param, 1));
    assert(pub2_len == 33);
    assert(memcmp(pub, pub2, 33) == 0);

    free(ec);
  }

  {
    wei_t *ec = malloc(sizeof(wei_t));
    wge_t g, p, q, r;
    jge_t jg, jp, jq, jr;
    unsigned char entropy[66];
    unsigned char p_raw[67];
    size_t p_size;

    const unsigned char g_raw[67] = {
      0x02, 0x00, 0xc6, 0x85, 0x8e, 0x06, 0xb7, 0x04,
      0x04, 0xe9, 0xcd, 0x9e, 0x3e, 0xcb, 0x66, 0x23,
      0x95, 0xb4, 0x42, 0x9c, 0x64, 0x81, 0x39, 0x05,
      0x3f, 0xb5, 0x21, 0xf8, 0x28, 0xaf, 0x60, 0x6b,
      0x4d, 0x3d, 0xba, 0xa1, 0x4b, 0x5e, 0x77, 0xef,
      0xe7, 0x59, 0x28, 0xfe, 0x1d, 0xc1, 0x27, 0xa2,
      0xff, 0xa8, 0xde, 0x33, 0x48, 0xb3, 0xc1, 0x85,
      0x6a, 0x42, 0x9b, 0xf9, 0x7e, 0x7e, 0x31, 0xc2,
      0xe5, 0xbd, 0x66
    };

    const unsigned char g2_raw[67] = {
      0x02, 0x00, 0x43, 0x3c, 0x21, 0x90, 0x24, 0x27,
      0x7e, 0x7e, 0x68, 0x2f, 0xcb, 0x28, 0x81, 0x48,
      0xc2, 0x82, 0x74, 0x74, 0x03, 0x27, 0x9b, 0x1c,
      0xcc, 0x06, 0x35, 0x2c, 0x6e, 0x55, 0x05, 0xd7,
      0x69, 0xbe, 0x97, 0xb3, 0xb2, 0x04, 0xda, 0x6e,
      0xf5, 0x55, 0x07, 0xaa, 0x10, 0x4a, 0x3a, 0x35,
      0xc5, 0xaf, 0x41, 0xcf, 0x2f, 0xa3, 0x64, 0xd6,
      0x0f, 0xd9, 0x67, 0xf4, 0x3e, 0x39, 0x33, 0xba,
      0x6d, 0x78, 0x3d
    };

    const unsigned char g3_raw[67] = {
      0x03, 0x01, 0xa7, 0x3d, 0x35, 0x24, 0x43, 0xde,
      0x29, 0x19, 0x5d, 0xd9, 0x1d, 0x6a, 0x64, 0xb5,
      0x95, 0x94, 0x79, 0xb5, 0x2a, 0x6e, 0x5b, 0x12,
      0x3d, 0x9a, 0xb9, 0xe5, 0xad, 0x7a, 0x11, 0x2d,
      0x7a, 0x8d, 0xd1, 0xad, 0x3f, 0x16, 0x4a, 0x3a,
      0x48, 0x32, 0x05, 0x1d, 0xa6, 0xbd, 0x16, 0xb5,
      0x9f, 0xe2, 0x1b, 0xae, 0xb4, 0x90, 0x86, 0x2c,
      0x32, 0xea, 0x05, 0xa5, 0x91, 0x9d, 0x2e, 0xde,
      0x37, 0xad, 0x7d
    };

    wei_init(ec, &curve_p521);

    random_bytes(entropy, sizeof(entropy));

    wei_randomize(ec, entropy);

    wge_set(ec, &g, &ec->g);
    wge_to_jge(ec, &jg, &ec->g);

    assert(wge_import(ec, &p, g_raw, 67));

    wge_to_jge(ec, &jp, &p);
    wge_to_jge(ec, &jq, &ec->g);

    assert(wge_validate(ec, &p));
    assert(jge_validate(ec, &jp));
    assert(jge_validate(ec, &jq));
    assert(wge_equal(ec, &p, &ec->g));
    assert(jge_equal(ec, &jp, &jq));

    assert(wge_import(ec, &q, g2_raw, 67));
    assert(wge_import(ec, &r, g3_raw, 67));

    wge_to_jge(ec, &jq, &q);
    wge_to_jge(ec, &jr, &r);

    wge_dbl(ec, &p, &ec->g);

    assert(wge_equal(ec, &p, &q));

    wge_add(ec, &p, &p, &ec->g);

    assert(wge_equal(ec, &p, &r));

    jge_dbl(ec, &jp, &jg);

    assert(jge_equal(ec, &jp, &jq));

    jge_add(ec, &jp, &jp, &jg);

    assert(jge_equal(ec, &jp, &jr));

    jge_sub(ec, &jp, &jp, &jg);

    assert(jge_equal(ec, &jp, &jq));

    jge_mixed_add(ec, &jp, &jp, &g);

    assert(jge_equal(ec, &jp, &jr));

    jge_mixed_sub(ec, &jp, &jp, &g);

    assert(jge_equal(ec, &jp, &jq));

    assert(jge_validate(ec, &jg));
    assert(jge_validate(ec, &jp));
    assert(jge_validate(ec, &jq));
    assert(jge_validate(ec, &jr));

    assert(!jge_is_zero(ec, &jg));
    assert(!jge_is_zero(ec, &jp));
    assert(!jge_is_zero(ec, &jq));
    assert(!jge_is_zero(ec, &jr));

    jge_to_wge(ec, &p, &jp);

    assert(wge_equal(ec, &p, &q));

    assert(wge_export(ec, p_raw, &p_size, &p, 1));
    assert(p_size == 67);

    assert(memcmp(p_raw, g2_raw, 67) == 0);

    free(ec);
  }

  {
    wei_t *ec = malloc(sizeof(wei_t));
    unsigned char entropy[66];
    unsigned char priv[66];
    unsigned char msg[66];
    unsigned char sig[66 * 2];
    unsigned char pub[66 + 1];
    size_t pub_len;
    unsigned char pub2[66 + 1];
    size_t pub2_len;
    unsigned int param;

    wei_init(ec, &curve_p521);

    random_bytes(entropy, sizeof(entropy));
    random_bytes(priv, sizeof(priv));
    random_bytes(msg, sizeof(msg));

    priv[0] = 0;

    wei_randomize(ec, entropy);

    assert(ecdsa_sign(ec, sig, &param, msg, 66, priv));
    assert(ecdsa_pubkey_create(ec, pub, &pub_len, priv, 1));
    assert(pub_len == 66 + 1);
    assert(ecdsa_verify(ec, msg, 66, sig, pub, 66 + 1));
    assert(ecdsa_recover(ec, pub2, &pub2_len, msg, 66, sig, param, 1));
    assert(pub2_len == 66 + 1);
    assert(memcmp(pub, pub2, 66 + 1) == 0);

    free(ec);
  }

  {
    edwards_t *ec = malloc(sizeof(edwards_t));
    ege_t g, p, q, r;
    xge_t jg, jp, jq, jr;
    unsigned char entropy[32];
    unsigned char p_raw[32];

    const unsigned char g_raw[32] = {
      0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
      0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
      0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
      0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
    };

    const unsigned char g2_raw[32] = {
      0xc9, 0xa3, 0xf8, 0x6a, 0xae, 0x46, 0x5f, 0x0e,
      0x56, 0x51, 0x38, 0x64, 0x51, 0x0f, 0x39, 0x97,
      0x56, 0x1f, 0xa2, 0xc9, 0xe8, 0x5e, 0xa2, 0x1d,
      0xc2, 0x29, 0x23, 0x09, 0xf3, 0xcd, 0x60, 0x22
    };

    const unsigned char g3_raw[32] = {
      0xd4, 0xb4, 0xf5, 0x78, 0x48, 0x68, 0xc3, 0x02,
      0x04, 0x03, 0x24, 0x67, 0x17, 0xec, 0x16, 0x9f,
      0xf7, 0x9e, 0x26, 0x60, 0x8e, 0xa1, 0x26, 0xa1,
      0xab, 0x69, 0xee, 0x77, 0xd1, 0xb1, 0x67, 0x12
    };

    edwards_init(ec, &curve_ed25519);

    random_bytes(entropy, sizeof(entropy));

    edwards_randomize(ec, entropy);

    ege_set(ec, &g, &ec->g);
    ege_to_xge(ec, &jg, &ec->g);

    assert(ege_import(ec, &p, g_raw));

    ege_to_xge(ec, &jp, &p);
    ege_to_xge(ec, &jq, &ec->g);

    assert(ege_validate(ec, &p));
    assert(xge_validate(ec, &jp));
    assert(xge_validate(ec, &jq));
    assert(ege_equal(ec, &p, &ec->g));
    assert(xge_equal(ec, &jp, &jq));

    assert(ege_import(ec, &q, g2_raw));
    assert(ege_import(ec, &r, g3_raw));

    ege_to_xge(ec, &jq, &q);
    ege_to_xge(ec, &jr, &r);

    ege_dbl(ec, &p, &ec->g);

    assert(ege_equal(ec, &p, &q));

    ege_add(ec, &p, &p, &ec->g);

    assert(ege_equal(ec, &p, &r));

    xge_dbl(ec, &jp, &jg);

    assert(xge_equal(ec, &jp, &jq));

    xge_add(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jr));

    xge_sub(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jq));

    xge_add(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jr));

    xge_sub(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jq));

    assert(xge_validate(ec, &jg));
    assert(xge_validate(ec, &jp));
    assert(xge_validate(ec, &jq));
    assert(xge_validate(ec, &jr));

    assert(!xge_is_zero(ec, &jg));
    assert(!xge_is_zero(ec, &jp));
    assert(!xge_is_zero(ec, &jq));
    assert(!xge_is_zero(ec, &jr));

    xge_to_ege(ec, &p, &jp);

    assert(ege_equal(ec, &p, &q));

    ege_export(ec, p_raw, &p);
    assert(memcmp(p_raw, g2_raw, 32) == 0);

    free(ec);
  }

  free(ec);

  return 0;
}

