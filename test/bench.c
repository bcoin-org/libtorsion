/*!
 * bench.c - benchmarks for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <torsion/drbg.h>
#include <torsion/ecc.h>
#include <torsion/hash.h>
#include <torsion/rsa.h>
#include <torsion/stream.h>
#include <torsion/util.h>

#include "utils.h"

/*
 * Helpers
 */

static void
drbg_init_rand(drbg_t *rng) {
  unsigned char entropy[ENTROPY_SIZE];
  size_t i;

  for (i = 0; i < ENTROPY_SIZE; i++)
    entropy[i] = rand();

  drbg_init(rng, HASH_SHA256, entropy, ENTROPY_SIZE);
}

/*
 * Benchmarks
 */

typedef uint64_t bench_t;

static void
bench_start(bench_t *start, const char *name) {
  printf("Benchmarking %s...\n", name);
  *start = torsion_hrtime();
}

static void
bench_puts64(const char *pre, uint64_t x) {
  /* Implement our own PRIu64. */
  char str[20 + 1];
  int len = 0;

  do {
    str[len++] = (x % 10) + '0';
    x /= 10;
  } while (x != 0);

  {
    int i = 0;
    int j = len - 1;
    int k = len >> 1;
    int ch;

    while (k--) {
      ch = str[i];
      str[i++] = str[j];
      str[j--] = ch;
    }
  }

  str[len] = '\0';

  printf("%s%s\n", pre, str);
}

static void
bench_end(bench_t *start, uint64_t ops) {
  bench_t nsec = torsion_hrtime() - *start;
  double sec = (double)nsec / 1000000000.0;

  bench_puts64("  Operations: ", ops);
  bench_puts64("  Nanoseconds: ", nsec);

  printf("  Seconds:     %f\n", sec);
  printf("  Ops/Sec:     %f\n", (double)ops / sec);
}

static void
bench_ecdsa_pubkey_create(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char pub[33];
  size_t pub_len;
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));

  ecdsa_privkey_generate(ec, priv, entropy);

  bench_start(&tv, "ecdsa_pubkey_create");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_pubkey_create(ec, pub, &pub_len, priv, 1));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdsa_pubkey_tweak_add(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char entropy1[ENTROPY_SIZE];
  unsigned char entropy2[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char tweak[32];
  unsigned char pub[33];
  unsigned char out[33];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy1, sizeof(entropy1));
  drbg_generate(rng, entropy2, sizeof(entropy2));

  ecdsa_privkey_generate(ec, priv, entropy1);
  ecdsa_privkey_generate(ec, tweak, entropy2);

  ASSERT(ecdsa_pubkey_create(ec, pub, NULL, priv, 1));

  bench_start(&tv, "ecdsa_pubkey_tweak_add");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_pubkey_tweak_add(ec, out, NULL, pub, 33, tweak, 1));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdsa_sign(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));
  drbg_generate(rng, msg, sizeof(msg));

  ecdsa_privkey_generate(ec, priv, entropy);

  bench_start(&tv, "ecdsa_sign");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_sign(ec, sig, NULL, msg, 32, priv));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdsa_verify(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char pub[33];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));
  drbg_generate(rng, msg, sizeof(msg));

  ecdsa_privkey_generate(ec, priv, entropy);

  ASSERT(ecdsa_sign(ec, sig, NULL, msg, 32, priv));
  ASSERT(ecdsa_pubkey_create(ec, pub, NULL, priv, 1));

  bench_start(&tv, "ecdsa_verify");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_verify(ec, msg, 32, sig, pub, 33));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdsa_derive(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char entropy1[ENTROPY_SIZE];
  unsigned char entropy2[ENTROPY_SIZE];
  unsigned char k0[32];
  unsigned char k1[32];
  unsigned char p0[65];
  unsigned char p1[33];
  size_t p0_len, p1_len;
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy1, sizeof(entropy1));
  drbg_generate(rng, entropy2, sizeof(entropy2));

  ecdsa_privkey_generate(ec, k0, entropy1);
  ecdsa_privkey_generate(ec, k1, entropy2);

  ASSERT(ecdsa_pubkey_create(ec, p0, &p0_len, k0, 0));

  bench_start(&tv, "ecdsa_derive");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_derive(ec, p1, &p1_len, p0, p0_len, k1, 1));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdh_derive(drbg_t *rng) {
  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char pub[32];
  unsigned char secret[32];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));

  ecdh_privkey_generate(ec, priv, entropy);
  ecdh_pubkey_create(ec, pub, priv);

  bench_start(&tv, "ecdh_derive");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdh_derive(ec, secret, pub, priv));

  bench_end(&tv, i);

  mont_curve_destroy(ec);
}

static void
bench_eddsa_sign(drbg_t *rng) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));
  drbg_generate(rng, msg, sizeof(msg));

  eddsa_privkey_generate(ec, priv, entropy);

  bench_start(&tv, "eddsa_sign");

  for (i = 0; i < 10000; i++)
    eddsa_sign(ec, sig, msg, 32, priv, -1, NULL, 0);

  bench_end(&tv, i);

  edwards_curve_destroy(ec);
}

static void
bench_eddsa_verify(drbg_t *rng) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char pub[32];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));
  drbg_generate(rng, msg, sizeof(msg));

  eddsa_privkey_generate(ec, priv, entropy);
  eddsa_sign(ec, sig, msg, 32, priv, -1, NULL, 0);
  eddsa_pubkey_create(ec, pub, priv);

  bench_start(&tv, "eddsa_verify");

  for (i = 0; i < 10000; i++)
    ASSERT(eddsa_verify(ec, msg, 32, sig, pub, -1, NULL, 0));

  bench_end(&tv, i);

  edwards_curve_destroy(ec);
}

static void
bench_eddsa_derive(drbg_t *rng) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char entropy1[ENTROPY_SIZE];
  unsigned char entropy2[ENTROPY_SIZE];
  unsigned char k0[32];
  unsigned char k1[32];
  unsigned char p0[32];
  unsigned char p1[32];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy1, sizeof(entropy1));
  drbg_generate(rng, entropy2, sizeof(entropy2));

  eddsa_privkey_generate(ec, k0, entropy1);
  eddsa_privkey_generate(ec, k1, entropy2);

  eddsa_pubkey_create(ec, p0, k0);

  bench_start(&tv, "eddsa_derive");

  for (i = 0; i < 10000; i++)
    ASSERT(eddsa_derive(ec, p1, p0, k1));

  bench_end(&tv, i);

  edwards_curve_destroy(ec);
}

typedef void mp_start_f(bench_t *start, const char *name);
typedef void mp_end_f(bench_t *start, uint64_t ops);
typedef void mp_rng_f(void *out, size_t size, void *arg);

void
torsion__bench_mpi_internal(mp_start_f *start, mp_end_f *end,
                            mp_rng_f *rng, void *arg);

static void
bench_mpi_internal(drbg_t *rng) {
  torsion__bench_mpi_internal(bench_start, bench_end, drbg_rng, rng);
}

static void
bench_rsa_generate(drbg_t *rng) {
  static const unsigned char entropy[ENTROPY_SIZE] = {0};
  unsigned char priv[RSA_MAX_PRIV_SIZE];
  size_t priv_len;
  bench_t tv;
  int i;

  (void)rng;

  bench_start(&tv, "rsa_generate");

  for (i = 0; i < 100; i++)
    ASSERT(rsa_privkey_generate(priv, &priv_len, 1024, 65537, entropy));

  bench_end(&tv, i);
}

static void
bench_rsa_sign(drbg_t *rng) {
  static const unsigned char priv[] = {
    0x30, 0x82, 0x04, 0xa5, 0x02, 0x01, 0x00, 0x02, 0x82, 0x01, 0x01, 0x00,
    0xc4, 0xc3, 0x93, 0x49, 0x7e, 0xca, 0xdd, 0x3d, 0xac, 0x9c, 0x46, 0xe4,
    0xab, 0x26, 0x04, 0xcb, 0xfe, 0x46, 0x79, 0x57, 0xa8, 0xd5, 0x6c, 0xee,
    0x4d, 0x6b, 0xe2, 0x6f, 0xd9, 0xe6, 0x23, 0x2c, 0x61, 0x3d, 0x09, 0x9d,
    0xc4, 0xa0, 0xb9, 0x5b, 0xd4, 0x37, 0x9e, 0x46, 0x3b, 0xbe, 0x62, 0xfa,
    0x26, 0xf2, 0x3c, 0x18, 0x73, 0x78, 0x2f, 0x7b, 0x21, 0xd0, 0x00, 0x9d,
    0xf9, 0x50, 0xb8, 0xe7, 0xb1, 0x06, 0x02, 0xc8, 0x0e, 0xd0, 0xda, 0x55,
    0x5c, 0x8a, 0x4b, 0x1c, 0x59, 0x49, 0x0e, 0x75, 0x03, 0x27, 0xe9, 0x64,
    0x45, 0x4d, 0x13, 0xfc, 0x6a, 0x64, 0xcf, 0xbc, 0xc7, 0x3f, 0x83, 0xc6,
    0x33, 0x0b, 0x3a, 0x77, 0xd7, 0xf7, 0x73, 0xaf, 0x36, 0x3a, 0xf3, 0xe7,
    0xcb, 0x09, 0x31, 0xc0, 0xbd, 0xcc, 0x09, 0x30, 0x07, 0x0a, 0x49, 0x39,
    0xe7, 0x54, 0x96, 0xe7, 0xc4, 0xd3, 0x8a, 0x30, 0xf4, 0xe3, 0x16, 0xd1,
    0x6d, 0x62, 0x3e, 0xc6, 0x54, 0xd6, 0xf7, 0x95, 0x63, 0x51, 0x94, 0xbb,
    0x29, 0x98, 0x54, 0x96, 0x1d, 0x89, 0x78, 0x47, 0x97, 0x40, 0xb6, 0x47,
    0xde, 0xda, 0xcb, 0x91, 0x9b, 0x22, 0xf9, 0x87, 0x0f, 0x2a, 0xef, 0xf0,
    0xbb, 0x69, 0x59, 0xf7, 0xd7, 0xa0, 0xb1, 0x06, 0xed, 0xa3, 0x1b, 0xbe,
    0x4b, 0x18, 0x70, 0xbb, 0xe1, 0x79, 0x72, 0xa7, 0x1d, 0x04, 0x9d, 0x06,
    0x34, 0x50, 0x51, 0xc2, 0xb1, 0x0f, 0x22, 0x67, 0x5e, 0xc8, 0x00, 0x70,
    0x71, 0x8b, 0x32, 0xc2, 0x61, 0x6f, 0xc1, 0x93, 0x75, 0x09, 0x75, 0x2c,
    0xfe, 0x21, 0x63, 0xc0, 0x49, 0x9c, 0x10, 0x76, 0x2b, 0x25, 0x51, 0x73,
    0x81, 0x1a, 0xdb, 0x85, 0x87, 0x7a, 0x43, 0xc5, 0x4a, 0xfb, 0x95, 0xa9,
    0x4e, 0x6f, 0xdd, 0x42, 0x8e, 0x05, 0x43, 0xd8, 0x68, 0x7d, 0x3a, 0x3d,
    0x8b, 0x0f, 0xc3, 0xbb, 0x02, 0x03, 0x01, 0x00, 0x01, 0x02, 0x82, 0x01,
    0x01, 0x00, 0xb7, 0xa6, 0x57, 0x65, 0xa4, 0xab, 0x51, 0xfe, 0x4b, 0x8a,
    0x7d, 0x7c, 0xd6, 0xe5, 0xa0, 0x5a, 0x8a, 0x15, 0x5b, 0x12, 0x5f, 0x69,
    0xfc, 0xc7, 0x1b, 0x8a, 0x13, 0x8a, 0x14, 0x56, 0x02, 0x04, 0x5d, 0x29,
    0xec, 0x3c, 0xce, 0x16, 0xb9, 0x8b, 0x25, 0x33, 0x58, 0x4d, 0xf3, 0x5b,
    0x4a, 0xe4, 0x72, 0xcf, 0x6b, 0x19, 0xe3, 0x44, 0x8d, 0x04, 0x9f, 0x55,
    0x96, 0x0d, 0xdc, 0x72, 0xe4, 0x72, 0x94, 0x3e, 0xa8, 0xff, 0xf3, 0x1c,
    0x2a, 0x7c, 0xbb, 0xe7, 0xf4, 0x1d, 0x1c, 0x94, 0xdc, 0xa2, 0x88, 0x74,
    0x8b, 0x19, 0x64, 0xb9, 0x81, 0x6a, 0xfa, 0x1e, 0xe7, 0xea, 0x2a, 0x0a,
    0x75, 0x42, 0xdb, 0xc7, 0xa2, 0x25, 0xd3, 0x74, 0x8a, 0x0c, 0x42, 0x50,
    0x99, 0xf0, 0x82, 0x08, 0x2d, 0xe0, 0xd9, 0x05, 0x84, 0x99, 0xc8, 0x28,
    0x68, 0x9a, 0x5b, 0xf9, 0x0e, 0xf4, 0x7a, 0x38, 0x7b, 0x64, 0x7a, 0x65,
    0xd6, 0xa6, 0x53, 0xc0, 0x47, 0xe4, 0x31, 0xe6, 0xca, 0xcf, 0xed, 0xba,
    0xdb, 0x27, 0xd2, 0x18, 0x2f, 0x41, 0x95, 0xa3, 0x03, 0xfd, 0x61, 0xbc,
    0xf8, 0x4d, 0xe8, 0x91, 0x06, 0x1d, 0x17, 0xa8, 0x0e, 0x8a, 0xfd, 0xa8,
    0x58, 0x9c, 0x80, 0xd2, 0x06, 0x1a, 0xb6, 0xb3, 0xf7, 0xcf, 0x7a, 0xde,
    0x44, 0x21, 0x1c, 0xbe, 0x1a, 0xa1, 0x2e, 0x32, 0xb0, 0x70, 0x7f, 0xcd,
    0xe0, 0x90, 0xb3, 0xcc, 0x38, 0x93, 0x0c, 0x2b, 0xe6, 0x9b, 0x4c, 0x27,
    0xd6, 0x93, 0xc7, 0x2c, 0x73, 0x89, 0x77, 0xf6, 0xeb, 0xde, 0xc8, 0x62,
    0xde, 0x95, 0xea, 0x6e, 0x67, 0x2a, 0xbe, 0xda, 0x81, 0x7e, 0xf3, 0x11,
    0x25, 0x88, 0x3c, 0xf3, 0x98, 0x1d, 0xb6, 0xeb, 0x8b, 0x77, 0xe5, 0xc8,
    0x02, 0x87, 0x9a, 0x19, 0x27, 0x86, 0xaa, 0x2a, 0x6a, 0x5e, 0x87, 0x93,
    0x27, 0x89, 0x4f, 0x76, 0x1b, 0xe9, 0x02, 0x81, 0x81, 0x00, 0xc7, 0x1e,
    0x56, 0x7b, 0x64, 0x8d, 0xa7, 0xa9, 0xa1, 0x18, 0xc5, 0x85, 0x39, 0x7a,
    0x07, 0x12, 0x36, 0xc7, 0x34, 0x02, 0xc4, 0x23, 0xcb, 0x50, 0xd0, 0x77,
    0x6b, 0xcc, 0x5d, 0x70, 0x77, 0x4e, 0x8b, 0xf9, 0x13, 0xf8, 0xe1, 0xca,
    0xd7, 0x1e, 0x8d, 0xad, 0xd7, 0xa6, 0xa7, 0xa3, 0x39, 0x23, 0x50, 0x74,
    0x7a, 0xed, 0x17, 0xfb, 0xa3, 0x51, 0x90, 0xad, 0xe3, 0xe5, 0x68, 0x27,
    0x53, 0xbe, 0x5e, 0x5f, 0x85, 0x1e, 0x35, 0x35, 0xb6, 0x87, 0xd5, 0x6d,
    0xc1, 0x8b, 0xb6, 0xd4, 0xe9, 0x73, 0x32, 0x21, 0x7e, 0x4e, 0x2f, 0x58,
    0x4b, 0xb3, 0xdf, 0x21, 0x6e, 0x86, 0xa9, 0xfb, 0x29, 0xe9, 0x67, 0x6d,
    0x44, 0xc1, 0x93, 0x66, 0x1b, 0x61, 0x75, 0x3b, 0x31, 0x0e, 0x66, 0x61,
    0x0a, 0xdd, 0xd7, 0x02, 0x98, 0xd9, 0x51, 0xeb, 0xbf, 0xf2, 0x0d, 0xea,
    0xc4, 0x18, 0x6c, 0x73, 0xdb, 0xf7, 0x02, 0x81, 0x81, 0x00, 0xfc, 0xf9,
    0x0c, 0x67, 0x7c, 0x1b, 0x71, 0x18, 0xbc, 0x48, 0xeb, 0xd0, 0x11, 0xf5,
    0x59, 0x06, 0x64, 0x99, 0x10, 0xe6, 0x29, 0x27, 0x17, 0xc3, 0x36, 0x73,
    0xce, 0x2d, 0x23, 0xbb, 0x48, 0x8c, 0x1b, 0x7e, 0x0f, 0x0a, 0xd6, 0xb9,
    0xa0, 0xd4, 0x73, 0x09, 0x33, 0xb9, 0xca, 0x93, 0xd4, 0x52, 0xfa, 0x34,
    0xbd, 0x0d, 0x7a, 0x83, 0x98, 0x8e, 0x53, 0x3d, 0x63, 0xf2, 0x3b, 0x65,
    0xdc, 0xb4, 0x15, 0x1d, 0x8c, 0xa4, 0x76, 0xc3, 0xbc, 0x27, 0xd8, 0x0d,
    0xa2, 0x45, 0x5a, 0x36, 0xf2, 0x5d, 0x18, 0x13, 0x9e, 0x06, 0xbb, 0x99,
    0xdf, 0x34, 0x6c, 0x13, 0x73, 0x21, 0x3b, 0x88, 0x23, 0xd1, 0x76, 0x05,
    0x8c, 0x18, 0xd7, 0xcc, 0x30, 0xef, 0xe8, 0x5f, 0x83, 0xbb, 0x5e, 0xb9,
    0xde, 0xe7, 0xa7, 0xa7, 0x2e, 0xbd, 0xc2, 0x45, 0x32, 0x85, 0x30, 0x1d,
    0x86, 0x40, 0x69, 0x29, 0x3d, 0x5d, 0x02, 0x81, 0x81, 0x00, 0x97, 0xfd,
    0x2f, 0x54, 0x46, 0xdf, 0xdd, 0xf0, 0x1c, 0x58, 0xd5, 0x44, 0xa9, 0x27,
    0xdd, 0x47, 0xe8, 0xea, 0x4b, 0x68, 0x25, 0x21, 0x91, 0x6b, 0x61, 0x85,
    0x16, 0x92, 0xcb, 0x6c, 0x32, 0x95, 0x91, 0x40, 0x92, 0x1f, 0x32, 0xf2,
    0xeb, 0x1b, 0x96, 0x57, 0xf1, 0x39, 0x73, 0xd2, 0xa2, 0xa5, 0xb3, 0x1f,
    0x06, 0x49, 0xfe, 0x39, 0x85, 0x63, 0x98, 0x45, 0x33, 0xa5, 0x03, 0xc8,
    0xa9, 0x22, 0xb1, 0xd4, 0xc5, 0xbe, 0xd6, 0x2c, 0xe6, 0xe4, 0x6e, 0x64,
    0xb6, 0x0d, 0x18, 0x85, 0x12, 0xa1, 0x6c, 0xcd, 0xa6, 0x24, 0xb5, 0xfc,
    0xf6, 0xe4, 0x18, 0xd8, 0xe3, 0x0e, 0x05, 0xa8, 0x03, 0x48, 0xf7, 0x3a,
    0xaf, 0xf5, 0xf6, 0xb6, 0x45, 0x06, 0x32, 0x3e, 0xf9, 0x66, 0x1d, 0x7d,
    0xcb, 0x96, 0xa4, 0x2d, 0x86, 0x50, 0xb4, 0x38, 0x78, 0xae, 0xa2, 0x32,
    0xe6, 0x76, 0x22, 0x2a, 0x99, 0xe7, 0x02, 0x81, 0x81, 0x00, 0xc1, 0x9c,
    0x22, 0x60, 0x39, 0x5e, 0x0f, 0x4a, 0xed, 0x1f, 0xaa, 0x4b, 0x0e, 0xd3,
    0x86, 0x15, 0x1c, 0x7d, 0x01, 0xb0, 0x05, 0xa3, 0x03, 0xce, 0xc6, 0x28,
    0x0f, 0x8e, 0x00, 0xa0, 0xdf, 0xbf, 0x4b, 0x73, 0x49, 0x33, 0xf4, 0x6f,
    0x11, 0xa6, 0x47, 0x7c, 0xad, 0x77, 0xee, 0x91, 0x01, 0x99, 0x98, 0x21,
    0x30, 0xe7, 0xd5, 0xf2, 0x4d, 0x99, 0xf0, 0x1f, 0x36, 0x15, 0x38, 0x5c,
    0x97, 0x73, 0xc4, 0x0d, 0x5f, 0x8c, 0xa7, 0xd0, 0xda, 0x7a, 0x6c, 0x22,
    0xd3, 0x24, 0xdd, 0x0c, 0xdc, 0xa5, 0x5f, 0x3d, 0xf4, 0x5e, 0x16, 0xca,
    0x87, 0x47, 0xe9, 0xc7, 0x60, 0xff, 0xf8, 0x3e, 0x13, 0x9b, 0xc6, 0x06,
    0x2c, 0xd8, 0xfe, 0xa0, 0x2a, 0x7c, 0x12, 0x8e, 0xb7, 0x95, 0x79, 0xc4,
    0x2b, 0xd3, 0x84, 0x3e, 0xb1, 0xc9, 0x4d, 0x9c, 0x04, 0x34, 0x67, 0x44,
    0xd1, 0x71, 0x0e, 0x8b, 0x1f, 0x89, 0x02, 0x81, 0x80, 0x08, 0xa1, 0x22,
    0x4b, 0xde, 0x11, 0xa4, 0x46, 0x57, 0xf4, 0x28, 0xd5, 0xa4, 0x93, 0xcd,
    0x50, 0xee, 0x22, 0xaf, 0xed, 0x0f, 0x44, 0x08, 0xda, 0xef, 0x3a, 0x3f,
    0x2b, 0xa9, 0x66, 0xbe, 0x9b, 0xdf, 0x8c, 0x0a, 0xb4, 0xfe, 0xb4, 0xff,
    0x6d, 0xd6, 0x1e, 0x44, 0x5d, 0xe2, 0x03, 0xe3, 0x5b, 0xf5, 0x77, 0x79,
    0xcf, 0x8c, 0x53, 0x26, 0x46, 0xd1, 0xd4, 0x4c, 0x9d, 0x32, 0xf4, 0xc0,
    0x54, 0xfa, 0x07, 0x96, 0x6d, 0x29, 0x8a, 0x6d, 0xd7, 0x35, 0x39, 0xc6,
    0xb5, 0x69, 0x86, 0xe8, 0x99, 0x09, 0xc5, 0x1d, 0xca, 0xc9, 0x40, 0xd6,
    0x2e, 0xd4, 0x87, 0xb3, 0x08, 0xbe, 0x44, 0xa4, 0x89, 0xa2, 0x1d, 0x59,
    0x85, 0xfa, 0x76, 0x6a, 0x86, 0x45, 0x6e, 0x10, 0x1a, 0xbe, 0x04, 0x4e,
    0xbd, 0x23, 0x61, 0x04, 0x80, 0x6e, 0xb0, 0x79, 0xac, 0x6f, 0x4c, 0x46,
    0x38, 0x4c, 0x5a, 0xd9, 0x71
  };

  static const unsigned char msg[] = {
    0x31, 0x26, 0x09, 0x86, 0xee, 0x94, 0x0f, 0xa7, 0x1d, 0x2c, 0x4c, 0xc7,
    0xc0, 0x0d, 0x4b, 0x1e, 0xc2, 0x13, 0x1b, 0x24, 0xf2, 0xb6, 0x24, 0x3f,
    0x48, 0xc2, 0xcb, 0xd3, 0xb7, 0xb8, 0x2e, 0xa3
  };

  static const unsigned char sig[] = {
    0xa2, 0xf8, 0xee, 0x7c, 0x03, 0x7a, 0xb3, 0xc5, 0x0f, 0x3a, 0x1c, 0x9a,
    0xe9, 0x3a, 0x66, 0xa4, 0xb8, 0x67, 0x49, 0x8f, 0xf6, 0xd9, 0x5f, 0xfe,
    0xc0, 0xa6, 0x71, 0x68, 0x5d, 0x5d, 0xa2, 0xb7, 0xde, 0xa2, 0x6e, 0x86,
    0x13, 0x9c, 0xcc, 0xdb, 0xdb, 0xc3, 0xa1, 0x2a, 0x47, 0x21, 0xe3, 0x9e,
    0xd4, 0x04, 0xfd, 0xae, 0x81, 0xd3, 0x9d, 0xe3, 0x2b, 0x79, 0x50, 0x5d,
    0x71, 0xba, 0x46, 0x5a, 0x99, 0x74, 0xa6, 0x75, 0xbe, 0x76, 0xe4, 0xeb,
    0x1d, 0x86, 0x65, 0xe6, 0xfe, 0x31, 0x31, 0x72, 0x3c, 0x0c, 0x84, 0x71,
    0xab, 0x22, 0x07, 0xeb, 0xff, 0x35, 0x4e, 0x93, 0xa1, 0x33, 0x64, 0x5a,
    0x8a, 0x30, 0xd3, 0x82, 0x3b, 0xe9, 0x71, 0x53, 0x12, 0x0f, 0x47, 0xd0,
    0x70, 0x49, 0xd2, 0x60, 0xe4, 0xee, 0xa9, 0xce, 0x5f, 0xce, 0xfe, 0x6f,
    0xc9, 0xf8, 0x39, 0x6d, 0x31, 0x36, 0x96, 0x47, 0xd1, 0x33, 0x59, 0x95,
    0x59, 0x52, 0x2c, 0xaa, 0x83, 0x51, 0xa0, 0x88, 0x9a, 0x8a, 0xe8, 0x9b,
    0xdd, 0x36, 0x5b, 0xc9, 0x79, 0x15, 0x8f, 0xda, 0xa7, 0x7f, 0xe2, 0xe0,
    0xa5, 0x70, 0xfe, 0xce, 0x6a, 0xff, 0x85, 0x63, 0x7c, 0xae, 0x78, 0xd1,
    0xcd, 0x9d, 0xde, 0xe8, 0xe9, 0xdf, 0xbd, 0x96, 0xdb, 0x74, 0x0b, 0xec,
    0xcb, 0x93, 0xc8, 0x3f, 0x27, 0x52, 0x42, 0xe6, 0x52, 0x65, 0xb0, 0xb7,
    0x5b, 0x6a, 0x52, 0x48, 0x61, 0x40, 0x17, 0x34, 0x72, 0xa5, 0x4b, 0x74,
    0x1e, 0x13, 0x46, 0x15, 0x5e, 0xf9, 0xaf, 0x0f, 0x32, 0xcd, 0x28, 0x2c,
    0x16, 0x41, 0x63, 0x54, 0xb5, 0x69, 0x30, 0xb6, 0xc4, 0x63, 0x3d, 0xc5,
    0x3a, 0x05, 0xfe, 0x59, 0xea, 0x21, 0x4e, 0x3f, 0x18, 0xad, 0xa8, 0x7b,
    0x3d, 0x1c, 0x0a, 0x5a, 0xf1, 0x60, 0x65, 0xe7, 0xec, 0x73, 0xb8, 0x46,
    0xe9, 0xba, 0x9d, 0xec
  };

  unsigned char entropy[ENTROPY_SIZE];
  unsigned char out[sizeof(sig) * 2];
  size_t len;
  bench_t tv;
  int i;

  drbg_generate(rng, entropy, sizeof(entropy));

  bench_start(&tv, "rsa_sign");

  for (i = 0; i < 1000; i++) {
    ASSERT(rsa_sign(out, &len, HASH_SHA256, msg, sizeof(msg),
                                            priv, sizeof(priv),
                                            entropy));
    ASSERT(len == sizeof(sig));
    ASSERT(torsion_memcmp_var(out, sig, sizeof(sig)) == 0);
  }

  bench_end(&tv, i);
}

static void
bench_rsa_verify(drbg_t *rng) {
  static const unsigned char pub[] = {
    0x30, 0x82, 0x01, 0x0a, 0x02, 0x82, 0x01, 0x01, 0x00, 0xc4, 0xc3, 0x93,
    0x49, 0x7e, 0xca, 0xdd, 0x3d, 0xac, 0x9c, 0x46, 0xe4, 0xab, 0x26, 0x04,
    0xcb, 0xfe, 0x46, 0x79, 0x57, 0xa8, 0xd5, 0x6c, 0xee, 0x4d, 0x6b, 0xe2,
    0x6f, 0xd9, 0xe6, 0x23, 0x2c, 0x61, 0x3d, 0x09, 0x9d, 0xc4, 0xa0, 0xb9,
    0x5b, 0xd4, 0x37, 0x9e, 0x46, 0x3b, 0xbe, 0x62, 0xfa, 0x26, 0xf2, 0x3c,
    0x18, 0x73, 0x78, 0x2f, 0x7b, 0x21, 0xd0, 0x00, 0x9d, 0xf9, 0x50, 0xb8,
    0xe7, 0xb1, 0x06, 0x02, 0xc8, 0x0e, 0xd0, 0xda, 0x55, 0x5c, 0x8a, 0x4b,
    0x1c, 0x59, 0x49, 0x0e, 0x75, 0x03, 0x27, 0xe9, 0x64, 0x45, 0x4d, 0x13,
    0xfc, 0x6a, 0x64, 0xcf, 0xbc, 0xc7, 0x3f, 0x83, 0xc6, 0x33, 0x0b, 0x3a,
    0x77, 0xd7, 0xf7, 0x73, 0xaf, 0x36, 0x3a, 0xf3, 0xe7, 0xcb, 0x09, 0x31,
    0xc0, 0xbd, 0xcc, 0x09, 0x30, 0x07, 0x0a, 0x49, 0x39, 0xe7, 0x54, 0x96,
    0xe7, 0xc4, 0xd3, 0x8a, 0x30, 0xf4, 0xe3, 0x16, 0xd1, 0x6d, 0x62, 0x3e,
    0xc6, 0x54, 0xd6, 0xf7, 0x95, 0x63, 0x51, 0x94, 0xbb, 0x29, 0x98, 0x54,
    0x96, 0x1d, 0x89, 0x78, 0x47, 0x97, 0x40, 0xb6, 0x47, 0xde, 0xda, 0xcb,
    0x91, 0x9b, 0x22, 0xf9, 0x87, 0x0f, 0x2a, 0xef, 0xf0, 0xbb, 0x69, 0x59,
    0xf7, 0xd7, 0xa0, 0xb1, 0x06, 0xed, 0xa3, 0x1b, 0xbe, 0x4b, 0x18, 0x70,
    0xbb, 0xe1, 0x79, 0x72, 0xa7, 0x1d, 0x04, 0x9d, 0x06, 0x34, 0x50, 0x51,
    0xc2, 0xb1, 0x0f, 0x22, 0x67, 0x5e, 0xc8, 0x00, 0x70, 0x71, 0x8b, 0x32,
    0xc2, 0x61, 0x6f, 0xc1, 0x93, 0x75, 0x09, 0x75, 0x2c, 0xfe, 0x21, 0x63,
    0xc0, 0x49, 0x9c, 0x10, 0x76, 0x2b, 0x25, 0x51, 0x73, 0x81, 0x1a, 0xdb,
    0x85, 0x87, 0x7a, 0x43, 0xc5, 0x4a, 0xfb, 0x95, 0xa9, 0x4e, 0x6f, 0xdd,
    0x42, 0x8e, 0x05, 0x43, 0xd8, 0x68, 0x7d, 0x3a, 0x3d, 0x8b, 0x0f, 0xc3,
    0xbb, 0x02, 0x03, 0x01, 0x00, 0x01
  };

  static const unsigned char msg[] = {
    0x31, 0x26, 0x09, 0x86, 0xee, 0x94, 0x0f, 0xa7, 0x1d, 0x2c, 0x4c, 0xc7,
    0xc0, 0x0d, 0x4b, 0x1e, 0xc2, 0x13, 0x1b, 0x24, 0xf2, 0xb6, 0x24, 0x3f,
    0x48, 0xc2, 0xcb, 0xd3, 0xb7, 0xb8, 0x2e, 0xa3
  };

  static const unsigned char sig[] = {
    0xa2, 0xf8, 0xee, 0x7c, 0x03, 0x7a, 0xb3, 0xc5, 0x0f, 0x3a, 0x1c, 0x9a,
    0xe9, 0x3a, 0x66, 0xa4, 0xb8, 0x67, 0x49, 0x8f, 0xf6, 0xd9, 0x5f, 0xfe,
    0xc0, 0xa6, 0x71, 0x68, 0x5d, 0x5d, 0xa2, 0xb7, 0xde, 0xa2, 0x6e, 0x86,
    0x13, 0x9c, 0xcc, 0xdb, 0xdb, 0xc3, 0xa1, 0x2a, 0x47, 0x21, 0xe3, 0x9e,
    0xd4, 0x04, 0xfd, 0xae, 0x81, 0xd3, 0x9d, 0xe3, 0x2b, 0x79, 0x50, 0x5d,
    0x71, 0xba, 0x46, 0x5a, 0x99, 0x74, 0xa6, 0x75, 0xbe, 0x76, 0xe4, 0xeb,
    0x1d, 0x86, 0x65, 0xe6, 0xfe, 0x31, 0x31, 0x72, 0x3c, 0x0c, 0x84, 0x71,
    0xab, 0x22, 0x07, 0xeb, 0xff, 0x35, 0x4e, 0x93, 0xa1, 0x33, 0x64, 0x5a,
    0x8a, 0x30, 0xd3, 0x82, 0x3b, 0xe9, 0x71, 0x53, 0x12, 0x0f, 0x47, 0xd0,
    0x70, 0x49, 0xd2, 0x60, 0xe4, 0xee, 0xa9, 0xce, 0x5f, 0xce, 0xfe, 0x6f,
    0xc9, 0xf8, 0x39, 0x6d, 0x31, 0x36, 0x96, 0x47, 0xd1, 0x33, 0x59, 0x95,
    0x59, 0x52, 0x2c, 0xaa, 0x83, 0x51, 0xa0, 0x88, 0x9a, 0x8a, 0xe8, 0x9b,
    0xdd, 0x36, 0x5b, 0xc9, 0x79, 0x15, 0x8f, 0xda, 0xa7, 0x7f, 0xe2, 0xe0,
    0xa5, 0x70, 0xfe, 0xce, 0x6a, 0xff, 0x85, 0x63, 0x7c, 0xae, 0x78, 0xd1,
    0xcd, 0x9d, 0xde, 0xe8, 0xe9, 0xdf, 0xbd, 0x96, 0xdb, 0x74, 0x0b, 0xec,
    0xcb, 0x93, 0xc8, 0x3f, 0x27, 0x52, 0x42, 0xe6, 0x52, 0x65, 0xb0, 0xb7,
    0x5b, 0x6a, 0x52, 0x48, 0x61, 0x40, 0x17, 0x34, 0x72, 0xa5, 0x4b, 0x74,
    0x1e, 0x13, 0x46, 0x15, 0x5e, 0xf9, 0xaf, 0x0f, 0x32, 0xcd, 0x28, 0x2c,
    0x16, 0x41, 0x63, 0x54, 0xb5, 0x69, 0x30, 0xb6, 0xc4, 0x63, 0x3d, 0xc5,
    0x3a, 0x05, 0xfe, 0x59, 0xea, 0x21, 0x4e, 0x3f, 0x18, 0xad, 0xa8, 0x7b,
    0x3d, 0x1c, 0x0a, 0x5a, 0xf1, 0x60, 0x65, 0xe7, 0xec, 0x73, 0xb8, 0x46,
    0xe9, 0xba, 0x9d, 0xec
  };

  bench_t tv;
  int i;

  (void)rng;

  bench_start(&tv, "rsa_verify");

  for (i = 0; i < 50000; i++) {
    ASSERT(rsa_verify(HASH_SHA256, msg, sizeof(msg),
                                   sig, sizeof(sig),
                                   pub, sizeof(pub)));
  }

  bench_end(&tv, i);
}

static void
bench_hash(drbg_t *rng) {
  unsigned char chain[32];
  bench_t tv;
  hash_t hash;
  size_t i;

  drbg_generate(rng, chain, sizeof(chain));

  bench_start(&tv, "hash");

  for (i = 0; i < 10000000; i++) {
    hash_init(&hash, HASH_SHA256);
    hash_update(&hash, chain, sizeof(chain));
    hash_final(&hash, chain, sizeof(chain));
  }

  bench_end(&tv, i);
}

static void
bench_sha256(drbg_t *rng) {
  unsigned char chain[32];
  bench_t tv;
  sha256_t sha;
  size_t i;

  drbg_generate(rng, chain, sizeof(chain));

  bench_start(&tv, "sha256");

  for (i = 0; i < 10000000; i++) {
    sha256_init(&sha);
    sha256_update(&sha, chain, sizeof(chain));
    sha256_final(&sha, chain);
  }

  bench_end(&tv, i);
}

static void
bench_sha3(drbg_t *rng) {
  unsigned char chain[32];
  bench_t tv;
  sha3_t sha;
  size_t i;

  drbg_generate(rng, chain, sizeof(chain));

  bench_start(&tv, "sha3");

  for (i = 0; i < 10000000; i++) {
    sha3_256_init(&sha);
    sha3_256_update(&sha, chain, sizeof(chain));
    sha3_256_final(&sha, chain);
  }

  bench_end(&tv, i);
}

static void
bench_aes_ctr(drbg_t *rng) {
  unsigned char raw[977];
  unsigned char key[32];
  unsigned char iv[16];
  cipher_t cipher;
  ctr_t mode;
  bench_t tv;
  size_t i;

  drbg_generate(rng, raw, sizeof(raw));
  drbg_generate(rng, key, sizeof(key));
  drbg_generate(rng, iv, sizeof(iv));

  bench_start(&tv, "aes_ctr");

  cipher_init(&cipher, CIPHER_AES256, key, sizeof(key));

  ctr_init(&mode, &cipher, iv);

  for (i = 0; i < 100000; i++)
    ctr_crypt(&mode, &cipher, raw, raw, sizeof(raw));

  bench_end(&tv, i);
}

static void
bench_aes_gcm(drbg_t *rng) {
  unsigned char raw[977];
  unsigned char key[32];
  unsigned char iv[12];
  cipher_t cipher;
  gcm_t mode;
  bench_t tv;
  size_t i;

  drbg_generate(rng, raw, sizeof(raw));
  drbg_generate(rng, key, sizeof(key));
  drbg_generate(rng, iv, sizeof(iv));

  bench_start(&tv, "aes_gcm");

  cipher_init(&cipher, CIPHER_AES256, key, sizeof(key));

  gcm_init(&mode, &cipher, iv, sizeof(iv));

  for (i = 0; i < 100000; i++)
    gcm_encrypt(&mode, &cipher, raw, raw, sizeof(raw));

  bench_end(&tv, i);
}

static void
bench_chacha20(drbg_t *rng) {
  unsigned char raw[977];
  unsigned char key[32];
  unsigned char iv[8];
  chacha20_t ctx;
  bench_t tv;
  size_t i;

  drbg_generate(rng, raw, sizeof(raw));
  drbg_generate(rng, key, sizeof(key));
  drbg_generate(rng, iv, sizeof(iv));

  bench_start(&tv, "chacha20");

  chacha20_init(&ctx, key, sizeof(key), iv, sizeof(iv), 0);

  for (i = 0; i < 1000000; i++)
    chacha20_crypt(&ctx, raw, raw, sizeof(raw));

  bench_end(&tv, i);
}

/*
 * Benchmark Registry
 */

static const struct {
  const char *name;
  void (*run)(drbg_t *);
} torsion_benches[] = {
#define B(name) { #name, bench_ ## name }
  B(ecdsa_pubkey_create),
  B(ecdsa_pubkey_tweak_add),
  B(ecdsa_sign),
  B(ecdsa_verify),
  B(ecdsa_derive),
  B(ecdh_derive),
  B(eddsa_sign),
  B(eddsa_verify),
  B(eddsa_derive),
  B(mpi_internal),
  B(rsa_generate),
  B(rsa_sign),
  B(rsa_verify),
  B(hash),
  B(sha256),
  B(sha3),
  B(aes_ctr),
  B(aes_gcm),
  B(chacha20)
#undef B
};

/*
 * Main
 */

int
main(int argc, char **argv) {
  drbg_t rng;
  size_t i, j;
  int found;

  drbg_init_rand(&rng);

  if (argc <= 1) {
    for (i = 0; i < ARRAY_SIZE(torsion_benches); i++)
      torsion_benches[i].run(&rng);

    return 0;
  }

  for (i = 1; i < (size_t)argc; i++) {
    found = 0;

    for (j = 0; j < ARRAY_SIZE(torsion_benches); j++) {
      if (strcmp(torsion_benches[j].name, argv[i]) == 0) {
        torsion_benches[j].run(&rng);
        found = 1;
        break;
      }
    }

    if (!found) {
      fprintf(stderr, "Unknown benchmark: %s.\n", argv[i]);
      return 1;
    }
  }

  return 0;
}
