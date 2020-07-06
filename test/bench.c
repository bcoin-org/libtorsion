/*!
 * bench.c - benchmarks for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <inttypes.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <torsion/drbg.h>
#include <torsion/ecc.h>
#include <torsion/hash.h>

#include "hrtime.h"
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
bench_end(bench_t *start, uint64_t ops) {
  bench_t nsec = torsion_hrtime() - *start;
  double sec = (double)nsec / 1000000000.0;

  printf("  Operations:  %" PRIu64 "\n", ops);
  printf("  Nanoseconds: %" PRIu64 "\n", nsec);
  printf("  Seconds:     %f\n", sec);
  printf("  Ops/Sec:     %f\n", (double)ops / sec);
}

static void
bench_ecdsa_pubkey_create(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char pub[33];
  size_t pub_len;
  unsigned char bytes[32];
  bench_t tv;
  size_t i;

  drbg_generate(rng, bytes, 32);

  bench_start(&tv, "ecdsa_pubkey_create");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_pubkey_create(ec, pub, &pub_len, bytes, 1));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdsa_derive(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char k0[32];
  unsigned char k1[32];
  unsigned char p0[65];
  unsigned char p1[33];
  size_t p0_len, p1_len;
  bench_t tv;
  size_t i;

  drbg_generate(rng, k0, 32);
  drbg_generate(rng, k1, 32);

  ASSERT(ecdsa_pubkey_create(ec, p0, &p0_len, k0, 0));

  bench_start(&tv, "ecdsa_derive");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_derive(ec, p1, &p1_len, p0, p0_len, k1, 1));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdsa(drbg_t *rng) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char pub[33];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));
  drbg_generate(rng, priv, sizeof(priv));
  drbg_generate(rng, msg, sizeof(msg));

  priv[0] = 0;

  wei_curve_randomize(ec, entropy);

  ASSERT(ecdsa_sign(ec, sig, NULL, msg, 32, priv));
  ASSERT(ecdsa_pubkey_create(ec, pub, NULL, priv, 1));

  bench_start(&tv, "ecdsa_verify");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_verify(ec, msg, 32, sig, pub, 33));

  bench_end(&tv, i);

  bench_start(&tv, "ecdsa_sign");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdsa_sign(ec, sig, NULL, msg, 32, priv));

  bench_end(&tv, i);

  wei_curve_destroy(ec);
}

static void
bench_ecdh(drbg_t *rng) {
  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char priv[MONT_MAX_SCALAR_SIZE];
  unsigned char pub[MONT_MAX_FIELD_SIZE];
  unsigned char secret[MONT_MAX_FIELD_SIZE];
  bench_t tv;
  size_t i;

  drbg_generate(rng, priv, sizeof(priv));

  ecdh_pubkey_create(ec, pub, priv);

  bench_start(&tv, "ecdh_derive");

  for (i = 0; i < 10000; i++)
    ASSERT(ecdh_derive(ec, secret, pub, priv));

  bench_end(&tv, i);

  mont_curve_destroy(ec);
}

static void
bench_eddsa(drbg_t *rng) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char pub[32];
  bench_t tv;
  size_t i;

  drbg_generate(rng, entropy, sizeof(entropy));
  drbg_generate(rng, priv, sizeof(priv));
  drbg_generate(rng, msg, sizeof(msg));

  edwards_curve_randomize(ec, entropy);

  eddsa_sign(ec, sig, msg, 32, priv, -1, NULL, 0);
  eddsa_pubkey_create(ec, pub, priv);

  bench_start(&tv, "eddsa_verify");

  for (i = 0; i < 10000; i++)
    ASSERT(eddsa_verify(ec, msg, 32, sig, pub, -1, NULL, 0));

  bench_end(&tv, i);

  edwards_curve_destroy(ec);
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

/*
 * Benchmark Registry
 */

static const struct {
  const char *name;
  void (*run)(drbg_t *);
} torsion_benches[] = {
#define B(name) { #name, bench_ ## name }
  B(ecdsa_pubkey_create),
  B(ecdsa_derive),
  B(ecdsa),
  B(ecdh),
  B(eddsa),
  B(hash),
  B(sha256)
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
