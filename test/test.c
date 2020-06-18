/*!
 * test.c - tests for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <torsion/aead.h>
#include <torsion/chacha20.h>
#include <torsion/drbg.h>
#include <torsion/dsa.h>
#include <torsion/ecc.h>
#include <torsion/encoding.h>
#include <torsion/hash.h>
#include <torsion/kdf.h>
#include <torsion/poly1305.h>
#ifdef TORSION_HAVE_RNG
#include <torsion/rand.h>
#endif
#include <torsion/rsa.h>
#include <torsion/salsa20.h>
#include <torsion/secretbox.h>
#include <torsion/siphash.h>
#include <torsion/util.h>

#include "testutil.h"

#include "data/aead-vectors.h"
#include "data/chacha20-vectors.h"
#include "data/hash-drbg-vectors.h"
#include "data/hmac-drbg-vectors.h"
#include "data/ctr-drbg-vectors.h"
#include "data/cipher-vectors.h"
#include "data/cipher-mode-vectors.h"
#include "data/cipher-aead-vectors.h"
#include "data/dsa-vectors.h"
#include "data/ecdsa-vectors.h"
#include "data/schnorr-legacy-vectors.h"
#include "data/schnorr-vectors.h"
#include "data/eddsa-vectors.h"
#include "data/hash-vectors.h"
#include "data/hmac-vectors.h"
#include "data/eb2k-vectors.h"
#include "data/hkdf-vectors.h"
#include "data/pbkdf2-vectors.h"
#include "data/poly1305-vectors.h"
#include "data/rsa-vectors.h"

/*
 * String Maps
 */

static const char *cipher_names[24] = {
  "AES128",
  "AES192",
  "AES256",
  "BLOWFISH",
  "CAMELLIA128",
  "CAMELLIA192",
  "CAMELLIA256",
  "CAST5",
  "DES",
  "DES_EDE",
  "DES_EDE3",
  "IDEA",
  "RC2",
  "RC2_GUTMANN",
  "RC2_40",
  "RC2_64",
  "RC2_128",
  "RC2_128_GUTMANN",
  "SERPENT128",
  "SERPENT192",
  "SERPENT256",
  "TWOFISH128",
  "TWOFISH192",
  "TWOFISH256",
};

static const char *cipher_modes[11] = {
  "RAW",
  "ECB",
  "CBC",
  "CTS",
  "XTS",
  "CTR",
  "CFB",
  "OFB",
  "GCM",
  "CCM",
  "EAX"
};

static const char *wei_curves[6] = {
  "P192",
  "P224",
  "P256",
  "P384",
  "P521",
  "SECP256K1"
};

static const char *mont_curves[2] = {
  "X25519",
  "X448"
};

static const char *edwards_curves[3] = {
  "ED25519",
  "ED448",
  "ED1174"
};

static const char *hash_names[32] = {
  "BLAKE2B_160",
  "BLAKE2B_256",
  "BLAKE2B_384",
  "BLAKE2B_512",
  "BLAKE2S_128",
  "BLAKE2S_160",
  "BLAKE2S_224",
  "BLAKE2S_256",
  "GOST94",
  "HASH160",
  "HASH256",
  "KECCAK224",
  "KECCAK256",
  "KECCAK384",
  "KECCAK512",
  "MD2",
  "MD4",
  "MD5",
  "MD5SHA1",
  "RIPEMD160",
  "SHA1",
  "SHA224",
  "SHA256",
  "SHA384",
  "SHA512",
  "SHA3_224",
  "SHA3_256",
  "SHA3_384",
  "SHA3_512",
  "SHAKE128",
  "SHAKE256",
  "WHIRLPOOL"
};

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

static unsigned int
drbg_uniform(drbg_t *rng, unsigned int mod) {
  unsigned int x;

  if (mod == 0)
    return 0;

  drbg_generate(rng, &x, sizeof(x));

  return x % mod;
}

static void
hex_decode(unsigned char *out, size_t *out_len, const char *str) {
  size_t size = strlen(str);

  ASSERT((size & 1) == 0);
  ASSERT(*out_len >= size / 2);
  ASSERT(base16_decode(out, out_len, str, size));
}

static void
hex_parse(unsigned char *out, size_t len, const char *str) {
  size_t out_len = len;

  hex_decode(out, &out_len, str);

  ASSERT(out_len == len);
}

/*
 * AEAD
 */

static void
test_aead(void) {
  static unsigned char input[8192];
  static unsigned char aad[128];
  static unsigned char key[32];
  static unsigned char nonce[32];
  static unsigned char raw[8192];
  static unsigned char data[8192];
  static unsigned char output[8192];
  static unsigned char tag[16];
  static unsigned char mac[16];
  aead_t ctx;
  size_t i;

  printf("Testing AEAD...\n");

  for (i = 0; i < ARRAY_SIZE(aead_vectors); i++) {
    size_t input_len = sizeof(input);
    size_t aad_len = sizeof(aad);
    size_t nonce_len = sizeof(nonce);
    size_t raw_len = sizeof(raw);
    size_t data_len, output_len;

    printf("  - AEAD vector #%lu\n", i + 1);

    hex_decode(input, &input_len, aead_vectors[i][0]);
    hex_decode(aad, &aad_len, aead_vectors[i][1]);
    hex_parse(key, 32, aead_vectors[i][2]);
    hex_decode(nonce, &nonce_len, aead_vectors[i][3]);
    hex_decode(raw, &raw_len, aead_vectors[i][4]);

    ASSERT(raw_len >= 16);

    data_len = input_len;
    output_len = raw_len - 16;

    memcpy(data, input, input_len);
    memcpy(output, raw, output_len);
    memcpy(tag, raw + output_len, 16);

    aead_init(&ctx, key, nonce, nonce_len);
    aead_aad(&ctx, aad, aad_len);
    aead_encrypt(&ctx, data, data, data_len);

    ASSERT(memcmp(data, output, output_len) == 0);

    aead_final(&ctx, mac);

    ASSERT(memcmp(mac, tag, 16) == 0);
    ASSERT(aead_verify(mac, tag));

    aead_init(&ctx, key, nonce, nonce_len);
    aead_aad(&ctx, aad, aad_len);
    aead_auth(&ctx, data, data_len);

    aead_final(&ctx, mac);

    ASSERT(memcmp(mac, tag, 16) == 0);
    ASSERT(aead_verify(mac, tag));

    aead_init(&ctx, key, nonce, nonce_len);
    aead_aad(&ctx, aad, aad_len);
    aead_decrypt(&ctx, data, data, data_len);

    ASSERT(memcmp(data, input, input_len) == 0);

    aead_final(&ctx, mac);

    ASSERT(memcmp(mac, tag, 16) == 0);
    ASSERT(aead_verify(mac, tag));
  }
}

/*
 * ChaCha20
 */

static void
test_chacha20(void) {
  unsigned char key[32];
  unsigned char nonce[24];
  unsigned char input[4096];
  unsigned char output[4096];
  unsigned char data[4096];
  chacha20_t ctx;
  size_t i;

  printf("Testing ChaCha20...\n");

  for (i = 0; i < ARRAY_SIZE(chacha20_vectors); i++) {
    size_t key_len = sizeof(key);
    size_t nonce_len = sizeof(nonce);
    unsigned int counter = chacha20_vectors[i].counter;
    size_t input_len = sizeof(input);
    size_t output_len = sizeof(output);
    size_t data_len;

    printf("  - ChaCha20 vector #%lu\n", i + 1);

    hex_decode(key, &key_len, chacha20_vectors[i].key);
    hex_decode(nonce, &nonce_len, chacha20_vectors[i].nonce);
    hex_decode(input, &input_len, chacha20_vectors[i].input);
    hex_decode(output, &output_len, chacha20_vectors[i].output);

    ASSERT(input_len == output_len);

    data_len = input_len;

    memcpy(data, input, input_len);

    chacha20_init(&ctx, key, key_len, nonce, nonce_len, counter);
    chacha20_encrypt(&ctx, data, data, data_len);

    ASSERT(memcmp(data, output, output_len) == 0);

    chacha20_init(&ctx, key, key_len, nonce, nonce_len, counter);
    chacha20_encrypt(&ctx, data, data, data_len);

    ASSERT(memcmp(data, input, input_len) == 0);
  }
}

/*
 * Cipher
 */

static void
test_ciphers(void) {
  unsigned char key[32];
  unsigned char iv[CIPHER_MAX_BLOCK_SIZE];
  unsigned char expect[CIPHER_MAX_BLOCK_SIZE];
  unsigned char data[CIPHER_MAX_BLOCK_SIZE];
  cipher_t ctx;
  size_t i, j;

  printf("Testing ciphers...\n");

  for (i = 0; i < ARRAY_SIZE(cipher_vectors); i++) {
    int type = cipher_vectors[i].type;
    size_t size = cipher_block_size(type);
    size_t key_len = sizeof(key);

    ASSERT(type >= 0 && type <= CIPHER_MAX);
    ASSERT(size != 0);

    printf("  - Cipher vector #%lu (%s)\n", i + 1, cipher_names[type]);

    hex_decode(key, &key_len, cipher_vectors[i].key);
    hex_parse(iv, size, cipher_vectors[i].iv);
    hex_parse(expect, size, cipher_vectors[i].expect);

    memcpy(data, iv, size);

    cipher_init(&ctx, type, key, key_len);

    for (j = 0; j < 1000; j++)
      cipher_encrypt(&ctx, data, data);

    ASSERT(memcmp(data, expect, size) == 0);

    for (j = 0; j < 1000; j++)
      cipher_decrypt(&ctx, data, data);

    ASSERT(memcmp(data, iv, size) == 0);
  }
}

static void
test_cipher_modes(void) {
  unsigned char key[64];
  unsigned char iv[CIPHER_MAX_BLOCK_SIZE];
  unsigned char input[64];
  unsigned char output[64];
  unsigned char data[64];
  size_t i;

  printf("Testing cipher modes...\n");

  for (i = 0; i < ARRAY_SIZE(cipher_mode_vectors); i++) {
    int type = cipher_mode_vectors[i].type;
    int mode = cipher_mode_vectors[i].mode;
    size_t key_len = sizeof(key);
    size_t iv_len = sizeof(iv);
    size_t input_len = sizeof(input);
    size_t output_len = sizeof(output);
    size_t len;

    ASSERT(type >= 0 && type <= CIPHER_MAX);
    ASSERT(mode >= 0 && mode <= CIPHER_MODE_MAX);

    printf("  - Cipher mode vector #%lu (%s-%s)\n",
           i + 1, cipher_names[type], cipher_modes[mode]);

    hex_decode(key, &key_len, cipher_mode_vectors[i].key);
    hex_decode(iv, &iv_len, cipher_mode_vectors[i].iv);
    hex_decode(input, &input_len, cipher_mode_vectors[i].input);
    hex_decode(output, &output_len, cipher_mode_vectors[i].output);

    ASSERT(sizeof(data) >= CIPHER_MAX_ENCRYPT_SIZE(input_len));
    ASSERT(sizeof(data) >= CIPHER_MAX_DECRYPT_SIZE(output_len));

    ASSERT(cipher_static_encrypt(data, &len, type, mode,
                                 key, key_len, iv, iv_len,
                                 input, input_len));

    ASSERT(len == output_len);
    ASSERT(memcmp(data, output, len) == 0);

    ASSERT(cipher_static_decrypt(data, &len, type, mode,
                                 key, key_len, iv, iv_len,
                                 output, output_len));

    ASSERT(len == input_len);
    ASSERT(memcmp(data, input, len) == 0);
  }
}

static void
test_cipher_aead(void) {
  unsigned char key[32];
  unsigned char iv[16];
  unsigned char aad[64];
  unsigned char input[64];
  unsigned char output[64];
  unsigned char tag[16];
  unsigned char data[64];
  unsigned char mac[16];
  cipher_stream_t ctx;
  size_t i;

  printf("Testing cipher AEADs...\n");

  for (i = 0; i < ARRAY_SIZE(cipher_aead_vectors); i++) {
    int type = cipher_aead_vectors[i].type;
    int mode = cipher_aead_vectors[i].mode;
    size_t key_len = sizeof(key);
    size_t iv_len = sizeof(iv);
    size_t aad_len = sizeof(aad);
    size_t input_len = sizeof(input);
    size_t output_len = sizeof(output);
    size_t tag_len = sizeof(tag);
    size_t len;

    ASSERT(type >= 0 && type <= CIPHER_MAX);
    ASSERT(mode >= 0 && mode <= CIPHER_MODE_MAX);

    printf("  - Cipher AEAD vector #%lu (%s-%s)\n",
           i + 1, cipher_names[type], cipher_modes[mode]);

    hex_decode(key, &key_len, cipher_aead_vectors[i].key);
    hex_decode(iv, &iv_len, cipher_aead_vectors[i].iv);
    hex_decode(aad, &aad_len, cipher_aead_vectors[i].aad);
    hex_decode(input, &input_len, cipher_aead_vectors[i].input);
    hex_decode(output, &output_len, cipher_aead_vectors[i].output);
    hex_decode(tag, &tag_len, cipher_aead_vectors[i].tag);

    ASSERT(input_len == output_len);
    ASSERT(cipher_stream_init(&ctx, type, mode, 1, key, key_len, iv, iv_len));

    if (mode == CIPHER_MODE_CCM)
      ASSERT(cipher_stream_set_ccm(&ctx, input_len, tag_len, aad, aad_len));
    else
      ASSERT(cipher_stream_set_aad(&ctx, aad, aad_len));

    ASSERT(sizeof(data) >= cipher_stream_update_size(&ctx, input_len));

    cipher_stream_update(&ctx, data, &len, input, input_len);

    ASSERT(len == output_len);

    ASSERT(cipher_stream_final(&ctx, data, &len));
    ASSERT(len == 0);
    ASSERT(memcmp(data, output, output_len) == 0);

    ASSERT(cipher_stream_get_tag(&ctx, mac, &len));
    ASSERT(len == tag_len && memcmp(mac, tag, tag_len) == 0);

    ASSERT(cipher_stream_init(&ctx, type, mode, 0, key, key_len, iv, iv_len));

    if (mode == CIPHER_MODE_CCM)
      ASSERT(cipher_stream_set_ccm(&ctx, output_len, tag_len, aad, aad_len));
    else
      ASSERT(cipher_stream_set_aad(&ctx, aad, aad_len));

    ASSERT(cipher_stream_set_tag(&ctx, tag, tag_len));
    ASSERT(sizeof(data) >= cipher_stream_update_size(&ctx, output_len));

    cipher_stream_update(&ctx, data, &len, output, output_len);

    ASSERT(len == input_len);

    ASSERT(cipher_stream_final(&ctx, data, &len));
    ASSERT(len == 0);
    ASSERT(memcmp(data, input, input_len) == 0);
  }
}

/*
 * DRBG
 */

static void
test_hash_drbg(void) {
  unsigned char entropy[256];
  unsigned char reseed[256];
  unsigned char add1[256];
  unsigned char add2[256];
  unsigned char expect[256];
  unsigned char data[256];
  hash_drbg_t drbg;
  size_t i;

  printf("Testing Hash-DRBG...\n");

  for (i = 0; i < ARRAY_SIZE(hash_drbg_vectors); i++) {
    int type = hash_drbg_vectors[i].type;
    size_t entropy_len = sizeof(entropy);
    size_t reseed_len = sizeof(reseed);
    size_t add1_len = sizeof(add1);
    size_t add2_len = sizeof(add2);
    size_t expect_len = sizeof(expect);
    size_t data_len;

    printf("  - Hash-DRBG vector #%lu\n", i + 1);

    hex_decode(entropy, &entropy_len, hash_drbg_vectors[i].entropy);
    hex_decode(reseed, &reseed_len, hash_drbg_vectors[i].reseed);
    hex_decode(add1, &add1_len, hash_drbg_vectors[i].add1);
    hex_decode(add2, &add2_len, hash_drbg_vectors[i].add2);
    hex_decode(expect, &expect_len, hash_drbg_vectors[i].expect);

    data_len = expect_len;

    hash_drbg_init(&drbg, type, entropy, entropy_len);
    hash_drbg_reseed(&drbg, reseed, reseed_len);
    hash_drbg_generate(&drbg, data, data_len, add1, add1_len);
    hash_drbg_generate(&drbg, data, data_len, add2, add2_len);

    ASSERT(memcmp(data, expect, expect_len) == 0);
  }
}

static void
test_hmac_drbg(void) {
  unsigned char entropy[256];
  unsigned char reseed[256];
  unsigned char add1[256];
  unsigned char add2[256];
  unsigned char expect[256];
  unsigned char data[256];
  hmac_drbg_t drbg;
  size_t i;

  printf("Testing HMAC-DRBG...\n");

  for (i = 0; i < ARRAY_SIZE(hmac_drbg_vectors); i++) {
    int type = hmac_drbg_vectors[i].type;
    size_t entropy_len = sizeof(entropy);
    size_t reseed_len = sizeof(reseed);
    size_t add1_len = sizeof(add1);
    size_t add2_len = sizeof(add2);
    size_t expect_len = sizeof(expect);
    size_t data_len;

    printf("  - HMAC-DRBG vector #%lu\n", i + 1);

    hex_decode(entropy, &entropy_len, hmac_drbg_vectors[i].entropy);
    hex_decode(reseed, &reseed_len, hmac_drbg_vectors[i].reseed);
    hex_decode(add1, &add1_len, hmac_drbg_vectors[i].add1);
    hex_decode(add2, &add2_len, hmac_drbg_vectors[i].add2);
    hex_decode(expect, &expect_len, hmac_drbg_vectors[i].expect);

    data_len = expect_len;

    hmac_drbg_init(&drbg, type, entropy, entropy_len);
    hmac_drbg_reseed(&drbg, reseed, reseed_len);
    hmac_drbg_generate(&drbg, data, data_len, add1, add1_len);
    hmac_drbg_generate(&drbg, data, data_len, add2, add2_len);

    ASSERT(memcmp(data, expect, expect_len) == 0);
  }
}

static void
test_ctr_drbg(void) {
  unsigned char entropy[256];
  unsigned char pers[256];
  unsigned char reseed[256];
  unsigned char add[256];
  unsigned char add1[256];
  unsigned char add2[256];
  unsigned char expect[256];
  unsigned char data[256];
  ctr_drbg_t drbg;
  size_t i;

  printf("Testing CTR-DRBG...\n");

  for (i = 0; i < ARRAY_SIZE(ctr_drbg_vectors); i++) {
    unsigned int bits = ctr_drbg_vectors[i].bits;
    int df = ctr_drbg_vectors[i].derivation;
    size_t entropy_len = sizeof(entropy);
    size_t pers_len = sizeof(pers);
    size_t reseed_len = sizeof(reseed);
    size_t add_len = sizeof(add);
    size_t add1_len = sizeof(add1);
    size_t add2_len = sizeof(add2);
    size_t expect_len = sizeof(expect);
    size_t data_len;

    printf("  - CTR-DRBG vector #%lu\n", i + 1);

    hex_decode(entropy, &entropy_len, ctr_drbg_vectors[i].entropy);
    hex_decode(pers, &pers_len, ctr_drbg_vectors[i].pers);
    hex_decode(reseed, &reseed_len, ctr_drbg_vectors[i].reseed);
    hex_decode(add, &add_len, ctr_drbg_vectors[i].add);
    hex_decode(add1, &add1_len, ctr_drbg_vectors[i].add1);
    hex_decode(add2, &add2_len, ctr_drbg_vectors[i].add2);
    hex_decode(expect, &expect_len, ctr_drbg_vectors[i].expect);

    data_len = expect_len;

    ctr_drbg_init(&drbg, bits, df, entropy, entropy_len, pers, pers_len);
    ctr_drbg_reseed(&drbg, reseed, reseed_len, add, add_len);
    ctr_drbg_generate(&drbg, data, data_len, add1, add1_len);
    ctr_drbg_generate(&drbg, data, data_len, add2, add2_len);

    ASSERT(memcmp(data, expect, expect_len) == 0);
  }
}

/*
 * DSA
 */

static void
test_dsa(void) {
  unsigned char params[DSA_MAX_PARAMS_SIZE];
  unsigned char priv[DSA_MAX_PRIV_SIZE];
  unsigned char pub[DSA_MAX_PUB_SIZE];
  unsigned char msg[64];
  unsigned char sig[DSA_MAX_SIG_SIZE];
  unsigned char der[DSA_MAX_DER_SIZE];
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char oparams[DSA_MAX_PARAMS_SIZE];
  unsigned char opub[DSA_MAX_PUB_SIZE];
  unsigned char osig[DSA_MAX_SIG_SIZE];
  unsigned char oder[DSA_MAX_DER_SIZE];
  size_t i;

  printf("Testing DSA...\n");

  for (i = 0; i < ENTROPY_SIZE; i++)
    entropy[i] = i + 1;

  for (i = 0; i < ARRAY_SIZE(dsa_vectors); i++) {
    size_t params_len = sizeof(params);
    size_t priv_len = sizeof(priv);
    size_t pub_len = sizeof(pub);
    size_t msg_len = sizeof(msg);
    size_t sig_len = sizeof(sig);
    size_t der_len = sizeof(der);
    size_t qsize, olen;

    printf("  - DSA vector #%lu\n", i + 1);

    hex_decode(params, &params_len, dsa_vectors[i][0]);
    hex_decode(priv, &priv_len, dsa_vectors[i][1]);
    hex_decode(pub, &pub_len, dsa_vectors[i][2]);
    hex_decode(msg, &msg_len, dsa_vectors[i][3]);
    hex_decode(sig, &sig_len, dsa_vectors[i][4]);
    hex_decode(der, &der_len, dsa_vectors[i][5]);

    ASSERT(dsa_params_verify(params, params_len));
    ASSERT(dsa_privkey_verify(priv, priv_len));
    ASSERT(dsa_pubkey_verify(pub, pub_len));

    ASSERT(dsa_params_create(oparams, &olen, priv, priv_len));
    ASSERT(dsa_params_verify(oparams, olen));
    ASSERT(olen == params_len);
    ASSERT(memcmp(oparams, params, params_len) == 0);

    ASSERT(dsa_params_create(oparams, &olen, pub, pub_len));
    ASSERT(dsa_params_verify(oparams, olen));
    ASSERT(olen == params_len);
    ASSERT(memcmp(oparams, params, params_len) == 0);

    ASSERT(dsa_pubkey_create(opub, &olen, priv, priv_len));
    ASSERT(dsa_pubkey_verify(opub, olen));
    ASSERT(olen == pub_len);
    ASSERT(memcmp(opub, pub, pub_len) == 0);

    qsize = dsa_pubkey_qbits(pub, pub_len) / 8;

    ASSERT(qsize != 0);

    ASSERT(dsa_sig_export(oder, &olen, sig, sig_len, qsize));
    ASSERT(olen == der_len);
    ASSERT(memcmp(oder, der, der_len) == 0);

    ASSERT(dsa_sig_export(oder, &olen, sig, sig_len, 0));
    ASSERT(olen == der_len);
    ASSERT(memcmp(oder, der, der_len) == 0);

    ASSERT(dsa_sig_import(osig, &olen, der, der_len, qsize));
    ASSERT(olen == sig_len);
    ASSERT(memcmp(osig, sig, sig_len) == 0);

    ASSERT(dsa_verify(msg, msg_len, sig, sig_len, pub, pub_len));
    ASSERT(dsa_sign(osig, &olen, msg, msg_len, priv, priv_len, entropy));
    ASSERT(dsa_verify(msg, msg_len, osig, olen, pub, pub_len));

    osig[4] ^= 1;

    ASSERT(!dsa_verify(msg, msg_len, osig, olen, pub, pub_len));

    osig[4] ^= 1;
    pub[4] ^= 1;

    ASSERT(!dsa_verify(msg, msg_len, osig, olen, pub, pub_len));

    pub[4] ^= 1;
    msg[4] ^= 1;

    ASSERT(!dsa_verify(msg, msg_len, osig, olen, pub, pub_len));

    msg[4] ^= 1;

    ASSERT(dsa_verify(msg, msg_len, osig, olen, pub, pub_len));
  }
}

static void
test_dsa_keygen(drbg_t *rng) {
  unsigned char params[DSA_MAX_PARAMS_SIZE];
  unsigned char alice_priv[DSA_MAX_PRIV_SIZE];
  unsigned char alice_pub[DSA_MAX_PUB_SIZE];
  unsigned char alice_sec[DSA_MAX_PUB_SIZE];
  unsigned char bob_priv[DSA_MAX_PRIV_SIZE];
  unsigned char bob_pub[DSA_MAX_PUB_SIZE];
  unsigned char bob_sec[DSA_MAX_PUB_SIZE];
  unsigned char entropy[ENTROPY_SIZE];
  size_t alice_priv_len, alice_pub_len, alice_sec_len;
  size_t bob_priv_len, bob_pub_len, bob_sec_len;
  size_t params_len;

  /* Params. */
  drbg_generate(rng, entropy, ENTROPY_SIZE);

  ASSERT(dsa_params_generate(params, &params_len, 1024, entropy));

  /* Alice key pair. */
  drbg_generate(rng, entropy, ENTROPY_SIZE);
  ASSERT(dsa_privkey_create(alice_priv, &alice_priv_len,
                            params, params_len, entropy));
  ASSERT(dsa_pubkey_create(alice_pub, &alice_pub_len,
                           alice_priv, alice_priv_len));

  /* Bob key pair. */
  drbg_generate(rng, entropy, ENTROPY_SIZE);
  ASSERT(dsa_privkey_create(bob_priv, &bob_priv_len,
                            params, params_len, entropy));
  ASSERT(dsa_pubkey_create(bob_pub, &bob_pub_len,
                           bob_priv, bob_priv_len));

  /* Verify. */
  ASSERT(dsa_params_verify(params, params_len));
  ASSERT(dsa_privkey_verify(alice_priv, alice_priv_len));
  ASSERT(dsa_privkey_verify(bob_priv, bob_priv_len));
  ASSERT(dsa_pubkey_verify(alice_pub, alice_pub_len));
  ASSERT(dsa_pubkey_verify(bob_pub, bob_pub_len));

  /* Diffie hellman. */
  ASSERT(dsa_derive(alice_sec, &alice_sec_len,
                    bob_pub, bob_pub_len,
                    alice_priv, alice_priv_len));
  ASSERT(dsa_derive(bob_sec, &bob_sec_len,
                    alice_pub, alice_pub_len,
                    bob_priv, bob_priv_len));

  ASSERT(alice_sec_len == bob_sec_len);
  ASSERT(memcmp(alice_sec, bob_sec, bob_sec_len) == 0);
}

/*
 * ECC
 */

#ifdef TORSION_TEST
void
__torsion_test_ecc(drbg_t *rng);
#endif

static void
test_ecdsa(void) {
  unsigned char priv[ECDSA_MAX_PRIV_SIZE];
  unsigned char pub[ECDSA_MAX_PUB_SIZE];
  unsigned char tweak[ECDSA_MAX_PRIV_SIZE];
  unsigned char privadd[ECDSA_MAX_PRIV_SIZE];
  unsigned char privmul[ECDSA_MAX_PRIV_SIZE];
  unsigned char privneg[ECDSA_MAX_PRIV_SIZE];
  unsigned char privinv[ECDSA_MAX_PRIV_SIZE];
  unsigned char pubadd[ECDSA_MAX_PUB_SIZE];
  unsigned char pubmul[ECDSA_MAX_PUB_SIZE];
  unsigned char pubneg[ECDSA_MAX_PUB_SIZE];
  unsigned char pubdbl[ECDSA_MAX_PUB_SIZE];
  unsigned char pubconv[ECDSA_MAX_PUB_SIZE];
  unsigned char pubhybrid[ECDSA_MAX_PUB_SIZE];
  unsigned char msg[128];
  unsigned char sig[ECDSA_MAX_SIG_SIZE];
  unsigned char der[ECDSA_MAX_DER_SIZE];
  unsigned char other[ECDSA_MAX_PRIV_SIZE];
  unsigned char secret[ECDSA_MAX_PUB_SIZE];
  unsigned char tweakneg[ECDSA_MAX_PRIV_SIZE];
  unsigned char tweakinv[ECDSA_MAX_PRIV_SIZE];
  unsigned char out[ECDSA_MAX_DER_SIZE];
  wei_curve_t *curves[WEI_CURVE_MAX + 1];
  const unsigned char *pubs[3];
  size_t publens[3];
  unsigned int flag;
  size_t i;

  for (i = 0; i < ARRAY_SIZE(curves); i++)
    curves[i] = wei_curve_create(i);

  printf("Testing ECDSA...\n");

  for (i = 0; i < ARRAY_SIZE(ecdsa_vectors); i++) {
    int type = ecdsa_vectors[i].type;
    wei_curve_t *ec = curves[type];
    size_t sc_size = wei_curve_scalar_size(ec);
    size_t sig_size = ecdsa_sig_size(ec);
    size_t pub_size = ecdsa_pubkey_size(ec, 1);
    size_t upub_size = ecdsa_pubkey_size(ec, 0);
    size_t msg_len = sizeof(msg);
    size_t der_len = sizeof(der);
    unsigned int param = ecdsa_vectors[i].param;
    size_t len;

    printf("  - ECDSA vector #%lu (%s)\n", i + 1, wei_curves[type]);

    hex_parse(priv, sc_size, ecdsa_vectors[i].priv);
    hex_parse(pub, pub_size, ecdsa_vectors[i].pub);
    hex_parse(tweak, sc_size, ecdsa_vectors[i].tweak);
    hex_parse(privadd, sc_size, ecdsa_vectors[i].privadd);
    hex_parse(privmul, sc_size, ecdsa_vectors[i].privmul);
    hex_parse(privneg, sc_size, ecdsa_vectors[i].privneg);
    hex_parse(privinv, sc_size, ecdsa_vectors[i].privinv);
    hex_parse(pubadd, pub_size, ecdsa_vectors[i].pubadd);
    hex_parse(pubmul, pub_size, ecdsa_vectors[i].pubmul);
    hex_parse(pubneg, pub_size, ecdsa_vectors[i].pubneg);
    hex_parse(pubdbl, pub_size, ecdsa_vectors[i].pubdbl);
    hex_parse(pubconv, upub_size, ecdsa_vectors[i].pubconv);
    hex_parse(pubhybrid, upub_size, ecdsa_vectors[i].pubhybrid);
    hex_decode(msg, &msg_len, ecdsa_vectors[i].msg);
    hex_parse(sig, sig_size, ecdsa_vectors[i].sig);
    hex_decode(der, &der_len, ecdsa_vectors[i].der);
    hex_parse(other, sc_size, ecdsa_vectors[i].other);
    hex_parse(secret, pub_size, ecdsa_vectors[i].secret);

    ASSERT(ecdsa_privkey_verify(ec, priv));
    ASSERT(ecdsa_pubkey_verify(ec, pub, pub_size));
    ASSERT(ecdsa_pubkey_verify(ec, pubconv, upub_size));
    ASSERT(ecdsa_pubkey_verify(ec, pubhybrid, upub_size));

    ASSERT(ecdsa_privkey_negate(ec, tweakneg, tweak));
    ASSERT(ecdsa_privkey_invert(ec, tweakinv, tweak));

    ASSERT(ecdsa_pubkey_create(ec, out, &len, priv, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_privkey_reduce(ec, out, priv, sc_size));
    ASSERT(memcmp(out, priv, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_add(ec, out, priv, tweak));
    ASSERT(memcmp(out, privadd, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_add(ec, out, out, tweakneg));
    ASSERT(memcmp(out, priv, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_mul(ec, out, priv, tweak));
    ASSERT(memcmp(out, privmul, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_mul(ec, out, out, tweakinv));
    ASSERT(memcmp(out, priv, sc_size) == 0);

    ASSERT(ecdsa_privkey_negate(ec, out, priv));
    ASSERT(memcmp(out, privneg, sc_size) == 0);

    ASSERT(ecdsa_privkey_invert(ec, out, priv));
    ASSERT(memcmp(out, privinv, sc_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_add(ec, out, &len, pub, pub_size, tweak, 1));
    ASSERT(len == pub_size && memcmp(out, pubadd, pub_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_add(ec, out, &len, pubadd, pub_size, tweakneg, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_mul(ec, out, &len, pub, pub_size, tweak, 1));
    ASSERT(len == pub_size && memcmp(out, pubmul, pub_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_mul(ec, out, &len, pubmul, pub_size, tweakinv, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_pubkey_negate(ec, out, &len, pub, pub_size, 1));
    ASSERT(len == pub_size && memcmp(out, pubneg, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pub;

    publens[0] = pub_size;
    publens[1] = pub_size;

    ASSERT(ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 2, 1));
    ASSERT(len == pub_size && memcmp(out, pubdbl, pub_size) == 0);

    pubs[0] = pubdbl;
    pubs[1] = pubneg;

    publens[0] = pub_size;
    publens[1] = pub_size;

    ASSERT(ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 2, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pubneg;
    pubs[2] = pubconv;

    publens[0] = pub_size;
    publens[1] = pub_size;
    publens[2] = upub_size;

    ASSERT(ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 3, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    ASSERT(!ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 2, 1));

    ASSERT(ecdsa_pubkey_create(ec, out, &len, priv, 0));
    ASSERT(len == upub_size && memcmp(out, pubconv, upub_size) == 0);

    ASSERT(ecdsa_pubkey_convert(ec, out, &len, pub, pub_size, 0));
    ASSERT(len == upub_size && memcmp(out, pubconv, upub_size) == 0);

    ASSERT(ecdsa_pubkey_convert(ec, out, &len, pubconv, upub_size, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_is_low_s(ec, sig));

    ASSERT(ecdsa_sig_export(ec, out, &len, sig));
    ASSERT(len == der_len && memcmp(out, der, der_len) == 0);

    ASSERT(ecdsa_sig_import(ec, out, der, der_len));
    ASSERT(memcmp(out, sig, sig_size) == 0);

    ASSERT(ecdsa_recover(ec, out, &len, msg, msg_len, sig, param, 1));
    ASSERT(len == pub_size && memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_recover(ec, out, &len, msg, msg_len, sig, param, 0));
    ASSERT(len == upub_size && memcmp(out, pubconv, upub_size) == 0);

    ASSERT(ecdsa_derive(ec, out, &len, pub, pub_size, other, 1));
    ASSERT(len == pub_size && memcmp(out, secret, pub_size) == 0);

    ASSERT(ecdsa_derive(ec, out, &len, pubconv, upub_size, other, 1));
    ASSERT(len == pub_size && memcmp(out, secret, pub_size) == 0);

    ASSERT(ecdsa_derive(ec, out, &len, pubhybrid, upub_size, other, 1));
    ASSERT(len == pub_size && memcmp(out, secret, pub_size) == 0);

    ASSERT(ecdsa_sign(ec, out, &flag, msg, msg_len, priv));
    ASSERT(memcmp(out, sig, sig_size) == 0);
    ASSERT(flag == param);

    ASSERT(ecdsa_verify(ec, msg, msg_len, sig, pub, pub_size));
    ASSERT(ecdsa_verify(ec, msg, msg_len, sig, pubconv, upub_size));
    ASSERT(ecdsa_verify(ec, msg, msg_len, sig, pubhybrid, upub_size));

    msg[2] ^= 1;

    ASSERT(!ecdsa_verify(ec, msg, msg_len, sig, pub, pub_size));

    msg[2] ^= 1;
    sig[2] ^= 1;

    ASSERT(!ecdsa_verify(ec, msg, msg_len, sig, pub, pub_size));

    sig[2] ^= 1;
    pub[2] ^= 1;

    ASSERT(!ecdsa_verify(ec, msg, msg_len, sig, pub, pub_size));

    pub[2] ^= 1;

    ASSERT(ecdsa_verify(ec, msg, msg_len, sig, pub, pub_size));
  }

  for (i = 0; i < ARRAY_SIZE(curves); i++)
    wei_curve_destroy(curves[i]);
}

static void
test_ecdsa_random(drbg_t *rng) {
  size_t i, j;

  printf("Randomized ECDSA testing...\n");

  for (i = 0; i < ARRAY_SIZE(wei_curves); i++) {
    const char *id = wei_curves[i];
    wei_curve_t *ec = wei_curve_create(i);
    size_t fe_size = wei_curve_field_size(ec);
    size_t fe_bits = wei_curve_field_bits(ec);
    size_t sc_size = wei_curve_scalar_size(ec);

    printf("  - %s\n", id);

    for (j = 0; j < 10; j++) {
      unsigned char entropy[ENTROPY_SIZE];
      unsigned char priv[ECDSA_MAX_PRIV_SIZE];
      unsigned char msg[WEI_MAX_SCALAR_SIZE];
      unsigned char sig[ECDSA_MAX_SIG_SIZE];
      unsigned char pub[ECDSA_MAX_PUB_SIZE];
      unsigned char rec[ECDSA_MAX_PUB_SIZE];
      size_t pub_len, rec_len;
      unsigned int param;
      size_t k;

      drbg_generate(rng, entropy, sizeof(entropy));
      drbg_generate(rng, priv, sizeof(priv));
      drbg_generate(rng, msg, sizeof(msg));

      priv[0] = 0;

      wei_curve_randomize(ec, entropy);

      ASSERT(ecdsa_sign(ec, sig, &param, msg, sc_size, priv));
      ASSERT(ecdsa_pubkey_create(ec, pub, &pub_len, priv, 1));
      ASSERT(pub_len == fe_size + 1);
      ASSERT(ecdsa_verify(ec, msg, sc_size, sig, pub, pub_len));
      ASSERT(ecdsa_recover(ec, rec, &rec_len, msg, sc_size, sig, param, 1));
      ASSERT(rec_len == pub_len);
      ASSERT(memcmp(pub, rec, pub_len) == 0);

      k = drbg_uniform(rng, sc_size);

      msg[k] ^= 1;

      if (fe_bits != 521 || k < 65)
        ASSERT(!ecdsa_verify(ec, msg, sc_size, sig, pub, pub_len));

      msg[k] ^= 1;
      pub[k] ^= 1;

      ASSERT(!ecdsa_verify(ec, msg, sc_size, sig, pub, pub_len));

      pub[k] ^= 1;
      sig[k] ^= 1;

      ASSERT(!ecdsa_verify(ec, msg, sc_size, sig, pub, pub_len));

      sig[k] ^= 1;
      sig[sc_size + k] ^= 1;

      ASSERT(!ecdsa_verify(ec, msg, sc_size, sig, pub, pub_len));

      sig[sc_size + k] ^= 1;

      ASSERT(ecdsa_verify(ec, msg, sc_size, sig, pub, pub_len));
    }

    wei_curve_destroy(ec);
  }
}

static void
test_ecdsa_sswu(void) {
  const unsigned char bytes[32] = {
    0x6a, 0x53, 0xb3, 0x44, 0xc1, 0x88, 0x8a, 0x62,
    0xf2, 0x81, 0xaf, 0xbd, 0xbe, 0x67, 0x6b, 0xc3,
    0x9e, 0xc6, 0xb1, 0x46, 0x6e, 0x1d, 0x97, 0x53,
    0x45, 0x01, 0xda, 0xf1, 0x63, 0xff, 0xc9, 0x98
  };

  const unsigned char expect[33] = {
    0x02, 0x61, 0x92, 0xb9, 0xcd, 0x61, 0x6b, 0x4c,
    0x06, 0x83, 0xdd, 0xc9, 0x40, 0x16, 0x1d, 0x60,
    0x8f, 0x68, 0xd3, 0xbd, 0xa9, 0x1f, 0x86, 0x7b,
    0x1f, 0x65, 0x34, 0xdf, 0xdd, 0x04, 0xbd, 0xb2,
    0x8a
  };

  wei_curve_t *ec = wei_curve_create(WEI_CURVE_P256);
  unsigned char out[33];
  size_t out_len;

  printf("Testing SSWU (P256).\n");

  ecdsa_pubkey_from_uniform(ec, out, &out_len, bytes, 1);

  ASSERT(out_len == 33);
  ASSERT(memcmp(out, expect, 33) == 0);
  ASSERT(ecdsa_pubkey_to_uniform(ec, out, expect, 33, 3));
  ASSERT(memcmp(out, bytes, 32) == 0);

  wei_curve_destroy(ec);
}

static void
test_ecdsa_svdw(void) {
  const unsigned char bytes[32] = {
    0xb0, 0xf0, 0xa9, 0x2d, 0x14, 0xa9, 0x82, 0xeb,
    0x12, 0x04, 0x78, 0x1a, 0x91, 0x6f, 0xdf, 0x38,
    0x2c, 0x4d, 0x84, 0x69, 0x38, 0xe6, 0x3f, 0x55,
    0xca, 0x59, 0x22, 0xb1, 0x0a, 0xb6, 0x82, 0xa0
  };

  const unsigned char expect[33] = {
    0x02, 0xa3, 0xb0, 0xbc, 0xa2, 0xaa, 0x06, 0xe3,
    0x78, 0x83, 0x14, 0xb8, 0x73, 0x54, 0xbd, 0x01,
    0x04, 0xf1, 0x10, 0x85, 0xa8, 0x67, 0xab, 0xeb,
    0x4f, 0x43, 0xd2, 0xf6, 0x22, 0xdb, 0xb3, 0x29,
    0x20
  };

  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char out[33];
  size_t out_len;

  printf("Testing SVDW (secp256k1).\n");

  ecdsa_pubkey_from_uniform(ec, out, &out_len, bytes, 1);

  ASSERT(out_len == 33);
  ASSERT(memcmp(out, expect, 33) == 0);
  ASSERT(ecdsa_pubkey_to_uniform(ec, out, expect, 33, 1));
  ASSERT(memcmp(out, bytes, 32) == 0);

  wei_curve_destroy(ec);
}

static void
test_schnorr_legacy(void) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  wei_scratch_t *scratch = wei_scratch_create(ec, 10);
  unsigned char priv[32];
  unsigned char pub[65];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char out[64];
  const unsigned char *msgs[1];
  const unsigned char *pubs[1];
  const unsigned char *sigs[1];
  size_t i, len;

  printf("Testing Schnorr-Legacy...\n");

  for (i = 0; i < ARRAY_SIZE(schnorr_legacy_vectors); i++) {
    size_t priv_len = sizeof(priv);
    size_t pub_len = sizeof(pub);
    size_t msg_len = sizeof(msg);
    size_t sig_len = sizeof(sig);
    int result = schnorr_legacy_vectors[i].result;
    const char *comment = schnorr_legacy_vectors[i].comment;

    printf("  - Schnorr-Legacy vector #%lu (%s)\n", i + 1, comment);

    hex_decode(priv, &priv_len, schnorr_legacy_vectors[i].priv);
    hex_decode(pub, &pub_len, schnorr_legacy_vectors[i].pub);
    hex_decode(msg, &msg_len, schnorr_legacy_vectors[i].msg);
    hex_decode(sig, &sig_len, schnorr_legacy_vectors[i].sig);

    ASSERT(priv_len == 0 || priv_len == 32);
    ASSERT(pub_len == 33);
    ASSERT(msg_len == 32);
    ASSERT(sig_len == 0 || sig_len == 64);

    if (sig_len == 0)
      memset(sig, 0, 64);

    if (priv_len > 0) {
      ASSERT(schnorr_legacy_privkey_verify(ec, priv));
      ASSERT(schnorr_legacy_pubkey_create(ec, out, &len, priv, 1));
      ASSERT(len == pub_len && memcmp(out, pub, pub_len) == 0);
      ASSERT(schnorr_legacy_sign(ec, out, msg, 32, priv));
      ASSERT(memcmp(out, sig, sig_len) == 0);
    }

    ASSERT(schnorr_legacy_verify(ec, msg, msg_len,
                                 sig, pub, pub_len) == result);

    msgs[0] = msg;
    sigs[0] = sig;
    pubs[0] = pub;

    ASSERT(schnorr_legacy_verify_batch(ec, msgs, &msg_len, sigs,
                                       pubs, &pub_len, 1, scratch) == result);
  }

  wei_scratch_destroy(ec, scratch);
  wei_curve_destroy(ec);
}

static void
test_schnorr_legacy_random(drbg_t *rng) {
  size_t i, j;

  printf("Randomized Schnorr-Legacy testing...\n");

  for (i = 0; i < ARRAY_SIZE(wei_curves); i++) {
    const char *id = wei_curves[i];
    wei_curve_t *ec;

    if (i == WEI_CURVE_P224)
      continue;

    ec = wei_curve_create(i);

    printf("  - %s\n", id);

    for (j = 0; j < 10; j++) {
      unsigned char entropy[ENTROPY_SIZE];
      unsigned char priv[SCHNORR_LEGACY_MAX_PRIV_SIZE];
      unsigned char sig[SCHNORR_LEGACY_MAX_SIG_SIZE];
      unsigned char pub[SCHNORR_LEGACY_MAX_PUB_SIZE];
      unsigned char msg[32];
      size_t pub_len;

      drbg_generate(rng, entropy, sizeof(entropy));
      drbg_generate(rng, priv, sizeof(priv));
      drbg_generate(rng, msg, sizeof(msg));

      priv[0] = 0;

      wei_curve_randomize(ec, entropy);

      ASSERT(schnorr_legacy_sign(ec, sig, msg, 32, priv));
      ASSERT(schnorr_legacy_pubkey_create(ec, pub, &pub_len, priv, 1));
      ASSERT(schnorr_legacy_verify(ec, msg, 32, sig, pub, pub_len));

      msg[0] ^= 1;

      ASSERT(!schnorr_legacy_verify(ec, msg, 32, sig, pub, pub_len));

      msg[0] ^= 1;
      pub[1] ^= 1;

      ASSERT(!schnorr_legacy_verify(ec, msg, 32, sig, pub, pub_len));

      pub[1] ^= 1;
      sig[0] ^= 1;

      ASSERT(!schnorr_legacy_verify(ec, msg, 32, sig, pub, pub_len));

      sig[0] ^= 1;

      ASSERT(schnorr_legacy_verify(ec, msg, 32, sig, pub, pub_len));
    }

    wei_curve_destroy(ec);
  }
}

static void
test_schnorr(void) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  wei_scratch_t *scratch = wei_scratch_create(ec, 10);
  unsigned char priv[32];
  unsigned char pub[32];
  unsigned char aux[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char out[64];
  const unsigned char *msgs[1];
  const unsigned char *pubs[1];
  const unsigned char *sigs[1];
  size_t i;

  printf("Testing Schnorr...\n");

  for (i = 0; i < ARRAY_SIZE(schnorr_vectors); i++) {
    size_t priv_len = sizeof(priv);
    size_t pub_len = sizeof(pub);
    size_t aux_len = sizeof(aux);
    size_t msg_len = sizeof(msg);
    size_t sig_len = sizeof(sig);
    int result = schnorr_vectors[i].result;
    const char *comment = schnorr_vectors[i].comment;

    printf("  - Schnorr vector #%lu (%s)\n", i + 1, comment);

    hex_decode(priv, &priv_len, schnorr_vectors[i].priv);
    hex_decode(pub, &pub_len, schnorr_vectors[i].pub);
    hex_decode(aux, &aux_len, schnorr_vectors[i].aux);
    hex_decode(msg, &msg_len, schnorr_vectors[i].msg);
    hex_decode(sig, &sig_len, schnorr_vectors[i].sig);

    ASSERT(priv_len == 0 || priv_len == 32);
    ASSERT(pub_len == 32);
    ASSERT(aux_len == 0 || aux_len == 32);
    ASSERT(msg_len == 32);
    ASSERT(sig_len == 64);

    if (aux_len == 0)
      memset(aux, 0, 32);

    if (priv_len > 0) {
      ASSERT(schnorr_privkey_verify(ec, priv));
      ASSERT(schnorr_pubkey_create(ec, out, priv));
      ASSERT(memcmp(out, pub, pub_len) == 0);
      ASSERT(schnorr_sign(ec, out, msg, 32, priv, aux));
      ASSERT(memcmp(out, sig, sig_len) == 0);
    }

    ASSERT(schnorr_verify(ec, msg, msg_len, sig, pub) == result);

    msgs[0] = msg;
    sigs[0] = sig;
    pubs[0] = pub;

    ASSERT(schnorr_verify_batch(ec, msgs, &msg_len, sigs,
                                pubs, 1, scratch) == result);
  }

  wei_scratch_destroy(ec, scratch);
  wei_curve_destroy(ec);
}

static void
test_schnorr_random(drbg_t *rng) {
  size_t i, j;

  printf("Randomized Schnorr-Legacy testing...\n");

  for (i = 0; i < ARRAY_SIZE(wei_curves); i++) {
    const char *id = wei_curves[i];
    wei_curve_t *ec;

    if (i == WEI_CURVE_P224)
      continue;

    ec = wei_curve_create(i);

    printf("  - %s\n", id);

    for (j = 0; j < 10; j++) {
      unsigned char entropy[ENTROPY_SIZE];
      unsigned char priv[SCHNORR_MAX_PRIV_SIZE];
      unsigned char sig[SCHNORR_MAX_SIG_SIZE];
      unsigned char pub[SCHNORR_MAX_PUB_SIZE];
      unsigned char msg[32];
      unsigned char aux[32];

      drbg_generate(rng, entropy, sizeof(entropy));
      drbg_generate(rng, priv, sizeof(priv));
      drbg_generate(rng, msg, sizeof(msg));
      drbg_generate(rng, aux, sizeof(aux));

      priv[0] = 0;

      wei_curve_randomize(ec, entropy);

      ASSERT(schnorr_sign(ec, sig, msg, 32, priv, aux));
      ASSERT(schnorr_pubkey_create(ec, pub, priv));
      ASSERT(schnorr_verify(ec, msg, 32, sig, pub));

      msg[0] ^= 1;

      ASSERT(!schnorr_verify(ec, msg, 32, sig, pub));

      msg[0] ^= 1;
      pub[1] ^= 1;

      ASSERT(!schnorr_verify(ec, msg, 32, sig, pub));

      pub[1] ^= 1;
      sig[0] ^= 1;

      ASSERT(!schnorr_verify(ec, msg, 32, sig, pub));

      sig[0] ^= 1;

      ASSERT(schnorr_verify(ec, msg, 32, sig, pub));
    }

    wei_curve_destroy(ec);
  }
}

static void
test_ecdh_x25519(void) {
  /* From RFC 7748 */
  const unsigned char intervals[3][32] = {
    {
      0x42, 0x2c, 0x8e, 0x7a, 0x62, 0x27, 0xd7, 0xbc,
      0xa1, 0x35, 0x0b, 0x3e, 0x2b, 0xb7, 0x27, 0x9f,
      0x78, 0x97, 0xb8, 0x7b, 0xb6, 0x85, 0x4b, 0x78,
      0x3c, 0x60, 0xe8, 0x03, 0x11, 0xae, 0x30, 0x79
    },
    {
      0x68, 0x4c, 0xf5, 0x9b, 0xa8, 0x33, 0x09, 0x55,
      0x28, 0x00, 0xef, 0x56, 0x6f, 0x2f, 0x4d, 0x3c,
      0x1c, 0x38, 0x87, 0xc4, 0x93, 0x60, 0xe3, 0x87,
      0x5f, 0x2e, 0xb9, 0x4d, 0x99, 0x53, 0x2c, 0x51
    },
    {
      0x7c, 0x39, 0x11, 0xe0, 0xab, 0x25, 0x86, 0xfd,
      0x86, 0x44, 0x97, 0x29, 0x7e, 0x57, 0x5e, 0x6f,
      0x3b, 0xc6, 0x01, 0xc0, 0x88, 0x3c, 0x30, 0xdf,
      0x5f, 0x4d, 0xd2, 0xd2, 0x4f, 0x66, 0x54, 0x24
    }
  };

  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char k[32];
  unsigned char u[32];
  unsigned char t[32];
  size_t i = 0;

  printf("Testing X25519 vectors...\n");

  memset(k, 0x00, 32);
  memset(u, 0x00, 32);
  memset(t, 0x00, 32);

  k[0] = 9;
  u[0] = 9;

  for (; i < 1; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 32);
    memcpy(k, t, 32);
  }

  ASSERT(memcmp(k, intervals[0], 32) == 0);

  for (; i < 1000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 32);
    memcpy(k, t, 32);
  }

  ASSERT(memcmp(k, intervals[1], 32) == 0);

#ifdef TORSION_TEST_SLOW
  for (; i < 1000000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 32);
    memcpy(k, t, 32);
  }

  ASSERT(memcmp(k, intervals[2], 32) == 0);
#endif

  ASSERT(ecdh_pubkey_verify(ec, t));

  mont_curve_destroy(ec);
}

static void
test_ecdh_x448(void) {
  /* From RFC 7748 */
  const unsigned char intervals[4][56] = {
    {
      0x3f, 0x48, 0x2c, 0x8a, 0x9f, 0x19, 0xb0, 0x1e,
      0x6c, 0x46, 0xee, 0x97, 0x11, 0xd9, 0xdc, 0x14,
      0xfd, 0x4b, 0xf6, 0x7a, 0xf3, 0x07, 0x65, 0xc2,
      0xae, 0x2b, 0x84, 0x6a, 0x4d, 0x23, 0xa8, 0xcd,
      0x0d, 0xb8, 0x97, 0x08, 0x62, 0x39, 0x49, 0x2c,
      0xaf, 0x35, 0x0b, 0x51, 0xf8, 0x33, 0x86, 0x8b,
      0x9b, 0xc2, 0xb3, 0xbc, 0xa9, 0xcf, 0x41, 0x13
    },
    {
      0xcc, 0xa0, 0x3d, 0x8e, 0xd3, 0xf5, 0x4b, 0xaf,
      0x8d, 0x1a, 0xa0, 0x88, 0xb1, 0xf2, 0x4b, 0xc6,
      0x8a, 0xed, 0x53, 0x8d, 0x06, 0x48, 0x5f, 0x02,
      0x5f, 0x17, 0xa5, 0x43, 0x1d, 0xed, 0x28, 0xf2,
      0x56, 0xd3, 0x4f, 0x6b, 0xdd, 0x3d, 0x63, 0xcc,
      0x5e, 0x04, 0x7c, 0x45, 0x8e, 0x81, 0x38, 0x55,
      0x19, 0xa9, 0x29, 0x99, 0xbd, 0xdc, 0x26, 0x53
    },
    {
      0xaa, 0x3b, 0x47, 0x49, 0xd5, 0x5b, 0x9d, 0xaf,
      0x1e, 0x5b, 0x00, 0x28, 0x88, 0x26, 0xc4, 0x67,
      0x27, 0x4c, 0xe3, 0xeb, 0xbd, 0xd5, 0xc1, 0x7b,
      0x97, 0x5e, 0x09, 0xd4, 0xaf, 0x6c, 0x67, 0xcf,
      0x10, 0xd0, 0x87, 0x20, 0x2d, 0xb8, 0x82, 0x86,
      0xe2, 0xb7, 0x9f, 0xce, 0xea, 0x3e, 0xc3, 0x53,
      0xef, 0x54, 0xfa, 0xa2, 0x6e, 0x21, 0x9f, 0x38
    },
    {
      0x07, 0x7f, 0x45, 0x36, 0x81, 0xca, 0xca, 0x36,
      0x93, 0x19, 0x84, 0x20, 0xbb, 0xe5, 0x15, 0xca,
      0xe0, 0x00, 0x24, 0x72, 0x51, 0x9b, 0x3e, 0x67,
      0x66, 0x1a, 0x7e, 0x89, 0xca, 0xb9, 0x46, 0x95,
      0xc8, 0xf4, 0xbc, 0xd6, 0x6e, 0x61, 0xb9, 0xb9,
      0xc9, 0x46, 0xda, 0x8d, 0x52, 0x4d, 0xe3, 0xd6,
      0x9b, 0xd9, 0xd9, 0xd6, 0x6b, 0x99, 0x7e, 0x37
    }
  };

  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X448);
  unsigned char k[56];
  unsigned char u[56];
  unsigned char t[56];
  size_t i = 0;

  printf("Testing X448 vectors...\n");

  memset(k, 0x00, 56);
  memset(u, 0x00, 56);
  memset(t, 0x00, 56);

  k[0] = 5;
  u[0] = 5;

  for (; i < 1; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(memcmp(k, intervals[0], 56) == 0);

  for (; i < 100; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(memcmp(k, intervals[1], 56) == 0);

  for (; i < 1000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(memcmp(k, intervals[2], 56) == 0);

#ifdef TORSION_TEST_SLOW
  for (; i < 1000000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(memcmp(k, intervals[3], 56) == 0);
#endif

  ASSERT(ecdh_pubkey_verify(ec, t));

  mont_curve_destroy(ec);
}

static void
test_ecdh_random(drbg_t *rng) {
  size_t i, j;

  printf("Randomized ECDH testing...\n");

  for (i = 0; i < ARRAY_SIZE(mont_curves); i++) {
    const char *id = mont_curves[i];
    mont_curve_t *ec = mont_curve_create(i);
    size_t fe_size = mont_curve_field_size(ec);

    printf("  - %s\n", id);

    for (j = 0; j < 10; j++) {
      unsigned char alice_priv[ECDH_MAX_PRIV_SIZE];
      unsigned char alice_pub[ECDH_MAX_PUB_SIZE];
      unsigned char alice_secret[ECDH_MAX_PUB_SIZE];
      unsigned char bob_priv[ECDH_MAX_PRIV_SIZE];
      unsigned char bob_pub[ECDH_MAX_PUB_SIZE];
      unsigned char bob_secret[ECDH_MAX_PUB_SIZE];

      drbg_generate(rng, alice_priv, sizeof(alice_priv));
      drbg_generate(rng, bob_priv, sizeof(bob_priv));

      ecdh_pubkey_create(ec, alice_pub, alice_priv);
      ecdh_pubkey_create(ec, bob_pub, bob_priv);

      ASSERT(ecdh_pubkey_verify(ec, alice_pub));
      ASSERT(ecdh_pubkey_verify(ec, bob_pub));

      ASSERT(ecdh_derive(ec, alice_secret, bob_pub, alice_priv));
      ASSERT(ecdh_derive(ec, bob_secret, alice_pub, bob_priv));

      ASSERT(memcmp(alice_secret, bob_secret, fe_size) == 0);
    }

    mont_curve_destroy(ec);
  }
}

static void
test_ecdh_elligator2(void) {
  const unsigned char bytes[32] = {
    0xf9, 0xd4, 0x08, 0xba, 0x8a, 0x4e, 0x6a, 0x04,
    0xb9, 0xeb, 0x8d, 0x10, 0x38, 0x9b, 0xc0, 0x13,
    0x12, 0x9f, 0x74, 0xed, 0xb7, 0x66, 0xa1, 0x96,
    0xd6, 0x15, 0x16, 0x2b, 0x62, 0xa1, 0xe6, 0x24
  };

  const unsigned char expect[32] = {
    0xe2, 0x0d, 0xc1, 0xb4, 0xd5, 0xd3, 0x27, 0xbe,
    0x28, 0xdd, 0x80, 0x3a, 0x91, 0xdc, 0x94, 0xa2,
    0xfc, 0xfa, 0x0b, 0x79, 0xe7, 0xc9, 0xe9, 0x09,
    0x52, 0x2a, 0x2f, 0xff, 0x35, 0x16, 0x0e, 0x04
  };

  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char out[32];

  printf("Testing Elligator 2 (x25519).\n");

  ecdh_pubkey_from_uniform(ec, out, bytes);

  ASSERT(memcmp(out, expect, 32) == 0);
  ASSERT(ecdh_pubkey_to_uniform(ec, out, expect, 0));
  ASSERT(memcmp(out, bytes, 32) == 0);

  mont_curve_destroy(ec);
}

static void
test_eddsa(void) {
  unsigned char priv[EDDSA_MAX_PRIV_SIZE];
  unsigned char scalar[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char prefix[EDDSA_MAX_PREFIX_SIZE];
  unsigned char reduced[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char pub[EDDSA_MAX_PUB_SIZE];
  unsigned char tweak[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char privadd[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char privmul[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char privneg[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char privinv[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char pubadd[EDDSA_MAX_PUB_SIZE];
  unsigned char pubmul[EDDSA_MAX_PUB_SIZE];
  unsigned char pubneg[EDDSA_MAX_PUB_SIZE];
  unsigned char pubdbl[EDDSA_MAX_PUB_SIZE];
  unsigned char pubconv[ECDH_MAX_PUB_SIZE];
  unsigned char msg[128];
  unsigned char sig[EDDSA_MAX_SIG_SIZE];
  unsigned char sigadd[EDDSA_MAX_SIG_SIZE];
  unsigned char sigmul[EDDSA_MAX_SIG_SIZE];
  unsigned char other[EDDSA_MAX_PRIV_SIZE];
  unsigned char secret[EDDSA_MAX_PUB_SIZE];
  unsigned char scalar_[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char prefix_[EDDSA_MAX_PREFIX_SIZE];
  unsigned char tweakneg[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char tweakinv[EDWARDS_MAX_SCALAR_SIZE];
  unsigned char inf[EDDSA_MAX_PUB_SIZE];
  unsigned char out[EDDSA_MAX_SIG_SIZE];
  edwards_curve_t *curves[EDWARDS_CURVE_MAX + 1];
  edwards_scratch_t *scratches[EDWARDS_CURVE_MAX + 1];
  const unsigned char *pubs[3];
  const unsigned char *msgs[1];
  const unsigned char *sigs[1];
  size_t i;

  printf("Testing EdDSA...\n");

  for (i = 0; i < ARRAY_SIZE(curves); i++) {
    curves[i] = edwards_curve_create(i);
    scratches[i] = edwards_scratch_create(curves[i], 10);
  }

  memset(inf, 0, sizeof(inf));

  inf[0] = 0x01;

  for (i = 0; i < ARRAY_SIZE(eddsa_vectors); i++) {
    int type = eddsa_vectors[i].type;
    edwards_curve_t *ec = curves[type];
    edwards_scratch_t *scratch = scratches[type];
    size_t fe_size = edwards_curve_field_size(ec);
    size_t sc_size = edwards_curve_scalar_size(ec);
    size_t priv_size = eddsa_privkey_size(ec);
    size_t pub_size = eddsa_pubkey_size(ec);
    size_t sig_size = eddsa_sig_size(ec);
    size_t msg_len = sizeof(msg);
    int ph = eddsa_vectors[i].ph;

    printf("  - EdDSA vector #%lu (%s)\n", i + 1, edwards_curves[type]);

    hex_parse(priv, priv_size, eddsa_vectors[i].priv);
    hex_parse(scalar, sc_size, eddsa_vectors[i].scalar);
    hex_parse(prefix, pub_size, eddsa_vectors[i].prefix);
    hex_parse(reduced, sc_size, eddsa_vectors[i].reduced);
    hex_parse(pub, pub_size, eddsa_vectors[i].pub);
    hex_parse(tweak, sc_size, eddsa_vectors[i].tweak);
    hex_parse(privadd, sc_size, eddsa_vectors[i].privadd);
    hex_parse(privmul, sc_size, eddsa_vectors[i].privmul);
    hex_parse(privneg, sc_size, eddsa_vectors[i].privneg);
    hex_parse(privinv, sc_size, eddsa_vectors[i].privinv);
    hex_parse(pubadd, pub_size, eddsa_vectors[i].pubadd);
    hex_parse(pubmul, pub_size, eddsa_vectors[i].pubmul);
    hex_parse(pubneg, pub_size, eddsa_vectors[i].pubneg);
    hex_parse(pubdbl, pub_size, eddsa_vectors[i].pubdbl);
    hex_parse(pubconv, fe_size, eddsa_vectors[i].pubconv);
    hex_decode(msg, &msg_len, eddsa_vectors[i].msg);
    hex_parse(sig, sig_size, eddsa_vectors[i].sig);
    hex_parse(sigadd, sig_size, eddsa_vectors[i].sigadd);
    hex_parse(sigmul, sig_size, eddsa_vectors[i].sigmul);
    hex_parse(other, priv_size, eddsa_vectors[i].other);
    hex_parse(secret, pub_size, eddsa_vectors[i].secret);

    ASSERT(eddsa_privkey_verify(ec, priv));
    ASSERT(eddsa_pubkey_verify(ec, pub));

    eddsa_scalar_negate(ec, tweakneg, tweak);
    eddsa_scalar_invert(ec, tweakinv, tweak);

    eddsa_pubkey_create(ec, out, priv);
    ASSERT(memcmp(out, pub, pub_size) == 0);

    eddsa_pubkey_from_scalar(ec, out, scalar);
    ASSERT(memcmp(out, pub, pub_size) == 0);

    eddsa_privkey_expand(ec, scalar_, prefix_, priv);
    ASSERT(memcmp(scalar_, scalar, sc_size) == 0);
    ASSERT(memcmp(prefix_, prefix, pub_size) == 0);

    eddsa_scalar_reduce(ec, out, scalar, sc_size);
    ASSERT(memcmp(out, reduced, sc_size) == 0);

    eddsa_scalar_tweak_add(ec, out, scalar, tweak);
    ASSERT(memcmp(out, privadd, sc_size) == 0);

    eddsa_scalar_tweak_add(ec, out, privadd, tweakneg);
    ASSERT(memcmp(out, reduced, sc_size) == 0);

    eddsa_scalar_tweak_mul(ec, out, scalar, tweak);
    ASSERT(memcmp(out, privmul, sc_size) == 0);

    eddsa_scalar_tweak_mul(ec, out, privmul, tweakinv);
    ASSERT(memcmp(out, reduced, sc_size) == 0);

    eddsa_scalar_negate(ec, out, scalar);
    ASSERT(memcmp(out, privneg, sc_size) == 0);

    eddsa_scalar_invert(ec, out, scalar);
    ASSERT(memcmp(out, privinv, sc_size) == 0);

    ASSERT(eddsa_pubkey_tweak_add(ec, out, pub, tweak));
    ASSERT(memcmp(out, pubadd, pub_size) == 0);

    ASSERT(eddsa_pubkey_tweak_add(ec, out, pubadd, tweakneg));
    ASSERT(memcmp(out, pub, pub_size) == 0);

    ASSERT(eddsa_pubkey_tweak_mul(ec, out, pub, tweak));
    ASSERT(memcmp(out, pubmul, pub_size) == 0);

    ASSERT(eddsa_pubkey_tweak_mul(ec, out, pubmul, tweakinv));
    ASSERT(memcmp(out, pub, pub_size) == 0);

    ASSERT(eddsa_pubkey_negate(ec, out, pub));
    ASSERT(memcmp(out, pubneg, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pub;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 2));
    ASSERT(memcmp(out, pubdbl, pub_size) == 0);

    pubs[0] = pubdbl;
    pubs[1] = pubneg;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 2));
    ASSERT(memcmp(out, pub, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pubneg;
    pubs[2] = pub;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 3));
    ASSERT(memcmp(out, pub, pub_size) == 0);

    eddsa_privkey_convert(ec, out, priv);
    ASSERT(memcmp(out, scalar, sc_size) == 0);

    ASSERT(eddsa_pubkey_convert(ec, out, pub));
    ASSERT(memcmp(out, pubconv, fe_size) == 0);

    ASSERT(eddsa_pubkey_combine(ec, out, NULL, 0));
    ASSERT(memcmp(out, inf, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pubneg;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 2));
    ASSERT(memcmp(out, inf, pub_size) == 0);

    ASSERT(eddsa_derive(ec, out, pub, other));
    ASSERT(memcmp(out, secret, sc_size) == 0);

    eddsa_privkey_convert(ec, out, other);
    ASSERT(eddsa_derive_with_scalar(ec, out, pub, out));
    ASSERT(memcmp(out, secret, sc_size) == 0);

    eddsa_sign(ec, out, msg, msg_len, priv, ph, NULL, 0);
    ASSERT(memcmp(out, sig, sig_size) == 0);

    eddsa_sign_with_scalar(ec, out, msg, msg_len, scalar, prefix, ph, NULL, 0);
    ASSERT(memcmp(out, sig, sig_size) == 0);

    msgs[0] = msg;
    pubs[0] = pub;
    sigs[0] = sig;

    ASSERT(eddsa_verify(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(eddsa_verify_single(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(eddsa_verify_batch(ec, msgs, &msg_len, sigs, pubs,
                              1, ph, NULL, 0, scratch));

    msg[0] ^= 1;

    ASSERT(!eddsa_verify(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(!eddsa_verify_single(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(!eddsa_verify_batch(ec, msgs, &msg_len, sigs, pubs,
                               1, ph, NULL, 0, scratch));

    msg[0] ^= 1;
    sig[0] ^= 1;

    ASSERT(!eddsa_verify(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(!eddsa_verify_single(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(!eddsa_verify_batch(ec, msgs, &msg_len, sigs, pubs,
                               1, ph, NULL, 0, scratch));

    sig[0] ^= 1;
    pub[0] ^= 1;

    ASSERT(!eddsa_verify(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(!eddsa_verify_single(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(!eddsa_verify_batch(ec, msgs, &msg_len, sigs, pubs,
                               1, ph, NULL, 0, scratch));

    pub[0] ^= 1;

    ASSERT(eddsa_verify(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(eddsa_verify_single(ec, msg, msg_len, sig, pub, ph, NULL, 0));
    ASSERT(eddsa_verify_batch(ec, msgs, &msg_len, sigs, pubs,
                              1, ph, NULL, 0, scratch));

    eddsa_sign_tweak_add(ec, out, msg, msg_len, priv, tweak, ph, NULL, 0);
    ASSERT(memcmp(out, sigadd, sig_size) == 0);

    ASSERT(eddsa_verify(ec, msg, msg_len, sigadd, pubadd, ph, NULL, 0));
    ASSERT(eddsa_verify_single(ec, msg, msg_len, sigadd, pubadd, ph, NULL, 0));

    eddsa_sign_tweak_mul(ec, out, msg, msg_len, priv, tweak, ph, NULL, 0);
    ASSERT(memcmp(out, sigmul, sig_size) == 0);

    ASSERT(eddsa_verify(ec, msg, msg_len, sigmul, pubmul, ph, NULL, 0));
    ASSERT(eddsa_verify_single(ec, msg, msg_len, sigmul, pubmul, ph, NULL, 0));
  }

  for (i = 0; i < ARRAY_SIZE(curves); i++) {
    edwards_scratch_destroy(curves[i], scratches[i]);
    edwards_curve_destroy(curves[i]);
  }
}

static void
test_eddsa_random(drbg_t *rng) {
  size_t i, j;

  printf("Randomized EdDSA testing...\n");

  for (i = 0; i < ARRAY_SIZE(edwards_curves); i++) {
    const char *id = edwards_curves[i];
    edwards_curve_t *ec = edwards_curve_create(i);
    size_t fe_size = edwards_curve_field_size(ec);
    size_t sc_size = edwards_curve_scalar_size(ec);

    printf("  - %s\n", id);

    for (j = 0; j < 10; j++) {
      unsigned char entropy[ENTROPY_SIZE];
      unsigned char priv[EDDSA_MAX_PRIV_SIZE];
      unsigned char msg[EDWARDS_MAX_SCALAR_SIZE];
      unsigned char sig[EDDSA_MAX_SIG_SIZE];
      unsigned char pub[EDDSA_MAX_PUB_SIZE];
      size_t k;

      drbg_generate(rng, entropy, sizeof(entropy));
      drbg_generate(rng, priv, sizeof(priv));
      drbg_generate(rng, msg, sizeof(msg));

      edwards_curve_randomize(ec, entropy);

      eddsa_sign(ec, sig, msg, sc_size, priv, -1, NULL, 0);
      eddsa_pubkey_create(ec, pub, priv);

      ASSERT(eddsa_verify(ec, msg, sc_size, sig, pub, -1, NULL, 0));

      k = drbg_uniform(rng, sc_size);

      msg[k] ^= 1;

      ASSERT(!eddsa_verify(ec, msg, sc_size, sig, pub, -1, NULL, 0));

      msg[k] ^= 1;
      pub[k] ^= 1;

      ASSERT(!eddsa_verify(ec, msg, sc_size, sig, pub, -1, NULL, 0));

      pub[k] ^= 1;
      sig[k] ^= 1;

      ASSERT(!eddsa_verify(ec, msg, sc_size, sig, pub, -1, NULL, 0));

      sig[k] ^= 1;
      sig[fe_size + k] ^= 1;

      ASSERT(!eddsa_verify(ec, msg, sc_size, sig, pub, -1, NULL, 0));

      sig[fe_size + k] ^= 1;

      ASSERT(eddsa_verify(ec, msg, sc_size, sig, pub, -1, NULL, 0));
    }

    edwards_curve_destroy(ec);
  }
}

static void
test_eddsa_elligator2(void) {
  const unsigned char bytes[32] = {
    0xd3, 0xef, 0xfb, 0x44, 0xc4, 0xc9, 0x8d, 0x69,
    0x33, 0xbf, 0xa6, 0x17, 0xfb, 0x88, 0x4f, 0x89,
    0xeb, 0x0a, 0x0c, 0x1c, 0x67, 0x7a, 0xff, 0x86,
    0x7c, 0xea, 0x5e, 0xe6, 0xde, 0xc4, 0x3f, 0x16
  };

  const unsigned char expect[32] = {
    0x91, 0xb8, 0xf0, 0x0a, 0x85, 0x44, 0x0a, 0x07,
    0x2a, 0xbd, 0xe9, 0x80, 0x8d, 0xb2, 0x69, 0x05,
    0x5b, 0xf7, 0x0b, 0x86, 0x83, 0xdb, 0x31, 0x5c,
    0x98, 0x5c, 0xc0, 0x9a, 0x34, 0x80, 0xd0, 0x0d
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char out[32];

  printf("Testing Elligator 2 (ed25519).\n");

  eddsa_pubkey_from_uniform(ec, out, bytes);

  ASSERT(memcmp(out, expect, 32) == 0);
  ASSERT(eddsa_pubkey_to_uniform(ec, out, expect, 0));
  ASSERT(memcmp(out, bytes, 32) == 0);

  edwards_curve_destroy(ec);
}

/*
 * Encoding
 */

static void
test_base16(void) {
  /* https://tools.ietf.org/html/rfc4648#section-10 */
  static const char *vectors[7][2] = {
    {"", ""},
    {"f", "66"},
    {"fo", "666f"},
    {"foo", "666f6f"},
    {"foob", "666f6f62"},
    {"fooba", "666f6f6261"},
    {"foobar", "666f6f626172"}
  };

  static const char *vectors_le[7][2] = {
    {"", ""},
    {"f", "66"},
    {"fo", "6f66"},
    {"foo", "6f6f66"},
    {"foob", "626f6f66"},
    {"fooba", "61626f6f66"},
    {"foobar", "7261626f6f66"}
  };

  static const char *invalid[6] = {
    "6",
    "6x",
    "x6",
    "66 ",
    " 66",
    "666fxa"
  };

  unsigned char buf[32];
  size_t i, len;

  printf("Testing Base16...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *hex = vectors[i][1];

    printf("  - Base16 vector #%lu\n", i + 1);

    ASSERT(base16_test(hex, strlen(hex)));
    ASSERT(base16_decode(buf, &len, hex, strlen(hex)));
    ASSERT(len == strlen(str) && memcmp(buf, str, len) == 0);

    base16_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(hex) && memcmp(buf, hex, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(vectors_le); i++) {
    const char *str = vectors_le[i][0];
    const char *hex = vectors_le[i][1];

    printf("  - Base16-LE vector #%lu\n", i + 1);

    ASSERT(base16le_test(hex, strlen(hex)));
    ASSERT(base16le_decode(buf, &len, hex, strlen(hex)));
    ASSERT(len == strlen(str) && memcmp(buf, str, len) == 0);

    base16le_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(hex) && memcmp(buf, hex, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *hex = invalid[i];

    printf("  - Base16 (invalid) vector #%lu\n", i + 1);

    ASSERT(!base16_test(hex, strlen(hex)));
    ASSERT(!base16le_test(hex, strlen(hex)));
    ASSERT(!base16_decode(buf, &len, hex, strlen(hex)));
    ASSERT(!base16le_decode(buf, &len, hex, strlen(hex)));
  }
}

static void
test_base32(void) {
  /* https://tools.ietf.org/html/rfc4648#section-10 */
  static const char *vectors[7][2] = {
    {"", ""},
    {"f", "my======"},
    {"fo", "mzxq===="},
    {"foo", "mzxw6==="},
    {"foob", "mzxw6yq="},
    {"fooba", "mzxw6ytb"},
    {"foobar", "mzxw6ytboi======"}
  };

  static const char *vectors_hex[7][2] = {
    {"", ""},
    {"f", "co======"},
    {"fo", "cpng===="},
    {"foo", "cpnmu==="},
    {"foob", "cpnmuog="},
    {"fooba", "cpnmuoj1"},
    {"foobar", "cpnmuoj1e8======"}
  };

  unsigned char buf[32];
  size_t i, len;

  printf("Testing Base32...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *b32 = vectors[i][1];

    printf("  - Base32 vector #%lu\n", i + 1);

    ASSERT(base32_test(b32, strlen(b32), 1));
    ASSERT(base32_decode(buf, &len, b32, strlen(b32), 1));
    ASSERT(len == strlen(str) && memcmp(buf, str, len) == 0);

    base32_encode((char *)buf, &len, (unsigned char *)str, strlen(str), 1);
    ASSERT(len == strlen(b32) && memcmp(buf, b32, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(vectors_hex); i++) {
    const char *str = vectors_hex[i][0];
    const char *b32 = vectors_hex[i][1];

    printf("  - Base32-Hex vector #%lu\n", i + 1);

    ASSERT(base32hex_test(b32, strlen(b32), 1));
    ASSERT(base32hex_decode(buf, &len, b32, strlen(b32), 1));
    ASSERT(len == strlen(str) && memcmp(buf, str, len) == 0);

    base32hex_encode((char *)buf, &len, (unsigned char *)str, strlen(str), 1);
    ASSERT(len == strlen(b32) && memcmp(buf, b32, len) == 0);
  }
}

static void
test_base58(void) {
  static const char *vectors[13][2] = {
    {"", ""},
    {"61", "2g"},
    {"626262", "a3gV"},
    {"636363", "aPEr"},
    {
      "73696d706c792061206c6f6e6720737472696e67",
      "2cFupjhnEsSn59qHXstmK2ffpLv2"
    },
    {
      "00eb15231dfceb60925886b67d065299925915aeb172c06647",
      "1NS17iag9jJgTHD1VXjvLCEnZuQ3rJDE9L"
    },
    {"516b6fcd0f", "ABnLTmg"},
    {"bf4f89001e670274dd", "3SEo3LWLoPntC"},
    {"572e4794", "3EFU7m"},
    {"ecac89cad93923c02321", "EJDM8drfXA6uyA"},
    {"10c8511e", "Rt5zm"},
    {"00000000000000000000", "1111111111"},
    {"000000deadbeef", "1116h8cQN"}
  };

  unsigned char buf[128];
  unsigned char data[128];
  size_t i, len, dlen;

  printf("Testing Base58...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *b58 = vectors[i][1];

    ASSERT(base16_decode(data, &dlen, str, strlen(str)));

    printf("  - Base58 vector #%lu\n", i + 1);

    ASSERT(base58_test(b58, strlen(b58)));
    ASSERT(base58_decode(buf, &len, b58, strlen(b58)));
    ASSERT(len == dlen && memcmp(buf, data, len) == 0);

    ASSERT(base58_encode((char *)buf, &len, data, dlen));
    ASSERT(len == strlen(b58) && memcmp(buf, b58, len) == 0);
  }
}

static void
test_base64(void) {
  /* https://tools.ietf.org/html/rfc4648#section-10 */
  static const char *vectors[8][2] = {
    {"", ""},
    {"f", "Zg=="},
    {"fo", "Zm8="},
    {"foo", "Zm9v"},
    {"foob", "Zm9vYg=="},
    {"fooba", "Zm9vYmE="},
    {"foobar", "Zm9vYmFy"},
    {"\x53\xe9\x36\x3b\x29\x62\xfc\xaf", "U+k2Oyli/K8="}
  };

  static const char *invalid[9] = {
    " Zg==",
    "Zg== ",
    "Z g==",
    "Zg ==",
    "Zg",
    "Zm8",
    "Zm9vYg",
    "Zm9vYmE",
    "U-k2Oyli_K8"
  };

  static const char *vectors_url[8][2] = {
    {"", ""},
    {"f", "Zg"},
    {"fo", "Zm8"},
    {"foo", "Zm9v"},
    {"foob", "Zm9vYg"},
    {"fooba", "Zm9vYmE"},
    {"foobar", "Zm9vYmFy"},
    {"\x53\xe9\x36\x3b\x29\x62\xfc\xaf", "U-k2Oyli_K8"}
  };

  static const char *invalid_url[5] = {
    "Zg==",
    "Zm8=",
    "Zm9vYg==",
    "Zm9vYmE=",
    "U+k2Oyli/K8="
  };

  unsigned char buf[32];
  size_t i, len;

  printf("Testing Base64...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *b64 = vectors[i][1];

    printf("  - Base64 vector #%lu\n", i + 1);

    ASSERT(base64_test(b64, strlen(b64)));
    ASSERT(base64_decode(buf, &len, b64, strlen(b64)));
    ASSERT(len == strlen(str) && memcmp(buf, str, len) == 0);

    base64_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(b64) && memcmp(buf, b64, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *b64 = invalid[i];

    printf("  - Base64 (invalid) vector #%lu\n", i + 1);

    ASSERT(!base64_test(b64, strlen(b64)));
    ASSERT(!base64_decode(buf, &len, b64, strlen(b64)));
  }

  for (i = 0; i < ARRAY_SIZE(vectors_url); i++) {
    const char *str = vectors_url[i][0];
    const char *b64 = vectors_url[i][1];

    printf("  - Base64-URL vector #%lu\n", i + 1);

    ASSERT(base64url_test(b64, strlen(b64)));
    ASSERT(base64url_decode(buf, &len, b64, strlen(b64)));
    ASSERT(len == strlen(str) && memcmp(buf, str, len) == 0);

    base64url_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(b64) && memcmp(buf, b64, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid_url); i++) {
    const char *b64 = invalid_url[i];

    printf("  - Base64-URL (invalid) vector #%lu\n", i + 1);

    ASSERT(!base64url_test(b64, strlen(b64)));
    ASSERT(!base64url_decode(buf, &len, b64, strlen(b64)));
  }
}

static void
test_bech32(void) {
  static const char *valid[6][2] = {
    {
      "BC1QW508D6QEJXTDG4Y5R3ZARVARY0C5XW7KV8F3T4",
      "0014751e76e8199196d454941c45d1b3a323f1433bd6"
    },
    {
      "tb1qrp33g0q5c5txsp9arysrx4k6zdkfs4nce4xj0gdcccefvpysxf3q0sl5k7",
      "00201863143c14c5166804bd19203356da136c985678cd4d27a1b8c6329604903262"
    },
    {
      "bc1pw508d6qejxtdg4y5r3zarvary0c5xw7kw508d6"
      "qejxtdg4y5r3zarvary0c5xw7k7grplx",
      "8128751e76e8199196d454941c45d1b3a323f1433b"
      "d6751e76e8199196d454941c45d1b3a323f1433bd6"
    },
    {
      "BC1SW50QA3JX3S",
      "9002751e"
    },
    {
      "bc1zw508d6qejxtdg4y5r3zarvaryvg6kdaj",
      "8210751e76e8199196d454941c45d1b3a323"
    },
    {
      "tb1qqqqqp399et2xygdj5xreqhjjvcmzhxw4aywxecjdzew6hylgvsesrxh6hy",
      "0020000000c4a5cad46221b2a187905e5266362b99d5e91c6ce24d165dab93e86433"
    }
  };

  /* To avoid strcasecmp(3) on windows. */
  static const char *expect[6] = {
    "bc1qw508d6qejxtdg4y5r3zarvary0c5xw7kv8f3t4",
    "tb1qrp33g0q5c5txsp9arysrx4k6zdkfs4nce4xj0gdcccefvpysxf3q0sl5k7",
    "bc1pw508d6qejxtdg4y5r3zarvary0c5xw7kw508d6"
    "qejxtdg4y5r3zarvary0c5xw7k7grplx",
    "bc1sw50qa3jx3s",
    "bc1zw508d6qejxtdg4y5r3zarvaryvg6kdaj",
    "tb1qqqqqp399et2xygdj5xreqhjjvcmzhxw4aywxecjdzew6hylgvsesrxh6hy"
  };

  static const char *invalid[6] = {
    "bc1qw508d6qejxtdg4y5r3zarvary0c5xw7kv8f3t5",
    "bc10w508d6qejxtdg4y5r3zarvary0c5xw7kw508d6"
    "qejxtdg4y5r3zarvary0c5xw7kw5rljs90",
    "bca0w508d6qejxtdg4y5r3zarvary0c5xw7kw508d6"
    "qejxtdg4y5r3zarvary0c5xw7kw5rljs90234567789035",
    "tb1qrp33g0q5c5txsp9arysrx4k6zdkfs4nce4xj0gdcccefvpysxf3q0sL5k7",
    "wtfbbqhelpnoshitwe2z5nuhllhu6z8pptu8m36clz"
    "ge37dnfsdquht73wsx4cmwcwql322x3gmmwq2gjuxp6eaaus",
    "bcfbbqhelpnoshitwe2z7anje5j3wvz8hw3rxadzcp"
    "pgghm0aec23ttfstphjegfx08hwk5uhmusa7j28yrk8cx4qj"
  };

  char hrp[BECH32_MAX_HRP_SIZE + 1];
  unsigned int version;
  uint8_t script[BECH32_MAX_DECODE_SIZE + 2];
  uint8_t data[BECH32_MAX_DECODE_SIZE];
  char out[BECH32_MAX_ENCODE_SIZE + 1];
  size_t i, script_len, data_len;

  printf("Testing Bech32...\n");

  for (i = 0; i < ARRAY_SIZE(valid); i++) {
    const char *addr = valid[i][0];
    const char *hex = valid[i][1];

    printf("  - Bech32 vector #%lu (%s)\n", i + 1, addr);

    ASSERT(base16_decode(script, &script_len, hex, strlen(hex)));

    ASSERT(bech32_test(addr));
    ASSERT(bech32_is(addr));
    ASSERT(bech32_decode(hrp, &version, data, &data_len, addr));
    ASSERT(strlen(hrp) >= 2 && memcmp(hrp, expect[i], 2) == 0);
    ASSERT(2 + data_len == script_len);
    ASSERT(script[0] == (version ? version + 0x80 : 0));
    ASSERT(script[1] == data_len);
    ASSERT(memcmp(data, script + 2, data_len) == 0);

    ASSERT(bech32_encode(out, hrp, version, data, data_len));
    ASSERT(strcmp(out, expect[i]) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *addr = invalid[i];

    printf("  - Bech32 (invalid) vector #%lu (%s)\n", i + 1, addr);

    ASSERT(!bech32_test(addr));
    ASSERT(!bech32_decode(hrp, &version, data, &data_len, addr));
  }
}

static void
test_cash32(void) {
  struct {
    unsigned int type;
    const char *addr;
    const char *prefix;
    const char *hash;
  } valid[6] = {
    {
      0,
      "bitcoincash:qr6m7j9njldwwzlg9v7v53unlr4jkmx6eylep8ekg2",
      "bitcoincash",
      "f5bf48b397dae70be82b3cca4793f8eb2b6cdac9"
    },
    {
      1,
      "bchtest:pr6m7j9njldwwzlg9v7v53unlr4jkmx6eyvwc0uz5t",
      "bchtest",
      "f5bf48b397dae70be82b3cca4793f8eb2b6cdac9"
    },
    {
      1,
      "pref:pr6m7j9njldwwzlg9v7v53unlr4jkmx6ey65nvtks5",
      "pref",
      "f5bf48b397dae70be82b3cca4793f8eb2b6cdac9"
    },
    {
      15,
      "prefix:0r6m7j9njldwwzlg9v7v53unlr4jkmx6ey3qnjwsrf",
      "prefix",
      "f5bf48b397dae70be82b3cca4793f8eb2b6cdac9"
    },
    {
      0,
      "bitcoincash:q9adhakpwzztepkpwp5z0dq62m6u5v5xtyj7j3h2ws4mr9g0",
      "bitcoincash",
      "7adbf6c17084bc86c1706827b41a56f5ca32865925e946ea"
    },
    {
      1,
      "bchtest:p9adhakpwzztepkpwp5z0dq62m6u5v5xtyj7j3h2u94tsynr",
      "bchtest",
      "7adbf6c17084bc86c1706827b41a56f5ca32865925e946ea"
    }
  };

  /* addr, prefix */
  static const char *invalid[6][2] = {
    {
      "prgq00r9x8r3q09aj0g75wmjjr7tkvreusdp3g5tw8",
      "notbitcoincash"
    },
    {
      "bitcoincash:prgq00r9x8r3q09aj0g75wmjjr7tkvreusdp3g5tw8",
      "notbitcoincash"
    },
    {
      "pr3f1x:prgq00r9x8r3q09aj0g75wmjjr7tkvreusdmpmxgtw",
      "pr3f1x"
    },
    {
      "bitcoincash:pr5vxqxg0xrwl2zvxlq9rxffqx00sm44ks62zuqyrr",
      "bitcoincash"
    },
    {
      "bitcoincash:qpg4nt2nwm9mm2a6s6gmcmjx2yr7c9kvta",
      "bitcoincash"
    },
    {
      "bitcoincash:pruptvpkmxamee0f72Sq40gm70wfr624zq0yyxtycm",
      "bitcoincash"
    }
  };

  unsigned int type;
  uint8_t data[CASH32_MAX_DECODE_SIZE];
  uint8_t hash[CASH32_MAX_DECODE_SIZE];
  char out[CASH32_MAX_ENCODE_SIZE + 1];
  size_t i, data_len, hash_len;

  printf("Testing Cash32...\n");

  for (i = 0; i < ARRAY_SIZE(valid); i++) {
    const char *addr = valid[i].addr;
    const char *prefix = valid[i].prefix;
    const char *hex = valid[i].hash;

    printf("  - Cash32 vector #%lu (%s)\n", i + 1, addr);

    ASSERT(base16_decode(hash, &hash_len, hex, strlen(hex)));

    ASSERT(cash32_test(addr, prefix));
    ASSERT(cash32_is(addr, prefix));
    ASSERT(cash32_decode(&type, data, &data_len, addr, prefix));
    ASSERT(type == valid[i].type);
    ASSERT(data_len == hash_len && memcmp(data, hash, hash_len) == 0);

    ASSERT(cash32_encode(out, prefix, type, data, data_len));
    ASSERT(strcmp(out, addr) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *addr = invalid[i][0];
    const char *prefix = invalid[i][1];

    printf("  - Cash32 (invalid) vector #%lu (%s)\n", i + 1, addr);

    ASSERT(!cash32_test(addr, prefix));
    ASSERT(!cash32_decode(&type, data, &data_len, addr, prefix));
  }
}

/*
 * Hash
 */

static void
test_hash(void) {
  unsigned char iv[256];
  unsigned char expect[HASH_MAX_OUTPUT_SIZE];
  unsigned char out[HASH_MAX_OUTPUT_SIZE];
  hash_t hash;
  size_t i, j;

  printf("Testing hashes...\n");

  for (i = 0; i < ARRAY_SIZE(hash_vectors); i++) {
    int type = hash_vectors[i].type;
    size_t size = hash_output_size(type);
    size_t iv_len = sizeof(iv);

    ASSERT(size != 0);

    printf("  - Hash vector #%lu (%s)\n", i + 1, hash_names[type]);

    hex_decode(iv, &iv_len, hash_vectors[i].iv);
    hex_parse(expect, size, hash_vectors[i].expect);

    hash_init(&hash, type);
    hash_update(&hash, iv, iv_len);
    hash_final(&hash, out, size);

    for (j = 0; j < 100; j++) {
      hash_init(&hash, type);
      hash_update(&hash, out, size);
      hash_final(&hash, out, size);
    }

    ASSERT(memcmp(out, expect, size) == 0);
  }
}

static void
test_hmac(void) {
  unsigned char data[256];
  unsigned char key[256];
  unsigned char expect[HASH_MAX_OUTPUT_SIZE];
  unsigned char out[HASH_MAX_OUTPUT_SIZE];
  hmac_t hmac;
  size_t i;

  printf("Testing HMACs...\n");

  for (i = 0; i < ARRAY_SIZE(hmac_vectors); i++) {
    int type = hmac_vectors[i].type;
    size_t size = hash_output_size(type);
    size_t data_len = sizeof(data);
    size_t key_len = sizeof(key);

    ASSERT(size != 0);

    printf("  - HMAC vector #%lu (%s)\n", i + 1, hash_names[type]);

    hex_decode(data, &data_len, hmac_vectors[i].data);
    hex_decode(key, &key_len, hmac_vectors[i].key);
    hex_parse(expect, size, hmac_vectors[i].expect);

    hmac_init(&hmac, type, key, key_len);
    hmac_update(&hmac, data, data_len);
    hmac_final(&hmac, out);

    ASSERT(memcmp(out, expect, size) == 0);
  }
}

/*
 * KDF
 */

static void
test_eb2k(void) {
  unsigned char salt[64];
  unsigned char key[64];
  unsigned char iv[64];
  unsigned char out1[64];
  unsigned char out2[64];
  size_t i;

  printf("Testing EB2K...\n");

  for (i = 0; i < ARRAY_SIZE(eb2k_vectors); i++) {
    const unsigned char *pass = (const unsigned char *)eb2k_vectors[i].pass;
    size_t pass_len = strlen(eb2k_vectors[i].pass);
    size_t salt_len = sizeof(salt);
    size_t key_len = sizeof(key);
    size_t iv_len = sizeof(iv);

    printf("  - EB2K vector #%lu (%s)\n", i + 1, pass);

    hex_decode(salt, &salt_len, eb2k_vectors[i].salt);
    hex_decode(key, &key_len, eb2k_vectors[i].key);
    hex_decode(iv, &iv_len, eb2k_vectors[i].iv);

    ASSERT(sizeof(out1) >= key_len);
    ASSERT(sizeof(out2) >= iv_len);

    ASSERT(eb2k_derive(out1, out2, HASH_MD5, pass, pass_len,
                       salt, salt_len, key_len, iv_len));

    ASSERT(memcmp(out1, key, key_len) == 0);
    ASSERT(memcmp(out2, iv, iv_len) == 0);
  }
}

static void
test_hkdf(void) {
  unsigned char ikm[128];
  unsigned char salt[128];
  unsigned char info[128];
  unsigned char prk[HASH_MAX_OUTPUT_SIZE];
  unsigned char okm[128];
  unsigned char out[128];
  size_t i;

  printf("Testing HKDF...\n");

  for (i = 0; i < ARRAY_SIZE(hkdf_vectors); i++) {
    int type = hkdf_vectors[i].type;
    size_t size = hash_output_size(type);
    size_t ikm_len = sizeof(ikm);
    size_t salt_len = sizeof(salt);
    size_t info_len = sizeof(info);
    size_t len = hkdf_vectors[i].len;

    ASSERT(size != 0);
    ASSERT(sizeof(okm) >= len);
    ASSERT(sizeof(out) >= len);

    printf("  - HKDF vector #%lu (%s)\n", i + 1, hkdf_vectors[i].okm);

    hex_decode(ikm, &ikm_len, hkdf_vectors[i].ikm);
    hex_decode(salt, &salt_len, hkdf_vectors[i].salt);
    hex_decode(info, &info_len, hkdf_vectors[i].info);
    hex_parse(prk, size, hkdf_vectors[i].prk);
    hex_parse(okm, len, hkdf_vectors[i].okm);

    ASSERT(hkdf_extract(out, type, ikm, ikm_len, salt, salt_len));
    ASSERT(memcmp(out, prk, size) == 0);

    ASSERT(hkdf_expand(out, type, prk, info, info_len, len));
    ASSERT(memcmp(out, okm, len) == 0);
  }
}

static void
test_pbkdf2(void) {
  unsigned char pass[256];
  unsigned char salt[256];
  unsigned char expect[256];
  unsigned char out[256];
  size_t i;

  printf("Testing HKDF...\n");

  for (i = 0; i < ARRAY_SIZE(pbkdf2_vectors); i++) {
    int type = pbkdf2_vectors[i].type;
    size_t pass_len = sizeof(pass);
    size_t salt_len = sizeof(salt);
    unsigned int iter = pbkdf2_vectors[i].iter;
    size_t len = pbkdf2_vectors[i].len;

    ASSERT(sizeof(expect) >= len);
    ASSERT(sizeof(out) >= len);

    printf("  - PBKDF2 vector #%lu (%s)\n", i + 1, pbkdf2_vectors[i].expect);

    hex_decode(pass, &pass_len, pbkdf2_vectors[i].pass);
    hex_decode(salt, &salt_len, pbkdf2_vectors[i].salt);
    hex_parse(expect, len, pbkdf2_vectors[i].expect);

    ASSERT(pbkdf2_derive(out, type, pass, pass_len, salt, salt_len, iter, len));
    ASSERT(memcmp(out, expect, len) == 0);
  }
}

static void
test_scrypt(void) {
  static const unsigned char expect1[64] = {
    0x77, 0xd6, 0x57, 0x62, 0x38, 0x65, 0x7b, 0x20, 0x3b, 0x19, 0xca, 0x42,
    0xc1, 0x8a, 0x04, 0x97, 0xf1, 0x6b, 0x48, 0x44, 0xe3, 0x07, 0x4a, 0xe8,
    0xdf, 0xdf, 0xfa, 0x3f, 0xed, 0xe2, 0x14, 0x42, 0xfc, 0xd0, 0x06, 0x9d,
    0xed, 0x09, 0x48, 0xf8, 0x32, 0x6a, 0x75, 0x3a, 0x0f, 0xc8, 0x1f, 0x17,
    0xe8, 0xd3, 0xe0, 0xfb, 0x2e, 0x0d, 0x36, 0x28, 0xcf, 0x35, 0xe2, 0x0c,
    0x38, 0xd1, 0x89, 0x06
  };

  static const unsigned char pass2[] = "password";
  static const unsigned char salt2[] = "NaCl";
  static const unsigned char expect2[64] = {
    0xfd, 0xba, 0xbe, 0x1c, 0x9d, 0x34, 0x72, 0x00, 0x78, 0x56, 0xe7, 0x19,
    0x0d, 0x01, 0xe9, 0xfe, 0x7c, 0x6a, 0xd7, 0xcb, 0xc8, 0x23, 0x78, 0x30,
    0xe7, 0x73, 0x76, 0x63, 0x4b, 0x37, 0x31, 0x62, 0x2e, 0xaf, 0x30, 0xd9,
    0x2e, 0x22, 0xa3, 0x88, 0x6f, 0xf1, 0x09, 0x27, 0x9d, 0x98, 0x30, 0xda,
    0xc7, 0x27, 0xaf, 0xb9, 0x4a, 0x83, 0xee, 0x6d, 0x83, 0x60, 0xcb, 0xdf,
    0xa2, 0xcc, 0x06, 0x40
  };

  static const unsigned char pass3[] = "pleaseletmein";
  static const unsigned char salt3[] = "SodiumChloride";
  static const unsigned char expect3[64] = {
    0x70, 0x23, 0xbd, 0xcb, 0x3a, 0xfd, 0x73, 0x48, 0x46, 0x1c, 0x06, 0xcd,
    0x81, 0xfd, 0x38, 0xeb, 0xfd, 0xa8, 0xfb, 0xba, 0x90, 0x4f, 0x8e, 0x3e,
    0xa9, 0xb5, 0x43, 0xf6, 0x54, 0x5d, 0xa1, 0xf2, 0xd5, 0x43, 0x29, 0x55,
    0x61, 0x3f, 0x0f, 0xcf, 0x62, 0xd4, 0x97, 0x05, 0x24, 0x2a, 0x9a, 0xf9,
    0xe6, 0x1e, 0x85, 0xdc, 0x0d, 0x65, 0x1e, 0x40, 0xdf, 0xcf, 0x01, 0x7b,
    0x45, 0x57, 0x58, 0x87
  };

  unsigned char out[64];

  printf("Testing Scrypt...\n");

  ASSERT(scrypt_derive(out, NULL, 0, NULL, 0, 16, 1, 1, 64));
  ASSERT(memcmp(out, expect1, 64) == 0);

  ASSERT(scrypt_derive(out, pass2, sizeof(pass2) - 1,
                            salt2, sizeof(salt2) - 1, 1024, 8, 16, 64));

  ASSERT(memcmp(out, expect2, 64) == 0);

  ASSERT(scrypt_derive(out, pass3, sizeof(pass3) - 1,
                            salt3, sizeof(salt3) - 1, 16384, 8, 1, 64));

  ASSERT(memcmp(out, expect3, 64) == 0);
}

static void
test_pgpdf(void) {
  static unsigned char expect1[32] = {
    0xc3, 0xab, 0x8f, 0xf1, 0x37, 0x20, 0xe8, 0xad, 0x90, 0x47, 0xdd, 0x39,
    0x46, 0x6b, 0x3c, 0x89, 0x74, 0xe5, 0x92, 0xc2, 0xfa, 0x38, 0x3d, 0x4a,
    0x39, 0x60, 0x71, 0x4c, 0xae, 0xf0, 0xc4, 0xf2
  };

  static unsigned char expect2[32] = {
    0xe6, 0x1e, 0x5b, 0xe4, 0x15, 0x6e, 0x2b, 0x48, 0xae, 0x4e, 0xec, 0xf3,
    0x80, 0xc3, 0x33, 0x5e, 0x64, 0xd5, 0x21, 0xb0, 0x71, 0x4c, 0x55, 0xdd,
    0x15, 0x75, 0xe6, 0x08, 0x24, 0x59, 0xc4, 0xf3
  };

  static unsigned char expect3[32] = {
    0x06, 0x61, 0x52, 0x30, 0xc1, 0x8f, 0x78, 0x5f, 0xf4, 0x95, 0x93, 0x43,
    0xae, 0xba, 0x16, 0xff, 0x26, 0x8d, 0x11, 0x18, 0x89, 0xbc, 0x65, 0xa2,
    0xbc, 0xae, 0xb8, 0xa9, 0xea, 0xad, 0x37, 0xfd
  };

  static const unsigned char pass[] = "foobar";
  static size_t pass_len = sizeof(pass) - 1;
  static const unsigned char salt[] = "salty!";
  static size_t salt_len = sizeof(salt) - 1;
  unsigned char out[32];

  printf("Testing PGPDF...\n");

  ASSERT(pgpdf_derive_simple(out, HASH_SHA256, pass, pass_len, 32));
  ASSERT(memcmp(out, expect1, 32) == 0);

  ASSERT(pgpdf_derive_salted(out, HASH_SHA256, pass, pass_len,
                             salt, salt_len, 32));
  ASSERT(memcmp(out, expect2, 32) == 0);

  ASSERT(pgpdf_derive_iterated(out, HASH_SHA256, pass, pass_len,
                               salt, salt_len, 100, 32));
  ASSERT(memcmp(out, expect3, 32) == 0);
}

static void
test_bcrypt(void) {
  static const struct {
    const char *pass;
    unsigned int rounds;
    const char *salt;
    const char *record;
  } vectors[5] = {
    {
      "<.S.2K(Zq'",
      4,
      "VYAclAMpaXY/oqAo9yUpku",
      "$2a$04$VYAclAMpaXY/oqAo9yUpkuWmoYywaPzyhu56HxXpVltnBIfmO9tgu"
    },
    {
      "5.rApO%5jA",
      5,
      "kVNDrnYKvbNr5AIcxNzeIu",
      "$2a$05$kVNDrnYKvbNr5AIcxNzeIuRcyIF5cZk6UrwHGxENbxP5dVv.WQM/G"
    },
    {
      "oW++kSrQW^",
      6,
      "QLKkRMH9Am6irtPeSKN5sO",
      "$2a$06$QLKkRMH9Am6irtPeSKN5sObJGr3j47cO6Pdf5JZ0AsJXuze0IbsNm"
    },
    {
      "ggJ\\KbTnDG",
      7,
      "4H896R09bzjhapgCPS/LYu",
      "$2a$07$4H896R09bzjhapgCPS/LYuMzAQluVgR5iu/ALF8L8Aln6lzzYXwbq"
    },
    {
      "49b0:;VkH/",
      8,
      "hfvO2retKrSrx5f2RXikWe",
      "$2a$08$hfvO2retKrSrx5f2RXikWeFWdtSesPlbj08t/uXxCeZoHRWDz/xFe"
    }
  };

  size_t i;

  printf("Testing Bcrypt...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *pass = vectors[i].pass;
    unsigned int rounds = vectors[i].rounds;
    const char *salt = vectors[i].salt;
    const char *record = vectors[i].record;
    char out[100];

    printf("  - Bcrypt vector #%lu (%s)\n", i + 1, record);

    ASSERT(bcrypt_generate_with_salt64(out, (const unsigned char *)pass,
                                       strlen(pass), salt, rounds, 'a'));

    ASSERT(strcmp(out, record) == 0);

    ASSERT(bcrypt_verify((const unsigned char *)pass, strlen(pass), record));
  }
}

static void
test_bcrypt_pbkdf(void) {
  static const unsigned char pass[] = "foo";

  static const unsigned char salt[] = {
    0xd8, 0xd5, 0x10, 0x52, 0x71, 0x00, 0x3f, 0x18, 0xaf, 0xb7, 0x51, 0x58,
    0x4a, 0xc9, 0xdf, 0x4d
  };

  static const unsigned char expect[] = {
    0xef, 0x31, 0x20, 0x85, 0x08, 0x6a, 0x78, 0x46, 0xef, 0x9d, 0x43, 0x64,
    0x4c, 0xa3, 0x36, 0x1d, 0x26, 0x76, 0x22, 0xef, 0xff, 0x8e, 0x5c, 0xeb,
    0xed, 0x1e, 0xa9, 0x51, 0x3a, 0x6a, 0xa0, 0xb5, 0x16, 0x0d, 0xe8, 0xbf,
    0x76, 0x11, 0x2e, 0x0c, 0x0e, 0xde, 0xc6, 0x85, 0x94, 0x77, 0x7f, 0x75
  };

  unsigned char out[48];

  printf("Testing Bcrypt-PBKDF...\n");

  ASSERT(bcrypt_pbkdf(out, pass, sizeof(pass) - 1, salt, sizeof(salt), 16, 48));
  ASSERT(memcmp(out, expect, 48) == 0);
}

/*
 * MPI
 */

static void
test_mpi_primes(drbg_t *rng) {
  /* TODO */
  (void)rng;
}

/*
 * Poly1305
 */

static void
test_poly1305(void) {
  unsigned char key[32];
  unsigned char msg[8192];
  unsigned char tag[16];
  unsigned char mac[16];
  poly1305_t ctx;
  size_t i;

  printf("Testing Poly1305...\n");

  for (i = 0; i < ARRAY_SIZE(poly1305_vectors); i++) {
    size_t msg_len = sizeof(msg);

    printf("  - Poly1305 vector #%lu (%s)\n", i + 1, poly1305_vectors[i][0]);

    hex_parse(key, 32, poly1305_vectors[i][0]);
    hex_decode(msg, &msg_len, poly1305_vectors[i][1]);
    hex_parse(tag, 16, poly1305_vectors[i][2]);

    poly1305_init(&ctx, key);
    poly1305_update(&ctx, msg, msg_len);
    poly1305_final(&ctx, mac);

    ASSERT(poly1305_verify(mac, tag));
    ASSERT(memcmp(mac, tag, 16) == 0);
  }
}

/*
 * Random
 */

#ifdef TORSION_HAVE_RNG
static int
looks_random(const void *data, size_t size) {
  const unsigned char *raw = data;
  uint64_t sum = 0;
  uint64_t avg;
  size_t i, j;

  for (i = 0; i < size; i++) {
    for (j = 0; j < 8; j++)
      sum += (raw[i] >> j) & 1;
  }

  avg = (sum * 100) / (size * 8);

  return avg >= 48 && avg <= 52;
}

static void
test_rand_rng(void) {
  uint32_t data[65536 / 4];
  rng_t rng;
  size_t i;

  printf("Testing RNG...\n");

  rng_init(&rng);
  rng_generate(&rng, data, sizeof(data));

  ASSERT(looks_random(data, sizeof(data)));

  for (i = 0; i < ARRAY_SIZE(data); i++)
    data[i] = rng_random(&rng);

  ASSERT(looks_random(data, sizeof(data)));

  for (i = 0; i < ARRAY_SIZE(data); i++)
    data[i] = rng_uniform(&rng, 0x7fffffff);

  ASSERT(looks_random(data, sizeof(data)));
}

static void
test_rand_entropy(void) {
  uint32_t data[65536 / 4];

  printf("Testing getentropy...\n");

  ASSERT(torsion_getentropy(data, sizeof(data)));
  ASSERT(looks_random(data, sizeof(data)));
}
#endif

/*
 * RC4
 */

static void
test_rc4(void) {
  static const unsigned char expect[32] = {
    0x42, 0x83, 0x06, 0x31, 0x1d, 0xd2, 0x34, 0x98,
    0x51, 0x3a, 0x18, 0x7c, 0x36, 0x4c, 0x03, 0xc0,
    0x56, 0x7b, 0x5c, 0x82, 0x94, 0x70, 0x29, 0xfa,
    0x1d, 0x26, 0x24, 0x9f, 0x86, 0x25, 0x1a, 0xa0
  };

  unsigned char key[32];
  unsigned char msg[32];
  unsigned char data[32];
  rc4_t ctx;
  size_t i;

  printf("Testing RC4...\n");

  for (i = 0; i < 32; i++) {
    key[i] = i + 10;
    msg[i] = i + 20;
  }

  memcpy(data, msg, 32);

  rc4_init(&ctx, key, 32);

  for (i = 0; i < 1000; i++)
    rc4_encrypt(&ctx, data, data, 32);

  ASSERT(memcmp(data, expect, 32) == 0);

  rc4_init(&ctx, key, 32);

  for (i = 0; i < 1000; i++)
    rc4_encrypt(&ctx, data, data, 32);

  ASSERT(memcmp(data, msg, 32) == 0);
}

/*
 * RSA
 */

static void
test_rsa(void) {
  static unsigned char priv[RSA_MAX_PRIV_SIZE];
  static unsigned char pub[RSA_MAX_PUB_SIZE];
  static unsigned char msg[HASH_MAX_OUTPUT_SIZE];
  static unsigned char sig1[RSA_MAX_MOD_SIZE];
  static unsigned char sig2[RSA_MAX_MOD_SIZE];
  static unsigned char ct1[RSA_MAX_MOD_SIZE];
  static unsigned char ct2[RSA_MAX_MOD_SIZE];
  static unsigned char out[RSA_MAX_PRIV_SIZE];
  static unsigned char tmp[RSA_MAX_PRIV_SIZE];

  static const unsigned char label[] = {
    0x62, 0x63, 0x72, 0x79, 0x70, 0x74, 0x6f /* "bcrypto" */
  };

  unsigned char entropy[ENTROPY_SIZE];
  size_t i, len;

  for (i = 0; i < ENTROPY_SIZE; i++)
    entropy[i] = i;

  printf("Testing RSA...\n");

  for (i = 0; i < ARRAY_SIZE(rsa_vectors); i++) {
    size_t priv_len = sizeof(priv);
    size_t pub_len = sizeof(pub);
    int hash = rsa_vectors[i].hash;
    int salt_len = rsa_vectors[i].salt_len;
    size_t msg_len = sizeof(msg);
    size_t sig1_len = sizeof(sig1);
    size_t sig2_len = sizeof(sig2);
    size_t ct1_len = sizeof(ct1);
    size_t ct2_len = sizeof(ct2);

    printf("  - RSA vector #%lu\n", i + 1);

    hex_decode(priv, &priv_len, rsa_vectors[i].priv);
    hex_decode(pub, &pub_len, rsa_vectors[i].pub);
    hex_decode(msg, &msg_len, rsa_vectors[i].msg);
    hex_decode(sig1, &sig1_len, rsa_vectors[i].sig1);
    hex_decode(sig2, &sig2_len, rsa_vectors[i].sig2);
    hex_decode(ct1, &ct1_len, rsa_vectors[i].ct1);
    hex_decode(ct2, &ct2_len, rsa_vectors[i].ct2);

    /* Key functions. */
    ASSERT(rsa_privkey_verify(priv, priv_len));
    ASSERT(rsa_pubkey_verify(pub, pub_len));
    ASSERT(rsa_pubkey_create(out, &len, priv, priv_len));
    ASSERT(len == pub_len && memcmp(out, pub, pub_len) == 0);

    /* RSA-PKCS1v1.5 type 1. */
    ASSERT(rsa_sign(out, &len, hash, msg, msg_len,
                    priv, priv_len, entropy));
    ASSERT(len == sig1_len && memcmp(out, sig1, sig1_len) == 0);

    ASSERT(rsa_verify(hash, msg, msg_len, sig1, sig1_len, pub, pub_len));

    msg[0] ^= 1;

    ASSERT(!rsa_verify(hash, msg, msg_len, sig1, sig1_len, pub, pub_len));

    msg[0] ^= 1;
    sig1[0] ^= 1;

    ASSERT(!rsa_verify(hash, msg, msg_len, sig1, sig1_len, pub, pub_len));

    sig1[0] ^= 1;
    pub[0] ^= 1;

    ASSERT(!rsa_verify(hash, msg, msg_len, sig1, sig1_len, pub, pub_len));

    pub[0] ^= 1;

    ASSERT(rsa_verify(hash, msg, msg_len, sig1, sig1_len, pub, pub_len));

    /* RSA-PSS. */
    ASSERT(rsa_sign_pss(out, &len, hash, msg, msg_len,
                        priv, priv_len, salt_len, entropy));

    ASSERT(rsa_verify_pss(hash, msg, msg_len, out, len,
                          pub, pub_len, salt_len));

    ASSERT(rsa_verify_pss(hash, msg, msg_len, sig2, sig2_len,
                          pub, pub_len, salt_len));

    msg[0] ^= 1;

    ASSERT(!rsa_verify_pss(hash, msg, msg_len, sig2, sig2_len,
                           pub, pub_len, salt_len));

    msg[0] ^= 1;
    sig2[0] ^= 1;

    ASSERT(!rsa_verify_pss(hash, msg, msg_len, sig2, sig2_len,
                           pub, pub_len, salt_len));

    sig2[0] ^= 1;
    pub[0] ^= 1;

    ASSERT(!rsa_verify_pss(hash, msg, msg_len, sig2, sig2_len,
                           pub, pub_len, salt_len));

    pub[0] ^= 1;

    ASSERT(rsa_verify_pss(hash, msg, msg_len, sig2, sig2_len,
                          pub, pub_len, salt_len));

    /* RSA-PKCS1v1.5 type 2. */
    ASSERT(rsa_decrypt(out, &len, ct1, ct1_len, priv, priv_len, entropy));
    ASSERT(len == msg_len && memcmp(out, msg, msg_len) == 0);

    ASSERT(rsa_encrypt(tmp, &len, msg, msg_len, pub, pub_len, entropy));
    ASSERT(rsa_decrypt(out, &len, tmp, len, priv, priv_len, entropy));
    ASSERT(len == msg_len && memcmp(out, msg, msg_len) == 0);

    /* RSA-OAEP. */
    ASSERT(rsa_decrypt_oaep(out, &len, hash, ct2, ct2_len,
                            priv, priv_len, label, sizeof(label), entropy));
    ASSERT(len == msg_len && memcmp(out, msg, msg_len) == 0);

    ASSERT(rsa_encrypt_oaep(tmp, &len, hash, msg, msg_len,
                            pub, pub_len, label, sizeof(label), entropy));
    ASSERT(rsa_decrypt_oaep(out, &len, hash, tmp, len,
                            priv, priv_len, label, sizeof(label), entropy));
    ASSERT(len == msg_len && memcmp(out, msg, msg_len) == 0);
  }
}

static void
test_rsa_random(drbg_t *rng) {
  static unsigned char priv[RSA_MAX_PRIV_SIZE];
  static unsigned char pub[RSA_MAX_PUB_SIZE];
  static unsigned char sig[RSA_MAX_MOD_SIZE];
  static unsigned char ct[RSA_MAX_MOD_SIZE];
  static unsigned char pt[RSA_MAX_MOD_SIZE];
  size_t priv_len, pub_len, sig_len, ct_len, pt_len;
  unsigned char msg[32];
  unsigned char entropy[ENTROPY_SIZE];
  size_t i, j;

  printf("Randomized RSA testing...\n");

  for (i = 0; i < 10; i++) {
    drbg_generate(rng, entropy, sizeof(entropy));

    ASSERT(rsa_privkey_generate(priv, &priv_len, 1024, 65537, entropy));
    ASSERT(rsa_privkey_verify(priv, priv_len));

    ASSERT(rsa_pubkey_create(pub, &pub_len, priv, priv_len));
    ASSERT(rsa_pubkey_verify(pub, pub_len));

    drbg_generate(rng, msg, 32);
    drbg_generate(rng, entropy, sizeof(entropy));

    j = drbg_uniform(rng, 128);

    ASSERT(rsa_sign(sig, &sig_len, HASH_SHA256, msg, 32,
                    priv, priv_len, entropy));

    ASSERT(rsa_verify(HASH_SHA256, msg, 32, sig, sig_len, pub, pub_len));

    sig[j] ^= 1;

    ASSERT(!rsa_verify(HASH_SHA256, msg, 32, sig, sig_len, pub, pub_len));

    drbg_generate(rng, msg, 32);
    drbg_generate(rng, entropy, sizeof(entropy));

    ASSERT(rsa_encrypt(ct, &ct_len, msg, 32, pub, pub_len, entropy));
    ASSERT(rsa_decrypt(pt, &pt_len, ct, ct_len, priv, priv_len, entropy));
    ASSERT(pt_len == 32 && memcmp(pt, msg, 32) == 0);

    ct[j] ^= 1;

    ASSERT(!rsa_decrypt(pt, &pt_len, ct, ct_len, priv, priv_len, entropy));

    drbg_generate(rng, msg, 32);
    drbg_generate(rng, entropy, sizeof(entropy));

    ASSERT(rsa_sign_pss(sig, &sig_len, HASH_SHA256, msg, 32,
                        priv, priv_len, 0, entropy));

    ASSERT(rsa_verify_pss(HASH_SHA256, msg, 32, sig, sig_len, pub, pub_len, 0));

    sig[j] ^= 1;

    ASSERT(!rsa_verify_pss(HASH_SHA256, msg, 32,
                           sig, sig_len, pub, pub_len, 0));

    ASSERT(rsa_sign_pss(sig, &sig_len, HASH_SHA256, msg, 32,
                        priv, priv_len, -1, entropy));

    ASSERT(rsa_verify_pss(HASH_SHA256, msg, 32,
                          sig, sig_len, pub, pub_len, -1));

    sig[j] ^= 1;

    ASSERT(!rsa_verify_pss(HASH_SHA256, msg, 32,
                           sig, sig_len, pub, pub_len, -1));

    drbg_generate(rng, msg, 32);
    drbg_generate(rng, entropy, sizeof(entropy));

    ASSERT(rsa_encrypt_oaep(ct, &ct_len, HASH_SHA256, msg, 32,
                            pub, pub_len, NULL, 0, entropy));

    ASSERT(rsa_decrypt_oaep(pt, &pt_len, HASH_SHA256, ct, ct_len,
                            priv, priv_len, NULL, 0, entropy));

    ASSERT(pt_len == 32 && memcmp(pt, msg, 32) == 0);

    ct[j] ^= 1;

    ASSERT(!rsa_decrypt_oaep(pt, &pt_len, HASH_SHA256, ct, ct_len,
                             priv, priv_len, NULL, 0, entropy));
  }
}

/*
 * Salsa20
 */

static void
test_salsa20(void) {
  /* https://github.com/golang/crypto/blob/master/salsa20/salsa20_test.go */
  struct {
    unsigned char key[32];
    unsigned char nonce[8];
    unsigned char expect[64];
  } vectors[4] = {
    {
      {
        0x00, 0x53, 0xa6, 0xf9, 0x4c, 0x9f, 0xf2, 0x45, 0x98, 0xeb, 0x3e, 0x91,
        0xe4, 0x37, 0x8a, 0xdd, 0x30, 0x83, 0xd6, 0x29, 0x7c, 0xcf, 0x22, 0x75,
        0xc8, 0x1b, 0x6e, 0xc1, 0x14, 0x67, 0xba, 0x0d
      },
      { 0x0d, 0x74, 0xdb, 0x42, 0xa9, 0x10, 0x77, 0xde },
      {
        0xc3, 0x49, 0xb6, 0xa5, 0x1a, 0x3e, 0xc9, 0xb7, 0x12, 0xea, 0xed, 0x3f,
        0x90, 0xd8, 0xbc, 0xee, 0x69, 0xb7, 0x62, 0x86, 0x45, 0xf2, 0x51, 0xa9,
        0x96, 0xf5, 0x52, 0x60, 0xc6, 0x2e, 0xf3, 0x1f, 0xd6, 0xc6, 0xb0, 0xae,
        0xa9, 0x4e, 0x13, 0x6c, 0x9d, 0x98, 0x4a, 0xd2, 0xdf, 0x35, 0x78, 0xf7,
        0x8e, 0x45, 0x75, 0x27, 0xb0, 0x3a, 0x04, 0x50, 0x58, 0x0d, 0xd8, 0x74,
        0xf6, 0x3b, 0x1a, 0xb9
      }
    },
    {
      {
        0x05, 0x58, 0xab, 0xfe, 0x51, 0xa4, 0xf7, 0x4a, 0x9d, 0xf0, 0x43, 0x96,
        0xe9, 0x3c, 0x8f, 0xe2, 0x35, 0x88, 0xdb, 0x2e, 0x81, 0xd4, 0x27, 0x7a,
        0xcd, 0x20, 0x73, 0xc6, 0x19, 0x6c, 0xbf, 0x12
      },
      { 0x16, 0x7d, 0xe4, 0x4b, 0xb2, 0x19, 0x80, 0xe7 },
      {
        0xc3, 0xea, 0xaf, 0x32, 0x83, 0x6b, 0xac, 0xe3, 0x2d, 0x04, 0xe1, 0x12,
        0x42, 0x31, 0xef, 0x47, 0xe1, 0x01, 0x36, 0x7d, 0x63, 0x05, 0x41, 0x3a,
        0x0e, 0xeb, 0x07, 0xc6, 0x06, 0x98, 0xa2, 0x87, 0x6e, 0x4d, 0x03, 0x18,
        0x70, 0xa7, 0x39, 0xd6, 0xff, 0xdd, 0xd2, 0x08, 0x59, 0x7a, 0xff, 0x0a,
        0x47, 0xac, 0x17, 0xed, 0xb0, 0x16, 0x7d, 0xd6, 0x7e, 0xba, 0x84, 0xf1,
        0x88, 0x3d, 0x4d, 0xfd
      }
    },
    {
      {
        0x0a, 0x5d, 0xb0, 0x03, 0x56, 0xa9, 0xfc, 0x4f, 0xa2, 0xf5, 0x48, 0x9b,
        0xee, 0x41, 0x94, 0xe7, 0x3a, 0x8d, 0xe0, 0x33, 0x86, 0xd9, 0x2c, 0x7f,
        0xd2, 0x25, 0x78, 0xcb, 0x1e, 0x71, 0xc4, 0x17
      },
      { 0x1f, 0x86, 0xed, 0x54, 0xbb, 0x22, 0x89, 0xf0 },
      {
        0x3c, 0xd2, 0x3c, 0x3d, 0xc9, 0x02, 0x01, 0xac, 0xc0, 0xcf, 0x49, 0xb4,
        0x40, 0xb6, 0xc4, 0x17, 0xf0, 0xdc, 0x8d, 0x84, 0x10, 0xa7, 0x16, 0xd5,
        0x31, 0x4c, 0x05, 0x9e, 0x14, 0xb1, 0xa8, 0xd9, 0xa9, 0xfb, 0x8e, 0xa3,
        0xd9, 0xc8, 0xda, 0xe1, 0x2b, 0x21, 0x40, 0x2f, 0x67, 0x4a, 0xa9, 0x5c,
        0x67, 0xb1, 0xfc, 0x51, 0x4e, 0x99, 0x4c, 0x9d, 0x3f, 0x3a, 0x6e, 0x41,
        0xdf, 0xf5, 0xbb, 0xa6
      }
    },
    {
      {
        0x0f, 0x62, 0xb5, 0x08, 0x5b, 0xae, 0x01, 0x54, 0xa7, 0xfa, 0x4d, 0xa0,
        0xf3, 0x46, 0x99, 0xec, 0x3f, 0x92, 0xe5, 0x38, 0x8b, 0xde, 0x31, 0x84,
        0xd7, 0x2a, 0x7d, 0xd0, 0x23, 0x76, 0xc9, 0x1c
      },
      { 0x28, 0x8f, 0xf6, 0x5d, 0xc4, 0x2b, 0x92, 0xf9 },
      {
        0xe0, 0x0e, 0xbc, 0xcd, 0x70, 0xd6, 0x91, 0x52, 0x72, 0x5f, 0x99, 0x87,
        0x98, 0x21, 0x78, 0xa2, 0xe2, 0xe1, 0x39, 0xc7, 0xbc, 0xbe, 0x04, 0xca,
        0x8a, 0x0e, 0x99, 0xe3, 0x18, 0xd9, 0xab, 0x76, 0xf9, 0x88, 0xc8, 0x54,
        0x9f, 0x75, 0xad, 0xd7, 0x90, 0xba, 0x4f, 0x81, 0xc1, 0x76, 0xda, 0x65,
        0x3c, 0x1a, 0x04, 0x3f, 0x11, 0xa9, 0x58, 0xe1, 0x69, 0xb6, 0xd2, 0x31,
        0x9f, 0x4e, 0xec, 0x1a
      }
    }
  };

  size_t size = 131072;
  unsigned char *out = malloc(size);
  unsigned char xor[64];
  size_t i, j, k;
  salsa20_t ctx;

  ASSERT(out != NULL);

  printf("Testing Salsa20...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const unsigned char *key = vectors[i].key;
    const unsigned char *nonce = vectors[i].nonce;
    const unsigned char *expect = vectors[i].expect;

    printf("  - Salsa20 vector #%lu\n", i + 1);

    memset(out, 0, size);
    memset(xor, 0, 64);

    salsa20_init(&ctx, key, 32, nonce, 8, 0);
    salsa20_encrypt(&ctx, out, out, size);

    for (j = 0; j < size; j += 64) {
      for (k = 0; k < 64; k++)
        xor[k] ^= out[j + k];
    }

    ASSERT(memcmp(xor, expect, 64) == 0);
  }

  free(out);
}

static void
test_xsalsa20(void) {
  /* https://github.com/golang/crypto/blob/master/salsa20/salsa20_test.go */
  static const unsigned char input[] = "Hello world!";
  static const unsigned char nonce[] = "24-byte nonce for xsalsa";
  static const unsigned char key[] = "this is 32-byte key for xsalsa20";
  static const unsigned char expect[12] = {
    0x00, 0x2d, 0x45, 0x13, 0x84, 0x3f,
    0xc2, 0x40, 0xc4, 0x01, 0xe5, 0x41
  };

  unsigned char data[12];
  salsa20_t ctx;

  printf("Testing XSalsa20...\n");

  salsa20_init(&ctx, key, 32, nonce, 24, 0);
  salsa20_encrypt(&ctx, data, input, 12);

  ASSERT(memcmp(data, expect, 12) == 0);
}

/*
 * Secretbox
 */

static void
test_secretbox(void) {
  static const unsigned char expect1[80] = {
    0x84, 0x42, 0xbc, 0x31, 0x3f, 0x46, 0x26, 0xf1, 0x35, 0x9e, 0x3b, 0x50,
    0x12, 0x2b, 0x6c, 0xe6, 0xfe, 0x66, 0xdd, 0xfe, 0x7d, 0x39, 0xd1, 0x4e,
    0x63, 0x7e, 0xb4, 0xfd, 0x5b, 0x45, 0xbe, 0xad, 0xab, 0x55, 0x19, 0x8d,
    0xf6, 0xab, 0x53, 0x68, 0x43, 0x97, 0x92, 0xa2, 0x3c, 0x87, 0xdb, 0x70,
    0xac, 0xb6, 0x15, 0x6d, 0xc5, 0xef, 0x95, 0x7a, 0xc0, 0x4f, 0x62, 0x76,
    0xcf, 0x60, 0x93, 0xb8, 0x4b, 0xe7, 0x7f, 0xf0, 0x84, 0x9c, 0xc3, 0x3e,
    0x34, 0xb7, 0x25, 0x4d, 0x5a, 0x8f, 0x65, 0xad
  };

  static const unsigned char expect2[80] = {
    0x78, 0xea, 0x30, 0xb1, 0x9d, 0x23, 0x41, 0xeb, 0xbd, 0xba, 0x54, 0x18,
    0x0f, 0x82, 0x1e, 0xec, 0x26, 0x5c, 0xf8, 0x63, 0x12, 0x54, 0x9b, 0xea,
    0x8a, 0x37, 0x65, 0x2a, 0x8b, 0xb9, 0x4f, 0x07, 0xb7, 0x8a, 0x73, 0xed,
    0x17, 0x08, 0x08, 0x5e, 0x6d, 0xdd, 0x0e, 0x94, 0x3b, 0xbd, 0xeb, 0x87,
    0x55, 0x07, 0x9a, 0x37, 0xeb, 0x31, 0xd8, 0x61, 0x63, 0xce, 0x24, 0x11,
    0x64, 0xa4, 0x76, 0x29, 0xc0, 0x53, 0x9f, 0x33, 0x0b, 0x49, 0x14, 0xcd,
    0x13, 0x5b, 0x38, 0x55, 0xbc, 0x2a, 0x2d, 0xfc
  };

  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char priv1[32];
  unsigned char priv2[32];
  unsigned char pub1[32];
  unsigned char secret[32];
  unsigned char key[32];
  unsigned char nonce[24];
  unsigned char msg[64];
  unsigned char sealed[SECRETBOX_SEAL_SIZE(64)];
  unsigned char opened[64];

  printf("Testing Secretbox...\n");

  ASSERT(sizeof(sealed) == 80);

  /* Raw. */
  memset(key, 1, sizeof(key));
  memset(nonce, 2, sizeof(nonce));
  memset(msg, 3, sizeof(msg));

  secretbox_seal(sealed, msg, sizeof(msg), key, nonce);

  ASSERT(secretbox_open(opened, sealed, sizeof(sealed), key, nonce));

  ASSERT(memcmp(sealed, expect1, sizeof(expect1)) == 0);
  ASSERT(memcmp(opened, msg, sizeof(msg)) == 0);

  /* crypto_secretbox_xsalsa20poly1305 */
  memset(priv1, 1, sizeof(priv1));
  memset(priv2, 2, sizeof(priv1));
  memset(msg, 3, sizeof(msg));
  memset(nonce, 4, sizeof(nonce));

  ecdh_pubkey_create(ec, pub1, priv1);

  ASSERT(ecdh_derive(ec, secret, pub1, priv2));

  secretbox_derive(key, secret);
  secretbox_seal(sealed, msg, sizeof(msg), key, nonce);

  ASSERT(secretbox_open(opened, sealed, sizeof(sealed), key, nonce));

  ASSERT(memcmp(sealed, expect2, sizeof(expect2)) == 0);
  ASSERT(memcmp(opened, msg, sizeof(msg)) == 0);

  mont_curve_destroy(ec);
}

static void
test_secretbox_random(drbg_t *rng) {
  unsigned char sealed[SECRETBOX_SEAL_SIZE(128)];
  unsigned char opened[128];
  unsigned char msg[128];
  unsigned char key[32];
  unsigned char nonce[24];
  size_t i, len, last;

  printf("Randomized Secretbox testing...\n");

  drbg_generate(rng, key, 32);
  drbg_generate(rng, nonce, 24);

  for (len = 0; len < 128; len += 17) {
    drbg_generate(rng, msg, len);

    secretbox_seal(sealed, msg, len, key, nonce);

    if (len > 0) {
      ASSERT(memcmp(sealed, msg, len) != 0);
      ASSERT(memcmp(sealed + 16, msg, len) != 0);
    }

    ASSERT(secretbox_open(opened, sealed, len + 16, key, nonce));

    last = len + 16;
  }

  for (i = 0; i < last; i++) {
    sealed[i] ^= 0x20;
    ASSERT(!secretbox_open(opened, sealed, last, key, nonce));
    sealed[i] ^= 0x20;
  }

  ASSERT(secretbox_open(opened, sealed, last, key, nonce));
}

/*
 * Siphash
 */

static void
test_siphash(void) {
  static const uint32_t v32 = 0x9dcb553a;
  static const uint64_t v64 = UINT64_C(0x73b4e2ae9316f6b2);
  unsigned char msg[32];
  unsigned char key[16];
  size_t i;

  printf("Testing Siphash...\n");

  for (i = 0; i < 32; i++)
    msg[i] = i;

  for (i = 0; i < 16; i++)
    key[i] = i + 32;

  ASSERT(siphash(msg, 32, key) == UINT64_C(10090947469682793545));
  ASSERT(siphash32(v32, key) == UINT32_C(828368916));
  ASSERT(siphash64(v64, key) == UINT64_C(620914895672640125));
  ASSERT(siphash32k256(v32, msg) == UINT32_C(3909845928));
  ASSERT(siphash64k256(v64, msg) == UINT64_C(16928650368383018294));
  ASSERT(sipmod(msg, 32, key, v64) == UINT64_C(4560894765423557143));
}

/*
 * Util
 */

static void
test_cleanse(drbg_t *rng) {
  static const unsigned char zero[32] = {0};
  unsigned char raw[32];

  printf("Testing Cleanse...\n");

  drbg_generate(rng, raw, 32);

  ASSERT(memcmp(raw, zero, 32) != 0);

  cleanse(raw, 32);

  ASSERT(memcmp(raw, zero, 32) == 0);
}

static void
test_murmur3(void) {
  static const struct {
    const char *str;
    uint32_t seed;
    uint32_t sum;
    uint32_t tweak;
  } vectors[10] = {
    {"", 0xf5afbb1c, 0x3bac103c, 0x75d142d3},
    {"5a", 0xcd8215c1, 0x6f381972, 0x723a6a82},
    {"1fd6", 0x3e1b1e33, 0xcdbb715d, 0x616e796a},
    {"b7b13c", 0xd3e83c4e, 0x894334d6, 0x9097817a},
    {"5fb6a79c", 0x1301f1c1, 0x8261fe1c, 0xb5b901f3},
    {"cbbe6fb48b", 0x7b98ed76, 0x19497f5b, 0xc0302c76},
    {"27529433dcfe", 0xb657be52, 0x7fc0568b, 0x4a48d5a4},
    {"430db045136914", 0x4fc10cb6, 0x045f5f16, 0x261a251b},
    {"fb0663b2ae9e3bc8", 0xee1c44a7, 0xac3e65dc, 0x95c15ad4},
    {"80bf2ed9f2b2cec461", 0x7c4cdae7, 0x7b10f7c8, 0xeb9760a2}
  };

  unsigned char data[32];
  size_t i, len;

  printf("Testing Murmur3...\n");

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i].str;
    uint32_t seed = vectors[i].seed;
    uint32_t sum = vectors[i].sum;
    uint32_t tweak = vectors[i].tweak;

    printf("  - Murmur3 vector #%lu (%u)\n", i + 1, seed);

    ASSERT(base16_decode(data, &len, str, strlen(str)));

    ASSERT(murmur3_sum(data, len, seed) == sum);
    ASSERT(murmur3_tweak(data, len, sum, seed) == tweak);
  }
}

/*
 * Benchmarks
 */

typedef uint64_t bench_t;

uint64_t torsion_hrtime(void);

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

static void
bench_all(drbg_t *rng) {
  bench_ecdsa_pubkey_create(rng);
  bench_ecdsa_derive(rng);
  bench_ecdsa(rng);
  bench_ecdh(rng);
  bench_eddsa(rng);
  bench_hash(rng);
  bench_sha256(rng);
}

/*
 * Main
 */

int
main(int argc, char **argv) {
  int ret = 0;
  drbg_t rng;

  drbg_init_rand(&rng);

  if (argc > 1 && strcmp(argv[1], "bench") == 0) {
    if (argc < 3) {
      bench_all(&rng);
    } else if (strcmp(argv[2], "ecdsa_pubkey_create") == 0) {
      bench_ecdsa_pubkey_create(&rng);
    } else if (strcmp(argv[2], "ecdsa_derive") == 0) {
      bench_ecdsa_derive(&rng);
    } else if (strcmp(argv[2], "ecdsa") == 0) {
      bench_ecdsa(&rng);
    } else if (strcmp(argv[2], "ecdh") == 0) {
      bench_ecdh(&rng);
    } else if (strcmp(argv[2], "eddsa") == 0) {
      bench_eddsa(&rng);
    } else if (strcmp(argv[2], "hash") == 0) {
      bench_hash(&rng);
    } else if (strcmp(argv[2], "sha256") == 0) {
      bench_sha256(&rng);
    } else {
      fprintf(stderr, "Unknown benchmark: %s\n", argv[2]);
      ret = 1;
    }
  } else {
    /* AEAD */
    test_aead();

    /* ChaCha20 */
    test_chacha20();

    /* Cipher */
    test_ciphers();
    test_cipher_modes();
    test_cipher_aead();

    /* DRBG */
    test_hash_drbg();
    test_hmac_drbg();
    test_ctr_drbg();

    /* DSA */
    test_dsa();
    test_dsa_keygen(&rng);

    /* ECC */
#ifdef TORSION_TEST
    __torsion_test_ecc(&rng);
#endif
    test_ecdsa();
    test_ecdsa_random(&rng);
    test_ecdsa_sswu();
    test_ecdsa_svdw();
    test_schnorr_legacy();
    test_schnorr_legacy_random(&rng);
    test_schnorr();
    test_schnorr_random(&rng);
    test_ecdh_x25519();
    test_ecdh_x448();
    test_ecdh_random(&rng);
    test_ecdh_elligator2();
    test_eddsa();
    test_eddsa_random(&rng);
    test_eddsa_elligator2();

    /* Encoding */
    test_base16();
    test_base32();
    test_base58();
    test_base64();
    test_bech32();
    test_cash32();

    /* Hash */
    test_hash();
    test_hmac();

    /* KDF */
    test_eb2k();
    test_hkdf();
    test_pbkdf2();
    test_scrypt();
    test_pgpdf();
    test_bcrypt();
    test_bcrypt_pbkdf();

    /* MPI */
    test_mpi_primes(&rng);

    /* Poly1305 */
    test_poly1305();

    /* Random */
#ifdef TORSION_HAVE_RNG
    test_rand_rng();
    test_rand_entropy();
#endif

    /* RC4 */
    test_rc4();

    /* RSA */
    test_rsa();
    test_rsa_random(&rng);

    /* Salsa20 */
    test_salsa20();
    test_xsalsa20();

    /* Secretbox */
    test_secretbox();
    test_secretbox_random(&rng);

    /* Siphash */
    test_siphash();

    /* Util */
    test_cleanse(&rng);
    test_murmur3();

    printf("All tests passed.\n");
  }

  return ret;
}
