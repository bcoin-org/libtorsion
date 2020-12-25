/*!
 * test.c - tests for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#ifdef TORSION_HAVE_FORK
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif
#ifdef TORSION_HAVE_ZLIB
#include <zlib.h>
#endif

#include <torsion/aead.h>
#include <torsion/drbg.h>
#include <torsion/dsa.h>
#include <torsion/ecc.h>
#include <torsion/encoding.h>
#include <torsion/hash.h>
#include <torsion/ies.h>
#include <torsion/kdf.h>
#include <torsion/mac.h>
#ifdef TORSION_HAVE_RNG
#include <torsion/rand.h>
#endif
#include <torsion/rsa.h>
#include <torsion/stream.h>
#include <torsion/util.h>

#ifdef TORSION_HAVE_THREADS
#include "thread.h"
#endif

#include "utils.h"

#include "data/chacha20_vectors.h"
#include "data/chachapoly_vectors.h"
#include "data/cipher_aead_vectors.h"
#include "data/cipher_mode_vectors.h"
#include "data/cipher_vectors.h"
#include "data/ctr_drbg_vectors.h"
#include "data/dsa_vectors.h"
#include "data/eb2k_vectors.h"
#include "data/ecdsa_vectors.h"
#include "data/eddsa_vectors.h"
#include "data/hash_drbg_vectors.h"
#include "data/hash_vectors.h"
#include "data/hkdf_vectors.h"
#include "data/hmac_drbg_vectors.h"
#include "data/hmac_vectors.h"
#include "data/pbkdf2_vectors.h"
#include "data/poly1305_vectors.h"
#include "data/rsa_vectors.h"
#include "data/schnorr_legacy_vectors.h"
#include "data/schnorr_vectors.h"

/*
 * String Maps
 */

static const char *cipher_names[24] = {
  "AES128",
  "AES192",
  "AES256",
  "ARC2",
  "ARC2_GUTMANN",
  "ARC2_40",
  "ARC2_64",
  "ARC2_128",
  "ARC2_128_GUTMANN",
  "BLOWFISH",
  "CAMELLIA128",
  "CAMELLIA192",
  "CAMELLIA256",
  "CAST5",
  "DES",
  "DES_EDE",
  "DES_EDE3",
  "IDEA",
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
#ifdef TORSION_HAVE_FORK
  pid_t pid = getpid();
#endif
  size_t i;

  for (i = 0; i < ENTROPY_SIZE; i++)
    entropy[i] = rand();

#ifdef TORSION_HAVE_FORK
  ASSERT(sizeof(pid) <= ENTROPY_SIZE);

  for (i = 0; i < sizeof(pid); i++)
    entropy[i] ^= ((unsigned char *)&pid)[i];
#endif

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
 * Memcmp
 */

static void
test_memcmp(drbg_t *unused) {
  static const unsigned char a[4] = {0, 1, 2, 3};
  static const unsigned char b[4] = {0, 1, 2, 3};
  static const unsigned char c[4] = {3, 2, 1, 0};
  static const unsigned char d[4] = {3, 2, 1, 0};

  (void)unused;

  ASSERT(torsion_memcmp(a, b, 4) == 0);
  ASSERT(torsion_memcmp(c, d, 4) == 0);
  ASSERT(torsion_memcmp(a, b, 4) >= 0);
  ASSERT(torsion_memcmp(c, d, 4) >= 0);
  ASSERT(torsion_memcmp(a, b, 4) <= 0);
  ASSERT(torsion_memcmp(c, d, 4) <= 0);
  ASSERT(torsion_memcmp(a, c, 4) != 0);
  ASSERT(torsion_memcmp(c, a, 4) != 0);
  ASSERT(torsion_memcmp(a, c, 4) < 0);
  ASSERT(torsion_memcmp(c, a, 4) > 0);
}

/*
 * AEAD
 */

static void
test_aead_chachapoly(drbg_t *unused) {
  static unsigned char input[8192];
  static unsigned char aad[128];
  static unsigned char key[32];
  static unsigned char nonce[32];
  static unsigned char raw[8192];
  static unsigned char data[8192];
  static unsigned char output[8192];
  static unsigned char tag[16];
  static unsigned char mac[16];
  chachapoly_t ctx;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(chachapoly_vectors); i++) {
    size_t input_len = sizeof(input);
    size_t aad_len = sizeof(aad);
    size_t nonce_len = sizeof(nonce);
    size_t raw_len = sizeof(raw);
    size_t data_len, output_len;

    printf("  - ChaCha20-Poly1305 vector #%u\n", i + 1);

    hex_decode(input, &input_len, chachapoly_vectors[i][0]);
    hex_decode(aad, &aad_len, chachapoly_vectors[i][1]);
    hex_parse(key, 32, chachapoly_vectors[i][2]);
    hex_decode(nonce, &nonce_len, chachapoly_vectors[i][3]);
    hex_decode(raw, &raw_len, chachapoly_vectors[i][4]);

    ASSERT(raw_len >= 16);

    data_len = input_len;
    output_len = raw_len - 16;

    memcpy(data, input, input_len);
    memcpy(output, raw, output_len);
    memcpy(tag, raw + output_len, 16);

    chachapoly_init(&ctx, key, nonce, nonce_len);
    chachapoly_aad(&ctx, aad, aad_len);
    chachapoly_encrypt(&ctx, data, data, data_len);

    ASSERT(torsion_memcmp(data, output, output_len) == 0);

    chachapoly_final(&ctx, mac);

    ASSERT(torsion_memcmp(mac, tag, 16) == 0);
    ASSERT(torsion_memequal(mac, tag, 16));

    chachapoly_init(&ctx, key, nonce, nonce_len);
    chachapoly_aad(&ctx, aad, aad_len);
    chachapoly_auth(&ctx, data, data_len);

    chachapoly_final(&ctx, mac);

    ASSERT(torsion_memcmp(mac, tag, 16) == 0);
    ASSERT(torsion_memequal(mac, tag, 16));

    chachapoly_init(&ctx, key, nonce, nonce_len);
    chachapoly_aad(&ctx, aad, aad_len);
    chachapoly_decrypt(&ctx, data, data, data_len);

    ASSERT(torsion_memcmp(data, input, input_len) == 0);

    chachapoly_final(&ctx, mac);

    ASSERT(torsion_memcmp(mac, tag, 16) == 0);
    ASSERT(torsion_memequal(mac, tag, 16));
  }
}

/*
 * Cipher
 */

static void
test_cipher_contexts(drbg_t *unused) {
  unsigned char key[32];
  unsigned char iv[CIPHER_MAX_BLOCK_SIZE];
  unsigned char expect[CIPHER_MAX_BLOCK_SIZE];
  unsigned char data[CIPHER_MAX_BLOCK_SIZE];
  cipher_t ctx;
  unsigned int i, j;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(cipher_vectors); i++) {
    int type = cipher_vectors[i].type;
    size_t size = cipher_block_size(type);
    size_t key_len = sizeof(key);

    ASSERT(type >= 0 && type <= CIPHER_MAX);
    ASSERT(size != 0);

    printf("  - Cipher vector #%u (%s)\n", i + 1, cipher_names[type]);

    hex_decode(key, &key_len, cipher_vectors[i].key);
    hex_parse(iv, size, cipher_vectors[i].iv);
    hex_parse(expect, size, cipher_vectors[i].expect);

    memcpy(data, iv, size);

    cipher_init(&ctx, type, key, key_len);

    for (j = 0; j < 1000; j++)
      cipher_encrypt(&ctx, data, data);

    ASSERT(torsion_memcmp(data, expect, size) == 0);

    for (j = 0; j < 1000; j++)
      cipher_decrypt(&ctx, data, data);

    ASSERT(torsion_memcmp(data, iv, size) == 0);
  }
}

static void
test_cipher_modes(drbg_t *unused) {
  unsigned char key[64];
  unsigned char iv[CIPHER_MAX_BLOCK_SIZE];
  unsigned char input[64];
  unsigned char output[64];
  unsigned char data[64];
  unsigned int i;

  (void)unused;

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

    printf("  - Cipher mode vector #%u (%s-%s)\n",
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
    ASSERT(torsion_memcmp(data, output, len) == 0);

    ASSERT(cipher_static_decrypt(data, &len, type, mode,
                                 key, key_len, iv, iv_len,
                                 output, output_len));

    ASSERT(len == input_len);
    ASSERT(torsion_memcmp(data, input, len) == 0);
  }
}

static void
test_cipher_aead(drbg_t *unused) {
  unsigned char key[32];
  unsigned char iv[16];
  unsigned char aad[64];
  unsigned char input[64];
  unsigned char output[64];
  unsigned char tag[16];
  unsigned char data[64];
  unsigned char mac[16];
  cipher_stream_t ctx;
  unsigned int i;

  (void)unused;

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

    printf("  - Cipher AEAD vector #%u (%s-%s)\n",
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
    ASSERT(torsion_memcmp(data, output, output_len) == 0);

    ASSERT(cipher_stream_get_tag(&ctx, mac, &len));
    ASSERT(len == tag_len && torsion_memcmp(mac, tag, tag_len) == 0);

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
    ASSERT(torsion_memcmp(data, input, input_len) == 0);
  }
}

/*
 * DRBG
 */

static void
test_drbg_hash(drbg_t *unused) {
  unsigned char entropy[256];
  unsigned char reseed[256];
  unsigned char add1[256];
  unsigned char add2[256];
  unsigned char expect[256];
  unsigned char data[256];
  hash_drbg_t drbg;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(hash_drbg_vectors); i++) {
    int type = hash_drbg_vectors[i].type;
    size_t entropy_len = sizeof(entropy);
    size_t reseed_len = sizeof(reseed);
    size_t add1_len = sizeof(add1);
    size_t add2_len = sizeof(add2);
    size_t expect_len = sizeof(expect);
    size_t data_len;

    printf("  - Hash-DRBG vector #%u\n", i + 1);

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

    ASSERT(torsion_memcmp(data, expect, expect_len) == 0);
  }
}

static void
test_drbg_hmac(drbg_t *unused) {
  unsigned char entropy[256];
  unsigned char reseed[256];
  unsigned char add1[256];
  unsigned char add2[256];
  unsigned char expect[256];
  unsigned char data[256];
  hmac_drbg_t drbg;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(hmac_drbg_vectors); i++) {
    int type = hmac_drbg_vectors[i].type;
    size_t entropy_len = sizeof(entropy);
    size_t reseed_len = sizeof(reseed);
    size_t add1_len = sizeof(add1);
    size_t add2_len = sizeof(add2);
    size_t expect_len = sizeof(expect);
    size_t data_len;

    printf("  - HMAC-DRBG vector #%u\n", i + 1);

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

    ASSERT(torsion_memcmp(data, expect, expect_len) == 0);
  }
}

static void
test_drbg_ctr(drbg_t *unused) {
  unsigned char entropy[256];
  unsigned char pers[256];
  unsigned char reseed[256];
  unsigned char add[256];
  unsigned char add1[256];
  unsigned char add2[256];
  unsigned char expect[256];
  unsigned char data[256];
  ctr_drbg_t drbg;
  unsigned int i;

  (void)unused;

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

    printf("  - CTR-DRBG vector #%u\n", i + 1);

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

    ASSERT(torsion_memcmp(data, expect, expect_len) == 0);
  }
}

/*
 * DSA
 */

static void
test_dsa_vectors(drbg_t *unused) {
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
  unsigned int i;

  (void)unused;

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

    printf("  - DSA vector #%u\n", i + 1);

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
    ASSERT(torsion_memcmp(oparams, params, params_len) == 0);

    ASSERT(dsa_params_create(oparams, &olen, pub, pub_len));
    ASSERT(dsa_params_verify(oparams, olen));
    ASSERT(olen == params_len);
    ASSERT(torsion_memcmp(oparams, params, params_len) == 0);

    ASSERT(dsa_pubkey_create(opub, &olen, priv, priv_len));
    ASSERT(dsa_pubkey_verify(opub, olen));
    ASSERT(olen == pub_len);
    ASSERT(torsion_memcmp(opub, pub, pub_len) == 0);

    qsize = dsa_pubkey_qbits(pub, pub_len) / 8;

    ASSERT(qsize != 0);

    ASSERT(dsa_sig_export(oder, &olen, sig, sig_len, qsize));
    ASSERT(olen == der_len);
    ASSERT(torsion_memcmp(oder, der, der_len) == 0);

    ASSERT(dsa_sig_export(oder, &olen, sig, sig_len, 0));
    ASSERT(olen == der_len);
    ASSERT(torsion_memcmp(oder, der, der_len) == 0);

    ASSERT(dsa_sig_import(osig, &olen, der, der_len, qsize));
    ASSERT(olen == sig_len);
    ASSERT(torsion_memcmp(osig, sig, sig_len) == 0);

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
  ASSERT(torsion_memcmp(alice_sec, bob_sec, bob_sec_len) == 0);
}

/*
 * ECC
 */

static void
test_ecdsa_vectors(drbg_t *unused) {
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
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(curves); i++)
    curves[i] = wei_curve_create(i);

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

    printf("  - ECDSA vector #%u (%s)\n", i + 1, wei_curves[type]);

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
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_privkey_tweak_add(ec, out, priv, tweak));
    ASSERT(torsion_memcmp(out, privadd, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_add(ec, out, out, tweakneg));
    ASSERT(torsion_memcmp(out, priv, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_mul(ec, out, priv, tweak));
    ASSERT(torsion_memcmp(out, privmul, sc_size) == 0);

    ASSERT(ecdsa_privkey_tweak_mul(ec, out, out, tweakinv));
    ASSERT(torsion_memcmp(out, priv, sc_size) == 0);

    ASSERT(ecdsa_privkey_negate(ec, out, priv));
    ASSERT(torsion_memcmp(out, privneg, sc_size) == 0);

    ASSERT(ecdsa_privkey_invert(ec, out, priv));
    ASSERT(torsion_memcmp(out, privinv, sc_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_add(ec, out, &len, pub, pub_size, tweak, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pubadd, pub_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_add(ec, out, &len, pubadd, pub_size, tweakneg, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_mul(ec, out, &len, pub, pub_size, tweak, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pubmul, pub_size) == 0);

    ASSERT(ecdsa_pubkey_tweak_mul(ec, out, &len, pubmul, pub_size, tweakinv, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_pubkey_negate(ec, out, &len, pub, pub_size, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pubneg, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pub;

    publens[0] = pub_size;
    publens[1] = pub_size;

    ASSERT(ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 2, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pubdbl, pub_size) == 0);

    pubs[0] = pubdbl;
    pubs[1] = pubneg;

    publens[0] = pub_size;
    publens[1] = pub_size;

    ASSERT(ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 2, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pubneg;
    pubs[2] = pubconv;

    publens[0] = pub_size;
    publens[1] = pub_size;
    publens[2] = upub_size;

    ASSERT(ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 3, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(!ecdsa_pubkey_combine(ec, out, &len, pubs, publens, 2, 1));

    ASSERT(ecdsa_pubkey_create(ec, out, &len, priv, 0));
    ASSERT(len == upub_size && torsion_memcmp(out, pubconv, upub_size) == 0);

    ASSERT(ecdsa_pubkey_convert(ec, out, &len, pub, pub_size, 0));
    ASSERT(len == upub_size && torsion_memcmp(out, pubconv, upub_size) == 0);

    ASSERT(ecdsa_pubkey_convert(ec, out, &len, pubconv, upub_size, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_is_low_s(ec, sig));

    ASSERT(ecdsa_sig_export(ec, out, &len, sig));
    ASSERT(len == der_len && torsion_memcmp(out, der, der_len) == 0);

    ASSERT(ecdsa_sig_import(ec, out, der, der_len));
    ASSERT(torsion_memcmp(out, sig, sig_size) == 0);

    ASSERT(ecdsa_recover(ec, out, &len, msg, msg_len, sig, param, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(ecdsa_recover(ec, out, &len, msg, msg_len, sig, param, 0));
    ASSERT(len == upub_size && torsion_memcmp(out, pubconv, upub_size) == 0);

    ASSERT(ecdsa_derive(ec, out, &len, pub, pub_size, other, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, secret, pub_size) == 0);

    ASSERT(ecdsa_derive(ec, out, &len, pubconv, upub_size, other, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, secret, pub_size) == 0);

    ASSERT(ecdsa_derive(ec, out, &len, pubhybrid, upub_size, other, 1));
    ASSERT(len == pub_size && torsion_memcmp(out, secret, pub_size) == 0);

    ASSERT(ecdsa_sign(ec, out, &flag, msg, msg_len, priv));
    ASSERT(torsion_memcmp(out, sig, sig_size) == 0);
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
      ASSERT(torsion_memcmp(pub, rec, pub_len) == 0);

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
test_ecdsa_sswu(drbg_t *unused) {
  static const unsigned char bytes[32] = {
    0x6a, 0x53, 0xb3, 0x44, 0xc1, 0x88, 0x8a, 0x62,
    0xf2, 0x81, 0xaf, 0xbd, 0xbe, 0x67, 0x6b, 0xc3,
    0x9e, 0xc6, 0xb1, 0x46, 0x6e, 0x1d, 0x97, 0x53,
    0x45, 0x01, 0xda, 0xf1, 0x63, 0xff, 0xc9, 0x98
  };

  static const unsigned char expect[33] = {
    0x02, 0x61, 0x92, 0xb9, 0xcd, 0x61, 0x6b, 0x4c,
    0x06, 0x83, 0xdd, 0xc9, 0x40, 0x16, 0x1d, 0x60,
    0x8f, 0x68, 0xd3, 0xbd, 0xa9, 0x1f, 0x86, 0x7b,
    0x1f, 0x65, 0x34, 0xdf, 0xdd, 0x04, 0xbd, 0xb2,
    0x8a
  };

  wei_curve_t *ec = wei_curve_create(WEI_CURVE_P256);
  unsigned char out[33];
  size_t out_len;

  (void)unused;

  ecdsa_pubkey_from_uniform(ec, out, &out_len, bytes, 1);

  ASSERT(out_len == 33);
  ASSERT(torsion_memcmp(out, expect, 33) == 0);
  ASSERT(ecdsa_pubkey_to_uniform(ec, out, expect, 33, 3));
  ASSERT(torsion_memcmp(out, bytes, 32) == 0);

  wei_curve_destroy(ec);
}

static void
test_ecdsa_svdw(drbg_t *unused) {
  static const unsigned char bytes[32] = {
    0xb0, 0xf0, 0xa9, 0x2d, 0x14, 0xa9, 0x82, 0xeb,
    0x12, 0x04, 0x78, 0x1a, 0x91, 0x6f, 0xdf, 0x38,
    0x2c, 0x4d, 0x84, 0x69, 0x38, 0xe6, 0x3f, 0x55,
    0xca, 0x59, 0x22, 0xb1, 0x0a, 0xb6, 0x82, 0xa0
  };

  static const unsigned char expect[33] = {
    0x02, 0xa3, 0xb0, 0xbc, 0xa2, 0xaa, 0x06, 0xe3,
    0x78, 0x83, 0x14, 0xb8, 0x73, 0x54, 0xbd, 0x01,
    0x04, 0xf1, 0x10, 0x85, 0xa8, 0x67, 0xab, 0xeb,
    0x4f, 0x43, 0xd2, 0xf6, 0x22, 0xdb, 0xb3, 0x29,
    0x20
  };

  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char out[33];
  size_t out_len;

  (void)unused;

  ecdsa_pubkey_from_uniform(ec, out, &out_len, bytes, 1);

  ASSERT(out_len == 33);
  ASSERT(torsion_memcmp(out, expect, 33) == 0);
  ASSERT(ecdsa_pubkey_to_uniform(ec, out, expect, 33, 1));
  ASSERT(torsion_memcmp(out, bytes, 32) == 0);

  wei_curve_destroy(ec);
}

static void
test_schnorr_legacy_vectors(drbg_t *unused) {
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
  unsigned int i;
  size_t len;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(schnorr_legacy_vectors); i++) {
    size_t priv_len = sizeof(priv);
    size_t pub_len = sizeof(pub);
    size_t msg_len = sizeof(msg);
    size_t sig_len = sizeof(sig);
    int result = schnorr_legacy_vectors[i].result;
    const char *comment = schnorr_legacy_vectors[i].comment;

    printf("  - Schnorr-Legacy vector #%u (%s)\n", i + 1, comment);

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
      ASSERT(len == pub_len && torsion_memcmp(out, pub, pub_len) == 0);
      ASSERT(schnorr_legacy_sign(ec, out, msg, 32, priv));
      ASSERT(torsion_memcmp(out, sig, sig_len) == 0);
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
test_schnorr_vectors(drbg_t *unused) {
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
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(schnorr_vectors); i++) {
    size_t priv_len = sizeof(priv);
    size_t pub_len = sizeof(pub);
    size_t aux_len = sizeof(aux);
    size_t msg_len = sizeof(msg);
    size_t sig_len = sizeof(sig);
    int result = schnorr_vectors[i].result;
    const char *comment = schnorr_vectors[i].comment;

    printf("  - Schnorr vector #%u (%s)\n", i + 1, comment);

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
      ASSERT(torsion_memcmp(out, pub, pub_len) == 0);
      ASSERT(schnorr_sign(ec, out, msg, 32, priv, aux));
      ASSERT(torsion_memcmp(out, sig, sig_len) == 0);
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

  for (i = 0; i < ARRAY_SIZE(wei_curves); i++) {
    const char *id = wei_curves[i];
    wei_curve_t *ec = wei_curve_create(i);

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
test_ecdh_x25519(drbg_t *unused) {
  /* From RFC 7748 */
  static const unsigned char intervals[3][32] = {
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

  (void)unused;

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

  ASSERT(torsion_memcmp(k, intervals[0], 32) == 0);

  for (; i < 1000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 32);
    memcpy(k, t, 32);
  }

  ASSERT(torsion_memcmp(k, intervals[1], 32) == 0);

#ifdef TORSION_TEST_SLOW
  for (; i < 1000000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 32);
    memcpy(k, t, 32);
  }

  ASSERT(torsion_memcmp(k, intervals[2], 32) == 0);
#endif

  ASSERT(ecdh_pubkey_verify(ec, t));

  mont_curve_destroy(ec);
}

static void
test_ecdh_x448(drbg_t *unused) {
  /* From RFC 7748 */
  static const unsigned char intervals[4][56] = {
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

  (void)unused;

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

  ASSERT(torsion_memcmp(k, intervals[0], 56) == 0);

  for (; i < 100; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(torsion_memcmp(k, intervals[1], 56) == 0);

  for (; i < 1000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(torsion_memcmp(k, intervals[2], 56) == 0);

#ifdef TORSION_TEST_SLOW
  for (; i < 1000000; i++) {
    ASSERT(ecdh_derive(ec, t, u, k));
    memcpy(u, k, 56);
    memcpy(k, t, 56);
  }

  ASSERT(torsion_memcmp(k, intervals[3], 56) == 0);
#endif

  ASSERT(ecdh_pubkey_verify(ec, t));

  mont_curve_destroy(ec);
}

static void
test_ecdh_random(drbg_t *rng) {
  size_t i, j;

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

      ASSERT(torsion_memcmp(alice_secret, bob_secret, fe_size) == 0);
    }

    mont_curve_destroy(ec);
  }
}

static void
test_ecdh_elligator2(drbg_t *unused) {
  static const unsigned char bytes[32] = {
    0xf9, 0xd4, 0x08, 0xba, 0x8a, 0x4e, 0x6a, 0x04,
    0xb9, 0xeb, 0x8d, 0x10, 0x38, 0x9b, 0xc0, 0x13,
    0x12, 0x9f, 0x74, 0xed, 0xb7, 0x66, 0xa1, 0x96,
    0xd6, 0x15, 0x16, 0x2b, 0x62, 0xa1, 0xe6, 0x24
  };

  static const unsigned char expect[32] = {
    0xe2, 0x0d, 0xc1, 0xb4, 0xd5, 0xd3, 0x27, 0xbe,
    0x28, 0xdd, 0x80, 0x3a, 0x91, 0xdc, 0x94, 0xa2,
    0xfc, 0xfa, 0x0b, 0x79, 0xe7, 0xc9, 0xe9, 0x09,
    0x52, 0x2a, 0x2f, 0xff, 0x35, 0x16, 0x0e, 0x04
  };

  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char out[32];

  (void)unused;

  ecdh_pubkey_from_uniform(ec, out, bytes);

  ASSERT(torsion_memcmp(out, expect, 32) == 0);
  ASSERT(ecdh_pubkey_to_uniform(ec, out, expect, 0));
  ASSERT(torsion_memcmp(out, bytes, 32) == 0);

  mont_curve_destroy(ec);
}

static void
test_eddsa_vectors(drbg_t *unused) {
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
  unsigned int i;

  (void)unused;

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

    printf("  - EdDSA vector #%u (%s)\n", i + 1, edwards_curves[type]);

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
    ASSERT(torsion_memcmp(out, pub, pub_size) == 0);

    eddsa_pubkey_from_scalar(ec, out, scalar);
    ASSERT(torsion_memcmp(out, pub, pub_size) == 0);

    eddsa_privkey_expand(ec, scalar_, prefix_, priv);
    ASSERT(torsion_memcmp(scalar_, scalar, sc_size) == 0);
    ASSERT(torsion_memcmp(prefix_, prefix, pub_size) == 0);

    eddsa_scalar_reduce(ec, out, scalar);
    ASSERT(torsion_memcmp(out, reduced, sc_size) == 0);

    eddsa_scalar_tweak_add(ec, out, scalar, tweak);
    ASSERT(torsion_memcmp(out, privadd, sc_size) == 0);

    eddsa_scalar_tweak_add(ec, out, privadd, tweakneg);
    ASSERT(torsion_memcmp(out, reduced, sc_size) == 0);

    eddsa_scalar_tweak_mul(ec, out, scalar, tweak);
    ASSERT(torsion_memcmp(out, privmul, sc_size) == 0);

    eddsa_scalar_tweak_mul(ec, out, privmul, tweakinv);
    ASSERT(torsion_memcmp(out, reduced, sc_size) == 0);

    eddsa_scalar_negate(ec, out, scalar);
    ASSERT(torsion_memcmp(out, privneg, sc_size) == 0);

    eddsa_scalar_invert(ec, out, scalar);
    ASSERT(torsion_memcmp(out, privinv, sc_size) == 0);

    ASSERT(eddsa_pubkey_tweak_add(ec, out, pub, tweak));
    ASSERT(torsion_memcmp(out, pubadd, pub_size) == 0);

    ASSERT(eddsa_pubkey_tweak_add(ec, out, pubadd, tweakneg));
    ASSERT(torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(eddsa_pubkey_tweak_mul(ec, out, pub, tweak));
    ASSERT(torsion_memcmp(out, pubmul, pub_size) == 0);

    ASSERT(eddsa_pubkey_tweak_mul(ec, out, pubmul, tweakinv));
    ASSERT(torsion_memcmp(out, pub, pub_size) == 0);

    ASSERT(eddsa_pubkey_negate(ec, out, pub));
    ASSERT(torsion_memcmp(out, pubneg, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pub;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 2));
    ASSERT(torsion_memcmp(out, pubdbl, pub_size) == 0);

    pubs[0] = pubdbl;
    pubs[1] = pubneg;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 2));
    ASSERT(torsion_memcmp(out, pub, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pubneg;
    pubs[2] = pub;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 3));
    ASSERT(torsion_memcmp(out, pub, pub_size) == 0);

    eddsa_privkey_convert(ec, out, priv);
    ASSERT(torsion_memcmp(out, scalar, sc_size) == 0);

    ASSERT(eddsa_pubkey_convert(ec, out, pub));
    ASSERT(torsion_memcmp(out, pubconv, fe_size) == 0);

    ASSERT(eddsa_pubkey_combine(ec, out, NULL, 0));
    ASSERT(torsion_memcmp(out, inf, pub_size) == 0);

    pubs[0] = pub;
    pubs[1] = pubneg;

    ASSERT(eddsa_pubkey_combine(ec, out, pubs, 2));
    ASSERT(eddsa_pubkey_is_infinity(ec, out));
    ASSERT(torsion_memcmp(out, inf, pub_size) == 0);

    ASSERT(eddsa_derive(ec, out, pub, other));
    ASSERT(torsion_memcmp(out, secret, sc_size) == 0);

    eddsa_privkey_convert(ec, out, other);
    ASSERT(eddsa_derive_with_scalar(ec, out, pub, out));
    ASSERT(torsion_memcmp(out, secret, sc_size) == 0);

    eddsa_sign(ec, out, msg, msg_len, priv, ph, NULL, 0);
    ASSERT(torsion_memcmp(out, sig, sig_size) == 0);

    eddsa_sign_with_scalar(ec, out, msg, msg_len, scalar, prefix, ph, NULL, 0);
    ASSERT(torsion_memcmp(out, sig, sig_size) == 0);

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
    ASSERT(torsion_memcmp(out, sigadd, sig_size) == 0);

    ASSERT(eddsa_verify(ec, msg, msg_len, sigadd, pubadd, ph, NULL, 0));
    ASSERT(eddsa_verify_single(ec, msg, msg_len, sigadd, pubadd, ph, NULL, 0));

    eddsa_sign_tweak_mul(ec, out, msg, msg_len, priv, tweak, ph, NULL, 0);
    ASSERT(torsion_memcmp(out, sigmul, sig_size) == 0);

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
test_eddsa_elligator2(drbg_t *unused) {
  static const unsigned char bytes[32] = {
    0xd3, 0xef, 0xfb, 0x44, 0xc4, 0xc9, 0x8d, 0x69,
    0x33, 0xbf, 0xa6, 0x17, 0xfb, 0x88, 0x4f, 0x89,
    0xeb, 0x0a, 0x0c, 0x1c, 0x67, 0x7a, 0xff, 0x86,
    0x7c, 0xea, 0x5e, 0xe6, 0xde, 0xc4, 0x3f, 0x16
  };

  static const unsigned char expect[32] = {
    0x91, 0xb8, 0xf0, 0x0a, 0x85, 0x44, 0x0a, 0x07,
    0x2a, 0xbd, 0xe9, 0x80, 0x8d, 0xb2, 0x69, 0x05,
    0x5b, 0xf7, 0x0b, 0x86, 0x83, 0xdb, 0x31, 0x5c,
    0x98, 0x5c, 0xc0, 0x9a, 0x34, 0x80, 0xd0, 0x0d
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char out[32];

  (void)unused;

  eddsa_pubkey_from_uniform(ec, out, bytes);

  ASSERT(torsion_memcmp(out, expect, 32) == 0);
  ASSERT(eddsa_pubkey_to_uniform(ec, out, expect, 0));
  ASSERT(torsion_memcmp(out, bytes, 32) == 0);

  edwards_curve_destroy(ec);
}

static void
test_ristretto_basepoint_multiples_ed25519(drbg_t *unused) {
  /* https://ristretto.group/test_vectors/ristretto255.html */
  static const unsigned char multiples[][32] = {
    {
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    },
    {
      0xe2, 0xf2, 0xae, 0x0a, 0x6a, 0xbc, 0x4e, 0x71,
      0xa8, 0x84, 0xa9, 0x61, 0xc5, 0x00, 0x51, 0x5f,
      0x58, 0xe3, 0x0b, 0x6a, 0xa5, 0x82, 0xdd, 0x8d,
      0xb6, 0xa6, 0x59, 0x45, 0xe0, 0x8d, 0x2d, 0x76
    },
    {
      0x6a, 0x49, 0x32, 0x10, 0xf7, 0x49, 0x9c, 0xd1,
      0x7f, 0xec, 0xb5, 0x10, 0xae, 0x0c, 0xea, 0x23,
      0xa1, 0x10, 0xe8, 0xd5, 0xb9, 0x01, 0xf8, 0xac,
      0xad, 0xd3, 0x09, 0x5c, 0x73, 0xa3, 0xb9, 0x19
    },
    {
      0x94, 0x74, 0x1f, 0x5d, 0x5d, 0x52, 0x75, 0x5e,
      0xce, 0x4f, 0x23, 0xf0, 0x44, 0xee, 0x27, 0xd5,
      0xd1, 0xea, 0x1e, 0x2b, 0xd1, 0x96, 0xb4, 0x62,
      0x16, 0x6b, 0x16, 0x15, 0x2a, 0x9d, 0x02, 0x59
    },
    {
      0xda, 0x80, 0x86, 0x27, 0x73, 0x35, 0x8b, 0x46,
      0x6f, 0xfa, 0xdf, 0xe0, 0xb3, 0x29, 0x3a, 0xb3,
      0xd9, 0xfd, 0x53, 0xc5, 0xea, 0x6c, 0x95, 0x53,
      0x58, 0xf5, 0x68, 0x32, 0x2d, 0xaf, 0x6a, 0x57
    },
    {
      0xe8, 0x82, 0xb1, 0x31, 0x01, 0x6b, 0x52, 0xc1,
      0xd3, 0x33, 0x70, 0x80, 0x18, 0x7c, 0xf7, 0x68,
      0x42, 0x3e, 0xfc, 0xcb, 0xb5, 0x17, 0xbb, 0x49,
      0x5a, 0xb8, 0x12, 0xc4, 0x16, 0x0f, 0xf4, 0x4e
    },
    {
      0xf6, 0x47, 0x46, 0xd3, 0xc9, 0x2b, 0x13, 0x05,
      0x0e, 0xd8, 0xd8, 0x02, 0x36, 0xa7, 0xf0, 0x00,
      0x7c, 0x3b, 0x3f, 0x96, 0x2f, 0x5b, 0xa7, 0x93,
      0xd1, 0x9a, 0x60, 0x1e, 0xbb, 0x1d, 0xf4, 0x03
    },
    {
      0x44, 0xf5, 0x35, 0x20, 0x92, 0x6e, 0xc8, 0x1f,
      0xbd, 0x5a, 0x38, 0x78, 0x45, 0xbe, 0xb7, 0xdf,
      0x85, 0xa9, 0x6a, 0x24, 0xec, 0xe1, 0x87, 0x38,
      0xbd, 0xcf, 0xa6, 0xa7, 0x82, 0x2a, 0x17, 0x6d
    },
    {
      0x90, 0x32, 0x93, 0xd8, 0xf2, 0x28, 0x7e, 0xbe,
      0x10, 0xe2, 0x37, 0x4d, 0xc1, 0xa5, 0x3e, 0x0b,
      0xc8, 0x87, 0xe5, 0x92, 0x69, 0x9f, 0x02, 0xd0,
      0x77, 0xd5, 0x26, 0x3c, 0xdd, 0x55, 0x60, 0x1c
    },
    {
      0x02, 0x62, 0x2a, 0xce, 0x8f, 0x73, 0x03, 0xa3,
      0x1c, 0xaf, 0xc6, 0x3f, 0x8f, 0xc4, 0x8f, 0xdc,
      0x16, 0xe1, 0xc8, 0xc8, 0xd2, 0x34, 0xb2, 0xf0,
      0xd6, 0x68, 0x52, 0x82, 0xa9, 0x07, 0x60, 0x31
    },
    {
      0x20, 0x70, 0x6f, 0xd7, 0x88, 0xb2, 0x72, 0x0a,
      0x1e, 0xd2, 0xa5, 0xda, 0xd4, 0x95, 0x2b, 0x01,
      0xf4, 0x13, 0xbc, 0xf0, 0xe7, 0x56, 0x4d, 0xe8,
      0xcd, 0xc8, 0x16, 0x68, 0x9e, 0x2d, 0xb9, 0x5f
    },
    {
      0xbc, 0xe8, 0x3f, 0x8b, 0xa5, 0xdd, 0x2f, 0xa5,
      0x72, 0x86, 0x4c, 0x24, 0xba, 0x18, 0x10, 0xf9,
      0x52, 0x2b, 0xc6, 0x00, 0x4a, 0xfe, 0x95, 0x87,
      0x7a, 0xc7, 0x32, 0x41, 0xca, 0xfd, 0xab, 0x42
    },
    {
      0xe4, 0x54, 0x9e, 0xe1, 0x6b, 0x9a, 0xa0, 0x30,
      0x99, 0xca, 0x20, 0x8c, 0x67, 0xad, 0xaf, 0xca,
      0xfa, 0x4c, 0x3f, 0x3e, 0x4e, 0x53, 0x03, 0xde,
      0x60, 0x26, 0xe3, 0xca, 0x8f, 0xf8, 0x44, 0x60
    },
    {
      0xaa, 0x52, 0xe0, 0x00, 0xdf, 0x2e, 0x16, 0xf5,
      0x5f, 0xb1, 0x03, 0x2f, 0xc3, 0x3b, 0xc4, 0x27,
      0x42, 0xda, 0xd6, 0xbd, 0x5a, 0x8f, 0xc0, 0xbe,
      0x01, 0x67, 0x43, 0x6c, 0x59, 0x48, 0x50, 0x1f
    },
    {
      0x46, 0x37, 0x6b, 0x80, 0xf4, 0x09, 0xb2, 0x9d,
      0xc2, 0xb5, 0xf6, 0xf0, 0xc5, 0x25, 0x91, 0x99,
      0x08, 0x96, 0xe5, 0x71, 0x6f, 0x41, 0x47, 0x7c,
      0xd3, 0x00, 0x85, 0xab, 0x7f, 0x10, 0x30, 0x1e
    },
    {
      0xe0, 0xc4, 0x18, 0xf7, 0xc8, 0xd9, 0xc4, 0xcd,
      0xd7, 0x39, 0x5b, 0x93, 0xea, 0x12, 0x4f, 0x3a,
      0xd9, 0x90, 0x21, 0xbb, 0x68, 0x1d, 0xfc, 0x33,
      0x02, 0xa9, 0xd9, 0x9a, 0x2e, 0x53, 0xe6, 0x4e
    }
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  const unsigned char *z = multiples[0];
  const unsigned char *g = multiples[1];
  const unsigned char *points[2];
  unsigned char p[32];
  size_t i;

  (void)unused;

  points[0] = z;

  ASSERT(ristretto_pubkey_combine(ec, p, points, 1));

  for (i = 0; i < ARRAY_SIZE(multiples); i++) {
    const unsigned char *raw = multiples[i];

    ASSERT(torsion_memcmp(p, raw, 32) == 0);

    points[0] = p;
    points[1] = g;

    ASSERT(ristretto_pubkey_combine(ec, p, points, 2));
  }

  edwards_curve_destroy(ec);
}

static void
test_ristretto_bad_points_ed25519(drbg_t *unused) {
  /* https://ristretto.group/test_vectors/ristretto255.html */
  static const unsigned char bad_points[][32] = {
    /* These are all bad because they're non-canonical field encodings. */
    {
      0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
    },
    {
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
    },
    {
      0xf3, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
    },
    {
      0xed, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
    },
    /* These are all bad because they're negative field elements. */
    {
      0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    },
    {
      0x01, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
    },
    {
      0xed, 0x57, 0xff, 0xd8, 0xc9, 0x14, 0xfb, 0x20,
      0x14, 0x71, 0xd1, 0xc3, 0xd2, 0x45, 0xce, 0x3c,
      0x74, 0x6f, 0xcb, 0xe6, 0x3a, 0x36, 0x79, 0xd5,
      0x1b, 0x6a, 0x51, 0x6e, 0xbe, 0xbe, 0x0e, 0x20
    },
    {
      0xc3, 0x4c, 0x4e, 0x18, 0x26, 0xe5, 0xd4, 0x03,
      0xb7, 0x8e, 0x24, 0x6e, 0x88, 0xaa, 0x05, 0x1c,
      0x36, 0xcc, 0xf0, 0xaa, 0xfe, 0xbf, 0xfe, 0x13,
      0x7d, 0x14, 0x8a, 0x2b, 0xf9, 0x10, 0x45, 0x62
    },
    {
      0xc9, 0x40, 0xe5, 0xa4, 0x40, 0x41, 0x57, 0xcf,
      0xb1, 0x62, 0x8b, 0x10, 0x8d, 0xb0, 0x51, 0xa8,
      0xd4, 0x39, 0xe1, 0xa4, 0x21, 0x39, 0x4e, 0xc4,
      0xeb, 0xcc, 0xb9, 0xec, 0x92, 0xa8, 0xac, 0x78
    },
    {
      0x47, 0xcf, 0xc5, 0x49, 0x7c, 0x53, 0xdc, 0x8e,
      0x61, 0xc9, 0x1d, 0x17, 0xfd, 0x62, 0x6f, 0xfb,
      0x1c, 0x49, 0xe2, 0xbc, 0xa9, 0x4e, 0xed, 0x05,
      0x22, 0x81, 0xb5, 0x10, 0xb1, 0x11, 0x7a, 0x24
    },
    {
      0xf1, 0xc6, 0x16, 0x5d, 0x33, 0x36, 0x73, 0x51,
      0xb0, 0xda, 0x8f, 0x6e, 0x45, 0x11, 0x01, 0x0c,
      0x68, 0x17, 0x4a, 0x03, 0xb6, 0x58, 0x12, 0x12,
      0xc7, 0x1c, 0x0e, 0x1d, 0x02, 0x6c, 0x3c, 0x72
    },
    {
      0x87, 0x26, 0x0f, 0x7a, 0x2f, 0x12, 0x49, 0x51,
      0x18, 0x36, 0x0f, 0x02, 0xc2, 0x6a, 0x47, 0x0f,
      0x45, 0x0d, 0xad, 0xf3, 0x4a, 0x41, 0x3d, 0x21,
      0x04, 0x2b, 0x43, 0xb9, 0xd9, 0x3e, 0x13, 0x09
    },
    /* These are all bad because they give a nonsquare x^2. */
    {
      0x26, 0x94, 0x8d, 0x35, 0xca, 0x62, 0xe6, 0x43,
      0xe2, 0x6a, 0x83, 0x17, 0x73, 0x32, 0xe6, 0xb6,
      0xaf, 0xeb, 0x9d, 0x08, 0xe4, 0x26, 0x8b, 0x65,
      0x0f, 0x1f, 0x5b, 0xbd, 0x8d, 0x81, 0xd3, 0x71
    },
    {
      0x4e, 0xac, 0x07, 0x7a, 0x71, 0x3c, 0x57, 0xb4,
      0xf4, 0x39, 0x76, 0x29, 0xa4, 0x14, 0x59, 0x82,
      0xc6, 0x61, 0xf4, 0x80, 0x44, 0xdd, 0x3f, 0x96,
      0x42, 0x7d, 0x40, 0xb1, 0x47, 0xd9, 0x74, 0x2f
    },
    {
      0xde, 0x6a, 0x7b, 0x00, 0xde, 0xad, 0xc7, 0x88,
      0xeb, 0x6b, 0x6c, 0x8d, 0x20, 0xc0, 0xae, 0x96,
      0xc2, 0xf2, 0x01, 0x90, 0x78, 0xfa, 0x60, 0x4f,
      0xee, 0x5b, 0x87, 0xd6, 0xe9, 0x89, 0xad, 0x7b
    },
    {
      0xbc, 0xab, 0x47, 0x7b, 0xe2, 0x08, 0x61, 0xe0,
      0x1e, 0x4a, 0x0e, 0x29, 0x52, 0x84, 0x14, 0x6a,
      0x51, 0x01, 0x50, 0xd9, 0x81, 0x77, 0x63, 0xca,
      0xf1, 0xa6, 0xf4, 0xb4, 0x22, 0xd6, 0x70, 0x42
    },
    {
      0x2a, 0x29, 0x2d, 0xf7, 0xe3, 0x2c, 0xab, 0xab,
      0xbd, 0x9d, 0xe0, 0x88, 0xd1, 0xd1, 0xab, 0xec,
      0x9f, 0xc0, 0x44, 0x0f, 0x63, 0x7e, 0xd2, 0xfb,
      0xa1, 0x45, 0x09, 0x4d, 0xc1, 0x4b, 0xea, 0x08
    },
    {
      0xf4, 0xa9, 0xe5, 0x34, 0xfc, 0x0d, 0x21, 0x6c,
      0x44, 0xb2, 0x18, 0xfa, 0x0c, 0x42, 0xd9, 0x96,
      0x35, 0xa0, 0x12, 0x7e, 0xe2, 0xe5, 0x3c, 0x71,
      0x2f, 0x70, 0x60, 0x96, 0x49, 0xfd, 0xff, 0x22
    },
    {
      0x82, 0x68, 0x43, 0x6f, 0x8c, 0x41, 0x26, 0x19,
      0x6c, 0xf6, 0x4b, 0x3c, 0x7d, 0xdb, 0xda, 0x90,
      0x74, 0x6a, 0x37, 0x86, 0x25, 0xf9, 0x81, 0x3d,
      0xd9, 0xb8, 0x45, 0x70, 0x77, 0x25, 0x67, 0x31
    },
    {
      0x28, 0x10, 0xe5, 0xcb, 0xc2, 0xcc, 0x4d, 0x4e,
      0xec, 0xe5, 0x4f, 0x61, 0xc6, 0xf6, 0x97, 0x58,
      0xe2, 0x89, 0xaa, 0x7a, 0xb4, 0x40, 0xb3, 0xcb,
      0xea, 0xa2, 0x19, 0x95, 0xc2, 0xf4, 0x23, 0x2b
    },
    /* These are all bad because they give a negative xy value. */
    {
      0x3e, 0xb8, 0x58, 0xe7, 0x8f, 0x5a, 0x72, 0x54,
      0xd8, 0xc9, 0x73, 0x11, 0x74, 0xa9, 0x4f, 0x76,
      0x75, 0x5f, 0xd3, 0x94, 0x1c, 0x0a, 0xc9, 0x37,
      0x35, 0xc0, 0x7b, 0xa1, 0x45, 0x79, 0x63, 0x0e
    },
    {
      0xa4, 0x5f, 0xdc, 0x55, 0xc7, 0x64, 0x48, 0xc0,
      0x49, 0xa1, 0xab, 0x33, 0xf1, 0x70, 0x23, 0xed,
      0xfb, 0x2b, 0xe3, 0x58, 0x1e, 0x9c, 0x7a, 0xad,
      0xe8, 0xa6, 0x12, 0x52, 0x15, 0xe0, 0x42, 0x20
    },
    {
      0xd4, 0x83, 0xfe, 0x81, 0x3c, 0x6b, 0xa6, 0x47,
      0xeb, 0xbf, 0xd3, 0xec, 0x41, 0xad, 0xca, 0x1c,
      0x61, 0x30, 0xc2, 0xbe, 0xee, 0xe9, 0xd9, 0xbf,
      0x06, 0x5c, 0x8d, 0x15, 0x1c, 0x5f, 0x39, 0x6e
    },
    {
      0x8a, 0x2e, 0x1d, 0x30, 0x05, 0x01, 0x98, 0xc6,
      0x5a, 0x54, 0x48, 0x31, 0x23, 0x96, 0x0c, 0xcc,
      0x38, 0xae, 0xf6, 0x84, 0x8e, 0x1e, 0xc8, 0xf5,
      0xf7, 0x80, 0xe8, 0x52, 0x37, 0x69, 0xba, 0x32
    },
    {
      0x32, 0x88, 0x84, 0x62, 0xf8, 0xb4, 0x86, 0xc6,
      0x8a, 0xd7, 0xdd, 0x96, 0x10, 0xbe, 0x51, 0x92,
      0xbb, 0xea, 0xf3, 0xb4, 0x43, 0x95, 0x1a, 0xc1,
      0xa8, 0x11, 0x84, 0x19, 0xd9, 0xfa, 0x09, 0x7b
    },
    {
      0x22, 0x71, 0x42, 0x50, 0x1b, 0x9d, 0x43, 0x55,
      0xcc, 0xba, 0x29, 0x04, 0x04, 0xbd, 0xe4, 0x15,
      0x75, 0xb0, 0x37, 0x69, 0x3c, 0xef, 0x1f, 0x43,
      0x8c, 0x47, 0xf8, 0xfb, 0xf3, 0x5d, 0x11, 0x65
    },
    {
      0x5c, 0x37, 0xcc, 0x49, 0x1d, 0xa8, 0x47, 0xcf,
      0xeb, 0x92, 0x81, 0xd4, 0x07, 0xef, 0xc4, 0x1e,
      0x15, 0x14, 0x4c, 0x87, 0x6e, 0x01, 0x70, 0xb4,
      0x99, 0xa9, 0x6a, 0x22, 0xed, 0x31, 0xe0, 0x1e
    },
    {
      0x44, 0x54, 0x25, 0x11, 0x7c, 0xb8, 0xc9, 0x0e,
      0xdc, 0xbc, 0x7c, 0x1c, 0xc0, 0xe7, 0x4f, 0x74,
      0x7f, 0x2c, 0x1e, 0xfa, 0x56, 0x30, 0xa9, 0x67,
      0xc6, 0x4f, 0x28, 0x77, 0x92, 0xa4, 0x8a, 0x4b
    },
    /* This is s = -1, which causes y = 0. */
    {
      0xec, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
      0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f
    }
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  size_t i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(bad_points); i++)
    ASSERT(!ristretto_pubkey_verify(ec, bad_points[i]));

  edwards_curve_destroy(ec);
}

static void
test_ristretto_basepoint_multiples_ed448(drbg_t *unused) {
  /* Self-generated. */
  static const unsigned char multiples[][56] = {
    {
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
      0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    },
    {
      0xf0, 0x76, 0x23, 0x5a, 0x63, 0x3f, 0x56, 0x30,
      0xd2, 0x52, 0x13, 0x34, 0x0b, 0x59, 0x1c, 0x51,
      0xd9, 0x73, 0x5a, 0x73, 0xda, 0x5d, 0x68, 0x6b,
      0xae, 0xb2, 0x11, 0x88, 0x2b, 0x5e, 0x77, 0x8c,
      0xa4, 0xf8, 0xde, 0x1a, 0xa8, 0xcb, 0x20, 0x59,
      0x02, 0xd3, 0xdf, 0xa9, 0xf6, 0xe9, 0x83, 0x7c,
      0x58, 0x8d, 0xe0, 0x5c, 0x0a, 0xd3, 0x04, 0x93
    },
    {
      0x56, 0xd7, 0xcc, 0xe1, 0x60, 0xec, 0x2b, 0x3a,
      0x40, 0x49, 0xed, 0xfc, 0x32, 0x37, 0x1b, 0xa7,
      0x32, 0xb8, 0xa2, 0x36, 0xb5, 0xcb, 0x62, 0x6d,
      0x01, 0x39, 0x7f, 0x66, 0x72, 0xf9, 0xdb, 0x63,
      0x1e, 0x07, 0x37, 0xff, 0x74, 0x3a, 0xae, 0xc6,
      0x72, 0x34, 0x20, 0xeb, 0xea, 0xd8, 0x04, 0x2e,
      0x09, 0xf4, 0x25, 0x3f, 0x7f, 0x0b, 0xea, 0xf6
    },
    {
      0x2a, 0xd6, 0x4e, 0x1f, 0x8a, 0x7d, 0xd1, 0x68,
      0xb8, 0x82, 0xf1, 0xc3, 0x49, 0x9d, 0x52, 0xd7,
      0x32, 0x45, 0xae, 0x8d, 0x01, 0x17, 0xa1, 0x87,
      0x19, 0x0a, 0x02, 0x7c, 0xa3, 0xa9, 0x31, 0x3a,
      0xd9, 0x16, 0x17, 0xf9, 0x2c, 0x21, 0xa2, 0xc7,
      0x9d, 0xe7, 0x1c, 0x43, 0xa2, 0x11, 0xa6, 0x63,
      0xec, 0x73, 0x17, 0xc6, 0xc2, 0x0a, 0x05, 0x9b
    },
    {
      0xd2, 0x2d, 0x2c, 0xde, 0x56, 0x48, 0xcc, 0xee,
      0xe9, 0xdf, 0xd9, 0x2e, 0xff, 0xd5, 0xb1, 0x78,
      0x6a, 0x9f, 0x99, 0x59, 0x04, 0x11, 0x3b, 0x2c,
      0x81, 0x33, 0x02, 0x65, 0x81, 0x6c, 0x1c, 0xe9,
      0x4a, 0x85, 0x5f, 0xc2, 0xf2, 0xbe, 0xa1, 0xec,
      0x08, 0xe9, 0x7a, 0x54, 0xf7, 0x25, 0x8d, 0x74,
      0x47, 0x98, 0x81, 0xf8, 0x53, 0xec, 0x21, 0x17
    },
    {
      0x1e, 0xb3, 0xf3, 0x1f, 0xb7, 0x03, 0x18, 0xd3,
      0x31, 0xab, 0x32, 0xc5, 0xd5, 0xb9, 0x83, 0x6c,
      0x0f, 0xff, 0x52, 0x60, 0x28, 0x23, 0x27, 0x64,
      0xd3, 0xcf, 0x58, 0x64, 0x3f, 0x53, 0x42, 0x62,
      0x41, 0xf8, 0x61, 0xd6, 0x97, 0x23, 0xb9, 0x5e,
      0x5f, 0x3d, 0xa5, 0xaf, 0x3f, 0xce, 0x50, 0x6d,
      0xcc, 0x38, 0xbd, 0x44, 0x84, 0x0f, 0x4d, 0x16
    },
    {
      0xa0, 0x44, 0xcf, 0xe0, 0x2c, 0x5a, 0x69, 0x80,
      0x04, 0x67, 0x97, 0xa9, 0x8b, 0xc4, 0x13, 0xba,
      0x10, 0x99, 0xc4, 0xff, 0x4c, 0x13, 0x29, 0xb6,
      0xaa, 0x8f, 0x2f, 0x42, 0xa1, 0xf3, 0x27, 0xdc,
      0xf1, 0xad, 0x8f, 0xb6, 0xbb, 0x9d, 0x88, 0xc7,
      0x28, 0xc5, 0xc0, 0x2b, 0x7e, 0x9a, 0x19, 0xd9,
      0x16, 0xf1, 0xc0, 0xce, 0xc2, 0x06, 0x35, 0x6b
    },
    {
      0xc6, 0x53, 0x19, 0xc2, 0x47, 0xca, 0x2c, 0x87,
      0x31, 0x9b, 0x69, 0x57, 0x84, 0x9f, 0x6a, 0xf9,
      0x04, 0xca, 0x0c, 0x6b, 0x17, 0x07, 0xcc, 0xf7,
      0x68, 0xef, 0x04, 0xae, 0x5a, 0x89, 0x6b, 0x24,
      0x1c, 0xe9, 0xc7, 0xde, 0x7a, 0x05, 0x5a, 0x2c,
      0x0a, 0x30, 0x07, 0x3f, 0x67, 0x31, 0x63, 0xde,
      0x96, 0x3a, 0x1a, 0xaa, 0x42, 0x7f, 0xe4, 0xc4
    },
    {
      0xd6, 0x14, 0x76, 0xc1, 0xd9, 0x65, 0xc9, 0xfd,
      0x26, 0x88, 0xfa, 0xc3, 0xc6, 0x1b, 0x7c, 0xd1,
      0x32, 0xab, 0x8f, 0x24, 0xf6, 0xfb, 0x7d, 0xb7,
      0x62, 0x14, 0x10, 0x2f, 0x9c, 0x24, 0xfa, 0xec,
      0x2a, 0xda, 0x82, 0x18, 0xb8, 0xc4, 0xb0, 0xd2,
      0xf4, 0x8a, 0x64, 0x73, 0x2a, 0xdc, 0x4f, 0xfb,
      0xd0, 0xea, 0xb3, 0xc3, 0x17, 0x39, 0x7c, 0x40
    },
    {
      0xb8, 0x62, 0xdc, 0x86, 0x35, 0x03, 0x76, 0x00,
      0x4c, 0xf1, 0x46, 0x12, 0xb8, 0x82, 0xb5, 0x27,
      0x35, 0x02, 0x9d, 0xc7, 0x1f, 0x75, 0x62, 0x37,
      0x17, 0xc7, 0xbd, 0x17, 0xae, 0x99, 0x63, 0xa1,
      0xc6, 0xd6, 0xb2, 0x3b, 0xfa, 0x16, 0x32, 0x07,
      0x0b, 0x36, 0x4b, 0x02, 0xbe, 0x84, 0x42, 0xe8,
      0xb6, 0x54, 0x6c, 0x5d, 0x49, 0x97, 0xfd, 0x33
    },
    {
      0xfc, 0x89, 0xc6, 0xca, 0x3b, 0xda, 0x12, 0xbe,
      0x8b, 0xe8, 0x7a, 0x8e, 0xf2, 0x15, 0x0b, 0x6f,
      0xc1, 0xae, 0x95, 0x53, 0x83, 0x25, 0x71, 0xdc,
      0x2d, 0x1f, 0x92, 0x0b, 0x56, 0x37, 0xba, 0x9b,
      0x47, 0x23, 0x6a, 0xaf, 0xf9, 0x41, 0xa0, 0x6a,
      0xf8, 0x8a, 0x8e, 0x40, 0xb9, 0xaa, 0xd8, 0xc0,
      0x4e, 0xc2, 0x13, 0xe5, 0x89, 0x40, 0xd5, 0xee
    },
    {
      0x98, 0x82, 0xd1, 0x00, 0x1e, 0xe0, 0x7c, 0x94,
      0x13, 0xe4, 0x9a, 0xb4, 0xfb, 0xa6, 0x52, 0x69,
      0xe0, 0x9d, 0xc7, 0x78, 0x95, 0xae, 0x52, 0x58,
      0x52, 0xde, 0x48, 0x25, 0x7d, 0xd7, 0xce, 0x61,
      0xec, 0x9d, 0x82, 0x1a, 0x1b, 0x1d, 0xb2, 0xa3,
      0x3f, 0x7c, 0xa1, 0xba, 0x1f, 0x1c, 0x9e, 0xb7,
      0xd1, 0x08, 0xca, 0x47, 0xf4, 0xa2, 0xf0, 0xec
    },
    {
      0x56, 0xd7, 0xc6, 0xfd, 0x8f, 0x00, 0x9a, 0x22,
      0x06, 0x2d, 0x34, 0x50, 0x44, 0x98, 0xca, 0x9d,
      0x21, 0xf8, 0x17, 0x4e, 0xb9, 0xd1, 0x3a, 0xb8,
      0x8a, 0xf6, 0x2a, 0x2d, 0x83, 0xc7, 0xb8, 0x98,
      0x28, 0x49, 0x3e, 0x77, 0xa9, 0xc6, 0xc5, 0x00,
      0x18, 0xb3, 0x38, 0xde, 0x60, 0x53, 0x17, 0x37,
      0xf7, 0xf4, 0x3b, 0xb6, 0xf8, 0xdd, 0x9b, 0x5d
    },
    {
      0x08, 0x26, 0x50, 0xf2, 0xee, 0x0e, 0x41, 0x45,
      0x0c, 0x18, 0x7a, 0x69, 0x1e, 0xb8, 0xa6, 0xe5,
      0x35, 0x33, 0xbc, 0x4b, 0xb5, 0x6b, 0xb6, 0x5d,
      0xfe, 0x2a, 0x8b, 0xee, 0x67, 0x1c, 0xe1, 0xc4,
      0xad, 0xbd, 0x9c, 0x7a, 0x64, 0x1c, 0x03, 0xc4,
      0xa9, 0x13, 0x25, 0x96, 0xd2, 0xee, 0x62, 0xbe,
      0x84, 0x5a, 0xa7, 0xe2, 0xda, 0xf5, 0x60, 0xf7
    },
    {
      0x76, 0xda, 0xc2, 0xa2, 0x49, 0xee, 0xe2, 0xf7,
      0xfa, 0xbe, 0x40, 0x18, 0xfb, 0x66, 0xe5, 0x8c,
      0x82, 0xd3, 0x71, 0x6d, 0xee, 0x1c, 0xb3, 0x56,
      0xbb, 0xf2, 0xd0, 0x29, 0x29, 0xff, 0xbc, 0x74,
      0xa4, 0x38, 0x54, 0x4c, 0xdc, 0xa5, 0x1a, 0x55,
      0x7f, 0x67, 0x06, 0xe2, 0x20, 0x41, 0x7b, 0xaa,
      0x97, 0x89, 0x1c, 0x0a, 0x4d, 0x39, 0x7e, 0xac
    },
    {
      0xe4, 0x78, 0x41, 0xc8, 0xa4, 0xf5, 0xb5, 0x9b,
      0x59, 0x90, 0x3d, 0xc0, 0x2f, 0xed, 0xb5, 0x94,
      0x4d, 0x7d, 0x84, 0x3d, 0x62, 0x5d, 0x7f, 0x30,
      0x15, 0xd0, 0xcc, 0xa0, 0x97, 0x76, 0x90, 0xf8,
      0x21, 0x7b, 0x98, 0xbb, 0x57, 0xda, 0xfd, 0xaa,
      0xc1, 0x5f, 0x27, 0x6c, 0xef, 0x4a, 0xac, 0x9b,
      0x4a, 0x3e, 0x16, 0x3b, 0x07, 0x21, 0x64, 0x68
    }
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED448);
  const unsigned char *z = multiples[0];
  const unsigned char *g = multiples[1];
  const unsigned char *points[2];
  unsigned char p[56];
  size_t i;

  (void)unused;

  points[0] = z;

  ASSERT(ristretto_pubkey_combine(ec, p, points, 1));

  for (i = 0; i < ARRAY_SIZE(multiples); i++) {
    const unsigned char *raw = multiples[i];

    ASSERT(torsion_memcmp(p, raw, 56) == 0);

    points[0] = p;
    points[1] = g;

    ASSERT(ristretto_pubkey_combine(ec, p, points, 2));
  }

  edwards_curve_destroy(ec);
}

static void
test_ristretto_elligator(drbg_t *unused) {
  /* https://ristretto.group/test_vectors/ristretto255.html */
  static const unsigned char bytes[][32] = {
    {
      0xb8, 0xf9, 0x87, 0x31, 0xfd, 0x7b, 0x59, 0x71,
      0x43, 0xa0, 0x06, 0xef, 0x07, 0x69, 0xd3, 0x29,
      0xc0, 0xf9, 0xb9, 0x39, 0x09, 0x66, 0x46, 0xc6,
      0x0f, 0x7f, 0x07, 0x1a, 0xa0, 0x66, 0x86, 0x47
    },
    {
      0xe5, 0x0e, 0xf1, 0xe3, 0x4b, 0x09, 0x76, 0x3c,
      0x80, 0x99, 0xe2, 0x15, 0xb7, 0xd9, 0x5b, 0x88,
      0x62, 0x00, 0xe7, 0x9c, 0x7c, 0x4d, 0x52, 0x8b,
      0x8e, 0x86, 0xa4, 0xa9, 0xa9, 0x3e, 0xfa, 0x34
    },
    {
      0x73, 0x6d, 0x24, 0xdc, 0xb4, 0xdf, 0x63, 0x06,
      0xcc, 0xa9, 0x13, 0x1d, 0xa9, 0x44, 0x54, 0x17,
      0x15, 0x6d, 0xbd, 0x95, 0x7f, 0xcd, 0x5b, 0x66,
      0xac, 0x23, 0x70, 0x23, 0x86, 0x45, 0xba, 0x22
    },
    {
      0x10, 0x31, 0x60, 0x6b, 0xab, 0xc7, 0xa4, 0x09,
      0x81, 0x10, 0x40, 0x3e, 0xf1, 0x3f, 0x84, 0xad,
      0xd1, 0xa0, 0x70, 0xd7, 0x69, 0x32, 0x9d, 0x51,
      0xfd, 0x69, 0x01, 0x9a, 0xe5, 0x19, 0x78, 0x53
    },
    {
      0x9c, 0x83, 0xa1, 0xa2, 0xec, 0xfb, 0x05, 0xbb,
      0xa7, 0xab, 0x11, 0xb2, 0x94, 0xd2, 0x5a, 0xcf,
      0x56, 0x15, 0x4f, 0xa1, 0xa7, 0xd7, 0xea, 0x01,
      0x88, 0xf2, 0xb6, 0xf8, 0x26, 0x55, 0x4f, 0x56
    },
    {
      0xfb, 0xb1, 0x7c, 0x36, 0x12, 0x65, 0x4b, 0xeb,
      0xf5, 0xba, 0x13, 0x2e, 0x85, 0x9d, 0xe5, 0x40,
      0x0a, 0x88, 0xb5, 0xb9, 0x4e, 0x90, 0xfe, 0xa7,
      0x89, 0x31, 0x6b, 0x0a, 0x3d, 0x0a, 0x15, 0x19
    },
    {
      0xe8, 0xc1, 0x14, 0x44, 0xf0, 0x4d, 0xba, 0x4d,
      0xb7, 0x28, 0x2c, 0x56, 0x96, 0x1f, 0xc6, 0xd4,
      0x4c, 0x51, 0x03, 0xd9, 0xc5, 0x08, 0x7e, 0x80,
      0x7e, 0x98, 0xa4, 0xd0, 0x99, 0x2c, 0xbd, 0x4d
    },
    {
      0xad, 0xe5, 0x95, 0xb1, 0x25, 0xe6, 0x1e, 0x45,
      0x3d, 0x38, 0xac, 0xbe, 0xdb, 0x73, 0xa7, 0xc2,
      0x47, 0x86, 0x3b, 0x4b, 0x1c, 0xf4, 0x76, 0x1a,
      0xa2, 0x61, 0x40, 0x10, 0x0f, 0xbd, 0x1e, 0x40
    },
    {
      0x6a, 0x47, 0x3d, 0x6b, 0xfa, 0x75, 0x2a, 0x97,
      0x5b, 0xca, 0xd4, 0x64, 0x34, 0xbc, 0xbe, 0x15,
      0x7d, 0xda, 0x1f, 0x12, 0xfd, 0xf1, 0xa0, 0x85,
      0x39, 0xf2, 0x03, 0xa4, 0xbd, 0x44, 0x6f, 0x4b
    },
    {
      0x70, 0xcc, 0xb6, 0x5a, 0xdc, 0xc6, 0x78, 0x49,
      0xad, 0x6b, 0xc1, 0x11, 0xe3, 0x28, 0xa2, 0x24,
      0x96, 0x8d, 0xeb, 0x37, 0xac, 0xb7, 0x0c, 0x27,
      0xc2, 0x88, 0x2b, 0x99, 0xf4, 0x76, 0x5b, 0x59
    },
    {
      0x6f, 0x18, 0xcb, 0x7b, 0xfe, 0xbd, 0x0b, 0xa2,
      0x33, 0xc4, 0xa3, 0x88, 0xcc, 0x8f, 0x0a, 0xde,
      0x21, 0x70, 0x51, 0xcd, 0x22, 0x23, 0x08, 0x42,
      0x5a, 0x06, 0xa4, 0x3a, 0xaa, 0xb1, 0x22, 0x19
    },
    {
      0xe1, 0xb7, 0x1e, 0x34, 0xec, 0x52, 0x06, 0xb7,
      0x6d, 0x19, 0xe3, 0xb5, 0x19, 0x52, 0x29, 0xc1,
      0x50, 0x4d, 0xa1, 0x50, 0xf2, 0xcb, 0x4f, 0xcc,
      0x88, 0xf5, 0x83, 0x6e, 0xed, 0x6a, 0x03, 0x3a
    },
    {
      0xcf, 0xf6, 0x26, 0x38, 0x1e, 0x56, 0xb0, 0x5a,
      0x1b, 0xc8, 0x3d, 0x2a, 0xdd, 0x1b, 0x38, 0xd2,
      0x4f, 0xb2, 0xbd, 0x78, 0x44, 0xc1, 0x78, 0xa7,
      0x4d, 0xb9, 0x35, 0xc5, 0x7c, 0x80, 0xbf, 0x7e
    },
    {
      0x01, 0x88, 0xd7, 0x50, 0xf0, 0x2e, 0x3f, 0x93,
      0x10, 0xf4, 0xe6, 0xcf, 0x52, 0xbd, 0x4a, 0x32,
      0x6a, 0xa9, 0x8a, 0x56, 0x1e, 0x83, 0xd6, 0xca,
      0xa6, 0x7d, 0xfb, 0xe4, 0x62, 0x18, 0x24, 0x15
    },
    {
      0xd2, 0xcf, 0xe4, 0x38, 0x9b, 0x74, 0xcf, 0x36,
      0x54, 0xc3, 0xfb, 0xd7, 0xf9, 0xc7, 0x74, 0x4b,
      0x6d, 0xef, 0xc4, 0xfb, 0xc2, 0xf6, 0xfc, 0xe4,
      0x46, 0x92, 0x9c, 0x23, 0x19, 0x27, 0xf1, 0x04
    },
    {
      0x22, 0x74, 0x7b, 0x09, 0x08, 0x28, 0x5d, 0xbd,
      0x09, 0x67, 0x39, 0x67, 0x42, 0xe3, 0x03, 0x02,
      0x9d, 0x6b, 0x86, 0xdb, 0xca, 0x4a, 0xe6, 0x9a,
      0x4e, 0x6b, 0xdb, 0xc3, 0xd6, 0x0e, 0x54, 0x50
    }
  };

  static const unsigned char images[][32] = {
    {
      0xb0, 0x9d, 0xed, 0x61, 0x42, 0x1d, 0x8c, 0xa6,
      0xa8, 0x5e, 0x1a, 0x9d, 0xd4, 0xd8, 0xe5, 0xa0,
      0xc3, 0xf6, 0xe8, 0xef, 0xa9, 0x70, 0x3f, 0xc1,
      0x40, 0x20, 0x98, 0x45, 0x0b, 0xbe, 0xf6, 0x56
    },
    {
      0xea, 0x8d, 0x4d, 0xcb, 0xb5, 0xe1, 0xfa, 0x4a,
      0xab, 0x3e, 0x0f, 0x76, 0x4e, 0xd4, 0x96, 0x13,
      0x83, 0x0e, 0xbc, 0xee, 0xc2, 0xf4, 0x8d, 0x8a,
      0xa6, 0xa2, 0x53, 0x7a, 0xe4, 0xc9, 0x13, 0x1a
    },
    {
      0xe8, 0xe7, 0x33, 0x5c, 0x05, 0xa8, 0x50, 0x24,
      0xad, 0xb3, 0x68, 0x44, 0xba, 0x95, 0x44, 0x28,
      0x8c, 0xaa, 0x1b, 0x67, 0x63, 0x8c, 0x15, 0xf2,
      0x2b, 0x3e, 0xfa, 0x86, 0xd0, 0xff, 0x3d, 0x59
    },
    {
      0xd0, 0x78, 0x8c, 0x81, 0xb1, 0xb3, 0xed, 0x9f,
      0xfc, 0xa0, 0x1c, 0x0d, 0xce, 0x05, 0xd3, 0xf1,
      0xc0, 0xda, 0x01, 0x61, 0x82, 0xf1, 0x14, 0xa9,
      0x77, 0x2e, 0xf6, 0x1d, 0x4f, 0x50, 0x4d, 0x54
    },
    {
      0xca, 0x0b, 0xec, 0x91, 0x3a, 0x0c, 0xb5, 0x9d,
      0xd1, 0x06, 0xd5, 0x58, 0x4b, 0x93, 0x0b, 0x77,
      0xbf, 0x8b, 0x2f, 0x8e, 0x21, 0x24, 0x99, 0xc1,
      0xdf, 0xb7, 0xb2, 0x08, 0xcd, 0x78, 0xf8, 0x6e
    },
    {
      0x1a, 0x42, 0xe7, 0x43, 0xcb, 0xaf, 0x74, 0x82,
      0x20, 0x88, 0x3e, 0xfd, 0xd7, 0x2e, 0x05, 0xd6,
      0xa6, 0xf8, 0x6c, 0xed, 0xd8, 0x47, 0xf4, 0xad,
      0x48, 0x85, 0x52, 0x06, 0x8f, 0xf0, 0x68, 0x29
    },
    {
      0x28, 0x9d, 0x66, 0x60, 0xc9, 0xdf, 0xc8, 0xc5,
      0x96, 0xb5, 0x6a, 0x53, 0x67, 0x7e, 0x8f, 0x21,
      0x91, 0xe6, 0x4e, 0x06, 0xab, 0x92, 0xd2, 0x8f,
      0x70, 0x05, 0xf5, 0x17, 0xb7, 0x8a, 0x12, 0x78
    },
    {
      0xdc, 0x25, 0x1b, 0xcb, 0xef, 0xc4, 0xb0, 0x83,
      0x25, 0x42, 0xbc, 0xf3, 0xb9, 0xfa, 0x71, 0x17,
      0xa7, 0xd3, 0x9a, 0xf3, 0xa8, 0xd7, 0x36, 0xab,
      0x9f, 0x24, 0xc3, 0x51, 0x0d, 0x96, 0x2b, 0x2b
    },
    {
      0xe8, 0x79, 0xb0, 0xde, 0xb7, 0xc4, 0x9f, 0x5a,
      0xee, 0xc1, 0x69, 0x34, 0x65, 0xa7, 0xf4, 0xaa,
      0x79, 0x72, 0xc4, 0x06, 0x43, 0x98, 0x50, 0xb9,
      0xdd, 0x07, 0x53, 0x69, 0xb0, 0xd0, 0xe0, 0x79
    },
    {
      0xe2, 0xb5, 0xb7, 0x34, 0xf1, 0xa3, 0x3d, 0xb3,
      0xdd, 0xcf, 0xdc, 0x49, 0xf5, 0xf2, 0x19, 0xec,
      0x43, 0x54, 0xb3, 0xde, 0xa7, 0x3e, 0xa7, 0xb6,
      0x20, 0x09, 0x5c, 0x1e, 0xa5, 0x7f, 0xcc, 0x44
    },
    {
      0xe2, 0x77, 0x10, 0xf2, 0xc8, 0x8b, 0xf0, 0x57,
      0x0b, 0xde, 0x5c, 0x92, 0x9c, 0xf3, 0x2e, 0x77,
      0x41, 0x3b, 0x01, 0xf8, 0x5c, 0xb7, 0x32, 0xaf,
      0x57, 0x28, 0xce, 0x35, 0xd0, 0xdc, 0x94, 0x0d
    },
    {
      0x46, 0xf0, 0x4f, 0x70, 0x36, 0x9d, 0xe4, 0x92,
      0x4a, 0x7a, 0xd8, 0x58, 0xe8, 0x3e, 0x9e, 0x0d,
      0x0e, 0x92, 0x73, 0x75, 0xb0, 0xde, 0x5a, 0xe1,
      0xf4, 0x17, 0x5e, 0xbe, 0x96, 0x07, 0x88, 0x60
    },
    {
      0x16, 0x47, 0xf1, 0x67, 0x2d, 0xc1, 0xc3, 0x90,
      0xb7, 0x65, 0x9a, 0x32, 0x27, 0x44, 0x31, 0x6e,
      0x33, 0x2c, 0x3e, 0x00, 0xe5, 0x71, 0x48, 0x51,
      0xa8, 0x1d, 0x49, 0x6a, 0x66, 0x28, 0x84, 0x18
    },
    {
      0xc4, 0x85, 0x6b, 0x0b, 0x82, 0x69, 0x4a, 0x21,
      0xcc, 0xab, 0x85, 0xdd, 0xae, 0xc1, 0xf1, 0x24,
      0x26, 0xb3, 0xc4, 0x6b, 0xdb, 0xb9, 0xb5, 0xfd,
      0xe4, 0x2f, 0x9b, 0x2a, 0xe7, 0x49, 0x29, 0x4e
    },
    {
      0x3a, 0xff, 0xe1, 0xc5, 0x73, 0xd0, 0xa0, 0x8f,
      0x27, 0xc5, 0x52, 0x45, 0x8f, 0xeb, 0x5c, 0xaa,
      0x4a, 0x28, 0x39, 0x0b, 0xab, 0xe3, 0x1a, 0xb9,
      0xd9, 0xcf, 0x5a, 0xb9, 0xc5, 0xbe, 0x23, 0x3c
    },
    {
      0x58, 0x2b, 0x5c, 0x76, 0xdf, 0x88, 0x69, 0x91,
      0xee, 0xba, 0x73, 0x08, 0xd6, 0x70, 0x99, 0xfd,
      0x26, 0x6c, 0xcd, 0xe6, 0x9d, 0x82, 0x0b, 0x42,
      0x65, 0x55, 0xfd, 0x6e, 0x6e, 0x0e, 0x94, 0x70
    }
  };

  static const unsigned int hints[] = {
    0,
    16,
    23,
    0,
    2,
    23,
    2,
    16,
    3,
    3,
    18,
    23,
    20,
    21,
    3,
    3
  };

  static const size_t totals[] = {
    3,
    5,
    2,
    4,
    4,
    5,
    5,
    4,
    7,
    5,
    4,
    5,
    6,
    4,
    5,
    6
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char p[32];
  unsigned char q[32];
  unsigned char r0[32];
  unsigned char r1[32];
  size_t i, j, total;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(images); i++) {
    ristretto_pubkey_from_uniform(ec, p, bytes[i]);

    ASSERT(torsion_memcmp(p, images[i], 32) == 0);

    ASSERT(ristretto_pubkey_to_uniform(ec, r0, p, hints[i]));
    ASSERT(torsion_memcmp(r0, bytes[i], 32) == 0);

    total = 0;

    for (j = 0; j < 8; j++) {
      if (ristretto_pubkey_to_uniform(ec, r1, p, j)) {
        ristretto_pubkey_from_uniform(ec, q, r1);

        ASSERT(torsion_memcmp(q, p, 32) == 0);

        total += 1;
      }
    }

    ASSERT(total == totals[i]);
  }

  edwards_curve_destroy(ec);
}

static void
test_ristretto_elligator_hash(drbg_t *rng) {
  /* https://ristretto.group/test_vectors/ristretto255.html */
  static const char *labels[] = {
    "Ristretto is traditionally a short shot of espresso coffee",
    "made with the normal amount of ground coffee but extracted with",
    "about half the amount of water in the same amount of time",
    "by using a finer grind.",
    "This produces a concentrated shot of coffee per volume.",
    "Just pulling a normal shot short will produce a weaker shot",
    "and is not a Ristretto as some believe."
  };

  static const unsigned char images[][32] = {
    {
      0x30, 0x66, 0xf8, 0x2a, 0x1a, 0x74, 0x7d, 0x45,
      0x12, 0x0d, 0x17, 0x40, 0xf1, 0x43, 0x58, 0x53,
      0x1a, 0x8f, 0x04, 0xbb, 0xff, 0xe6, 0xa8, 0x19,
      0xf8, 0x6d, 0xfe, 0x50, 0xf4, 0x4a, 0x0a, 0x46
    },
    {
      0xf2, 0x6e, 0x5b, 0x6f, 0x7d, 0x36, 0x2d, 0x2d,
      0x2a, 0x94, 0xc5, 0xd0, 0xe7, 0x60, 0x2c, 0xb4,
      0x77, 0x3c, 0x95, 0xa2, 0xe5, 0xc3, 0x1a, 0x64,
      0xf1, 0x33, 0x18, 0x9f, 0xa7, 0x6e, 0xd6, 0x1b
    },
    {
      0x00, 0x6c, 0xcd, 0x2a, 0x9e, 0x68, 0x67, 0xe6,
      0xa2, 0xc5, 0xce, 0xa8, 0x3d, 0x33, 0x02, 0xcc,
      0x9d, 0xe1, 0x28, 0xdd, 0x2a, 0x9a, 0x57, 0xdd,
      0x8e, 0xe7, 0xb9, 0xd7, 0xff, 0xe0, 0x28, 0x26
    },
    {
      0xf8, 0xf0, 0xc8, 0x7c, 0xf2, 0x37, 0x95, 0x3c,
      0x58, 0x90, 0xae, 0xc3, 0x99, 0x81, 0x69, 0x00,
      0x5d, 0xae, 0x3e, 0xca, 0x1f, 0xbb, 0x04, 0x54,
      0x8c, 0x63, 0x59, 0x53, 0xc8, 0x17, 0xf9, 0x2a
    },
    {
      0xae, 0x81, 0xe7, 0xde, 0xdf, 0x20, 0xa4, 0x97,
      0xe1, 0x0c, 0x30, 0x4a, 0x76, 0x5c, 0x17, 0x67,
      0xa4, 0x2d, 0x6e, 0x06, 0x02, 0x97, 0x58, 0xd2,
      0xd7, 0xe8, 0xef, 0x7c, 0xc4, 0xc4, 0x11, 0x79
    },
    {
      0xe2, 0x70, 0x56, 0x52, 0xff, 0x9f, 0x5e, 0x44,
      0xd3, 0xe8, 0x41, 0xbf, 0x1c, 0x25, 0x1c, 0xf7,
      0xdd, 0xdb, 0x77, 0xd1, 0x40, 0x87, 0x0d, 0x1a,
      0xb2, 0xed, 0x64, 0xf1, 0xa9, 0xce, 0x86, 0x28
    },
    {
      0x80, 0xbd, 0x07, 0x26, 0x25, 0x11, 0xcd, 0xde,
      0x48, 0x63, 0xf8, 0xa7, 0x43, 0x4c, 0xef, 0x69,
      0x67, 0x50, 0x68, 0x1c, 0xb9, 0x51, 0x0e, 0xea,
      0x55, 0x70, 0x88, 0xf7, 0x6d, 0x9e, 0x50, 0x65
    }
  };

  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char entropy[ENTROPY_SIZE];
  unsigned char bytes[64];
  unsigned char p[32];
  unsigned char q[32];
  sha512_t hash;
  size_t i;

  drbg_generate(rng, entropy, ENTROPY_SIZE);

  for (i = 0; i < ARRAY_SIZE(images); i++) {
    sha512_init(&hash);
    sha512_update(&hash, labels[i], strlen(labels[i]));
    sha512_final(&hash, bytes);

    ristretto_pubkey_from_hash(ec, p, bytes);

    ASSERT(torsion_memcmp(p, images[i], 32) == 0);

    ASSERT(ristretto_pubkey_to_hash(ec, bytes, p, entropy));

    ristretto_pubkey_from_hash(ec, q, bytes);

    ASSERT(torsion_memcmp(q, p, 32) == 0);
  }

  edwards_curve_destroy(ec);
}

static void
test_ristretto_tweak(drbg_t *rng) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char k1[32], k2[32], k3[32];
  unsigned char p1[32], p2[32], p3[32];
  unsigned char tweak[32], entropy[ENTROPY_SIZE];
  const unsigned char *points[2];

  drbg_generate(rng, entropy, ENTROPY_SIZE);

  ristretto_privkey_generate(ec, k1, entropy);

  drbg_generate(rng, entropy, ENTROPY_SIZE);

  ristretto_privkey_generate(ec, tweak, entropy);

  ASSERT(ristretto_pubkey_create(ec, p1, k1));
  ASSERT(ristretto_privkey_tweak_add(ec, k2, k1, tweak));
  ASSERT(ristretto_pubkey_tweak_add(ec, p2, p1, tweak));
  ASSERT(ristretto_pubkey_create(ec, p3, k2));
  ASSERT(torsion_memcmp(p3, p2, 32) == 0);

  ASSERT(ristretto_privkey_negate(ec, tweak, tweak));
  ASSERT(ristretto_privkey_tweak_add(ec, k3, k2, tweak));
  ASSERT(ristretto_pubkey_tweak_add(ec, p3, p2, tweak));
  ASSERT(torsion_memcmp(k3, k1, 32) == 0);
  ASSERT(torsion_memcmp(p3, p1, 32) == 0);

  ASSERT(ristretto_pubkey_create(ec, p1, k1));
  ASSERT(ristretto_privkey_tweak_mul(ec, k2, k1, tweak));
  ASSERT(ristretto_pubkey_tweak_mul(ec, p2, p1, tweak));
  ASSERT(ristretto_pubkey_create(ec, p3, k2));
  ASSERT(torsion_memcmp(p3, p2, 32) == 0);

  ASSERT(ristretto_privkey_invert(ec, tweak, tweak));
  ASSERT(ristretto_privkey_tweak_mul(ec, k3, k2, tweak));
  ASSERT(ristretto_pubkey_tweak_mul(ec, p3, p2, tweak));
  ASSERT(torsion_memcmp(k3, k1, 32) == 0);
  ASSERT(torsion_memcmp(p3, p1, 32) == 0);

  ASSERT(ristretto_privkey_negate(ec, k3, k2));
  ASSERT(ristretto_privkey_tweak_add(ec, k2, k2, k3));
  ASSERT(ristretto_privkey_is_zero(ec, k2));

  ASSERT(ristretto_pubkey_negate(ec, p3, p2));

  points[0] = p2;
  points[1] = p3;

  ASSERT(ristretto_pubkey_combine(ec, p2, points, 2));
  ASSERT(ristretto_pubkey_is_infinity(ec, p2));

  edwards_curve_destroy(ec);
}

static void
test_ristretto_derive(drbg_t *rng) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char alice_priv[32], alice_pub[32], alice_secret[32];
  unsigned char bob_priv[32], bob_pub[32], bob_secret[32];
  unsigned char alice_raw[32], entropy[64];

  drbg_generate(rng, entropy, ENTROPY_SIZE);

  ASSERT(ristretto_privkey_size(ec) == 32);
  ASSERT(ristretto_pubkey_size(ec) == 32);

  ristretto_privkey_generate(ec, alice_priv, entropy);

  drbg_generate(rng, entropy, 64);

  ristretto_privkey_from_uniform(ec, bob_priv, entropy);

  ASSERT(ristretto_privkey_verify(ec, alice_priv));
  ASSERT(ristretto_privkey_verify(ec, bob_priv));

  ASSERT(!ristretto_privkey_is_zero(ec, alice_priv));
  ASSERT(!ristretto_privkey_is_zero(ec, bob_priv));

  ASSERT(ristretto_pubkey_create(ec, alice_pub, alice_priv));
  ASSERT(ristretto_pubkey_create(ec, bob_pub, bob_priv));

  ASSERT(ristretto_pubkey_verify(ec, alice_pub));
  ASSERT(ristretto_pubkey_verify(ec, bob_pub));

  ASSERT(!ristretto_pubkey_is_infinity(ec, alice_pub));
  ASSERT(!ristretto_pubkey_is_infinity(ec, bob_pub));

  ASSERT(ristretto_derive(ec, alice_secret, bob_pub, alice_priv));
  ASSERT(ristretto_derive(ec, bob_secret, alice_pub, bob_priv));

  ASSERT(torsion_memcmp(alice_secret, bob_secret, 32) == 0);

  ASSERT(ristretto_privkey_export(ec, alice_raw, alice_priv));

  ASSERT(torsion_memcmp(alice_raw, alice_priv, 32) == 0);

  ASSERT(ristretto_privkey_import(ec, alice_priv, alice_raw, 32));

  ASSERT(torsion_memcmp(alice_priv, alice_raw, 32) == 0);

  ASSERT(!ristretto_privkey_import(ec, alice_priv, entropy, 64));

  edwards_curve_destroy(ec);
}

/*
 * Encoding
 */

static void
test_encoding_base16(drbg_t *unused) {
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
  unsigned int i;
  size_t len;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *hex = vectors[i][1];

    printf("  - Base16 vector #%u\n", i + 1);

    ASSERT(base16_test(hex, strlen(hex)));
    ASSERT(base16_decode(buf, &len, hex, strlen(hex)));
    ASSERT(len == strlen(str) && torsion_memcmp(buf, str, len) == 0);

    base16_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(hex) && torsion_memcmp(buf, hex, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(vectors_le); i++) {
    const char *str = vectors_le[i][0];
    const char *hex = vectors_le[i][1];

    printf("  - Base16-LE vector #%u\n", i + 1);

    ASSERT(base16le_test(hex, strlen(hex)));
    ASSERT(base16le_decode(buf, &len, hex, strlen(hex)));
    ASSERT(len == strlen(str) && torsion_memcmp(buf, str, len) == 0);

    base16le_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(hex) && torsion_memcmp(buf, hex, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *hex = invalid[i];

    printf("  - Base16 (invalid) vector #%u\n", i + 1);

    ASSERT(!base16_test(hex, strlen(hex)));
    ASSERT(!base16le_test(hex, strlen(hex)));
    ASSERT(!base16_decode(buf, &len, hex, strlen(hex)));
    ASSERT(!base16le_decode(buf, &len, hex, strlen(hex)));
  }
}

static void
test_encoding_base32(drbg_t *unused) {
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
  unsigned int i;
  size_t len;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *b32 = vectors[i][1];

    printf("  - Base32 vector #%u\n", i + 1);

    ASSERT(base32_test(b32, strlen(b32), 1));
    ASSERT(base32_decode(buf, &len, b32, strlen(b32), 1));
    ASSERT(len == strlen(str) && torsion_memcmp(buf, str, len) == 0);

    base32_encode((char *)buf, &len, (unsigned char *)str, strlen(str), 1);
    ASSERT(len == strlen(b32) && torsion_memcmp(buf, b32, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(vectors_hex); i++) {
    const char *str = vectors_hex[i][0];
    const char *b32 = vectors_hex[i][1];

    printf("  - Base32-Hex vector #%u\n", i + 1);

    ASSERT(base32hex_test(b32, strlen(b32), 1));
    ASSERT(base32hex_decode(buf, &len, b32, strlen(b32), 1));
    ASSERT(len == strlen(str) && torsion_memcmp(buf, str, len) == 0);

    base32hex_encode((char *)buf, &len, (unsigned char *)str, strlen(str), 1);
    ASSERT(len == strlen(b32) && torsion_memcmp(buf, b32, len) == 0);
  }
}

static void
test_encoding_base58(drbg_t *unused) {
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
  unsigned int i;
  size_t len, dlen;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *b58 = vectors[i][1];

    ASSERT(base16_decode(data, &dlen, str, strlen(str)));

    printf("  - Base58 vector #%u\n", i + 1);

    ASSERT(base58_test(b58, strlen(b58)));
    ASSERT(base58_decode(buf, &len, b58, strlen(b58)));
    ASSERT(len == dlen && torsion_memcmp(buf, data, len) == 0);

    ASSERT(base58_encode((char *)buf, &len, data, dlen));
    ASSERT(len == strlen(b58) && torsion_memcmp(buf, b58, len) == 0);
  }
}

static void
test_encoding_base64(drbg_t *unused) {
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
  unsigned int i;
  size_t len;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i][0];
    const char *b64 = vectors[i][1];

    printf("  - Base64 vector #%u\n", i + 1);

    ASSERT(base64_test(b64, strlen(b64)));
    ASSERT(base64_decode(buf, &len, b64, strlen(b64)));
    ASSERT(len == strlen(str) && torsion_memcmp(buf, str, len) == 0);

    base64_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(b64) && torsion_memcmp(buf, b64, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *b64 = invalid[i];

    printf("  - Base64 (invalid) vector #%u\n", i + 1);

    ASSERT(!base64_test(b64, strlen(b64)));
    ASSERT(!base64_decode(buf, &len, b64, strlen(b64)));
  }

  for (i = 0; i < ARRAY_SIZE(vectors_url); i++) {
    const char *str = vectors_url[i][0];
    const char *b64 = vectors_url[i][1];

    printf("  - Base64-URL vector #%u\n", i + 1);

    ASSERT(base64url_test(b64, strlen(b64)));
    ASSERT(base64url_decode(buf, &len, b64, strlen(b64)));
    ASSERT(len == strlen(str) && torsion_memcmp(buf, str, len) == 0);

    base64url_encode((char *)buf, &len, (unsigned char *)str, strlen(str));
    ASSERT(len == strlen(b64) && torsion_memcmp(buf, b64, len) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid_url); i++) {
    const char *b64 = invalid_url[i];

    printf("  - Base64-URL (invalid) vector #%u\n", i + 1);

    ASSERT(!base64url_test(b64, strlen(b64)));
    ASSERT(!base64url_decode(buf, &len, b64, strlen(b64)));
  }
}

static void
test_encoding_bech32(drbg_t *unused) {
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
  size_t script_len, data_len;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(valid); i++) {
    const char *addr = valid[i][0];
    const char *hex = valid[i][1];

    printf("  - Bech32 vector #%u (%s)\n", i + 1, addr);

    ASSERT(base16_decode(script, &script_len, hex, strlen(hex)));

    ASSERT(bech32_test(addr));
    ASSERT(bech32_is(addr));
    ASSERT(bech32_decode(hrp, &version, data, &data_len, addr));
    ASSERT(strlen(hrp) >= 2 && torsion_memcmp(hrp, expect[i], 2) == 0);
    ASSERT(2 + data_len == script_len);
    ASSERT(script[0] == (version ? version + 0x80 : 0));
    ASSERT(script[1] == data_len);
    ASSERT(torsion_memcmp(data, script + 2, data_len) == 0);

    ASSERT(bech32_encode(out, hrp, version, data, data_len));
    ASSERT(strcmp(out, expect[i]) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *addr = invalid[i];

    printf("  - Bech32 (invalid) vector #%u (%s)\n", i + 1, addr);

    ASSERT(!bech32_test(addr));
    ASSERT(!bech32_decode(hrp, &version, data, &data_len, addr));
  }
}

static void
test_encoding_cash32(drbg_t *unused) {
  static const struct {
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
  size_t data_len, hash_len;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(valid); i++) {
    const char *addr = valid[i].addr;
    const char *prefix = valid[i].prefix;
    const char *hex = valid[i].hash;

    printf("  - Cash32 vector #%u (%s)\n", i + 1, addr);

    ASSERT(base16_decode(hash, &hash_len, hex, strlen(hex)));

    ASSERT(cash32_test(addr, prefix));
    ASSERT(cash32_is(addr, prefix));
    ASSERT(cash32_decode(&type, data, &data_len, addr, prefix));
    ASSERT(type == valid[i].type);
    ASSERT(data_len == hash_len && torsion_memcmp(data, hash, hash_len) == 0);

    ASSERT(cash32_encode(out, prefix, type, data, data_len));
    ASSERT(strcmp(out, addr) == 0);
  }

  for (i = 0; i < ARRAY_SIZE(invalid); i++) {
    const char *addr = invalid[i][0];
    const char *prefix = invalid[i][1];

    printf("  - Cash32 (invalid) vector #%u (%s)\n", i + 1, addr);

    ASSERT(!cash32_test(addr, prefix));
    ASSERT(!cash32_decode(&type, data, &data_len, addr, prefix));
  }
}

/*
 * Hash
 */

static void
test_hash_digest(drbg_t *unused) {
  unsigned char iv[256];
  unsigned char expect[HASH_MAX_OUTPUT_SIZE];
  unsigned char out[HASH_MAX_OUTPUT_SIZE];
  hash_t hash;
  unsigned int i, j;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(hash_vectors); i++) {
    int type = hash_vectors[i].type;
    size_t size = hash_output_size(type);
    size_t iv_len = sizeof(iv);

    ASSERT(size != 0);

    printf("  - Hash vector #%u (%s)\n", i + 1, hash_names[type]);

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

    ASSERT(torsion_memcmp(out, expect, size) == 0);
  }
}

static void
test_hash_hmac(drbg_t *unused) {
  unsigned char data[256];
  unsigned char key[256];
  unsigned char expect[HASH_MAX_OUTPUT_SIZE];
  unsigned char out[HASH_MAX_OUTPUT_SIZE];
  hmac_t hmac;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(hmac_vectors); i++) {
    int type = hmac_vectors[i].type;
    size_t size = hash_output_size(type);
    size_t data_len = sizeof(data);
    size_t key_len = sizeof(key);

    ASSERT(size != 0);

    printf("  - HMAC vector #%u (%s)\n", i + 1, hash_names[type]);

    hex_decode(data, &data_len, hmac_vectors[i].data);
    hex_decode(key, &key_len, hmac_vectors[i].key);
    hex_parse(expect, size, hmac_vectors[i].expect);

    hmac_init(&hmac, type, key, key_len);
    hmac_update(&hmac, data, data_len);
    hmac_final(&hmac, out);

    ASSERT(torsion_memcmp(out, expect, size) == 0);
  }
}

/*
 * IES
 */

static void
test_ies_secretbox(drbg_t *unused) {
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

  (void)unused;

  ASSERT(sizeof(sealed) == 80);

  /* Raw. */
  memset(key, 1, sizeof(key));
  memset(nonce, 2, sizeof(nonce));
  memset(msg, 3, sizeof(msg));

  secretbox_seal(sealed, msg, sizeof(msg), key, nonce);

  ASSERT(secretbox_open(opened, sealed, sizeof(sealed), key, nonce));

  ASSERT(torsion_memcmp(sealed, expect1, sizeof(expect1)) == 0);
  ASSERT(torsion_memcmp(opened, msg, sizeof(msg)) == 0);

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

  ASSERT(torsion_memcmp(sealed, expect2, sizeof(expect2)) == 0);
  ASSERT(torsion_memcmp(opened, msg, sizeof(msg)) == 0);

  mont_curve_destroy(ec);
}

static void
test_ies_secretbox_random(drbg_t *rng) {
  unsigned char sealed[SECRETBOX_SEAL_SIZE(128)];
  unsigned char opened[128];
  unsigned char msg[128];
  unsigned char key[32];
  unsigned char nonce[24];
  size_t i, len, last;

  drbg_generate(rng, key, 32);
  drbg_generate(rng, nonce, 24);

  for (len = 0; len < 128; len += 17) {
    drbg_generate(rng, msg, len);

    secretbox_seal(sealed, msg, len, key, nonce);

    if (len > 0) {
      ASSERT(torsion_memcmp(sealed, msg, len) != 0);
      ASSERT(torsion_memcmp(sealed + 16, msg, len) != 0);
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
 * KDF
 */

static void
test_kdf_bcrypt(drbg_t *unused) {
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

  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *pass = vectors[i].pass;
    unsigned int rounds = vectors[i].rounds;
    const char *salt = vectors[i].salt;
    const char *record = vectors[i].record;
    char out[100];

    printf("  - Bcrypt vector #%u (%s)\n", i + 1, record);

    ASSERT(bcrypt_generate_with_salt64(out, (const unsigned char *)pass,
                                       strlen(pass), salt, rounds, 'a'));

    ASSERT(strcmp(out, record) == 0);

    ASSERT(bcrypt_verify((const unsigned char *)pass, strlen(pass), record));
  }
}

static void
test_kdf_bcrypt_pbkdf(drbg_t *unused) {
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

  (void)unused;

  ASSERT(bcrypt_pbkdf(out, pass, sizeof(pass) - 1, salt, sizeof(salt), 16, 48));
  ASSERT(torsion_memcmp(out, expect, 48) == 0);
}

static void
test_kdf_eb2k(drbg_t *unused) {
  unsigned char salt[64];
  unsigned char key[64];
  unsigned char iv[64];
  unsigned char out1[64];
  unsigned char out2[64];
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(eb2k_vectors); i++) {
    const unsigned char *pass = (const unsigned char *)eb2k_vectors[i].pass;
    size_t pass_len = strlen(eb2k_vectors[i].pass);
    size_t salt_len = sizeof(salt);
    size_t key_len = sizeof(key);
    size_t iv_len = sizeof(iv);

    printf("  - EB2K vector #%u (%s)\n", i + 1, pass);

    hex_decode(salt, &salt_len, eb2k_vectors[i].salt);
    hex_decode(key, &key_len, eb2k_vectors[i].key);
    hex_decode(iv, &iv_len, eb2k_vectors[i].iv);

    ASSERT(sizeof(out1) >= key_len);
    ASSERT(sizeof(out2) >= iv_len);

    ASSERT(eb2k_derive(out1, out2, HASH_MD5, pass, pass_len,
                       salt, salt_len, key_len, iv_len));

    ASSERT(torsion_memcmp(out1, key, key_len) == 0);
    ASSERT(torsion_memcmp(out2, iv, iv_len) == 0);
  }
}

static void
test_kdf_hkdf(drbg_t *unused) {
  unsigned char ikm[128];
  unsigned char salt[128];
  unsigned char info[128];
  unsigned char prk[HASH_MAX_OUTPUT_SIZE];
  unsigned char okm[128];
  unsigned char out[128];
  unsigned int i;

  (void)unused;

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

    printf("  - HKDF vector #%u (%s)\n", i + 1, hkdf_vectors[i].okm);

    hex_decode(ikm, &ikm_len, hkdf_vectors[i].ikm);
    hex_decode(salt, &salt_len, hkdf_vectors[i].salt);
    hex_decode(info, &info_len, hkdf_vectors[i].info);
    hex_parse(prk, size, hkdf_vectors[i].prk);
    hex_parse(okm, len, hkdf_vectors[i].okm);

    ASSERT(hkdf_extract(out, type, ikm, ikm_len, salt, salt_len));
    ASSERT(torsion_memcmp(out, prk, size) == 0);

    ASSERT(hkdf_expand(out, type, prk, info, info_len, len));
    ASSERT(torsion_memcmp(out, okm, len) == 0);
  }
}

static void
test_kdf_pbkdf2(drbg_t *unused) {
  unsigned char pass[256];
  unsigned char salt[256];
  unsigned char expect[256];
  unsigned char out[256];
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(pbkdf2_vectors); i++) {
    int type = pbkdf2_vectors[i].type;
    size_t pass_len = sizeof(pass);
    size_t salt_len = sizeof(salt);
    unsigned int iter = pbkdf2_vectors[i].iter;
    size_t len = pbkdf2_vectors[i].len;

    ASSERT(sizeof(expect) >= len);
    ASSERT(sizeof(out) >= len);

    printf("  - PBKDF2 vector #%u (%s)\n", i + 1, pbkdf2_vectors[i].expect);

    hex_decode(pass, &pass_len, pbkdf2_vectors[i].pass);
    hex_decode(salt, &salt_len, pbkdf2_vectors[i].salt);
    hex_parse(expect, len, pbkdf2_vectors[i].expect);

    ASSERT(pbkdf2_derive(out, type, pass, pass_len, salt, salt_len, iter, len));
    ASSERT(torsion_memcmp(out, expect, len) == 0);
  }
}

static void
test_kdf_pgpdf(drbg_t *unused) {
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

  (void)unused;

  ASSERT(pgpdf_derive_simple(out, HASH_SHA256, pass, pass_len, 32));
  ASSERT(torsion_memcmp(out, expect1, 32) == 0);

  ASSERT(pgpdf_derive_salted(out, HASH_SHA256, pass, pass_len,
                             salt, salt_len, 32));
  ASSERT(torsion_memcmp(out, expect2, 32) == 0);

  ASSERT(pgpdf_derive_iterated(out, HASH_SHA256, pass, pass_len,
                               salt, salt_len, 100, 32));
  ASSERT(torsion_memcmp(out, expect3, 32) == 0);
}

static void
test_kdf_scrypt(drbg_t *unused) {
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

  (void)unused;

  ASSERT(scrypt_derive(out, NULL, 0, NULL, 0, 16, 1, 1, 64));
  ASSERT(torsion_memcmp(out, expect1, 64) == 0);

  ASSERT(scrypt_derive(out, pass2, sizeof(pass2) - 1,
                            salt2, sizeof(salt2) - 1, 1024, 8, 16, 64));

  ASSERT(torsion_memcmp(out, expect2, 64) == 0);

  ASSERT(scrypt_derive(out, pass3, sizeof(pass3) - 1,
                            salt3, sizeof(salt3) - 1, 16384, 8, 1, 64));

  ASSERT(torsion_memcmp(out, expect3, 64) == 0);
}

/*
 * MAC
 */

static void
test_mac_poly1305(drbg_t *unused) {
  unsigned char key[32];
  unsigned char msg[8192];
  unsigned char tag[16];
  unsigned char mac[16];
  poly1305_t ctx;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(poly1305_vectors); i++) {
    size_t msg_len = sizeof(msg);

    printf("  - Poly1305 vector #%u (%s)\n", i + 1, poly1305_vectors[i][0]);

    hex_parse(key, 32, poly1305_vectors[i][0]);
    hex_decode(msg, &msg_len, poly1305_vectors[i][1]);
    hex_parse(tag, 16, poly1305_vectors[i][2]);

    poly1305_init(&ctx, key);
    poly1305_update(&ctx, msg, msg_len);
    poly1305_final(&ctx, mac);

    ASSERT(torsion_memequal(mac, tag, 16));
    ASSERT(torsion_memcmp(mac, tag, 16) == 0);
  }
}

static void
test_mac_siphash(drbg_t *unused) {
  static const uint32_t v32 = UINT32_C(0x9dcb553a);
  static const uint64_t v64 = UINT64_C(0x73b4e2ae9316f6b2);
  unsigned char msg[32];
  unsigned char key[16];
  size_t i;

  (void)unused;

  for (i = 0; i < 32; i++)
    msg[i] = i;

  for (i = 0; i < 16; i++)
    key[i] = i + 32;

  ASSERT(siphash_sum(msg, 32, key) == UINT64_C(10090947469682793545));
  ASSERT(siphash_mod(msg, 32, key, v64) == UINT64_C(4560894765423557143));
  ASSERT((uint32_t)siphash128_sum(v32, key) == UINT32_C(828368916));
  ASSERT(siphash128_sum(v64, key) == UINT64_C(620914895672640125));
  ASSERT((uint32_t)siphash256_sum(v32, msg) == UINT32_C(3909845928));
  ASSERT(siphash256_sum(v64, msg) == UINT64_C(16928650368383018294));
}

/*
 * MPI
 */

typedef void mp_rng_f(void *out, size_t size, void *arg);

void
__torsion_test_mpi_internal(mp_rng_f *rng, void *arg);

static void
test_mpi_internal(drbg_t *rng) {
  __torsion_test_mpi_internal(drbg_rng, rng);
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

  ASSERT(size != 0);

  for (i = 0; i < size; i++) {
    for (j = 0; j < 8; j++)
      sum += (raw[i] >> j) & 1;
  }

  avg = (sum * 100) / (size * 8);

  return avg >= 48 && avg <= 52;
}

static void
test_rand_getentropy(drbg_t *unused) {
  unsigned char data[16384];
  size_t i;

  (void)unused;

  for (i = 0; i < 10; i++) {
    ASSERT(torsion_getentropy(data, sizeof(data)));
    ASSERT(looks_random(data, sizeof(data)));
  }
}

static void
test_rand_getrandom(drbg_t *unused) {
  unsigned char data[16384];
  size_t i;

  (void)unused;

  for (i = 0; i < 10; i++) {
    ASSERT(torsion_getrandom(data, sizeof(data)));
    ASSERT(looks_random(data, sizeof(data)));
  }
}

#ifdef TORSION_HAVE_ZLIB
static size_t
rand_deflate_perc(const void *data, size_t size) {
  size_t len = compressBound(size);
  unsigned char *buf = malloc(len);

  ASSERT(size != 0);
  ASSERT(buf != NULL);
  ASSERT(compress2(buf, &len, data, size, 5) == Z_OK);

  free(buf);

  return (len * 100) / size;
}

static void
test_rand_deflate_sanity(drbg_t *unused) {
  size_t size = (size_t)4 << 20;
  unsigned char *data = malloc(size);

  (void)unused;

  ASSERT(data != NULL);

  memset(data, 0xaa, size);

  ASSERT(rand_deflate_perc(data, size) <= 1);

  free(data);
}

static void
test_rand_getentropy_deflate(drbg_t *unused) {
  size_t size = (size_t)4 << 20;
  unsigned char *data = malloc(size);

  (void)unused;

  ASSERT(data != NULL);
  ASSERT(torsion_getentropy(data, size));
  ASSERT(rand_deflate_perc(data, size) >= 99);

  free(data);
}

static void
test_rand_getrandom_deflate(drbg_t *unused) {
  size_t size = (size_t)4 << 20;
  unsigned char *data = malloc(size);

  (void)unused;

  ASSERT(data != NULL);
  ASSERT(torsion_getrandom(data, size));
  ASSERT(rand_deflate_perc(data, size) >= 99);

  free(data);
}
#endif /* TORSION_HAVE_ZLIB */

static void
test_rand_random(drbg_t *unused) {
  uint32_t data[65536 / 4];
  size_t i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(data); i++)
    ASSERT(torsion_random(&data[i]));

  ASSERT(looks_random(data, sizeof(data)));
}

static void
test_rand_uniform(drbg_t *unused) {
  uint32_t data[65536 / 4];
  uint32_t hi, lo;
  size_t i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(data); i++) {
    ASSERT(torsion_uniform(&hi, 0x10000));
    ASSERT(torsion_uniform(&lo, 0x10000));
    ASSERT(hi < 0x10000 && lo < 0x10000);

    data[i] = (hi << 16) | lo;
  }

  ASSERT(looks_random(data, sizeof(data)));
}

#ifdef TORSION_HAVE_THREADS
typedef struct rng_res_s {
  uint64_t ptr;
  unsigned char data[32];
} rng_res_t;

static void *
thread_random(void *ptr) {
  rng_res_t *obj = ptr;

  obj->ptr = torsion_randomaddr();

  ASSERT(torsion_getrandom(obj->data, 32));

  return NULL;
}

static void
test_rand_thread_safety(drbg_t *unused) {
  const char *name = NULL;
  torsion_thread_t *t1, *t2;
  rng_res_t x0, x1, x2;

  (void)unused;

  memset(&x0, 0, sizeof(x0));
  memset(&x1, 0, sizeof(x1));
  memset(&x2, 0, sizeof(x2));

  switch (torsion_threadsafety()) {
    case TORSION_THREAD_SAFETY_NONE:
      printf("  - Skipping RNG test (not thread safe).\n");
      return;
    case TORSION_THREAD_SAFETY_TLS:
      name = "TLS";
      break;
    case TORSION_THREAD_SAFETY_MUTEX:
      name = "MUTEX";
      break;
    default:
      ASSERT(0);
      break;
  }

  printf("  - Testing thread safety (%s)\n", name);

  t1 = torsion_thread_alloc();
  t2 = torsion_thread_alloc();

  ASSERT(torsion_thread_create(t1, NULL, thread_random, (void *)&x1) == 0);
  ASSERT(torsion_thread_create(t2, NULL, thread_random, (void *)&x2) == 0);

  ASSERT(torsion_thread_join(t1, NULL) == 0);
  ASSERT(torsion_thread_join(t2, NULL) == 0);

  torsion_thread_free(t1);
  torsion_thread_free(t2);

  thread_random(&x0);

  ASSERT(x0.ptr && x1.ptr && x2.ptr);

  if (torsion_threadsafety() == TORSION_THREAD_SAFETY_TLS) {
    ASSERT(x0.ptr != x1.ptr);
    ASSERT(x0.ptr != x2.ptr);

    ASSERT(x1.ptr != x0.ptr);
    ASSERT(x1.ptr != x2.ptr);

    ASSERT(x2.ptr != x0.ptr);
    ASSERT(x2.ptr != x1.ptr);
  }

  ASSERT(torsion_memcmp(x0.data, x1.data, 32) != 0);
  ASSERT(torsion_memcmp(x0.data, x2.data, 32) != 0);

  ASSERT(torsion_memcmp(x1.data, x0.data, 32) != 0);
  ASSERT(torsion_memcmp(x1.data, x2.data, 32) != 0);

  ASSERT(torsion_memcmp(x2.data, x0.data, 32) != 0);
  ASSERT(torsion_memcmp(x2.data, x1.data, 32) != 0);
}
#endif /* TORSION_HAVE_THREADS */

#ifdef TORSION_HAVE_FORK
static void
test_rand_fork_safety(drbg_t *unused) {
  unsigned char ours[32];
  unsigned char theirs[32];
  int pfds[2];
  int status;
  pid_t pid;

  (void)unused;

  memset(ours, 0, 32);
  memset(theirs, 0, 32);

  ASSERT(pipe(pfds) == 0);

  pid = fork();

  ASSERT(pid != -1);
  ASSERT(torsion_getrandom(ours, 32));

  if (pid) {
    ASSERT(close(pfds[1]) == 0);
    ASSERT(read(pfds[0], theirs, 32) == 32);
    ASSERT(close(pfds[0]) == 0);

    ASSERT(waitpid(pid, &status, 0) == pid);
    ASSERT(WIFEXITED(status));
    ASSERT(WEXITSTATUS(status) == 0);
  } else {
    ASSERT(close(pfds[0]) == 0);
    ASSERT(write(pfds[1], ours, 32) == 32);
    ASSERT(close(pfds[1]) == 0);

    exit(0);
  }

  ASSERT(torsion_memcmp(ours, theirs, 32) != 0);
}
#endif /* TORSION_HAVE_FORK */
#endif /* TORSION_HAVE_RNG */

/*
 * RSA
 */

static void
test_rsa_vectors(drbg_t *unused) {
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
  unsigned int i;
  size_t len;

  (void)unused;

  for (i = 0; i < ENTROPY_SIZE; i++)
    entropy[i] = i;

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

    printf("  - RSA vector #%u\n", i + 1);

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
    ASSERT(len == pub_len && torsion_memcmp(out, pub, pub_len) == 0);

    /* RSA-PKCS1v1.5 type 1. */
    ASSERT(rsa_sign(out, &len, hash, msg, msg_len,
                    priv, priv_len, entropy));
    ASSERT(len == sig1_len && torsion_memcmp(out, sig1, sig1_len) == 0);

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
    ASSERT(len == msg_len && torsion_memcmp(out, msg, msg_len) == 0);

    ASSERT(rsa_encrypt(tmp, &len, msg, msg_len, pub, pub_len, entropy));
    ASSERT(rsa_decrypt(out, &len, tmp, len, priv, priv_len, entropy));
    ASSERT(len == msg_len && torsion_memcmp(out, msg, msg_len) == 0);

    /* RSA-OAEP. */
    ASSERT(rsa_decrypt_oaep(out, &len, hash, ct2, ct2_len,
                            priv, priv_len, label, sizeof(label), entropy));
    ASSERT(len == msg_len && torsion_memcmp(out, msg, msg_len) == 0);

    ASSERT(rsa_encrypt_oaep(tmp, &len, hash, msg, msg_len,
                            pub, pub_len, label, sizeof(label), entropy));
    ASSERT(rsa_decrypt_oaep(out, &len, hash, tmp, len,
                            priv, priv_len, label, sizeof(label), entropy));
    ASSERT(len == msg_len && torsion_memcmp(out, msg, msg_len) == 0);
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
    ASSERT(pt_len == 32 && torsion_memcmp(pt, msg, 32) == 0);

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

    ASSERT(pt_len == 32 && torsion_memcmp(pt, msg, 32) == 0);

    ct[j] ^= 1;

    ASSERT(!rsa_decrypt_oaep(pt, &pt_len, HASH_SHA256, ct, ct_len,
                             priv, priv_len, NULL, 0, entropy));
  }
}

/*
 * Stream
 */

static void
test_stream_arc4(drbg_t *unused) {
  static const unsigned char expect[32] = {
    0x42, 0x83, 0x06, 0x31, 0x1d, 0xd2, 0x34, 0x98,
    0x51, 0x3a, 0x18, 0x7c, 0x36, 0x4c, 0x03, 0xc0,
    0x56, 0x7b, 0x5c, 0x82, 0x94, 0x70, 0x29, 0xfa,
    0x1d, 0x26, 0x24, 0x9f, 0x86, 0x25, 0x1a, 0xa0
  };

  unsigned char key[32];
  unsigned char msg[32];
  unsigned char data[32];
  arc4_t ctx;
  size_t i;

  (void)unused;

  for (i = 0; i < 32; i++) {
    key[i] = i + 10;
    msg[i] = i + 20;
  }

  memcpy(data, msg, 32);

  arc4_init(&ctx, key, 32);

  for (i = 0; i < 1000; i++)
    arc4_crypt(&ctx, data, data, 32);

  ASSERT(torsion_memcmp(data, expect, 32) == 0);

  arc4_init(&ctx, key, 32);

  for (i = 0; i < 1000; i++)
    arc4_crypt(&ctx, data, data, 32);

  ASSERT(torsion_memcmp(data, msg, 32) == 0);
}

static void
test_stream_chacha20(drbg_t *unused) {
  unsigned char key[32];
  unsigned char nonce[24];
  unsigned char input[4096];
  unsigned char output[4096];
  unsigned char data[4096];
  chacha20_t ctx;
  unsigned int i;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(chacha20_vectors); i++) {
    size_t key_len = sizeof(key);
    size_t nonce_len = sizeof(nonce);
    unsigned int counter = chacha20_vectors[i].counter;
    size_t input_len = sizeof(input);
    size_t output_len = sizeof(output);
    size_t data_len;

    printf("  - ChaCha20 vector #%u\n", i + 1);

    hex_decode(key, &key_len, chacha20_vectors[i].key);
    hex_decode(nonce, &nonce_len, chacha20_vectors[i].nonce);
    hex_decode(input, &input_len, chacha20_vectors[i].input);
    hex_decode(output, &output_len, chacha20_vectors[i].output);

    ASSERT(input_len == output_len);

    data_len = input_len;

    memcpy(data, input, input_len);

    chacha20_init(&ctx, key, key_len, nonce, nonce_len, counter);
    chacha20_crypt(&ctx, data, data, data_len);

    ASSERT(torsion_memcmp(data, output, output_len) == 0);

    chacha20_init(&ctx, key, key_len, nonce, nonce_len, counter);
    chacha20_crypt(&ctx, data, data, data_len);

    ASSERT(torsion_memcmp(data, input, input_len) == 0);
  }
}

static void
test_stream_salsa20(drbg_t *unused) {
  /* https://github.com/golang/crypto/blob/master/salsa20/salsa20_test.go */
  static const struct {
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
  unsigned int i, j, k;
  salsa20_t ctx;

  (void)unused;

  ASSERT(out != NULL);

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const unsigned char *key = vectors[i].key;
    const unsigned char *nonce = vectors[i].nonce;
    const unsigned char *expect = vectors[i].expect;

    printf("  - Salsa20 vector #%u\n", i + 1);

    memset(out, 0, size);
    memset(xor, 0, 64);

    salsa20_init(&ctx, key, 32, nonce, 8, 0);
    salsa20_crypt(&ctx, out, out, size);

    for (j = 0; j < size; j += 64) {
      for (k = 0; k < 64; k++)
        xor[k] ^= out[j + k];
    }

    ASSERT(torsion_memcmp(xor, expect, 64) == 0);
  }

  free(out);
}

static void
test_stream_xsalsa20(drbg_t *unused) {
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

  (void)unused;

  salsa20_init(&ctx, key, 32, nonce, 24, 0);
  salsa20_crypt(&ctx, data, input, 12);

  ASSERT(torsion_memcmp(data, expect, 12) == 0);
}

/*
 * Util
 */

static void
test_util_cleanse(drbg_t *rng) {
  static const unsigned char zero[32] = {0};
  unsigned char raw[32];

  drbg_generate(rng, raw, 32);

  ASSERT(torsion_memcmp(raw, zero, 32) != 0);

  torsion_cleanse(raw, 32);

  ASSERT(torsion_memcmp(raw, zero, 32) == 0);
}

static void
test_util_memequal(drbg_t *rng) {
  unsigned char s1[32];
  unsigned char s2[32];

  drbg_generate(rng, s1, 32);

  memcpy(s2, s1, 32);

  ASSERT(torsion_memequal(s1, s2, 32));

  s2[drbg_uniform(rng, 32)] ^= 1;

  ASSERT(!torsion_memequal(s1, s2, 32));
}

static void
test_util_murmur3(drbg_t *unused) {
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
  unsigned int i;
  size_t len;

  (void)unused;

  for (i = 0; i < ARRAY_SIZE(vectors); i++) {
    const char *str = vectors[i].str;
    uint32_t seed = vectors[i].seed;
    uint32_t sum = vectors[i].sum;
    uint32_t tweak = vectors[i].tweak;

    printf("  - Murmur3 vector #%u (%lu)\n", i + 1, (unsigned long)seed);

    ASSERT(base16_decode(data, &len, str, strlen(str)));

    ASSERT(murmur3_sum(data, len, seed) == sum);
    ASSERT(murmur3_tweak(data, len, sum, seed) == tweak);
  }
}

/*
 * Test Registry
 */

typedef struct torsion_test_s {
  const char *name;
  void (*run)(drbg_t *);
} torsion_test_t;

static const torsion_test_t torsion_tests[] = {
#define T(name) { #name, test_ ## name }
  /* Memcmp */
  T(memcmp),

  /* AEAD */
  T(aead_chachapoly),

  /* Cipher */
  T(cipher_contexts),
  T(cipher_modes),
  T(cipher_aead),

  /* DRBG */
  T(drbg_hash),
  T(drbg_hmac),
  T(drbg_ctr),

  /* DSA */
  T(dsa_vectors),
  T(dsa_keygen),

  /* ECC */
  T(ecc_internal),
  T(ecdsa_vectors),
  T(ecdsa_random),
  T(ecdsa_sswu),
  T(ecdsa_svdw),
  T(schnorr_legacy_vectors),
  T(schnorr_legacy_random),
  T(schnorr_vectors),
  T(schnorr_random),
  T(ecdh_x25519),
  T(ecdh_x448),
  T(ecdh_random),
  T(ecdh_elligator2),
  T(eddsa_vectors),
  T(eddsa_random),
  T(eddsa_elligator2),
  T(ristretto_basepoint_multiples_ed25519),
  T(ristretto_bad_points_ed25519),
  T(ristretto_basepoint_multiples_ed448),
  T(ristretto_elligator),
  T(ristretto_elligator_hash),
  T(ristretto_tweak),
  T(ristretto_derive),

  /* Encoding */
  T(encoding_base16),
  T(encoding_base32),
  T(encoding_base58),
  T(encoding_base64),
  T(encoding_bech32),
  T(encoding_cash32),

  /* Hash */
  T(hash_digest),
  T(hash_hmac),

  /* IES */
  T(ies_secretbox),
  T(ies_secretbox_random),

  /* KDF */
  T(kdf_bcrypt),
  T(kdf_bcrypt_pbkdf),
  T(kdf_eb2k),
  T(kdf_hkdf),
  T(kdf_pbkdf2),
  T(kdf_pgpdf),
  T(kdf_scrypt),

  /* MAC */
  T(mac_poly1305),
  T(mac_siphash),

  /* MPI */
  T(mpi_internal),

  /* Random */
#ifdef TORSION_HAVE_RNG
  T(rand_getentropy),
  T(rand_getrandom),
#ifdef TORSION_HAVE_ZLIB
  T(rand_deflate_sanity),
  T(rand_getentropy_deflate),
  T(rand_getrandom_deflate),
#endif
  T(rand_random),
  T(rand_uniform),
#ifdef TORSION_HAVE_THREADS
  T(rand_thread_safety),
#endif
#ifdef TORSION_HAVE_FORK
  T(rand_fork_safety),
#endif
#endif /* TORSION_HAVE_RNG */

  /* RSA */
  T(rsa_vectors),
  T(rsa_random),

  /* Stream */
  T(stream_arc4),
  T(stream_chacha20),
  T(stream_salsa20),
  T(stream_xsalsa20),

  /* Util */
  T(util_cleanse),
  T(util_memequal),
  T(util_murmur3)
#undef T
};

/*
 * Test Runner
 */

static int
run_test(const torsion_test_t *test, drbg_t *rng, int child) {
  printf("Testing %s...\n", test->name);

#ifdef TORSION_HAVE_FORK
  if (child) {
    pid_t pid = fork();
    int status;

    ASSERT(pid != -1);

    if (pid) {
      ASSERT(waitpid(pid, &status, 0) == pid);

      if (WIFEXITED(status) == 0)
        return 0;

      if (WEXITSTATUS(status) != 0)
        return 0;

      return 1;
    }

    drbg_init_rand(rng);

    test->run(rng);

    exit(0);

    return 0;
  }
#endif

  (void)child;

  test->run(rng);

  return 1;
}

/*
 * Main
 */

int
main(int argc, char **argv) {
  size_t args_len = 0;
  char *args[256];
  drbg_t rng;
  int grep = 0;
  int child = 0;
  int total = 0;
  int passing = 0;
  int found;
  size_t i, j;

  if ((size_t)argc > ARRAY_SIZE(args)) {
    fprintf(stderr, "Too many arguments.\n");
    return 1;
  }

  for (i = 1; i < (size_t)argc; i++) {
    if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--grep") == 0) {
      grep = 1;
      continue;
    }

    if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--fork") == 0) {
      child = 1;
      continue;
    }

    if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--list") == 0) {
      for (j = 0; j < ARRAY_SIZE(torsion_tests); j++)
        printf("%s\n", torsion_tests[j].name);
      return 0;
    }

    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      printf("%s",
        "\n"
        "  Usage: ./torsion_test [options] [test-name]\n"
        "\n"
        "  Options:\n"
        "\n"
        "    -g, --grep   treat [test-name] as a substring\n"
        "    -f, --fork   run each test in a subprocess\n"
        "    -l, --list   list all test names\n"
        "    -h, --help   output usage information\n"
        "\n"
      );
      return 0;
    }

    if (argv[i][0] == '-') {
      fprintf(stderr, "Unknown argument: %s.\n", argv[i]);
      return 1;
    }

    args[args_len++] = argv[i];
  }

  drbg_init_rand(&rng);

  if (grep) {
    for (i = 0; i < args_len; i++) {
      for (j = 0; j < ARRAY_SIZE(torsion_tests); j++) {
        if (strstr(torsion_tests[j].name, args[i]) != NULL) {
          passing += run_test(&torsion_tests[j], &rng, child);
          total += 1;
        }
      }
    }

    if (total == 0) {
      fprintf(stderr, "No tests found.\n");
      return 1;
    }

    goto done;
  }

  if (args_len == 0) {
    for (i = 0; i < ARRAY_SIZE(torsion_tests); i++) {
      passing += run_test(&torsion_tests[i], &rng, child);
      total += 1;
    }

    goto done;
  }

  for (i = 0; i < args_len; i++) {
    found = 0;

    for (j = 0; j < ARRAY_SIZE(torsion_tests); j++) {
      if (strcmp(torsion_tests[j].name, args[i]) == 0) {
        passing += run_test(&torsion_tests[j], &rng, child);
        total += 1;
        found = 1;
        break;
      }
    }

    if (!found) {
      fprintf(stderr, "Unknown test: %s.\n", args[i]);
      return 1;
    }
  }

done:
  if (passing != total) {
    printf("Tests failed (%d/%d).\n", passing, total);
    return 1;
  }

  printf("Tests passed (%d/%d).\n", passing, total);

  return 0;
}
