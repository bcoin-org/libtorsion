#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <openssl/sha.h>
#include "keccak.h"

#define HASH_SHA256 0
#define HASH_SHA384 1
#define HASH_SHA512 2
#define HASH_SHAKE256 3

#define MAX_HASH_SIZE 64
#define MAX_BLOCK_SIZE 128

/*
 * Structs
 */

typedef struct hash_s {
  int type;
  union {
    SHA256_CTX sha256;
    SHA512_CTX sha384;
    SHA512_CTX sha512;
    keccak_t keccak;
  } ctx;
} hash_t;

typedef struct hmac_s {
  int type;
  hash_t inner;
  hash_t outer;
} hmac_t;

typedef struct drbg_s {
  int type;
  hmac_t kmac;
  unsigned char K[MAX_HASH_SIZE];
  unsigned char V[MAX_HASH_SIZE];
} drbg_t;

/*
 * Constants
 */

static const unsigned char ZERO[1] = {0x00};
static const unsigned char ONE[1] = {0x01};

/*
 * Hash
 */

static void
hash_init(hash_t *hash, int type) {
  hash->type = type;
  switch (hash->type) {
    case HASH_SHA256:
      SHA256_Init(&hash->ctx.sha256);
      break;
    case HASH_SHA384:
      SHA384_Init(&hash->ctx.sha384);
      break;
    case HASH_SHA512:
      SHA512_Init(&hash->ctx.sha512);
      break;
    case HASH_SHAKE256:
      keccak_init(&hash->ctx.keccak, 256);
      break;
    default:
      assert(0);
      break;
  }
}

static void
hash_update(hash_t *hash, const void *data, size_t len) {
  switch (hash->type) {
    case HASH_SHA256:
      SHA256_Update(&hash->ctx.sha256, data, len);
      break;
    case HASH_SHA384:
      SHA384_Update(&hash->ctx.sha384, data, len);
      break;
    case HASH_SHA512:
      SHA512_Update(&hash->ctx.sha512, data, len);
      break;
    case HASH_SHAKE256:
      keccak_update(&hash->ctx.keccak, data, len);
      break;
    default:
      assert(0);
      break;
  }
}

static void
hash_final(hash_t *hash, unsigned char *out, size_t len) {
  switch (hash->type) {
    case HASH_SHA256:
      SHA256_Final(out, &hash->ctx.sha256);
      break;
    case HASH_SHA384:
      SHA384_Final(out, &hash->ctx.sha384);
      break;
    case HASH_SHA512:
      SHA512_Final(out, &hash->ctx.sha512);
      break;
    case HASH_SHAKE256:
      keccak_final(&hash->ctx.keccak, out, 0x1f, len);
      break;
    default:
      assert(0);
      break;
  }
}

static size_t
hash_output_size(int type) {
  switch (type) {
    case HASH_SHA256:
      return 32;
    case HASH_SHA384:
      return 48;
    case HASH_SHA512:
      return 64;
    case HASH_SHAKE256:
      return 32;
    default:
      assert(0);
      break;
  }
}

static size_t
hash_block_size(int type) {
  switch (type) {
    case HASH_SHA256:
      return 64;
    case HASH_SHA384:
      return 128;
    case HASH_SHA512:
      return 128;
    case HASH_SHAKE256:
      return 136;
    default:
      assert(0);
      break;
  }
}

/*
 * HMAC
 */

static void
hmac_init(hmac_t *hmac, int type, const unsigned char *key, size_t len) {
  size_t block_size = hash_block_size(type);
  unsigned char pad[MAX_BLOCK_SIZE];
  size_t i;

  hmac->type = type;

  assert(len <= block_size);

  for (i = 0; i < len; i++)
    pad[i] = key[i] ^ 0x36;

  for (i = len; i < block_size; i++)
    pad[i] = 0x36;

  hash_init(&hmac->inner, type);
  hash_update(&hmac->inner, pad, block_size);

  for (i = 0; i < len; i++)
    pad[i] = key[i] ^ 0x5c;

  for (i = len; i < block_size; i++)
    pad[i] = 0x5c;

  hash_init(&hmac->outer, type);
  hash_update(&hmac->outer, pad, block_size);

  memset(pad, 0x00, block_size);
}

static void
hmac_update(hmac_t *hmac, const void *data, size_t len) {
  hash_update(&hmac->inner, data, len);
}

static void
hmac_final(hmac_t *hmac, unsigned char *out) {
  size_t hash_size = hash_output_size(hmac->type);

  hash_final(&hmac->inner, out, hash_size);
  hash_update(&hmac->outer, out, hash_size);
  hash_final(&hmac->outer, out, hash_size);
}

/*
 * DRBG
 */

static void
drbg_update(drbg_t *drbg, const unsigned char *seed, size_t seed_len);

static void
drbg_init(drbg_t *drbg, int type, const unsigned char *seed, size_t seed_len) {
  size_t hash_size = hash_output_size(type);

  assert(seed != NULL);
  assert(seed_len >= 24);

  drbg->type = type;

  memset(drbg->K, 0x00, hash_size);
  memset(drbg->V, 0x01, hash_size);

  drbg_update(drbg, seed, seed_len);
}

static void
drbg_update(drbg_t *drbg, const unsigned char *seed, size_t seed_len) {
  size_t hash_size = hash_output_size(drbg->type);

  hmac_init(&drbg->kmac, drbg->type, drbg->K, hash_size);
  hmac_update(&drbg->kmac, drbg->V, hash_size);
  hmac_update(&drbg->kmac, ZERO, 1);

  if (seed_len != 0)
    hmac_update(&drbg->kmac, seed, seed_len);

  hmac_final(&drbg->kmac, drbg->K);

  hmac_init(&drbg->kmac, drbg->type, drbg->K, hash_size);
  hmac_update(&drbg->kmac, drbg->V, hash_size);
  hmac_final(&drbg->kmac, drbg->V);

  if (seed_len != 0) {
    hmac_init(&drbg->kmac, drbg->type, drbg->K, hash_size);
    hmac_update(&drbg->kmac, drbg->V, hash_size);
    hmac_update(&drbg->kmac, ONE, 1);
    hmac_update(&drbg->kmac, seed, seed_len);
    hmac_final(&drbg->kmac, drbg->K);

    hmac_init(&drbg->kmac, drbg->type, drbg->K, hash_size);
    hmac_update(&drbg->kmac, drbg->V, hash_size);
    hmac_final(&drbg->kmac, drbg->V);
  }
}

static void
drbg_generate(drbg_t *drbg, void *out, size_t len) {
  size_t hash_size = hash_output_size(drbg->type);
  unsigned char *bytes = (unsigned char *)out;
  size_t pos = 0;
  size_t left = len;
  size_t outlen = hash_size;

  while (pos < len) {
    hmac_init(&drbg->kmac, drbg->type, drbg->K, hash_size);
    hmac_update(&drbg->kmac, drbg->V, hash_size);
    hmac_final(&drbg->kmac, drbg->V);

    if (outlen > left)
      outlen = left;

    memcpy(bytes + pos, drbg->V, outlen);

    pos += outlen;
    left -= outlen;
  }

  assert(pos == len);
  assert(left == 0);

  drbg_update(drbg, NULL, 0);
}
