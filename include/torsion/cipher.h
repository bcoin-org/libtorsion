/*!
 * cipher.h - ciphers for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef _TORSION_CIPHER_H
#define _TORSION_CIPHER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>

/*
 * Symbol Aliases
 */

#define aes_init_encrypt torsion_aes_init_encrypt
#define aes_init_decrypt torsion_aes_init_decrypt
#define aes_encrypt torsion_aes_encrypt
#define aes_decrypt torsion_aes_decrypt
#define blowfish_init torsion_blowfish_init
#define blowfish_stream2word torsion_blowfish_stream2word
#define blowfish_expand0state torsion_blowfish_expand0state
#define blowfish_expandstate torsion_blowfish_expandstate
#define blowfish_enc torsion_blowfish_enc
#define blowfish_dec torsion_blowfish_dec
#define blowfish_encrypt torsion_blowfish_encrypt
#define blowfish_decrypt torsion_blowfish_decrypt
#define camellia_init torsion_camellia_init
#define camellia_encrypt torsion_camellia_encrypt
#define camellia_decrypt torsion_camellia_decrypt
#define cipher_init torsion_cipher_init
#define cipher_update torsion_cipher_update
#define cipher_final torsion_cipher_final

/*
 * Definitions
 */

#define CIPHER_AES128 0
#define CIPHER_AES192 1
#define CIPHER_AES256 2
#define CIPHER_BLOWFISH 3
#define CIPHER_CAMELLIA128 4
#define CIPHER_CAMELLIA192 5
#define CIPHER_CAMELLIA256 6
#define CIPHER_MAX 6

#define CIPHER_MODE_ECB 0
#define CIPHER_MODE_CBC 1
#define CIPHER_MODE_CTR 2
#define CIPHER_MODE_CFB 3
#define CIPHER_MODE_OFB 4
#define CIPHER_MODE_GCM 5
#define CIPHER_MODE_MAX 5

#define CIPHER_MAX_BLOCK_SIZE 16

/*
 * Structs
 */

typedef struct _aes_s {
  unsigned int rounds;
  uint32_t key[60];
} aes_t;

typedef struct _blowfish_s {
  uint32_t S[4][256];
  uint32_t P[18];
} blowfish_t;

typedef struct _camellia_s {
  unsigned int bits;
  uint32_t key[64];
} camellia_t;

typedef struct _cipher_s {
  int type;
  int mode;
  int encrypt;
  size_t block_size;
  size_t block_pos;
  size_t last_size;
  union {
    aes_t aes;
    blowfish_t blowfish;
    camellia_t camellia;
  } ctx;
  unsigned char block[CIPHER_MAX_BLOCK_SIZE];
  unsigned char last[CIPHER_MAX_BLOCK_SIZE];
  unsigned char prev[CIPHER_MAX_BLOCK_SIZE];
  unsigned char state[CIPHER_MAX_BLOCK_SIZE];
  unsigned char ctr[CIPHER_MAX_BLOCK_SIZE];
} cipher_t;

/*
 * AES
 */

void
aes_init_encrypt(aes_t *ctx, unsigned int bits, const unsigned char *key);

void
aes_init_decrypt(aes_t *ctx, unsigned int bits, const unsigned char *key);

void
aes_encrypt(const aes_t *ctx, unsigned char *dst, const unsigned char *src);

void
aes_decrypt(const aes_t *ctx, unsigned char *dst, const unsigned char *src);

/*
 * Blowfish
 */

void
blowfish_init(blowfish_t *ctx,
              const unsigned char *key, size_t key_len,
              const unsigned char *salt, size_t salt_len);

uint32_t
blowfish_stream2word(const unsigned char *data, size_t len, size_t *off);

void
blowfish_expand0state(blowfish_t *ctx,
                      const unsigned char *key,
                      size_t key_len);

void
blowfish_expandstate(blowfish_t *ctx,
                     const unsigned char *key, size_t key_len,
                     const unsigned char *data, size_t data_len);

void
blowfish_enc(blowfish_t *ctx, uint32_t *data, size_t len);

void
blowfish_dec(blowfish_t *ctx, uint32_t *data, size_t len);

void
blowfish_encrypt(blowfish_t *ctx, unsigned char *dst, const unsigned char *src);

void
blowfish_decrypt(blowfish_t *ctx, unsigned char *dst, const unsigned char *src);

/*
 * Camellia
 */

void
camellia_init(camellia_t *ctx, unsigned int bits, const unsigned char *key);

void
camellia_encrypt(camellia_t *ctx, unsigned char *dst, const unsigned char *src);

void
camellia_decrypt(camellia_t *ctx, unsigned char *dst, const unsigned char *src);

/*
 * Cipher
 */

int
cipher_init(cipher_t *ctx, int type, int mode, int encrypt,
            const unsigned char *key, size_t key_len,
            const unsigned char *iv, size_t iv_len);

void
cipher_update(cipher_t *ctx,
              unsigned char *output, size_t *output_len,
              const unsigned char *input, size_t input_len);

int
cipher_final(cipher_t *ctx, unsigned char *output, size_t *output_len);

#ifdef __cplusplus
}
#endif

#endif /* _TORSION_CIPHER_H */
