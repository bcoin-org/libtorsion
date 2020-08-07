/*!
 * stream.c - stream ciphers for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <torsion/stream.h>
#include <torsion/util.h>
#include "bio.h"
#include "internal.h"

/*
 * ARC4
 *
 * Resources:
 *   https://en.wikipedia.org/wiki/RC4
 *   http://cypherpunks.venona.com/archive/1994/09/msg00304.html
 *   https://web.archive.org/web/20080207125928/http://cypherpunks.venona.com/archive/1994/09/msg00304.html
 *   https://tools.ietf.org/html/rfc4345
 *   https://tools.ietf.org/html/rfc6229
 *   https://github.com/golang/go/blob/master/src/crypto/rc4/rc4.go
 */

void
arc4_init(arc4_t *ctx, const unsigned char *key, size_t key_len) {
  size_t k = key_len;
  uint8_t *s = ctx->s;
  uint8_t j, si;
  size_t i;

  CHECK(k >= 1 && k <= 256);

  for (i = 0; i < 256; i++)
    s[i] = i;

  j = 0;

  for (i = 0; i < 256; i++) {
    j += s[i] + key[i % k];
    j &= 0xff;

    si = s[i];
    s[i] = s[j];
    s[j] = si;
  }

  ctx->i = 0;
  ctx->j = 0;
}

void
arc4_crypt(arc4_t *ctx,
           unsigned char *dst,
           const unsigned char *src,
           size_t len) {
  uint8_t *s = ctx->s;
  uint8_t i = ctx->i;
  uint8_t j = ctx->j;
  uint8_t x, y;
  size_t k;

  for (k = 0; k < len; k++) {
    i = (i + 1) & 0xff;
    x = s[i];

    j = (j + x) & 0xff;
    y = s[j];

    s[i] = y;
    s[j] = x;

    dst[k] = src[k] ^ s[(x + y) & 0xff];
  }

  ctx->i = i;
  ctx->j = j;
}

/*
 * ChaCha20
 *
 * Resources:
 *   https://en.wikipedia.org/wiki/Chacha20
 *   https://tools.ietf.org/html/rfc7539#section-2
 *   https://cr.yp.to/chacha.html
 */

#define ROTL32(x, y) ((x) << (y)) | ((x) >> (32 - (y)))

#define QROUND(x, a, b, c, d)                   \
  x[a] += x[b]; x[d] = ROTL32(x[d] ^ x[a], 16); \
  x[c] += x[d]; x[b] = ROTL32(x[b] ^ x[c], 12); \
  x[a] += x[b]; x[d] = ROTL32(x[d] ^ x[a], 8);  \
  x[c] += x[d]; x[b] = ROTL32(x[b] ^ x[c], 7)

void
chacha20_init(chacha20_t *ctx,
              const unsigned char *key,
              size_t key_len,
              const unsigned char *nonce,
              size_t nonce_len,
              uint64_t counter) {
  uint8_t tmp[32];

  CHECK(key_len == 16 || key_len == 32);

  if (nonce_len >= 24) {
    chacha20_derive(tmp, key, key_len, nonce);
    key = tmp;
    key_len = 32;
    nonce += 16;
    nonce_len -= 16;
  }

  ctx->state[0] = 0x61707865;
  ctx->state[1] = key_len < 32 ? 0x3120646e : 0x3320646e;
  ctx->state[2] = key_len < 32 ? 0x79622d36 : 0x79622d32;
  ctx->state[3] = 0x6b206574;
  ctx->state[4] = read32le(key + 0);
  ctx->state[5] = read32le(key + 4);
  ctx->state[6] = read32le(key + 8);
  ctx->state[7] = read32le(key + 12);
  ctx->state[8] = read32le(key + 16 % key_len);
  ctx->state[9] = read32le(key + 20 % key_len);
  ctx->state[10] = read32le(key + 24 % key_len);
  ctx->state[11] = read32le(key + 28 % key_len);
  ctx->state[12] = counter;

  if (nonce_len == 8) {
    ctx->state[13] = counter >> 32;
    ctx->state[14] = read32le(nonce + 0);
    ctx->state[15] = read32le(nonce + 4);
  } else if (nonce_len == 12) {
    ctx->state[13] = read32le(nonce + 0);
    ctx->state[14] = read32le(nonce + 4);
    ctx->state[15] = read32le(nonce + 8);
  } else if (nonce_len == 16) {
    ctx->state[12] = read32le(nonce + 0);
    ctx->state[13] = read32le(nonce + 4);
    ctx->state[14] = read32le(nonce + 8);
    ctx->state[15] = read32le(nonce + 12);
  } else {
    torsion_abort(); /* LCOV_EXCL_LINE */
  }

  ctx->pos = 0;

  torsion_cleanse(tmp, sizeof(tmp));
}

static void
chacha20_block(chacha20_t *ctx, uint32_t *stream) {
#if defined(TORSION_HAVE_ASM_X64)
  /* Borrowed from:
   * https://github.com/gnutls/nettle/blob/master/x86_64/chacha-core-internal.asm
   *
   * Registers:
   *
   *   %rdi = dst pointer (stream)
   *   %rsi = src pointer (ctx->state)
   *   %edx = rounds integer (nettle does `20 >> 1`)
   *
   * For reference, our full range of clobbered registers:
   *
   *   %rsi, %rdi, %edx, %xmm[0-5]
   */
  __asm__ __volatile__(
    "movl $10, %%edx\n"

    "movups (%%rsi), %%xmm0\n"
    "movups 16(%%rsi), %%xmm1\n"
    "movups 32(%%rsi), %%xmm2\n"
    "movups 48(%%rsi), %%xmm3\n"

    ".align 16\n"

    "1:\n"

    "paddd %%xmm1, %%xmm0\n"
    "pxor %%xmm0, %%xmm3\n"
    "movaps %%xmm3, %%xmm4\n"

    "pshufhw $0xb1, %%xmm3, %%xmm3\n"
    "pshuflw $0xb1, %%xmm3, %%xmm3\n"

    "paddd %%xmm3, %%xmm2\n"
    "pxor %%xmm2, %%xmm1\n"
    "movaps %%xmm1, %%xmm4\n"
    "pslld $12, %%xmm1\n"
    "psrld $20, %%xmm4\n"
    "por %%xmm4, %%xmm1\n"

    "paddd %%xmm1, %%xmm0\n"
    "pxor %%xmm0, %%xmm3\n"
    "movaps %%xmm3, %%xmm4\n"
    "pslld $8, %%xmm3\n"
    "psrld $24, %%xmm4\n"
    "por %%xmm4, %%xmm3\n"

    "paddd %%xmm3, %%xmm2\n"
    "pxor %%xmm2, %%xmm1\n"
    "movaps %%xmm1, %%xmm4\n"
    "pslld $7, %%xmm1\n"
    "psrld $25, %%xmm4\n"
    "por %%xmm4, %%xmm1\n"

    "pshufd $0x39, %%xmm1, %%xmm1\n"
    "pshufd $0x4e, %%xmm2, %%xmm2\n"
    "pshufd $0x93, %%xmm3, %%xmm3\n"

    "paddd %%xmm1, %%xmm0\n"
    "pxor %%xmm0, %%xmm3\n"
    "movaps %%xmm3, %%xmm4\n"

    "pshufhw $0xb1, %%xmm3, %%xmm3\n"
    "pshuflw $0xb1, %%xmm3, %%xmm3\n"

    "paddd %%xmm3, %%xmm2\n"
    "pxor %%xmm2, %%xmm1\n"
    "movaps %%xmm1, %%xmm4\n"
    "pslld $12, %%xmm1\n"
    "psrld $20, %%xmm4\n"
    "por %%xmm4, %%xmm1\n"

    "paddd %%xmm1, %%xmm0\n"
    "pxor %%xmm0, %%xmm3\n"
    "movaps %%xmm3, %%xmm4\n"
    "pslld $8, %%xmm3\n"
    "psrld $24, %%xmm4\n"
    "por %%xmm4, %%xmm3\n"

    "paddd %%xmm3, %%xmm2\n"
    "pxor %%xmm2, %%xmm1\n"
    "movaps %%xmm1, %%xmm4\n"
    "pslld $7, %%xmm1\n"
    "psrld $25, %%xmm4\n"
    "por %%xmm4, %%xmm1\n"

    "pshufd $0x93, %%xmm1, %%xmm1\n"
    "pshufd $0x4e, %%xmm2, %%xmm2\n"
    "pshufd $0x39, %%xmm3, %%xmm3\n"

    "decl %%edx\n"
    "jnz 1b\n"

    "movups (%%rsi), %%xmm4\n"
    "movups 16(%%rsi), %%xmm5\n"
    "paddd %%xmm4, %%xmm0\n"
    "paddd %%xmm5, %%xmm1\n"
    "movups %%xmm0,(%%rdi)\n"
    "movups %%xmm1,16(%%rdi)\n"
    "movups 32(%%rsi), %%xmm4\n"
    "movups 48(%%rsi), %%xmm5\n"
    "paddd %%xmm4, %%xmm2\n"
    "paddd %%xmm5, %%xmm3\n"
    "movups %%xmm2,32(%%rdi)\n"
    "movups %%xmm3,48(%%rdi)\n"

    "incq 48(%%rsi)\n"
    :
    : "D" (stream), "S" (ctx->state)
    : "edx", "xmm0", "xmm1", "xmm2",
      "xmm3", "xmm4", "xmm5", "cc",
      "memory"
  );
#else
  uint64_t c;
  size_t i;

  for (i = 0; i < 16; i++)
    stream[i] = ctx->state[i];

  for (i = 0; i < 10; i++) {
    QROUND(stream, 0, 4,  8, 12);
    QROUND(stream, 1, 5,  9, 13);
    QROUND(stream, 2, 6, 10, 14);
    QROUND(stream, 3, 7, 11, 15);
    QROUND(stream, 0, 5, 10, 15);
    QROUND(stream, 1, 6, 11, 12);
    QROUND(stream, 2, 7,  8, 13);
    QROUND(stream, 3, 4,  9, 14);
  }

  for (i = 0; i < 16; i++)
    stream[i] += ctx->state[i];

  if (TORSION_BIGENDIAN) {
    for (i = 0; i < 16; i++)
      stream[i] = torsion_bswap32(stream[i]);
  }

  c = (uint64_t)ctx->state[12] + 1;

  ctx->state[12] = c;
  ctx->state[13] += (uint32_t)(c >> 32);
#endif
}

void
chacha20_crypt(chacha20_t *ctx,
               unsigned char *out,
               const unsigned char *data,
               size_t len) {
  unsigned char *bytes = (unsigned char *)ctx->stream;
  size_t i;

  for (i = 0; i < len; i++) {
    if ((ctx->pos & 63) == 0) {
      chacha20_block(ctx, ctx->stream);
      ctx->pos = 0;
    }

    out[i] = data[i] ^ bytes[ctx->pos++];
  }
}

void
chacha20_derive(unsigned char *out,
                const unsigned char *key,
                size_t key_len,
                const unsigned char *nonce16) {
  uint32_t state[16];
  size_t i;

  CHECK(key_len == 16 || key_len == 32);

  state[0] = 0x61707865;
  state[1] = key_len < 32 ? 0x3120646e : 0x3320646e;
  state[2] = key_len < 32 ? 0x79622d36 : 0x79622d32;
  state[3] = 0x6b206574;
  state[4] = read32le(key + 0);
  state[5] = read32le(key + 4);
  state[6] = read32le(key + 8);
  state[7] = read32le(key + 12);
  state[8] = read32le(key + 16 % key_len);
  state[9] = read32le(key + 20 % key_len);
  state[10] = read32le(key + 24 % key_len);
  state[11] = read32le(key + 28 % key_len);
  state[12] = read32le(nonce16 + 0);
  state[13] = read32le(nonce16 + 4);
  state[14] = read32le(nonce16 + 8);
  state[15] = read32le(nonce16 + 12);

  for (i = 0; i < 10; i++) {
    QROUND(state, 0, 4,  8, 12);
    QROUND(state, 1, 5,  9, 13);
    QROUND(state, 2, 6, 10, 14);
    QROUND(state, 3, 7, 11, 15);
    QROUND(state, 0, 5, 10, 15);
    QROUND(state, 1, 6, 11, 12);
    QROUND(state, 2, 7,  8, 13);
    QROUND(state, 3, 4,  9, 14);
  }

  write32le(out +  0, state[0]);
  write32le(out +  4, state[1]);
  write32le(out +  8, state[2]);
  write32le(out + 12, state[3]);
  write32le(out + 16, state[12]);
  write32le(out + 20, state[13]);
  write32le(out + 24, state[14]);
  write32le(out + 28, state[15]);

  torsion_cleanse(state, sizeof(state));
}

#undef ROTL32
#undef QROUND

/*
 * Salsa20
 *
 * Resources
 *   https://en.wikipedia.org/wiki/Salsa20
 *   https://cr.yp.to/snuffle.html
 *   https://cr.yp.to/snuffle/spec.pdf
 *   https://cr.yp.to/snuffle/812.pdf
 *   http://www.ecrypt.eu.org/stream/salsa20pf.html
 */

#define ROTL32(v, n) ((v) << (n)) | ((v) >> (32 - (n)))

#define QROUND(x, a, b, c, d)      \
  x[b] ^= ROTL32(x[a] + x[d], 7);  \
  x[c] ^= ROTL32(x[b] + x[a], 9);  \
  x[d] ^= ROTL32(x[c] + x[b], 13); \
  x[a] ^= ROTL32(x[d] + x[c], 18)

void
salsa20_init(salsa20_t *ctx,
             const unsigned char *key,
             size_t key_len,
             const unsigned char *nonce,
             size_t nonce_len,
             uint64_t counter) {
  uint8_t tmp[32];

  CHECK(key_len == 16 || key_len == 32);

  if (nonce_len >= 24) {
    salsa20_derive(tmp, key, key_len, nonce);
    key = tmp;
    key_len = 32;
    nonce += 16;
    nonce_len -= 16;
  }

  ctx->state[0] = 0x61707865;
  ctx->state[1] = read32le(key + 0);
  ctx->state[2] = read32le(key + 4);
  ctx->state[3] = read32le(key + 8);
  ctx->state[4] = read32le(key + 12);
  ctx->state[5] = key_len < 32 ? 0x3120646e : 0x3320646e;

  if (nonce_len == 8) {
    ctx->state[6] = read32le(nonce + 0);
    ctx->state[7] = read32le(nonce + 4);
    ctx->state[8] = counter;
    ctx->state[9] = counter >> 32;
  } else if (nonce_len == 12) {
    ctx->state[6] = read32le(nonce + 0);
    ctx->state[7] = read32le(nonce + 4);
    ctx->state[8] = read32le(nonce + 8);
    ctx->state[9] = counter;
  } else if (nonce_len == 16) {
    ctx->state[6] = read32le(nonce + 0);
    ctx->state[7] = read32le(nonce + 4);
    ctx->state[8] = read32le(nonce + 8);
    ctx->state[9] = read32le(nonce + 12);
  } else {
    torsion_abort(); /* LCOV_EXCL_LINE */
  }

  ctx->state[10] = key_len < 32 ? 0x79622d36 : 0x79622d32;
  ctx->state[11] = read32le(key + 16 % key_len);
  ctx->state[12] = read32le(key + 20 % key_len);
  ctx->state[13] = read32le(key + 24 % key_len);
  ctx->state[14] = read32le(key + 28 % key_len);
  ctx->state[15] = 0x6b206574;

  ctx->pos = 0;

  torsion_cleanse(tmp, sizeof(tmp));
}

static void
salsa20_block(salsa20_t *ctx, uint32_t *stream) {
#if defined(TORSION_HAVE_ASM_X64)
  /* Borrowed from:
   * https://github.com/gnutls/nettle/blob/master/x86_64/salsa20-core-internal.asm
   *
   * Layout:
   *
   *   %rdi = dst pointer (stream)
   *   %rsi = src pointer (ctx->state)
   *   %edx = rounds integer (nettle does `20 >> 1`)
   *
   * For reference, our full range of clobbered registers:
   *
   *   %rsi, %rdi, %edx, %xmm[0-8]
   */
  __asm__ __volatile__(
    "movl $-1, %%edx\n"
    "movd %%edx, %%xmm6\n"

    "movl $10, %%edx\n"

    "pshufd $0x09, %%xmm6, %%xmm8\n"
    "pshufd $0x41, %%xmm6, %%xmm7\n"
    "pshufd $0x22, %%xmm6, %%xmm6\n"

    "movups (%%rsi), %%xmm0\n"
    "movups 16(%%rsi), %%xmm1\n"
    "movups 32(%%rsi), %%xmm2\n"
    "movups 48(%%rsi), %%xmm3\n"

    "movaps %%xmm0, %%xmm4\n"
    "pxor %%xmm1, %%xmm0\n"
    "pand %%xmm6, %%xmm0\n"
    "pxor %%xmm0, %%xmm1\n"
    "pxor %%xmm4, %%xmm0\n"

    "movaps %%xmm2, %%xmm4\n"
    "pxor %%xmm3, %%xmm2\n"
    "pand %%xmm6, %%xmm2\n"
    "pxor %%xmm2, %%xmm3\n"
    "pxor %%xmm4, %%xmm2\n"

    "movaps %%xmm1, %%xmm4\n"
    "pxor %%xmm3, %%xmm1\n"
    "pand %%xmm7, %%xmm1\n"
    "pxor %%xmm1, %%xmm3\n"
    "pxor %%xmm4, %%xmm1\n"

    "movaps %%xmm0, %%xmm4\n"
    "pxor %%xmm2, %%xmm0\n"
    "pand %%xmm8, %%xmm0\n"
    "pxor %%xmm0, %%xmm2\n"
    "pxor %%xmm4, %%xmm0\n"

    ".align 16\n"

    "1:\n"

    "movaps %%xmm3, %%xmm4\n"
    "paddd %%xmm0, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $7, %%xmm4\n"
    "psrld $25, %%xmm5\n"
    "pxor %%xmm4, %%xmm1\n"
    "pxor %%xmm5, %%xmm1\n"

    "movaps %%xmm0, %%xmm4\n"
    "paddd %%xmm1, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $9, %%xmm4\n"
    "psrld $23, %%xmm5\n"
    "pxor %%xmm4, %%xmm2\n"
    "pxor %%xmm5, %%xmm2\n"

    "movaps %%xmm1, %%xmm4\n"
    "paddd %%xmm2, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $13, %%xmm4\n"
    "psrld $19, %%xmm5\n"
    "pxor %%xmm4, %%xmm3\n"
    "pxor %%xmm5, %%xmm3\n"

    "movaps %%xmm2, %%xmm4\n"
    "paddd %%xmm3, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $18, %%xmm4\n"
    "psrld $14, %%xmm5\n"
    "pxor %%xmm4, %%xmm0\n"
    "pxor %%xmm5, %%xmm0\n"

    "pshufd $0x93, %%xmm1, %%xmm1\n"
    "pshufd $0x4e, %%xmm2, %%xmm2\n"
    "pshufd $0x39, %%xmm3, %%xmm3\n"

    "movaps %%xmm1, %%xmm4\n"
    "paddd %%xmm0, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $7, %%xmm4\n"
    "psrld $25, %%xmm5\n"
    "pxor %%xmm4, %%xmm3\n"
    "pxor %%xmm5, %%xmm3\n"

    "movaps %%xmm0, %%xmm4\n"
    "paddd %%xmm3, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $9, %%xmm4\n"
    "psrld $23, %%xmm5\n"
    "pxor %%xmm4, %%xmm2\n"
    "pxor %%xmm5, %%xmm2\n"

    "movaps %%xmm3, %%xmm4\n"
    "paddd %%xmm2, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $13, %%xmm4\n"
    "psrld $19, %%xmm5\n"
    "pxor %%xmm4, %%xmm1\n"
    "pxor %%xmm5, %%xmm1\n"

    "movaps %%xmm2, %%xmm4\n"
    "paddd %%xmm1, %%xmm4\n"
    "movaps %%xmm4, %%xmm5\n"
    "pslld $18, %%xmm4\n"
    "psrld $14, %%xmm5\n"
    "pxor %%xmm4, %%xmm0\n"
    "pxor %%xmm5, %%xmm0\n"

    "pshufd $0x39, %%xmm1, %%xmm1\n"
    "pshufd $0x4e, %%xmm2, %%xmm2\n"
    "pshufd $0x93, %%xmm3, %%xmm3\n"

    "decl %%edx\n"
    "jnz 1b\n"

    "movaps %%xmm0, %%xmm4\n"
    "pxor %%xmm2, %%xmm0\n"
    "pand %%xmm8, %%xmm0\n"
    "pxor %%xmm0, %%xmm2\n"
    "pxor %%xmm4, %%xmm0\n"

    "movaps %%xmm1, %%xmm4\n"
    "pxor %%xmm3, %%xmm1\n"
    "pand %%xmm7, %%xmm1\n"
    "pxor %%xmm1, %%xmm3\n"
    "pxor %%xmm4, %%xmm1\n"

    "movaps %%xmm0, %%xmm4\n"
    "pxor %%xmm1, %%xmm0\n"
    "pand %%xmm6, %%xmm0\n"
    "pxor %%xmm0, %%xmm1\n"
    "pxor %%xmm4, %%xmm0\n"

    "movaps %%xmm2, %%xmm4\n"
    "pxor %%xmm3, %%xmm2\n"
    "pand %%xmm6, %%xmm2\n"
    "pxor %%xmm2, %%xmm3\n"
    "pxor %%xmm4, %%xmm2\n"

    "movups (%%rsi), %%xmm4\n"
    "movups 16(%%rsi), %%xmm5\n"
    "paddd %%xmm4, %%xmm0\n"
    "paddd %%xmm5, %%xmm1\n"
    "movups %%xmm0,(%%rdi)\n"
    "movups %%xmm1,16(%%rdi)\n"
    "movups 32(%%rsi), %%xmm4\n"
    "movups 48(%%rsi), %%xmm5\n"
    "paddd %%xmm4, %%xmm2\n"
    "paddd %%xmm5, %%xmm3\n"
    "movups %%xmm2,32(%%rdi)\n"
    "movups %%xmm3,48(%%rdi)\n"

    "incq 32(%%rsi)\n"
    :
    : "D" (stream), "S" (ctx->state)
    : "edx", "xmm0", "xmm1", "xmm2",
      "xmm3", "xmm4", "xmm5", "xmm6",
      "xmm7", "xmm8", "cc", "memory"
  );
#else
  uint64_t c;
  size_t i;

  for (i = 0; i < 16; i++)
    stream[i] = ctx->state[i];

  for (i = 0; i < 10; i++) {
    QROUND(stream,  0,  4,  8, 12);
    QROUND(stream,  5,  9, 13,  1);
    QROUND(stream, 10, 14,  2,  6);
    QROUND(stream, 15,  3,  7, 11);
    QROUND(stream,  0,  1,  2,  3);
    QROUND(stream,  5,  6,  7,  4);
    QROUND(stream, 10, 11,  8,  9);
    QROUND(stream, 15, 12, 13, 14);
  }

  for (i = 0; i < 16; i++)
    stream[i] += ctx->state[i];

  if (TORSION_BIGENDIAN) {
    for (i = 0; i < 16; i++)
      stream[i] = torsion_bswap32(stream[i]);
  }

  c = (uint64_t)ctx->state[8] + 1;

  ctx->state[8] = c;
  ctx->state[9] += (uint32_t)(c >> 32);
#endif
}

void
salsa20_crypt(salsa20_t *ctx,
              unsigned char *out,
              const unsigned char *data,
              size_t len) {
  unsigned char *bytes = (unsigned char *)ctx->stream;
  size_t i;

  for (i = 0; i < len; i++) {
    if ((ctx->pos & 63) == 0) {
      salsa20_block(ctx, ctx->stream);
      ctx->pos = 0;
    }

    out[i] = data[i] ^ bytes[ctx->pos++];
  }
}

void
salsa20_derive(unsigned char *out,
               const unsigned char *key,
               size_t key_len,
               const unsigned char *nonce16) {
  uint32_t state[16];
  size_t i;

  CHECK(key_len == 16 || key_len == 32);

  state[0] = 0x61707865;
  state[1] = read32le(key + 0);
  state[2] = read32le(key + 4);
  state[3] = read32le(key + 8);
  state[4] = read32le(key + 12);
  state[5] = key_len < 32 ? 0x3120646e : 0x3320646e;
  state[6] = read32le(nonce16 + 0);
  state[7] = read32le(nonce16 + 4);
  state[8] = read32le(nonce16 + 8);
  state[9] = read32le(nonce16 + 12);
  state[10] = key_len < 32 ? 0x79622d36 : 0x79622d32;
  state[11] = read32le(key + 16 % key_len);
  state[12] = read32le(key + 20 % key_len);
  state[13] = read32le(key + 24 % key_len);
  state[14] = read32le(key + 28 % key_len);
  state[15] = 0x6b206574;

  for (i = 0; i < 10; i++) {
    QROUND(state,  0,  4,  8, 12);
    QROUND(state,  5,  9, 13,  1);
    QROUND(state, 10, 14,  2,  6);
    QROUND(state, 15,  3,  7, 11);
    QROUND(state,  0,  1,  2,  3);
    QROUND(state,  5,  6,  7,  4);
    QROUND(state, 10, 11,  8,  9);
    QROUND(state, 15, 12, 13, 14);
  }

  write32le(out +  0, state[0]);
  write32le(out +  4, state[5]);
  write32le(out +  8, state[10]);
  write32le(out + 12, state[15]);
  write32le(out + 16, state[6]);
  write32le(out + 20, state[7]);
  write32le(out + 24, state[8]);
  write32le(out + 28, state[9]);

  torsion_cleanse(state, sizeof(state));
}

#undef ROTL32
#undef QROUND
