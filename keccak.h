#include <stdlib.h>
#include <stdint.h>

#define ROTL64(n,x) (((x)<<(n)) | ((x)>>((-(n))&63)))

static uint64_t
read64(const void *src);

typedef struct keccak_s {
  size_t bs;
  uint64_t state[25];
  uint8_t block[168];
  size_t pos;
} keccak_t;

static void
keccak_init(keccak_t *ctx, size_t bits) {
  size_t rate = 1600 - bits * 2;

  assert(bits >= 128);
  assert(bits <= 512);
  assert((rate & 63) == 0);

  ctx->bs = rate >> 3;
  ctx->pos = 0;

  memset(ctx->state, 0, sizeof(ctx->state));
}

static void
keccak_permute(keccak_t *ctx) {
  /* borrowed from nettle */
  static const uint64_t rc[24] = {
    0x0000000000000001ULL, 0X0000000000008082ULL,
    0X800000000000808AULL, 0X8000000080008000ULL,
    0X000000000000808BULL, 0X0000000080000001ULL,
    0X8000000080008081ULL, 0X8000000000008009ULL,
    0X000000000000008AULL, 0X0000000000000088ULL,
    0X0000000080008009ULL, 0X000000008000000AULL,
    0X000000008000808BULL, 0X800000000000008BULL,
    0X8000000000008089ULL, 0X8000000000008003ULL,
    0X8000000000008002ULL, 0X8000000000000080ULL,
    0X000000000000800AULL, 0X800000008000000AULL,
    0X8000000080008081ULL, 0X8000000000008080ULL,
    0X0000000080000001ULL, 0X8000000080008008ULL,
  };

  uint64_t C[5], D[5], T, X;
  unsigned i, y;

#define A ctx->state

  C[0] = A[0] ^ A[5+0] ^ A[10+0] ^ A[15+0] ^ A[20+0];
  C[1] = A[1] ^ A[5+1] ^ A[10+1] ^ A[15+1] ^ A[20+1];
  C[2] = A[2] ^ A[5+2] ^ A[10+2] ^ A[15+2] ^ A[20+2];
  C[3] = A[3] ^ A[5+3] ^ A[10+3] ^ A[15+3] ^ A[20+3];
  C[4] = A[4] ^ A[5+4] ^ A[10+4] ^ A[15+4] ^ A[20+4];

  for (i = 0; i < 24; i++) {
    D[0] = C[4] ^ ROTL64(1, C[1]);
    D[1] = C[0] ^ ROTL64(1, C[2]);
    D[2] = C[1] ^ ROTL64(1, C[3]);
    D[3] = C[2] ^ ROTL64(1, C[4]);
    D[4] = C[3] ^ ROTL64(1, C[0]);

    A[0] ^= D[0];
    X = A[ 1] ^ D[1];     T = ROTL64(1, X);
    X = A[ 6] ^ D[1]; A[ 1] = ROTL64 (44, X);
    X = A[ 9] ^ D[4]; A[ 6] = ROTL64 (20, X);
    X = A[22] ^ D[2]; A[ 9] = ROTL64 (61, X);
    X = A[14] ^ D[4]; A[22] = ROTL64 (39, X);
    X = A[20] ^ D[0]; A[14] = ROTL64 (18, X);
    X = A[ 2] ^ D[2]; A[20] = ROTL64 (62, X);
    X = A[12] ^ D[2]; A[ 2] = ROTL64 (43, X);
    X = A[13] ^ D[3]; A[12] = ROTL64 (25, X);
    X = A[19] ^ D[4]; A[13] = ROTL64 ( 8, X);
    X = A[23] ^ D[3]; A[19] = ROTL64 (56, X);
    X = A[15] ^ D[0]; A[23] = ROTL64 (41, X);
    X = A[ 4] ^ D[4]; A[15] = ROTL64 (27, X);
    X = A[24] ^ D[4]; A[ 4] = ROTL64 (14, X);
    X = A[21] ^ D[1]; A[24] = ROTL64 ( 2, X);
    X = A[ 8] ^ D[3]; A[21] = ROTL64 (55, X);
    X = A[16] ^ D[1]; A[ 8] = ROTL64 (45, X);
    X = A[ 5] ^ D[0]; A[16] = ROTL64 (36, X);
    X = A[ 3] ^ D[3]; A[ 5] = ROTL64 (28, X);
    X = A[18] ^ D[3]; A[ 3] = ROTL64 (21, X);
    X = A[17] ^ D[2]; A[18] = ROTL64 (15, X);
    X = A[11] ^ D[1]; A[17] = ROTL64 (10, X);
    X = A[ 7] ^ D[2]; A[11] = ROTL64 ( 6, X);
    X = A[10] ^ D[0]; A[ 7] = ROTL64 ( 3, X);
    A[10] = T;

    D[0] = ~A[1] & A[2];
    D[1] = ~A[2] & A[3];
    D[2] = ~A[3] & A[4];
    D[3] = ~A[4] & A[0];
    D[4] = ~A[0] & A[1];

    A[0] ^= D[0] ^ rc[i]; C[0] = A[0];
    A[1] ^= D[1]; C[1] = A[1];
    A[2] ^= D[2]; C[2] = A[2];
    A[3] ^= D[3]; C[3] = A[3];
    A[4] ^= D[4]; C[4] = A[4];

    for (y = 5; y < 25; y+= 5) {
      D[0] = ~A[y+1] & A[y+2];
      D[1] = ~A[y+2] & A[y+3];
      D[2] = ~A[y+3] & A[y+4];
      D[3] = ~A[y+4] & A[y+0];
      D[4] = ~A[y+0] & A[y+1];

      A[y+0] ^= D[0]; C[0] ^= A[y+0];
      A[y+1] ^= D[1]; C[1] ^= A[y+1];
      A[y+2] ^= D[2]; C[2] ^= A[y+2];
      A[y+3] ^= D[3]; C[3] ^= A[y+3];
      A[y+4] ^= D[4]; C[4] ^= A[y+4];
    }
  }
#undef A
}

static void
keccak_transform(keccak_t *ctx, const unsigned char *chunk) {
  size_t count = ctx->bs >> 3;
  size_t i;

  for (i = 0; i < count; i++)
    ctx->state[i] ^= read64(chunk + i * 8);

  keccak_permute(ctx);
}

static void
keccak_update(keccak_t *ctx, const void *data, size_t len) {
  const unsigned char *bytes = (const unsigned char *)data;
  size_t pos = ctx->pos;
  size_t off = 0;

  ctx->pos = (ctx->pos + len) % ctx->bs;

  if (pos > 0) {
    size_t want = ctx->bs - pos;

    if (want > len)
      want = len;

    memcpy(ctx->block + pos, bytes + off, want);

    pos += want;
    len -= want;
    off += want;

    if (pos < ctx->bs)
      return;

    keccak_transform(ctx, ctx->block);
  }

  while (len >= ctx->bs) {
    keccak_transform(ctx, bytes + off);
    off += ctx->bs;
    len -= ctx->bs;
  }

  if (len > 0)
    memcpy(ctx->block, bytes + off, len);
}

static void
keccak_final(keccak_t *ctx, unsigned char *out, int pad, size_t len) {
  size_t i;

  if (pad == 0)
    pad = 0x01;

  if (len == 0)
    len = 100 - (ctx->bs >> 1);

  memset(ctx->block + ctx->pos, 0x00, ctx->bs - ctx->pos);
  ctx->block[ctx->pos] |= pad;
  ctx->block[ctx->bs - 1] |= 0x80;
  keccak_transform(ctx, ctx->block);

  for (i = 0; i < len; i++)
    out[i] = ctx->state[i >> 3] >> (8 * (i & 7));

  memset(ctx->state, 0x00, sizeof(ctx->state));
  memset(ctx->block, 0x00, sizeof(ctx->block));

  ctx->bs = 0;
  ctx->pos = 0;
}

static uint64_t
read64(const void *src) {
#ifndef WORDS_BIGENDIAN
  uint64_t w;
  memcpy(&w, src, sizeof(w));
  return w;
#else
  const uint8_t *p = (const uint8_t *)src;
  return ((uint64_t)p[7] << 56)
       | ((uint64_t)p[6] << 48)
       | ((uint64_t)p[5] << 40)
       | ((uint64_t)p[4] << 32)
       | ((uint64_t)p[3] << 24)
       | ((uint64_t)p[2] << 16)
       | ((uint64_t)p[1] << 8)
       | ((uint64_t)p[0] << 0);
#endif
}
