/* permutation/compression functions from nettle */

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <torsion/hash.h>

#define ROTL32(n, x) (((x) << (n)) | ((x) >> ((-(n) & 31))))
#define ROTL64(n, x) (((x) << (n)) | ((x) >> ((-(n)) & 63)))

/*
 * Helpers
 */

static uint32_t
read32be(const void *src) {
#ifdef WORDS_BIGENDIAN
  uint32_t w;
  memcpy(&w, src, sizeof(w));
  return w;
#else
  const uint8_t *p = (const uint8_t *)src;
  return ((uint32_t)p[0] << 24)
       | ((uint32_t)p[1] << 16)
       | ((uint32_t)p[2] << 8)
       | ((uint32_t)p[3] << 0);
#endif
}

static void
write32be(void *dst, uint32_t w) {
#ifdef WORDS_BIGENDIAN
  memcpy(dst, &w, sizeof(w));
#else
  uint8_t *p = (uint8_t *)dst;
  p[0] = w >> 24;
  p[1] = w >> 16;
  p[2] = w >> 8;
  p[3] = w >> 0;
#endif
}

static uint64_t
read64be(const void *src) {
#ifdef WORDS_BIGENDIAN
  uint64_t w;
  memcpy(&w, src, sizeof(w));
  return w;
#else
  const uint8_t *p = (const uint8_t *)src;
  return ((uint64_t)p[0] << 56)
       | ((uint64_t)p[1] << 48)
       | ((uint64_t)p[2] << 40)
       | ((uint64_t)p[3] << 32)
       | ((uint64_t)p[4] << 24)
       | ((uint64_t)p[5] << 16)
       | ((uint64_t)p[6] << 8)
       | ((uint64_t)p[7] << 0);
#endif
}

static void
write64be(void *dst, uint64_t w) {
#ifdef WORDS_BIGENDIAN
  memcpy(dst, &w, sizeof(w));
#else
  uint8_t *p = (uint8_t *)dst;
  p[0] = w >> 56;
  p[1] = w >> 48;
  p[2] = w >> 40;
  p[3] = w >> 32;
  p[4] = w >> 24;
  p[5] = w >> 16;
  p[6] = w >> 8;
  p[7] = w >> 0;
#endif
}

static uint64_t
read64le(const void *src) {
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

/*
 * SHA256
 */

static const uint32_t sha256_K[64] = {
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
  0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
  0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
  0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
  0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
  0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
  0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
  0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
  0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

static const unsigned char sha256_P[64] = {
  0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

void
sha256_init(sha256_t *ctx) {
  ctx->state[0] = 0x6a09e667;
  ctx->state[1] = 0xbb67ae85;
  ctx->state[2] = 0x3c6ef372;
  ctx->state[3] = 0xa54ff53a;
  ctx->state[4] = 0x510e527f;
  ctx->state[5] = 0x9b05688c;
  ctx->state[6] = 0x1f83d9ab;
  ctx->state[7] = 0x5be0cd19;
  ctx->size = 0;
}

static void
sha256_transform(sha256_t *ctx, const unsigned char *chunk) {
  const uint32_t *k = sha256_K;
  uint32_t data[16];
  uint32_t A, B, C, D, E, F, G, H;
  unsigned i;
  uint32_t *d;

  for (i = 0; i < 16; i++, chunk += 4)
    data[i] = read32be(chunk);

  A = ctx->state[0];
  B = ctx->state[1];
  C = ctx->state[2];
  D = ctx->state[3];
  E = ctx->state[4];
  F = ctx->state[5];
  G = ctx->state[6];
  H = ctx->state[7];

#define Ch(x, y, z) ((z) ^ ((x) & ((y) ^ (z))))
#define Maj(x, y, z) (((x) & (y)) ^ ((z) & ((x) ^ (y))))

#define S0(x) (ROTL32(30, (x)) ^ ROTL32(19, (x)) ^ ROTL32(10, (x)))
#define S1(x) (ROTL32(26, (x)) ^ ROTL32(21, (x)) ^ ROTL32(7, (x)))

#define s0(x) (ROTL32(25, (x)) ^ ROTL32(14, (x)) ^ ((x) >> 3))
#define s1(x) (ROTL32(15, (x)) ^ ROTL32(13, (x)) ^ ((x) >> 10))

#define EXPAND(W,i) \
(W[(i) & 15] += (s1(W[((i)-2) & 15]) + W[((i)-7) & 15] + s0(W[((i)-15) & 15])))

#define ROUND(a, b, c, d, e, f, g, h, k, data) do { \
  h += S1(e) + Ch(e,f,g) + k + data;                \
  d += h;                                           \
  h += S0(a) + Maj(a,b,c);                          \
} while (0)

  for (i = 0, d = data; i < 16; i += 8, k += 8, d += 8) {
    ROUND(A, B, C, D, E, F, G, H, k[0], d[0]);
    ROUND(H, A, B, C, D, E, F, G, k[1], d[1]);
    ROUND(G, H, A, B, C, D, E, F, k[2], d[2]);
    ROUND(F, G, H, A, B, C, D, E, k[3], d[3]);
    ROUND(E, F, G, H, A, B, C, D, k[4], d[4]);
    ROUND(D, E, F, G, H, A, B, C, k[5], d[5]);
    ROUND(C, D, E, F, G, H, A, B, k[6], d[6]);
    ROUND(B, C, D, E, F, G, H, A, k[7], d[7]);
  }

  for (; i < 64; i += 16, k += 16) {
    ROUND(A, B, C, D, E, F, G, H, k[ 0], EXPAND(data,  0));
    ROUND(H, A, B, C, D, E, F, G, k[ 1], EXPAND(data,  1));
    ROUND(G, H, A, B, C, D, E, F, k[ 2], EXPAND(data,  2));
    ROUND(F, G, H, A, B, C, D, E, k[ 3], EXPAND(data,  3));
    ROUND(E, F, G, H, A, B, C, D, k[ 4], EXPAND(data,  4));
    ROUND(D, E, F, G, H, A, B, C, k[ 5], EXPAND(data,  5));
    ROUND(C, D, E, F, G, H, A, B, k[ 6], EXPAND(data,  6));
    ROUND(B, C, D, E, F, G, H, A, k[ 7], EXPAND(data,  7));
    ROUND(A, B, C, D, E, F, G, H, k[ 8], EXPAND(data,  8));
    ROUND(H, A, B, C, D, E, F, G, k[ 9], EXPAND(data,  9));
    ROUND(G, H, A, B, C, D, E, F, k[10], EXPAND(data, 10));
    ROUND(F, G, H, A, B, C, D, E, k[11], EXPAND(data, 11));
    ROUND(E, F, G, H, A, B, C, D, k[12], EXPAND(data, 12));
    ROUND(D, E, F, G, H, A, B, C, k[13], EXPAND(data, 13));
    ROUND(C, D, E, F, G, H, A, B, k[14], EXPAND(data, 14));
    ROUND(B, C, D, E, F, G, H, A, k[15], EXPAND(data, 15));
  }

#undef Ch
#undef Maj
#undef S0
#undef S1
#undef s0
#undef s1
#undef EXPAND
#undef ROUND

  ctx->state[0] += A;
  ctx->state[1] += B;
  ctx->state[2] += C;
  ctx->state[3] += D;
  ctx->state[4] += E;
  ctx->state[5] += F;
  ctx->state[6] += G;
  ctx->state[7] += H;
}

void
sha256_update(sha256_t *ctx, const void *data, size_t len) {
  const unsigned char *bytes = (const unsigned char *)data;
  size_t pos = ctx->size & 63;
  size_t off = 0;

  ctx->size += len;

  if (pos > 0) {
    size_t want = 64 - pos;

    if (want > len)
      want = len;

    memcpy(ctx->block + pos, bytes + off, want);

    pos += want;
    len -= want;
    off += want;

    if (pos < 64)
      return;

    sha256_transform(ctx, ctx->block);
  }

  while (len >= 64) {
    sha256_transform(ctx, bytes + off);
    off += 64;
    len -= 64;
  }

  if (len > 0)
    memcpy(ctx->block, bytes + off, len);
}

void
sha256_final(sha256_t *ctx, unsigned char *out) {
  size_t pos = ctx->size & 63;
  uint64_t len = ctx->size << 3;
  unsigned char D[8];
  size_t i;

  write64be(D, len);

  sha256_update(ctx, sha256_P, 1 + ((119 - pos) & 63));
  sha256_update(ctx, D, 8);

  for (i = 0; i < 8; i++)
    write32be(out + i * 4, ctx->state[i]);

  memset(ctx->state, 0x00, sizeof(ctx->state));
  memset(ctx->block, 0x00, sizeof(ctx->block));

  ctx->size = 0;
}

/*
 * SHA512
 */

static const uint64_t sha512_K[80] = {
  0x428a2f98d728ae22ull, 0x7137449123ef65cdull,
  0xb5c0fbcfec4d3b2full, 0xe9b5dba58189dbbcull,
  0x3956c25bf348b538ull, 0x59f111f1b605d019ull,
  0x923f82a4af194f9bull, 0xab1c5ed5da6d8118ull,
  0xd807aa98a3030242ull, 0x12835b0145706fbeull,
  0x243185be4ee4b28cull, 0x550c7dc3d5ffb4e2ull,
  0x72be5d74f27b896full, 0x80deb1fe3b1696b1ull,
  0x9bdc06a725c71235ull, 0xc19bf174cf692694ull,
  0xe49b69c19ef14ad2ull, 0xefbe4786384f25e3ull,
  0x0fc19dc68b8cd5b5ull, 0x240ca1cc77ac9c65ull,
  0x2de92c6f592b0275ull, 0x4a7484aa6ea6e483ull,
  0x5cb0a9dcbd41fbd4ull, 0x76f988da831153b5ull,
  0x983e5152ee66dfabull, 0xa831c66d2db43210ull,
  0xb00327c898fb213full, 0xbf597fc7beef0ee4ull,
  0xc6e00bf33da88fc2ull, 0xd5a79147930aa725ull,
  0x06ca6351e003826full, 0x142929670a0e6e70ull,
  0x27b70a8546d22ffcull, 0x2e1b21385c26c926ull,
  0x4d2c6dfc5ac42aedull, 0x53380d139d95b3dfull,
  0x650a73548baf63deull, 0x766a0abb3c77b2a8ull,
  0x81c2c92e47edaee6ull, 0x92722c851482353bull,
  0xa2bfe8a14cf10364ull, 0xa81a664bbc423001ull,
  0xc24b8b70d0f89791ull, 0xc76c51a30654be30ull,
  0xd192e819d6ef5218ull, 0xd69906245565a910ull,
  0xf40e35855771202aull, 0x106aa07032bbd1b8ull,
  0x19a4c116b8d2d0c8ull, 0x1e376c085141ab53ull,
  0x2748774cdf8eeb99ull, 0x34b0bcb5e19b48a8ull,
  0x391c0cb3c5c95a63ull, 0x4ed8aa4ae3418acbull,
  0x5b9cca4f7763e373ull, 0x682e6ff3d6b2b8a3ull,
  0x748f82ee5defb2fcull, 0x78a5636f43172f60ull,
  0x84c87814a1f0ab72ull, 0x8cc702081a6439ecull,
  0x90befffa23631e28ull, 0xa4506cebde82bde9ull,
  0xbef9a3f7b2c67915ull, 0xc67178f2e372532bull,
  0xca273eceea26619cull, 0xd186b8c721c0c207ull,
  0xeada7dd6cde0eb1eull, 0xf57d4f7fee6ed178ull,
  0x06f067aa72176fbaull, 0x0a637dc5a2c898a6ull,
  0x113f9804bef90daeull, 0x1b710b35131c471bull,
  0x28db77f523047d84ull, 0x32caab7b40c72493ull,
  0x3c9ebe0a15c9bebcull, 0x431d67c49c100d4cull,
  0x4cc5d4becb3e42b6ull, 0x597f299cfc657e2aull,
  0x5fcb6fab3ad6faecull, 0x6c44198c4a475817ull
};

static const unsigned char sha512_P[128] = {
  0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

void
sha512_init(sha512_t *ctx) {
  ctx->state[0] = 0x6a09e667f3bcc908ull;
  ctx->state[1] = 0xbb67ae8584caa73bull;
  ctx->state[2] = 0x3c6ef372fe94f82bull;
  ctx->state[3] = 0xa54ff53a5f1d36f1ull;
  ctx->state[4] = 0x510e527fade682d1ull;
  ctx->state[5] = 0x9b05688c2b3e6c1full;
  ctx->state[6] = 0x1f83d9abfb41bd6bull;
  ctx->state[7] = 0x5be0cd19137e2179ull;
  ctx->size = 0;
}

static void
sha512_transform(sha512_t *ctx, const unsigned char *chunk) {
  const uint64_t *k = sha512_K;
  uint64_t data[16];
  uint64_t A, B, C, D, E, F, G, H;
  unsigned i;
  uint64_t *d;

  for (i = 0; i < 16; i++, chunk += 8)
    data[i] = read64be(chunk);

  A = ctx->state[0];
  B = ctx->state[1];
  C = ctx->state[2];
  D = ctx->state[3];
  E = ctx->state[4];
  F = ctx->state[5];
  G = ctx->state[6];
  H = ctx->state[7];

#define Ch(x, y, z) ((z) ^ ((x) & ((y) ^ (z))))
#define Maj(x, y, z) (((x) & (y)) ^ ((z) & ((x) ^ (y))))

#define S0(x) (ROTL64(36, (x)) ^ ROTL64(30, (x)) ^ ROTL64(25, (x)))
#define S1(x) (ROTL64(50, (x)) ^ ROTL64(46, (x)) ^ ROTL64(23, (x)))

#define s0(x) (ROTL64(63, (x)) ^ ROTL64(56, (x)) ^ ((x) >> 7))
#define s1(x) (ROTL64(45, (x)) ^ ROTL64(3, (x)) ^ ((x) >> 6))

#define EXPAND(W,i) \
(W[(i) & 15] += (s1(W[((i)-2) & 15]) + W[((i)-7) & 15] + s0(W[((i)-15) & 15])))

#define ROUND(a, b, c, d, e, f, g, h, k, data) do { \
  h += S1(e) + Ch(e,f,g) + k + data;                \
  d += h;                                           \
  h += S0(a) + Maj(a,b,c);                          \
} while (0)

  for (i = 0, d = data; i < 16; i += 8, k += 8, d += 8) {
    ROUND(A, B, C, D, E, F, G, H, k[0], d[0]);
    ROUND(H, A, B, C, D, E, F, G, k[1], d[1]);
    ROUND(G, H, A, B, C, D, E, F, k[2], d[2]);
    ROUND(F, G, H, A, B, C, D, E, k[3], d[3]);
    ROUND(E, F, G, H, A, B, C, D, k[4], d[4]);
    ROUND(D, E, F, G, H, A, B, C, k[5], d[5]);
    ROUND(C, D, E, F, G, H, A, B, k[6], d[6]);
    ROUND(B, C, D, E, F, G, H, A, k[7], d[7]);
  }

  for (; i < 80; i += 16, k += 16) {
    ROUND(A, B, C, D, E, F, G, H, k[ 0], EXPAND(data,  0));
    ROUND(H, A, B, C, D, E, F, G, k[ 1], EXPAND(data,  1));
    ROUND(G, H, A, B, C, D, E, F, k[ 2], EXPAND(data,  2));
    ROUND(F, G, H, A, B, C, D, E, k[ 3], EXPAND(data,  3));
    ROUND(E, F, G, H, A, B, C, D, k[ 4], EXPAND(data,  4));
    ROUND(D, E, F, G, H, A, B, C, k[ 5], EXPAND(data,  5));
    ROUND(C, D, E, F, G, H, A, B, k[ 6], EXPAND(data,  6));
    ROUND(B, C, D, E, F, G, H, A, k[ 7], EXPAND(data,  7));
    ROUND(A, B, C, D, E, F, G, H, k[ 8], EXPAND(data,  8));
    ROUND(H, A, B, C, D, E, F, G, k[ 9], EXPAND(data,  9));
    ROUND(G, H, A, B, C, D, E, F, k[10], EXPAND(data, 10));
    ROUND(F, G, H, A, B, C, D, E, k[11], EXPAND(data, 11));
    ROUND(E, F, G, H, A, B, C, D, k[12], EXPAND(data, 12));
    ROUND(D, E, F, G, H, A, B, C, k[13], EXPAND(data, 13));
    ROUND(C, D, E, F, G, H, A, B, k[14], EXPAND(data, 14));
    ROUND(B, C, D, E, F, G, H, A, k[15], EXPAND(data, 15));
  }

#undef Ch
#undef Maj
#undef S0
#undef S1
#undef s0
#undef s1
#undef EXPAND
#undef ROUND

  ctx->state[0] += A;
  ctx->state[1] += B;
  ctx->state[2] += C;
  ctx->state[3] += D;
  ctx->state[4] += E;
  ctx->state[5] += F;
  ctx->state[6] += G;
  ctx->state[7] += H;
}

void
sha512_update(sha512_t *ctx, const void *data, size_t len) {
  const unsigned char *bytes = (const unsigned char *)data;
  size_t pos = ctx->size & 127;
  size_t off = 0;

  ctx->size += len;

  if (pos > 0) {
    size_t want = 128 - pos;

    if (want > len)
      want = len;

    memcpy(ctx->block + pos, bytes + off, want);

    pos += want;
    len -= want;
    off += want;

    if (pos < 128)
      return;

    sha512_transform(ctx, ctx->block);
  }

  while (len >= 128) {
    sha512_transform(ctx, bytes + off);
    off += 128;
    len -= 128;
  }

  if (len > 0)
    memcpy(ctx->block, bytes + off, len);
}

void
sha512_final(sha512_t *ctx, unsigned char *out) {
  size_t pos = ctx->size & 127;
  uint64_t len = ctx->size << 3;
  unsigned char D[16];
  size_t i;

  write64be(D + 0, 0);
  write64be(D + 8, len);

  sha512_update(ctx, sha512_P, 1 + ((239 - pos) & 127));
  sha512_update(ctx, D, 16);

  for (i = 0; i < 8; i++)
    write64be(out + i * 8, ctx->state[i]);

  memset(ctx->state, 0x00, sizeof(ctx->state));
  memset(ctx->block, 0x00, sizeof(ctx->block));

  ctx->size = 0;
}

/*
 * SHA384
 */

void
sha384_init(sha384_t *ctx) {
  ctx->state[0] = 0xcbbb9d5dc1059ed8ull;
  ctx->state[1] = 0x629a292a367cd507ull;
  ctx->state[2] = 0x9159015a3070dd17ull;
  ctx->state[3] = 0x152fecd8f70e5939ull;
  ctx->state[4] = 0x67332667ffc00b31ull;
  ctx->state[5] = 0x8eb44a8768581511ull;
  ctx->state[6] = 0xdb0c2e0d64f98fa7ull;
  ctx->state[7] = 0x47b5481dbefa4fa4ull;
  ctx->size = 0;
}

void
sha384_update(sha384_t *ctx, const void *data, size_t len) {
  sha512_update(ctx, data, len);
}

void
sha384_final(sha384_t *ctx, unsigned char *out) {
  uint8_t buf[64];

  sha512_final(ctx, buf);

  memcpy(out, buf, 48);
  memset(buf, 0x00, 48);
}

/*
 * Keccak
 */

void
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
  static const uint64_t rc[24] = {
    0x0000000000000001ull, 0x0000000000008082ull,
    0x800000000000808aull, 0x8000000080008000ull,
    0x000000000000808bull, 0x0000000080000001ull,
    0x8000000080008081ull, 0x8000000000008009ull,
    0x000000000000008aull, 0x0000000000000088ull,
    0x0000000080008009ull, 0x000000008000000aull,
    0x000000008000808bull, 0x800000000000008bull,
    0x8000000000008089ull, 0x8000000000008003ull,
    0x8000000000008002ull, 0x8000000000000080ull,
    0x000000000000800aull, 0x800000008000000aull,
    0x8000000080008081ull, 0x8000000000008080ull,
    0x0000000080000001ull, 0x8000000080008008ull
  };

  uint64_t C[5], D[5], T, X;
  unsigned i, y;

#define A ctx->state

  C[0] = A[0] ^ A[5 + 0] ^ A[10 + 0] ^ A[15 + 0] ^ A[20 + 0];
  C[1] = A[1] ^ A[5 + 1] ^ A[10 + 1] ^ A[15 + 1] ^ A[20 + 1];
  C[2] = A[2] ^ A[5 + 2] ^ A[10 + 2] ^ A[15 + 2] ^ A[20 + 2];
  C[3] = A[3] ^ A[5 + 3] ^ A[10 + 3] ^ A[15 + 3] ^ A[20 + 3];
  C[4] = A[4] ^ A[5 + 4] ^ A[10 + 4] ^ A[15 + 4] ^ A[20 + 4];

  for (i = 0; i < 24; i++) {
    D[0] = C[4] ^ ROTL64(1, C[1]);
    D[1] = C[0] ^ ROTL64(1, C[2]);
    D[2] = C[1] ^ ROTL64(1, C[3]);
    D[3] = C[2] ^ ROTL64(1, C[4]);
    D[4] = C[3] ^ ROTL64(1, C[0]);

    A[0] ^= D[0];
    X = A[ 1] ^ D[1];     T = ROTL64( 1, X);
    X = A[ 6] ^ D[1]; A[ 1] = ROTL64(44, X);
    X = A[ 9] ^ D[4]; A[ 6] = ROTL64(20, X);
    X = A[22] ^ D[2]; A[ 9] = ROTL64(61, X);
    X = A[14] ^ D[4]; A[22] = ROTL64(39, X);
    X = A[20] ^ D[0]; A[14] = ROTL64(18, X);
    X = A[ 2] ^ D[2]; A[20] = ROTL64(62, X);
    X = A[12] ^ D[2]; A[ 2] = ROTL64(43, X);
    X = A[13] ^ D[3]; A[12] = ROTL64(25, X);
    X = A[19] ^ D[4]; A[13] = ROTL64( 8, X);
    X = A[23] ^ D[3]; A[19] = ROTL64(56, X);
    X = A[15] ^ D[0]; A[23] = ROTL64(41, X);
    X = A[ 4] ^ D[4]; A[15] = ROTL64(27, X);
    X = A[24] ^ D[4]; A[ 4] = ROTL64(14, X);
    X = A[21] ^ D[1]; A[24] = ROTL64( 2, X);
    X = A[ 8] ^ D[3]; A[21] = ROTL64(55, X);
    X = A[16] ^ D[1]; A[ 8] = ROTL64(45, X);
    X = A[ 5] ^ D[0]; A[16] = ROTL64(36, X);
    X = A[ 3] ^ D[3]; A[ 5] = ROTL64(28, X);
    X = A[18] ^ D[3]; A[ 3] = ROTL64(21, X);
    X = A[17] ^ D[2]; A[18] = ROTL64(15, X);
    X = A[11] ^ D[1]; A[17] = ROTL64(10, X);
    X = A[ 7] ^ D[2]; A[11] = ROTL64( 6, X);
    X = A[10] ^ D[0]; A[ 7] = ROTL64( 3, X);
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
      D[0] = ~A[y + 1] & A[y + 2];
      D[1] = ~A[y + 2] & A[y + 3];
      D[2] = ~A[y + 3] & A[y + 4];
      D[3] = ~A[y + 4] & A[y + 0];
      D[4] = ~A[y + 0] & A[y + 1];

      A[y + 0] ^= D[0]; C[0] ^= A[y + 0];
      A[y + 1] ^= D[1]; C[1] ^= A[y + 1];
      A[y + 2] ^= D[2]; C[2] ^= A[y + 2];
      A[y + 3] ^= D[3]; C[3] ^= A[y + 3];
      A[y + 4] ^= D[4]; C[4] ^= A[y + 4];
    }
  }
#undef A
}

static void
keccak_transform(keccak_t *ctx, const unsigned char *chunk) {
  size_t count = ctx->bs >> 3;
  size_t i;

  for (i = 0; i < count; i++)
    ctx->state[i] ^= read64le(chunk + i * 8);

  keccak_permute(ctx);
}

void
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

void
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

/*
 * Hash
 */

void
hash_init(hash_t *hash, int type) {
  hash->type = type;
  switch (hash->type) {
    case HASH_SHA256:
      sha256_init(&hash->ctx.sha256);
      break;
    case HASH_SHA384:
      sha384_init(&hash->ctx.sha384);
      break;
    case HASH_SHA512:
      sha512_init(&hash->ctx.sha512);
      break;
    case HASH_SHAKE256:
      keccak_init(&hash->ctx.keccak, 256);
      break;
    default:
      assert(0);
      break;
  }
}

void
hash_update(hash_t *hash, const void *data, size_t len) {
  switch (hash->type) {
    case HASH_SHA256:
      sha256_update(&hash->ctx.sha256, data, len);
      break;
    case HASH_SHA384:
      sha384_update(&hash->ctx.sha384, data, len);
      break;
    case HASH_SHA512:
      sha512_update(&hash->ctx.sha512, data, len);
      break;
    case HASH_SHAKE256:
      keccak_update(&hash->ctx.keccak, data, len);
      break;
    default:
      assert(0);
      break;
  }
}

void
hash_final(hash_t *hash, unsigned char *out, size_t len) {
  switch (hash->type) {
    case HASH_SHA256:
      sha256_final(&hash->ctx.sha256, out);
      break;
    case HASH_SHA384:
      sha384_final(&hash->ctx.sha384, out);
      break;
    case HASH_SHA512:
      sha512_final(&hash->ctx.sha512, out);
      break;
    case HASH_SHAKE256:
      keccak_final(&hash->ctx.keccak, out, 0x1f, len);
      break;
    default:
      assert(0);
      break;
  }
}

size_t
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

size_t
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

void
hmac_init(hmac_t *hmac, int type, const unsigned char *key, size_t len) {
  size_t hash_size = hash_output_size(type);
  size_t block_size = hash_block_size(type);
  unsigned char tmp[MAX_HASH_SIZE];
  unsigned char pad[MAX_BLOCK_SIZE];
  size_t i;

  hmac->type = type;

  if (len > block_size) {
    hash_init(&hmac->inner, type);
    hash_update(&hmac->inner, key, len);
    hash_final(&hmac->inner, tmp, hash_size);
    key = tmp;
    len = hash_size;
  }

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

  memset(tmp, 0x00, hash_size);
  memset(pad, 0x00, block_size);
}

void
hmac_update(hmac_t *hmac, const void *data, size_t len) {
  hash_update(&hmac->inner, data, len);
}

void
hmac_final(hmac_t *hmac, unsigned char *out) {
  size_t hash_size = hash_output_size(hmac->type);

  hash_final(&hmac->inner, out, hash_size);
  hash_update(&hmac->outer, out, hash_size);
  hash_final(&hmac->outer, out, hash_size);
}
