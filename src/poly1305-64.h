/**
 * Parts of this software are based on poly1305-donna:
 * https://github.com/floodyberry/poly1305-donna
 *
 * MIT License
 * http://www.opensource.org/licenses/mit-license.php
 */

#if defined(_MSC_VER)
#include <intrin.h>

typedef struct _uint128_s {
  uint64_t lo;
  uint64_t hi;
} uint128_t;

#define MUL(out, x, y) \
  out.lo = _umul128((x), (y), &out.hi)

#define ADD(out, in) do {         \
  uint64_t t = out.lo;            \
  out.lo += in.lo;                \
  out.hi += (out.lo < t) + in.hi; \
} while (0)

#define ADDLO(out, in) do { \
  uint64_t t = out.lo;      \
  out.lo += in;             \
  out.hi += (out.lo < t);   \
} while (0)

#define SHR(in, shift) (__shiftright128(in.lo, in.hi, (shift)))
#define LO(in) (in.lo)

#elif defined(__GNUC__)

#if defined(__SIZEOF_INT128__)
typedef unsigned __int128 uint128_t;
#else
typedef unsigned uint128_t __attribute__((mode(TI)));
#endif

#define MUL(out, x, y) out = ((uint128_t)x * y)
#define ADD(out, in) out += in
#define ADDLO(out, in) out += in
#define SHR(in, shift) (uint64_t)(in >> (shift))
#define LO(in) (uint64_t)(in)

#endif

#define POLY1305_BLOCK_SIZE 16

/*
 * State
 */

typedef struct _poly1305_internal_s {
  uint64_t r[3];
  uint64_t h[3];
  uint64_t pad[2];
  size_t leftover;
  unsigned char buffer[POLY1305_BLOCK_SIZE];
  unsigned char final;
} poly1305_internal_t;

/*
 * Helpers
 */

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

static void
write64le(void *dst, uint64_t w) {
#ifndef WORDS_BIGENDIAN
  memcpy(dst, &w, sizeof(w));
#else
  uint8_t *p = (uint8_t *)dst;
  p[7] = w >> 56;
  p[6] = w >> 48;
  p[5] = w >> 40;
  p[4] = w >> 32;
  p[3] = w >> 24;
  p[2] = w >> 16;
  p[1] = w >> 8;
  p[0] = w >> 0;
#endif
}

/*
 * Poly1305
 */

void
poly1305_init(poly1305_t *ctx, const unsigned char *key) {
  poly1305_internal_t *st = (poly1305_internal_t *)ctx;
  uint64_t t0, t1;

  /* r &= 0xffffffc0ffffffc0ffffffc0fffffff */
  t0 = read64le(key + 0);
  t1 = read64le(key + 8);

  st->r[0] = t0 & 0xffc0fffffff;
  st->r[1] = ((t0 >> 44) | (t1 << 20)) & 0xfffffc0ffff;
  st->r[2] = (t1 >> 24) & 0x00ffffffc0f;

  /* h = 0 */
  st->h[0] = 0;
  st->h[1] = 0;
  st->h[2] = 0;

  /* save pad for later */
  st->pad[0] = read64le(key + 16);
  st->pad[1] = read64le(key + 24);

  st->leftover = 0;
  st->final = 0;
}

static void
poly1305_blocks(poly1305_internal_t *st, const unsigned char *m, size_t bytes) {
  uint64_t hibit = st->final ? 0 : ((uint64_t)1 << 40); /* 1 << 128 */
  uint64_t r0, r1, r2;
  uint64_t s1, s2;
  uint64_t h0, h1, h2;
  uint64_t c;
  uint128_t d0, d1, d2, d;

  r0 = st->r[0];
  r1 = st->r[1];
  r2 = st->r[2];

  h0 = st->h[0];
  h1 = st->h[1];
  h2 = st->h[2];

  s1 = r1 * (5 << 2);
  s2 = r2 * (5 << 2);

  while (bytes >= POLY1305_BLOCK_SIZE) {
    uint64_t t0, t1;

    /* h += m[i] */
    t0 = read64le(m + 0);
    t1 = read64le(m + 8);

    h0 += t0 & 0xfffffffffff;
    h1 += ((t0 >> 44) | (t1 << 20)) & 0xfffffffffff;
    h2 += (((t1 >> 24)) & 0x3ffffffffff) | hibit;

    /* h *= r */
    MUL(d0, h0, r0);
    MUL(d, h1, s2);
    ADD(d0, d);
    MUL(d, h2, s1);
    ADD(d0, d);

    MUL(d1, h0, r1);
    MUL(d, h1, r0);
    ADD(d1, d);
    MUL(d, h2, s2);
    ADD(d1, d);

    MUL(d2, h0, r2);
    MUL(d, h1, r1);
    ADD(d2, d);
    MUL(d, h2, r0);
    ADD(d2, d);

    /* (partial) h %= p */
    c = SHR(d0, 44);
    h0 = LO(d0) & 0xfffffffffff;

    ADDLO(d1, c);
    c = SHR(d1, 44);
    h1 = LO(d1) & 0xfffffffffff;

    ADDLO(d2, c);
    c = SHR(d2, 42);
    h2 = LO(d2) & 0x3ffffffffff;

    h0 += c * 5;
    c = (h0 >> 44);
    h0 = h0 & 0xfffffffffff;

    h1 += c;

    m += POLY1305_BLOCK_SIZE;
    bytes -= POLY1305_BLOCK_SIZE;
  }

  st->h[0] = h0;
  st->h[1] = h1;
  st->h[2] = h2;
}

void
poly1305_final(poly1305_t *ctx, unsigned char *mac) {
  poly1305_internal_t *st = (poly1305_internal_t *)ctx;
  uint64_t h0, h1, h2, c;
  uint64_t g0, g1, g2;
  uint64_t t0, t1;

  /* process the remaining block */
  if (st->leftover > 0) {
    size_t i = st->leftover;

    st->buffer[i] = 1;

    for (i = i + 1; i < POLY1305_BLOCK_SIZE; i++)
      st->buffer[i] = 0;

    st->final = 1;

    poly1305_blocks(st, st->buffer, POLY1305_BLOCK_SIZE);
  }

  /* fully carry h */
  h0 = st->h[0];
  h1 = st->h[1];
  h2 = st->h[2];

  c = (h1 >> 44);
  h1 &= 0xfffffffffff;

  h2 += c;
  c = (h2 >> 42);
  h2 &= 0x3ffffffffff;

  h0 += c * 5;
  c = (h0 >> 44);
  h0 &= 0xfffffffffff;

  h1 += c;
  c = (h1 >> 44);
  h1 &= 0xfffffffffff;

  h2 += c;
  c = (h2 >> 42);
  h2 &= 0x3ffffffffff;

  h0 += c * 5;
  c = (h0 >> 44);
  h0 &= 0xfffffffffff;
  h1 += c;

  /* compute h + -p */
  g0 = h0 + 5;
  c = (g0 >> 44);
  g0 &= 0xfffffffffff;

  g1 = h1 + c;
  c = (g1 >> 44);
  g1 &= 0xfffffffffff;
  g2 = h2 + c - ((uint64_t)1 << 42);

  /* select h if h < p, or h + -p if h >= p */
  c = (g2 >> 63) - 1;
  g0 &= c;
  g1 &= c;
  g2 &= c;
  c = ~c;
  h0 = (h0 & c) | g0;
  h1 = (h1 & c) | g1;
  h2 = (h2 & c) | g2;

  /* h = (h + pad) */
  t0 = st->pad[0];
  t1 = st->pad[1];

  h0 += (t0 & 0xfffffffffff);
  c = (h0 >> 44);
  h0 &= 0xfffffffffff;

  h1 += (((t0 >> 44) | (t1 << 20)) & 0xfffffffffff) + c;
  c = (h1 >> 44);
  h1 &= 0xfffffffffff;

  h2 += (((t1 >> 24)) & 0x3ffffffffff) + c;
  h2 &= 0x3ffffffffff;

  /* mac = h % (2^128) */
  h0 = (h0 | (h1 << 44));
  h1 = ((h1 >> 20) | (h2 << 24));

  write64le(mac + 0, h0);
  write64le(mac + 8, h1);

  /* zero out the state */
  st->h[0] = 0;
  st->h[1] = 0;
  st->h[2] = 0;
  st->r[0] = 0;
  st->r[1] = 0;
  st->r[2] = 0;
  st->pad[0] = 0;
  st->pad[1] = 0;
}
