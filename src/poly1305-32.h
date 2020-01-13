/**
 * Parts of this software are based on poly1305-donna:
 * https://github.com/floodyberry/poly1305-donna
 *
 * MIT License
 * http://www.opensource.org/licenses/mit-license.php
 */

#define POLY1305_BLOCK_SIZE 16

/*
 * State
 */

typedef struct _poly1305_internal_s {
  uint32_t r[5];
  uint32_t h[5];
  uint32_t pad[4];
  size_t leftover;
  unsigned char buffer[POLY1305_BLOCK_SIZE];
  unsigned char final;
} poly1305_internal_t;

/*
 * Helpers
 */

static uint32_t
read32le(const void *src) {
#ifndef WORDS_BIGENDIAN
  uint32_t w;
  memcpy(&w, src, sizeof(w));
  return w;
#else
  const uint8_t *p = (const uint8_t *)src;
  return ((uint32_t)p[3] << 24)
       | ((uint32_t)p[2] << 16)
       | ((uint32_t)p[1] << 8)
       | ((uint32_t)p[0] << 0);
#endif
}

static void
write32le(void *dst, uint32_t w) {
#ifndef WORDS_BIGENDIAN
  memcpy(dst, &w, sizeof(w));
#else
  uint8_t *p = (uint8_t *)dst;
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

  /* r &= 0xffffffc0ffffffc0ffffffc0fffffff */
  st->r[0] = read32le(key + 0) & 0x3ffffff;
  st->r[1] = (read32le(key + 3) >> 2) & 0x3ffff03;
  st->r[2] = (read32le(key + 6) >> 4) & 0x3ffc0ff;
  st->r[3] = (read32le(key + 9) >> 6) & 0x3f03fff;
  st->r[4] = (read32le(key + 12) >> 8) & 0x00fffff;

  /* h = 0 */
  st->h[0] = 0;
  st->h[1] = 0;
  st->h[2] = 0;
  st->h[3] = 0;
  st->h[4] = 0;

  /* save pad for later */
  st->pad[0] = read32le(key + 16);
  st->pad[1] = read32le(key + 20);
  st->pad[2] = read32le(key + 24);
  st->pad[3] = read32le(key + 28);

  st->leftover = 0;
  st->final = 0;
}

static void
poly1305_blocks(poly1305_internal_t *st, const unsigned char *m, size_t bytes) {
  uint32_t hibit = st->final ? 0 : (1UL << 24); /* 1 << 128 */
  uint32_t r0, r1, r2, r3, r4;
  uint32_t s1, s2, s3, s4;
  uint32_t h0, h1, h2, h3, h4;
  uint64_t d0, d1, d2, d3, d4;
  uint32_t c;

  r0 = st->r[0];
  r1 = st->r[1];
  r2 = st->r[2];
  r3 = st->r[3];
  r4 = st->r[4];

  s1 = r1 * 5;
  s2 = r2 * 5;
  s3 = r3 * 5;
  s4 = r4 * 5;

  h0 = st->h[0];
  h1 = st->h[1];
  h2 = st->h[2];
  h3 = st->h[3];
  h4 = st->h[4];

  while (bytes >= POLY1305_BLOCK_SIZE) {
    /* h += m[i] */
    h0 += (read32le(m + 0)) & 0x3ffffff;
    h1 += (read32le(m + 3) >> 2) & 0x3ffffff;
    h2 += (read32le(m + 6) >> 4) & 0x3ffffff;
    h3 += (read32le(m + 9) >> 6) & 0x3ffffff;
    h4 += (read32le(m + 12) >> 8) | hibit;

    /* h *= r */
    d0 = ((uint64_t)h0 * r0)
       + ((uint64_t)h1 * s4)
       + ((uint64_t)h2 * s3)
       + ((uint64_t)h3 * s2)
       + ((uint64_t)h4 * s1);

    d1 = ((uint64_t)h0 * r1)
       + ((uint64_t)h1 * r0)
       + ((uint64_t)h2 * s4)
       + ((uint64_t)h3 * s3)
       + ((uint64_t)h4 * s2);

    d2 = ((uint64_t)h0 * r2)
       + ((uint64_t)h1 * r1)
       + ((uint64_t)h2 * r0)
       + ((uint64_t)h3 * s4)
       + ((uint64_t)h4 * s3);

    d3 = ((uint64_t)h0 * r3)
       + ((uint64_t)h1 * r2)
       + ((uint64_t)h2 * r1)
       + ((uint64_t)h3 * r0)
       + ((uint64_t)h4 * s4);

    d4 = ((uint64_t)h0 * r4)
       + ((uint64_t)h1 * r3)
       + ((uint64_t)h2 * r2)
       + ((uint64_t)h3 * r1)
       + ((uint64_t)h4 * r0);

    /* (partial) h %= p */
    c = (uint32_t)(d0 >> 26);
    h0 = (uint32_t)d0 & 0x3ffffff;
    d1 += c;

    c = (uint32_t)(d1 >> 26);
    h1 = (uint32_t)d1 & 0x3ffffff;
    d2 += c;

    c = (uint32_t)(d2 >> 26);
    h2 = (uint32_t)d2 & 0x3ffffff;
    d3 += c;

    c = (uint32_t)(d3 >> 26);
    h3 = (uint32_t)d3 & 0x3ffffff;
    d4 += c;

    c = (uint32_t)(d4 >> 26);
    h4 = (uint32_t)d4 & 0x3ffffff;
    h0 += c * 5;

    c = (h0 >> 26);
    h0 = h0 & 0x3ffffff;
    h1 += c;

    m += POLY1305_BLOCK_SIZE;
    bytes -= POLY1305_BLOCK_SIZE;
  }

  st->h[0] = h0;
  st->h[1] = h1;
  st->h[2] = h2;
  st->h[3] = h3;
  st->h[4] = h4;
}

void
poly1305_final(poly1305_t *ctx, unsigned char *mac) {
  poly1305_internal_t *st = (poly1305_internal_t *)ctx;
  uint32_t h0, h1, h2, h3, h4, c;
  uint32_t g0, g1, g2, g3, g4;
  uint64_t f;
  uint32_t mask;

  /* process the remaining block */
  if (st->leftover > 0) {
    size_t i = st->leftover;

    st->buffer[i++] = 1;

    for (; i < POLY1305_BLOCK_SIZE; i++)
      st->buffer[i] = 0;

    st->final = 1;

    poly1305_blocks(st, st->buffer, POLY1305_BLOCK_SIZE);
  }

  /* fully carry h */
  h0 = st->h[0];
  h1 = st->h[1];
  h2 = st->h[2];
  h3 = st->h[3];
  h4 = st->h[4];

  c = h1 >> 26;
  h1 = h1 & 0x3ffffff;
  h2 += c;

  c = h2 >> 26;
  h2 = h2 & 0x3ffffff;
  h3 += c;

  c = h3 >> 26;
  h3 = h3 & 0x3ffffff;
  h4 += c;

  c = h4 >> 26;
  h4 = h4 & 0x3ffffff;
  h0 += c * 5;

  c = h0 >> 26;
  h0 = h0 & 0x3ffffff;

  h1 += c;

  /* compute h + -p */
  g0 = h0 + 5;
  c = g0 >> 26;
  g0 &= 0x3ffffff;

  g1 = h1 + c;
  c = g1 >> 26;
  g1 &= 0x3ffffff;

  g2 = h2 + c;
  c = g2 >> 26;
  g2 &= 0x3ffffff;

  g3 = h3 + c;
  c = g3 >> 26;
  g3 &= 0x3ffffff;
  g4 = h4 + c - (1UL << 26);

  /* select h if h < p, or h + -p if h >= p */
  mask = (g4 >> 31) - 1;
  g0 &= mask;
  g1 &= mask;
  g2 &= mask;
  g3 &= mask;
  g4 &= mask;
  mask = ~mask;
  h0 = (h0 & mask) | g0;
  h1 = (h1 & mask) | g1;
  h2 = (h2 & mask) | g2;
  h3 = (h3 & mask) | g3;
  h4 = (h4 & mask) | g4;

  /* h = h % (2^128) */
  h0 = (h0 | (h1 << 26)) & 0xffffffff;
  h1 = ((h1 >> 6) | (h2 << 20)) & 0xffffffff;
  h2 = ((h2 >> 12) | (h3 << 14)) & 0xffffffff;
  h3 = ((h3 >> 18) | (h4 <<  8)) & 0xffffffff;

  /* mac = (h + pad) % (2^128) */
  f = (uint64_t)h0 + st->pad[0];
  h0 = (uint32_t)f;

  f = (uint64_t)h1 + st->pad[1] + (f >> 32);
  h1 = (uint32_t)f;

  f = (uint64_t)h2 + st->pad[2] + (f >> 32);
  h2 = (uint32_t)f;

  f = (uint64_t)h3 + st->pad[3] + (f >> 32);
  h3 = (uint32_t)f;

  write32le(mac + 0, h0);
  write32le(mac + 4, h1);
  write32le(mac + 8, h2);
  write32le(mac + 12, h3);

  /* zero out the state */
  st->h[0] = 0;
  st->h[1] = 0;
  st->h[2] = 0;
  st->h[3] = 0;
  st->h[4] = 0;
  st->r[0] = 0;
  st->r[1] = 0;
  st->r[2] = 0;
  st->r[3] = 0;
  st->r[4] = 0;
  st->pad[0] = 0;
  st->pad[1] = 0;
  st->pad[2] = 0;
  st->pad[3] = 0;
}
