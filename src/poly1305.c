/**
 * Parts of this software are based on poly1305-donna:
 * https://github.com/floodyberry/poly1305-donna
 *
 * MIT License
 * http://www.opensource.org/licenses/mit-license.php
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <torsion/poly1305.h>

#ifdef TORSION_64BIT
#include "poly1305-64.h"
#else
#include "poly1305-32.h"
#endif

void
poly1305_update(poly1305_t *ctx, const unsigned char *m, size_t bytes) {
  poly1305_internal_t *st = (poly1305_internal_t *)ctx;
  size_t i;

  /* handle leftover */
  if (st->leftover > 0) {
    size_t want = POLY1305_BLOCK_SIZE - st->leftover;

    if (want > bytes)
      want = bytes;

    for (i = 0; i < want; i++)
      st->buffer[st->leftover + i] = m[i];

    bytes -= want;
    m += want;
    st->leftover += want;

    if (st->leftover < POLY1305_BLOCK_SIZE)
      return;

    poly1305_blocks(st, st->buffer, POLY1305_BLOCK_SIZE);

    st->leftover = 0;
  }

  /* process full blocks */
  if (bytes >= POLY1305_BLOCK_SIZE) {
    size_t want = bytes & ~(POLY1305_BLOCK_SIZE - 1);

    poly1305_blocks(st, m, want);

    m += want;
    bytes -= want;
  }

  /* store leftover */
  if (bytes > 0) {
    for (i = 0; i < bytes; i++)
      st->buffer[st->leftover + i] = m[i];

    st->leftover += bytes;
  }
}

void
poly1305_auth(unsigned char *mac,
              const unsigned char *m,
              size_t bytes,
              const unsigned char *key) {
  poly1305_t ctx;
  poly1305_init(&ctx, key);
  poly1305_update(&ctx, m, bytes);
  poly1305_final(&ctx, mac);
}

int
poly1305_verify(const unsigned char *mac1, const unsigned char *mac2) {
  uint32_t z = 0;
  size_t i;

  for (i = 0; i < 16; i++)
    z |= (uint32_t)mac1[i] ^ (uint32_t)mac2[i];

  return ((z - 1) >> 31) & 1;
}
