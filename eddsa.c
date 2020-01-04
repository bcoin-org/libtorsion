#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include "field.h"

/*
 * ECDSA
 */

int RAND_status(void);
int RAND_poll();
int RAND_bytes(unsigned char *buf, int num);

static int
ecdsa_random_bytes(unsigned char *dst, size_t len) {
  memset(dst, 0x00, len);

  if (len > (size_t)INT_MAX)
    return 0;

  for (;;) {
    int status = RAND_status();

    assert(status >= 0);

    if (status != 0)
      break;

    if (RAND_poll() == 0)
      break;
  }

  return RAND_bytes(dst, (int)len) == 1;
}

int
main(void) {
  {
    edwards_t *ec = malloc(sizeof(edwards_t));
    ege_t g, p, q, r;
    xge_t jg, jp, jq, jr;
    unsigned char entropy[32];
    unsigned char p_raw[32];

    const unsigned char g_raw[32] = {
      0x58, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
      0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
      0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66,
      0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66, 0x66
    };

    const unsigned char g2_raw[32] = {
      0xc9, 0xa3, 0xf8, 0x6a, 0xae, 0x46, 0x5f, 0x0e,
      0x56, 0x51, 0x38, 0x64, 0x51, 0x0f, 0x39, 0x97,
      0x56, 0x1f, 0xa2, 0xc9, 0xe8, 0x5e, 0xa2, 0x1d,
      0xc2, 0x29, 0x23, 0x09, 0xf3, 0xcd, 0x60, 0x22
    };

    const unsigned char g3_raw[32] = {
      0xd4, 0xb4, 0xf5, 0x78, 0x48, 0x68, 0xc3, 0x02,
      0x04, 0x03, 0x24, 0x67, 0x17, 0xec, 0x16, 0x9f,
      0xf7, 0x9e, 0x26, 0x60, 0x8e, 0xa1, 0x26, 0xa1,
      0xab, 0x69, 0xee, 0x77, 0xd1, 0xb1, 0x67, 0x12
    };

    edwards_init(ec, &curve_ed25519);

    ecdsa_random_bytes(entropy, sizeof(entropy));

    edwards_randomize(ec, entropy);

    ege_set(ec, &g, &ec->g);
    ege_to_xge(ec, &jg, &ec->g);

    assert(ege_import(ec, &p, g_raw));

    ege_to_xge(ec, &jp, &p);
    ege_to_xge(ec, &jq, &ec->g);

    assert(ege_validate(ec, &p));
    assert(xge_validate(ec, &jp));
    assert(xge_validate(ec, &jq));
    assert(ege_equal(ec, &p, &ec->g));
    assert(xge_equal(ec, &jp, &jq));

    assert(ege_import(ec, &q, g2_raw));
    assert(ege_import(ec, &r, g3_raw));

    ege_to_xge(ec, &jq, &q);
    ege_to_xge(ec, &jr, &r);

    ege_dbl(ec, &p, &ec->g);

    assert(ege_equal(ec, &p, &q));

    ege_add(ec, &p, &p, &ec->g);

    assert(ege_equal(ec, &p, &r));

    xge_dbl(ec, &jp, &jg);

    assert(xge_equal(ec, &jp, &jq));

    xge_add(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jr));

    xge_sub(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jq));

    xge_add(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jr));

    xge_sub(ec, &jp, &jp, &jg);

    assert(xge_equal(ec, &jp, &jq));

    assert(xge_validate(ec, &jg));
    assert(xge_validate(ec, &jp));
    assert(xge_validate(ec, &jq));
    assert(xge_validate(ec, &jr));

    assert(!xge_is_zero(ec, &jg));
    assert(!xge_is_zero(ec, &jp));
    assert(!xge_is_zero(ec, &jq));
    assert(!xge_is_zero(ec, &jr));

    xge_to_ege(ec, &p, &jp);

    assert(ege_equal(ec, &p, &q));

    assert(ege_export(ec, p_raw, &p));
    assert(memcmp(p_raw, g2_raw, 32) == 0);

    free(ec);
  }

  return 0;
}
