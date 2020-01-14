#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <torsion/drbg.h>
#include <torsion/hash.h>
#include <torsion/rsa.h>
#include <torsion/util.h>
#include "asn1.h"
#include "mpz.h"

/*
 * Constants
 */

static const unsigned char digest_info[32][24] = {
  { /* BLAKE2B160 */
    0x15, 0x30, 0x27, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x01, 0x05, 0x05, 0x00, 0x04, 0x14, 0x00, 0x00
  },
  { /* BLAKE2B256 */
    0x15, 0x30, 0x33, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x01, 0x08, 0x05, 0x00, 0x04, 0x20, 0x00, 0x00
  },
  { /* BLAKE2B384 */
    0x15, 0x30, 0x43, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x01, 0x0c, 0x05, 0x00, 0x04, 0x30, 0x00, 0x00
  },
  { /* BLAKE2B512 */
    0x15, 0x30, 0x53, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x01, 0x10, 0x05, 0x00, 0x04, 0x40, 0x00, 0x00
  },
  { /* BLAKE2S128 */
    0x15, 0x30, 0x23, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x02, 0x04, 0x05, 0x00, 0x04, 0x10, 0x00, 0x00
  },
  { /* BLAKE2S160 */
    0x15, 0x30, 0x27, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x02, 0x05, 0x05, 0x00, 0x04, 0x14, 0x00, 0x00
  },
  { /* BLAKE2S224 */
    0x15, 0x30, 0x2f, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x02, 0x07, 0x05, 0x00, 0x04, 0x1c, 0x00, 0x00
  },
  { /* BLAKE2S256 */
    0x15, 0x30, 0x33, 0x30, 0x0f, 0x06, 0x0b, 0x2b,
    0x06, 0x01, 0x04, 0x01, 0x8d, 0x3a, 0x0c, 0x02,
    0x02, 0x08, 0x05, 0x00, 0x04, 0x20, 0x00, 0x00
  },
  { /* GOST94 */
    0x10, 0x30, 0x2e, 0x30, 0x0a, 0x06, 0x06, 0x2a,
    0x85, 0x03, 0x02, 0x02, 0x14, 0x05, 0x00, 0x04,
    0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* HASH160 */
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* HASH256 */
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* KECCAK224 */
    0x13, 0x30, 0x2d, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x07,
    0x05, 0x00, 0x04, 0x1c, 0x00, 0x00, 0x00, 0x00
  },
  { /* KECCAK256 */
    0x13, 0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x08,
    0x05, 0x00, 0x04, 0x20, 0x00, 0x00, 0x00, 0x00
  },
  { /* KECCAK384 */
    0x13, 0x30, 0x41, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x09,
    0x05, 0x00, 0x04, 0x30, 0x00, 0x00, 0x00, 0x00
  },
  { /* KECCAK512 */
    0x13, 0x30, 0x51, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x0a,
    0x05, 0x00, 0x04, 0x40, 0x00, 0x00, 0x00, 0x00
  },
  { /* MD2 */
    0x12, 0x30, 0x20, 0x30, 0x0c, 0x06, 0x08, 0x2a,
    0x86, 0x48, 0x86, 0xf7, 0x0d, 0x02, 0x02, 0x05,
    0x00, 0x04, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* MD4 */
    0x12, 0x30, 0x20, 0x30, 0x0c, 0x06, 0x08, 0x2a,
    0x86, 0x48, 0x86, 0xf7, 0x0d, 0x02, 0x04, 0x05,
    0x00, 0x04, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* MD5 */
    0x12, 0x30, 0x20, 0x30, 0x0c, 0x06, 0x08, 0x2a,
    0x86, 0x48, 0x86, 0xf7, 0x0d, 0x02, 0x05, 0x05,
    0x00, 0x04, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* MD5SHA1 */
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* RIPEMD160 */
    0x10, 0x30, 0x22, 0x30, 0x0a, 0x06, 0x06, 0x28,
    0xcf, 0x06, 0x03, 0x00, 0x31, 0x05, 0x00, 0x04,
    0x14, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA1 */
    0x0f, 0x30, 0x21, 0x30, 0x09, 0x06, 0x05, 0x2b,
    0x0e, 0x03, 0x02, 0x1a, 0x05, 0x00, 0x04, 0x14,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA224 */
    0x13, 0x30, 0x2d, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x04,
    0x05, 0x00, 0x04, 0x1c, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA256 */
    0x13, 0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x01,
    0x05, 0x00, 0x04, 0x20, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA384 */
    0x13, 0x30, 0x41, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x02,
    0x05, 0x00, 0x04, 0x30, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA512 */
    0x13, 0x30, 0x51, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x03,
    0x05, 0x00, 0x04, 0x40, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA3_224 */
    0x13, 0x30, 0x2d, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x07,
    0x05, 0x00, 0x04, 0x1c, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA3_256 */
    0x13, 0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x08,
    0x05, 0x00, 0x04, 0x20, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA3_384 */
    0x13, 0x30, 0x41, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x09,
    0x05, 0x00, 0x04, 0x30, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHA3_512 */
    0x13, 0x30, 0x51, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x0a,
    0x05, 0x00, 0x04, 0x40, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHAKE128 */
    0x13, 0x30, 0x21, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x0b,
    0x05, 0x00, 0x04, 0x10, 0x00, 0x00, 0x00, 0x00
  },
  { /* SHAKE256 */
    0x13, 0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60,
    0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, 0x0c,
    0x05, 0x00, 0x04, 0x20, 0x00, 0x00, 0x00, 0x00
  },
  { /* WHIRLPOOL */
    0x10, 0x30, 0x4e, 0x30, 0x0a, 0x06, 0x06, 0x28,
    0xcf, 0x06, 0x03, 0x00, 0x37, 0x05, 0x00, 0x04,
    0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
  }
};

static const unsigned char pss_prefix[8] = {0, 0, 0, 0, 0, 0, 0, 0};

/*
 * Structs
 */

typedef struct _rsa_pub_s {
  mpz_t n;
  mpz_t e;
} rsa_pub_t;

typedef struct _rsa_priv_s {
  mpz_t n;
  mpz_t e;
  mpz_t d;
  mpz_t p;
  mpz_t q;
  mpz_t dp;
  mpz_t dq;
  mpz_t qi;
} rsa_priv_t;

/*
 * Helpers
 */

static uint32_t
safe_equal(uint32_t x, uint32_t y) {
  return ((x ^ y) - 1) >> 31;
}

static uint32_t
safe_select(uint32_t x, uint32_t y, uint32_t v) {
  return (x & (v - 1)) | (y & ~(v - 1));
}

static uint32_t
safe_lte(uint32_t x, uint32_t y) {
  return (x - y - 1) >> 31;
}

static uint32_t
safe_cmp(const unsigned char *x, const unsigned char *y, size_t len) {
  uint32_t v = 0;
  size_t i;

  for (i = 0; i < len; i++)
    v |= (uint32_t)x[i] ^ (uint32_t)y[i];

  return (v - 1) >> 31;
}

static void *
safe_malloc(size_t size) {
  if (size == 0)
    return NULL;

  return calloc(size, sizeof(unsigned char));
}

static void
safe_free(void *ptr, size_t size) {
  if (ptr != NULL) {
    cleanse(ptr, size);
    free(ptr);
  }
}

/*
 * Private Key
 */

static void
rsa_priv_init(rsa_priv_t *k) {
  mpz_init(k->n);
  mpz_init(k->e);
  mpz_init(k->d);
  mpz_init(k->p);
  mpz_init(k->q);
  mpz_init(k->dp);
  mpz_init(k->dq);
  mpz_init(k->qi);
}

static void
rsa_priv_clear(rsa_priv_t *k) {
  mpz_cleanse(k->n);
  mpz_cleanse(k->e);
  mpz_cleanse(k->d);
  mpz_cleanse(k->p);
  mpz_cleanse(k->q);
  mpz_cleanse(k->dp);
  mpz_cleanse(k->dq);
  mpz_cleanse(k->qi);
}

static int
rsa_priv_import(rsa_priv_t *k,
                const unsigned char *data,
                size_t len,
                int strict) {
  if (!asn1_read_seq(&data, &len, strict))
    return 0;

  /* Read version first. */
  if (!asn1_read_int(k->n, &data, &len, strict))
    return 0;

  /* Should be zero. */
  if (strict && mpz_sgn(k->n) != 0)
    return 0;

  if (!asn1_read_int(k->n, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->e, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->d, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->p, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->q, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->dp, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->dq, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->qi, &data, &len, strict))
    return 0;

  if (strict && len != 0)
    return 0;

  return 1;
}

static void
rsa_priv_export(unsigned char *out, size_t *out_len, const rsa_priv_t *k) {
  size_t size = 0;
  size_t pos = 0;

  size += 3; /* version */
  size += asn1_size_int(k->n);
  size += asn1_size_int(k->e);
  size += asn1_size_int(k->d);
  size += asn1_size_int(k->p);
  size += asn1_size_int(k->q);
  size += asn1_size_int(k->dp);
  size += asn1_size_int(k->dq);
  size += asn1_size_int(k->qi);

  /* 0x30 [size] [body] */
  out[pos++] = 0x30;
  pos += asn1_write_size(out, pos, size);
  out[pos++] = 0x02;
  out[pos++] = 0x01;
  out[pos++] = 0x00; /* version */
  pos += asn1_write_int(out, pos, k->n);
  pos += asn1_write_int(out, pos, k->e);
  pos += asn1_write_int(out, pos, k->d);
  pos += asn1_write_int(out, pos, k->p);
  pos += asn1_write_int(out, pos, k->q);
  pos += asn1_write_int(out, pos, k->dp);
  pos += asn1_write_int(out, pos, k->dq);
  pos += asn1_write_int(out, pos, k->qi);

  *out_len = pos;
}

static int
rsa_priv_generate(rsa_priv_t *k,
                  size_t bits,
                  uint64_t exp,
                  const unsigned char *entropy) {
  /* [RFC8017] Page 9, Section 3.2.
   * [FIPS186] Page 51, Appendix B.3.1
   *           Page 55, Appendix B.3.3
   *
   * There are two methods for choosing `d`.
   * Implementations differ on whether they
   * use Euler's totient or the Carmichael
   * function.
   *
   * The best explanation of Euler's phi vs.
   * Carmichael's lambda I've seen comes from
   * the crypto stackexchange[1].
   *
   * Note that both functions are _equivalent_
   * when used with RSA, however, Carmichael's
   * may lend itself to some perf benefits.
   *
   * We currently use Euler's totient in order
   * to maintain compatibility with OpenSSL.
   *
   * [1] https://crypto.stackexchange.com/a/29595
   */
  mpz_t pm1, qm1, phi;
  drbg_t rng;

  if (bits < RSA_MIN_MOD_BITS
      || bits > RSA_MAX_MOD_BITS
      || exp < RSA_MIN_EXP
      || exp > RSA_MAX_EXP
      || (exp & 1ull) == 0) {
    return 0;
  }

  drbg_init(&rng, HASH_SHA256, entropy, 32);

  mpz_init(pm1);
  mpz_init(qm1);
  mpz_init(phi);

  if ((exp >> 32) == 0) {
    mpz_set_ui(k->e, exp);
  } else {
    mpz_set_ui(k->e, exp >> 32);
    mpz_mul_2exp(k->e, k->e, 32);
    mpz_add_ui(k->e, k->e, exp & 0xffffffff);
  }

  for (;;) {
    mpz_random_prime(k->p, (bits >> 1) + (bits & 1), &rng);
    mpz_random_prime(k->q, bits >> 1, &rng);

    if (mpz_cmp(k->p, k->q) == 0)
      continue;

    if (mpz_cmp(k->p, k->q) < 0)
      mpz_swap(k->p, k->q);

    mpz_sub(pm1, k->p, k->q);

    if (mpz_bitlen(pm1) <= (bits >> 1) - 99)
      continue;

    mpz_mul(k->n, k->p, k->q);

    if (mpz_bitlen(k->n) != bits)
      continue;

    mpz_sub_ui(pm1, k->p, 1);
    mpz_sub_ui(qm1, k->q, 1);
    mpz_mul(phi, pm1, qm1);

    if (!mpz_invert(k->d, k->e, phi))
      continue;

    mpz_mod(k->dp, k->d, pm1);
    mpz_mod(k->dq, k->d, qm1);

    assert(mpz_invert(k->qi, k->q, k->p));

    break;
  }

  cleanse(&rng, sizeof(rng));

  mpz_cleanse(pm1);
  mpz_cleanse(qm1);
  mpz_cleanse(phi);

  return 1;
}

static int
rsa_priv_verify(const rsa_priv_t *k) {
  /* [RFC8017] Page 9, Section 3.2. */
  size_t nb = mpz_bitlen(k->n);
  size_t eb = mpz_bitlen(k->e);
  size_t db = mpz_bitlen(k->d);
  size_t pb = mpz_bitlen(k->p);
  size_t qb = mpz_bitlen(k->q);
  size_t dpb = mpz_bitlen(k->dp);
  size_t dqb = mpz_bitlen(k->dq);
  size_t qib = mpz_bitlen(k->qi);
  mpz_t pm1, qm1, lam, t;
  int r = 0;

  mpz_init(pm1);
  mpz_init(qm1);
  mpz_init(lam);
  mpz_init(t);

  if (nb < RSA_MIN_MOD_BITS || nb > RSA_MAX_MOD_BITS)
    goto fail;

  if (eb < RSA_MIN_EXP_BITS || eb > RSA_MAX_EXP_BITS)
    goto fail;

  /* d < (p - 1) * (q - 1) */
  if (db == 0 || db > nb)
    goto fail;

  /* p < n */
  if (pb <= 1 || pb > nb)
    goto fail;

  /* q < n */
  if (qb <= 1 || qb > nb)
    goto fail;

  /* dp < p - 1 */
  if (dpb == 0 || dpb > pb)
    goto fail;

  /* dq < q - 1 */
  if (dqb == 0 || dqb > qb)
    goto fail;

  /* qi < p */
  if (qib == 0 || qib > pb)
    goto fail;

  if (!mpz_odd_p(k->n))
    goto fail;

  if (!mpz_odd_p(k->e))
    goto fail;

  if (!mpz_odd_p(k->p))
    goto fail;

  if (!mpz_odd_p(k->q))
    goto fail;

  /* n == p * q */
  mpz_mul(t, k->p, k->q);

  if (mpz_cmp(t, k->n) != 0)
    goto fail;

  /* e * d mod lcm(p - 1, q - 1) == 1 */
  mpz_sub_ui(pm1, k->p, 1);
  mpz_sub_ui(qm1, k->q, 1);
  mpz_lcm(lam, pm1, qm1);
  mpz_mul(t, k->e, k->d);
  mpz_mod(t, t, lam);

  if (mpz_cmp_ui(t, 1) != 0)
    goto fail;

  /* dp == d mod (p - 1) */
  mpz_mod(t, k->d, pm1);

  if (mpz_cmp(t, k->dp) != 0)
    goto fail;

  /* dq == d mod (q - 1) */
  mpz_mod(t, k->d, qm1);

  if (mpz_cmp(t, k->dq) != 0)
    goto fail;

  /* q * qi mod p == 1 */
  mpz_mul(t, k->q, k->qi);
  mpz_mod(t, t, k->p);

  if (mpz_cmp_ui(t, 1) != 0)
    goto fail;

  r = 1;
fail:
  mpz_cleanse(pm1);
  mpz_cleanse(qm1);
  mpz_cleanse(lam);
  mpz_cleanse(t);
  return r;
}

static int
rsa_priv_recover_pqe(rsa_priv_t *k) {
  /* Recover from (p, q, e) or (p, q, d). */
  mpz_t pm1, qm1, phi;
  int r = 0;

  mpz_init(pm1);
  mpz_init(qm1);
  mpz_init(phi);
  mpz_init(lam);

  /* Sanity check. */
  if (mpz_sgn(k->e) == 0 && mpz_sgn(k->d) == 0)
    goto fail;

  mpz_mul(k->n, k->p, k->q);
  mpz_sub_ui(pm1, k->p, 1);
  mpz_sub_ui(qm1, k->q, 1);

  /* Another sanity check. */
  if (mpz_sgn(pm1) <= 0 || mpz_sgn(qm1) <= 0)
    goto fail;

  mpz_mul(phi, pm1, qm1);

  if (mpz_sgn(k->e) == 0) {
    mpz_lcm(lam, pm1, qm1);

    if (!mpz_invert(k->e, k->d, phi)
        || !mpz_invert(k->e, k->d, lam)) {
      goto fail;
    }
  } else {
    if (!mpz_invert(k->d, k->e, phi))
      goto fail;
  }

  mpz_mod(k->dp, k->d, pm1);
  mpz_mod(k->dq, k->d, qm1);

  if (!mpz_invert(k->qi, k->q, k->p))
    goto fail;

  if (!rsa_priv_verify(k))
    goto fail;

  r = 1;
fail:
  mpz_cleanse(pm1);
  mpz_cleanse(qm1);
  mpz_cleanse(phi);
  mpz_cleanse(lam);
  return r;
}

static int
rsa_priv_recover_ned(rsa_priv_t *k, const unsigned char *entropy) {
  /* Factor an RSA modulus given (n, e, d).
   *
   * This is basically the same logic as the
   * Miller-Rabin primality test[1][2].
   *
   * [1] https://crypto.stackexchange.com/questions/11509
   * [2] https://crypto.stackexchange.com/questions/22374
   */
  mpz_t f, nm1, nm3, g, a, b, c;
  size_t i, j, s;
  drbg_t rng;
  int r = 0;

  mpz_init(f);
  mpz_init(nm1);
  mpz_init(nm3);
  mpz_init(g);
  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  /* f = e * d - 1 */
  mpz_mul(f, k->e, k->d);
  mpz_sub_ui(f, f, 1);

  /* nm1 = n - 1 */
  mpz_sub_ui(nm1, k->n, 1);

  /* nm3 = nm1 - 2 */
  mpz_sub_ui(nm3, nm1, 2);

  /* Sanity check. */
  if (mpz_sgn(f) <= 0 || mpz_sgn(nm3) <= 0)
    goto fail;

  /* s = f factors of 2 */
  s = mpz_scan1(f, 0);

  /* g = f >> s */
  mpz_tdiv_q_2exp(g, f, s);

  /* Seed RNG. */
  drbg_init(&rng, HASH_SHA256, entropy, 32);

  for (i = 0; i < 128; i++) {
    /* a = random int in [2,n-1] */
    mpz_random_int(a, nm3, &rng);
    mpz_add_ui(a, a, 2);

    /* b = a^g mod n */
    mpz_powm(b, a, g, n);

    if (mpz_cmp_ui(b, 1) == 0 || mpz_cmp(b, nm1) == 0)
      continue;

    for (j = 1; j < s; j++) {
      /* c = b^2 mod n */
      mpz_mul(c, b, b);
      mpz_mod(c, c, n);

      if (mpz_cmp_ui(c, 1) == 0) {
        /* p = gcd(n, b - 1) */
        mpz_sub_ui(c, b, 1);
        mpz_gcd(k->p, k->n, c);

        /* q = gcd(n, b + 1) */
        mpz_add_ui(c, b, 1);
        mpz_gcd(k->q, k->n, c);

        if (mpz_cmp(k->p, k->q) < 0)
          mpz_swap(k->p, k->q);

        mpz_sub_ui(b, p, 1);
        mpz_sub_ui(c, q, 1);
        mpz_mod(k->dp, k->d, b);
        mpz_mod(k->dq, k->d, c);

        if (!mpz_invert(k->qi, k->q, k->p))
          goto fail;

        if (!rsa_priv_verify(k))
          goto fail;

        goto succeed;
      }

      if (mpz_cmp(c, nm1) == 0)
        break;

      mpz_set(b, c);
    }
  }

  goto fail;
succeed:
  r = 1;
fail:
  cleanse(&rng, sizeof(rng));
  mpz_cleanse(f);
  mpz_cleanse(nm1);
  mpz_cleanse(nm3);
  mpz_cleanse(g);
  mpz_cleanse(a);
  mpz_cleanse(b);
  mpz_cleanse(c);
  return r;
}

static int
rsa_priv_recover(rsa_priv_t *k, const unsigned char *entropy) {
  if (mpz_sgn(k->p) > 0 && mpz_sgn(k->q) > 0)
    return rsa_priv_recover_pq(k);

  return rsa_priv_recover_ned(k, entropy);
}

static int
rsa_priv_decrypt(const rsa_priv_t *k,
                 unsigned char *out,
                 const unsigned char *msg,
                 size_t msg_len,
                 const unsigned char *entropy) {
  /* [RFC8017] Page 13, Section 5.1.2.
   *           Page 15, Section 5.2.1.
   */
  mpz_t t, s, b, bi, m;
  drbg_t rng;
  int r = 0;

  drbg_init(&rng, HASH_SHA256, entropy, 32);

  mpz_init(t);
  mpz_init(s);
  mpz_init(b);
  mpz_init(bi);
  mpz_init(m);

  if (mpz_sgn(k->n) <= 0 || mpz_sgn(k->d) <= 0)
    goto fail;

  mpz_import(m, msg_len, 1, 1, 0, 0, msg);

  if (mpz_cmp(m, k->n) >= 0)
    goto fail;

  /* t = n - 1 */
  mpz_sub_ui(t, k->n, 1);

  /* Generate blinding factor. */
  for (;;) {
    /* s = random integer in [1,n-1] */
    mpz_random_int(s, t, &rng);
    mpz_add_ui(s, s, 1);

    /* bi = s^-1 mod n */
    if (!mpz_invert(bi, s, k->n))
      continue;

    /* b = s^e mod n */
    mpz_powm(b, s, k->e, k->n);

    break;
  }

  /* c' = c * b mod n (blind) */
  mpz_mul(m, m, b);
  mpz_mod(m, m, k->n);

  /* m' = c'^d mod n */
  if (mpz_sgn(k->d) > 0 && mpz_odd_p(k->n))
    mpz_powm_sec(m, m, k->d, k->n);
  else
    mpz_powm(m, m, k->d, k->n);

  /* m = m' * bi mod n (unblind) */
  mpz_mul(m, m, bi);
  mpz_mod(m, m, k->n);
  mpz_export_pad(out, m, mpz_bytelen(k->n), 1);

  r = 1;
fail:
  mpz_cleanse(t);
  mpz_cleanse(s);
  mpz_cleanse(b);
  mpz_cleanse(bi);
  mpz_cleanse(m);
  return r;
}

/*
 * Public Key
 */

static void
rsa_pub_init(rsa_pub_t *k) {
  mpz_init(k->n);
  mpz_init(k->e);
}

static void
rsa_pub_clear(rsa_pub_t *k) {
  mpz_cleanse(k->n);
  mpz_cleanse(k->e);
}

static int
rsa_pub_import(rsa_pub_t *k,
               const unsigned char *data,
               size_t len, int strict) {
  if (!asn1_read_seq(&data, &len, strict))
    return 0;

  if (!asn1_read_int(k->n, &data, &len, strict))
    return 0;

  if (!asn1_read_int(k->e, &data, &len, strict))
    return 0;

  if (strict && len != 0)
    return 0;

  return 1;
}

static void
rsa_pub_export(unsigned char *out, size_t *out_len, const rsa_pub_t *k) {
  size_t size = 0;
  size_t pos = 0;

  size += asn1_size_int(k->n);
  size += asn1_size_int(k->e);

  /* 0x30 [size] [body] */
  out[pos++] = 0x30;
  pos += asn1_write_size(out, pos, size);
  pos += asn1_write_int(out, pos, k->n);
  pos += asn1_write_int(out, pos, k->e);

  *out_len = pos;
}

static int
rsa_pub_verify(const rsa_pub_t *k) {
  size_t nb = mpz_bitlen(k->n);
  size_t eb = mpz_bitlen(k->e);

  if (nb < RSA_MIN_MOD_BITS || nb > RSA_MAX_MOD_BITS)
    return 0;

  if (eb < RSA_MIN_EXP_BITS || eb > RSA_MAX_EXP_BITS)
    return 0;

  if (!mpz_odd_p(k->n))
    return 0;

  if (!mpz_odd_p(k->e))
    return 0;

  return 1;
}

static int
rsa_pub_encrypt(const rsa_pub_t *k,
            unsigned char *out,
            const unsigned char *msg,
            size_t msg_len) {
  /* [RFC8017] Page 13, Section 5.1.1.
   *           Page 16, Section 5.2.2.
   */
  mpz_t m;
  int r = 0;

  mpz_init(m);

  if (mpz_sgn(k->n) <= 0 || mpz_sgn(k->e) <= 0)
    goto fail;

  mpz_import(m, msg_len, 1, 1, 0, 0, msg);

  if (mpz_cmp(m, k->n) >= 0)
    goto fail;

  /* c = m^e mod n */
  mpz_powm(m, m, k->e, k->n);
  mpz_export_pad(out, m, mpz_bytelen(k->n), 1);

  r = 1;
fail:
  mpz_cleanse(m);
  return r;
}

/*
 * Digest Info
 */

static int
get_digest_info(const unsigned char **data, size_t *len, int type) {
  const unsigned char *info;

  if (type == -1) {
    *data = NULL;
    *len = 0;
    return 1;
  }

  if (type < 0 || type > HASH_MAX)
    return 0;

  info = digest_info[type];

  *data = &info[1];
  *len = info[0];

  return 1;
}

/*
 * MGF1
 */

static void
mgf1xor(int type,
        unsigned char *out,
        size_t out_len,
        const unsigned char *seed,
        size_t seed_len) {
  /* [RFC8017] Page 67, Section B.2.1. */
  size_t hash_size = hash_output_size(type);
  unsigned char ctr[4] = {0, 0, 0, 0};
  unsigned char digest[HASH_MAX_OUTPUT_SIZE];
  hash_t ctx, hash;
  size_t i = 0;
  size_t j;
  int k;

  hash_init(&ctx, type);
  hash_update(&ctx, seed, seed_len);

  while (i < out_len) {
    memcpy(&hash, &ctx, sizeof(hash_t));
    hash_update(&hash, ctr, sizeof(ctr));
    hash_final(&hash, digest, hash_size);

    j = 0;

    while (i < out_len && j < hash_size)
      out[i++] ^= digest[j++];

    for (k = 3; k >= 0; k--) {
      ctr[k] += 1;

      if (ctr[k] != 0x00)
        break;
    }
  }

  cleanse(ctr, sizeof(ctr));
  cleanse(digest, sizeof(digest));
  cleanse(&ctx, sizeof(ctx));
  cleanse(&hash, sizeof(hash));
}

/*
 * PSS
 */

static int
pss_encode(unsigned char *out,
           size_t *out_len,
           int type,
           const unsigned char *msg,
           size_t msg_len,
           size_t embits,
           const unsigned char *salt,
           size_t salt_len) {
  /* [RFC8017] Page 42, Section 9.1.1. */
  unsigned char *em = out; /* maybe need to truncate outside? */
  size_t hlen = hash_output_size(type);
  size_t slen = salt_len;
  size_t emlen = (embits + 7) >> 3;
  size_t dlen = emlen - hlen - 1;
  unsigned char mask = 0xff >> (8 * emlen - embits);
  unsigned char h0[HASH_MAX_OUTPUT_SIZE];
  unsigned char *db, *h;
  hash_t hash;

  if (msg_len != hlen)
    return 0;

  if (emlen < hlen + slen + 2)
    return 0;

  /* EM = (PS || 0x01 || salt) || H || 0xbc */
  db = &em[0];
  h = &em[emlen - hlen - 1];

  hash_init(&hash, type);
  hash_update(&hash, pss_prefix, sizeof(pss_prefix));
  hash_update(&hash, msg, msg_len);
  hash_update(&hash, salt, salt_len);
  hash_final(&hash, h0, hlen);

  memset(db, 0x00, emlen - slen - hlen - 2);
  db[emlen - slen - hlen - 2] = 0x01;
  memcpy(db + emlen - slen - hlen - 1, salt, salt_len);
  memcpy(h, h0, hlen);
  em[emlen - 1] = 0xbc;

  mgf1xor(type, db, dlen, h, hlen);

  db[0] &= mask;

  *out_len = emlen;

  return 1;
}

static int
pss_verify(int type,
           const unsigned char *msg,
           size_t msg_len,
           unsigned char *em,
           size_t embits,
           size_t slen) {
  /* [RFC8017] Page 44, Section 9.1.2. */
  size_t hlen = hash_output_size(type);
  size_t emlen = (embits + 7) >> 3;
  size_t dlen = emlen - hlen - 1;
  unsigned char mask = 0xff >> (8 * emlen - embits);
  unsigned char h0[HASH_MAX_OUTPUT_SIZE];
  unsigned char *db, *h, *salt;
  size_t i;
  hash_t hash;

  if (msg_len != hlen)
    return 0;

  if (emlen < hlen + slen + 2)
    return 0;

  if (em[emlen - 1] != 0xbc)
    return 0;

  /* EM = (PS || 0x01 || salt) || H || 0xbc */
  db = &em[0];
  h = &em[emlen - hlen - 1];

  if (em[0] & ~mask)
    return 0;

  mgf1xor(type, db, dlen, h, hlen);

  db[0] &= mask;

  if (slen == 0) { /* Auto */
    slen = ~((size_t)0);

    for (i = 0; i < dlen; i++) {
      if (db[i] == 0x00)
        continue;

      if (db[i] == 0x01) {
        slen = dlen - (i + 1);
        break;
      }

      return 0;
    }

    if (slen == ~((size_t)0))
      return 0;
  } else {
    size_t len = dlen - slen - 1;

    for (i = 0; i < len; i++) {
      if (db[i] != 0x00)
        return 0;
    }

    if (db[len] != 0x01)
      return 0;
  }

  salt = &db[dlen - slen];

  hash_init(&hash, type);
  hash_update(&hash, pss_prefix, sizeof(pss_prefix));
  hash_update(&hash, msg, msg_len);
  hash_update(&hash, salt, slen);
  hash_final(&hash, h0, hlen);

  return safe_cmp(h0, h, hlen);
}

/*
 * RSA
 */

int
rsa_privkey_generate(unsigned char *out,
                     size_t *out_len,
                     unsigned long bits,
                     unsigned long long exp,
                     const unsigned char *entropy) {
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  if (!rsa_priv_generate(&k, bits, exp, entropy))
    goto fail;

  rsa_priv_export(out, out_len, &k);
  r = 1;
fail:
  rsa_priv_clear(&k);
  return r;
}

size_t
rsa_privkey_size(const unsigned char *key, size_t key_len) {
  rsa_priv_t k;
  size_t r = 0;

  rsa_priv_init(&k);

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  r = mpz_bitlen(k.n);
fail:
  rsa_priv_clear(&k);
  return r;
}

int
rsa_privkey_verify(const unsigned char *key, size_t key_len) {
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  r = 1;
fail:
  rsa_priv_clear(&k);
  return r;
}

int
rsa_privkey_recover(unsigned char *out,
                    size_t *out_len,
                    const unsigned char *key,
                    size_t key_len,
                    const unsigned char *entropy) {
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k)) {
    if (!rsa_priv_recover(&k, entropy))
      goto fail;
  }

  rsa_priv_export(out, out_len, &k);
  r = 1;
fail:
  rsa_priv_clear(&k);
  return r;
}

int
rsa_privkey_normalize(unsigned char *out,
                      size_t *out_len,
                      const unsigned char *key,
                      size_t key_len) {
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  /* Parse with the most lax rules possible. */
  if (!rsa_priv_import(&k, key, key_len, 0))
    goto fail;

  rsa_priv_export(out, out_len, &k);
  r = 1;
fail:
  rsa_priv_clear(&k);
  return r;
}

int
rsa_pubkey_create(unsigned char *out,
                  size_t *out_len,
                  const unsigned char *key,
                  size_t key_len,
                  const unsigned char *entropy) {
  rsa_priv_t k;
  rsa_pub_t k;
  int r = 0;

  rsa_priv_init(&k);
  rsa_pub_init(&p);

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  mpz_set(p->n, k->n);
  mpz_set(p->e, k->e);

  rsa_pub_export(out, out_len, &p);
  r = 1;
fail:
  rsa_priv_clear(&k);
  rsa_pub_clear(&p);
  return r;
}

size_t
rsa_pubkey_size(const unsigned char *key, size_t key_len) {
  rsa_pub_t k;
  size_t r = 0;

  rsa_pub_init(&k);

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  r = mpz_bitlen(k.n);
fail:
  rsa_pub_clear(&k);
  return r;
}

int
rsa_pubkey_verify(const unsigned char *key, size_t key_len) {
  rsa_pub_t k;
  int r = 0;

  rsa_pub_init(&k);

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  r = 1;
fail:
  rsa_pub_clear(&k);
  return r;
}

int
rsa_pubkey_normalize(unsigned char *out,
                     size_t *out_len,
                     const unsigned char *key,
                     size_t key_len) {
  rsa_pub_t k;
  int r = 0;

  rsa_pub_init(&k);

  /* Parse with the most lax rules possible. */
  if (!rsa_pub_import(&k, key, key_len, 0))
    goto fail;

  rsa_pub_export(out, out_len, &k);
  r = 1;
fail:
  rsa_pub_clear(&k);
  return r;
}

int
rsa_sign(unsigned char *out,
         size_t *out_len,
         int type,
         const unsigned char *msg,
         size_t msg_len,
         const unsigned char *key,
         size_t key_len,
         const unsigned char *entropy) {
  /* [RFC8017] Page 36, Section 8.2.1.
   *           Page 45, Section 9.2.
   */
  size_t hlen = hash_output_size(type);
  size_t i, prefix_len, tlen, klen;
  const unsigned char *prefix;
  unsigned char *em = out;
  int r = 0;
  rsa_priv_t k;

  rsa_priv_init(&k);

  if (!get_digest_info(&prefix, &prefix_len, type))
    goto fail;

  if (type == -1)
    hlen = msg_len;

  if (msg_len != hlen)
    goto fail;

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  tlen = prefix_len + hlen;
  klen = mpz_bytelen(k.n);

  if (klen < tlen + 11)
    goto fail;

  /* EM = 0x00 || 0x01 || PS || 0x00 || T */
  em[0] = 0x00;
  em[1] = 0x01;

  for (i = 2; i < klen - tlen - 1; i++)
    em[i] = 0xff;

  em[klen - tlen - 1] = 0x00;

  memcpy(em + klen - tlen, prefix, prefix_len);
  memcpy(em + klen - hlen, msg, msg_len);

  if (!rsa_priv_decrypt(&k, out, em, klen, entropy))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  rsa_priv_clear(&k);
  return r;
}

int
rsa_verify(int type,
           const unsigned char *msg,
           size_t msg_len,
           const unsigned char *sig,
           size_t sig_len,
           const unsigned char *key,
           size_t key_len,
           int strict) {
  /* [RFC8017] Page 37, Section 8.2.2.
   *           Page 45, Section 9.2.
   */
  size_t hlen = hash_output_size(type);
  size_t klen = 0;
  size_t i, prefix_len, tlen;
  const unsigned char *prefix;
  unsigned char *em = NULL;
  uint32_t ok;
  int r = 0;
  rsa_pub_t k;

  rsa_pub_init(&k);

  if (!get_digest_info(&prefix, &prefix_len, type))
    goto fail;

  if (type == -1)
    hlen = msg_len;

  if (msg_len != hlen)
    goto fail;

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  tlen = prefix_len + hlen;
  klen = mpz_bytelen(k.n);

  if (strict && sig_len != klen)
    goto fail;

  if (klen < tlen + 11)
    goto fail;

  em = safe_malloc(klen);

  if (em == NULL)
    goto fail;

  if (!rsa_pub_encrypt(&k, em, sig, sig_len))
    goto fail;

  /* EM = 0x00 || 0x01 || PS || 0x00 || T */
  ok = 1;

  ok &= safe_equal(em[0], 0x00);
  ok &= safe_equal(em[1], 0x01);

  for (i = 2; i < klen - tlen - 1; i++)
    ok &= safe_equal(em[i], 0xff);

  ok &= safe_equal(em[klen - tlen - 1], 0x00);
  ok &= safe_cmp(em + klen - tlen, prefix, prefix_len);
  ok &= safe_cmp(em + klen - hlen, msg, msg_len);

  r = (ok == 1);
fail:
  rsa_pub_clear(&k);
  safe_free(em, klen);
  return r;
}

int
rsa_encrypt(unsigned char *out,
            size_t *out_len,
            const unsigned char *msg,
            size_t msg_len,
            const unsigned char *key,
            size_t key_len,
            const unsigned char *entropy) {
  /* [RFC8017] Page 28, Section 7.2.1. */
  unsigned char *em = out;
  size_t i, mlen, plen;
  size_t klen = 0;
  drbg_t rng;
  rsa_pub_t k;
  int r = 0;

  drbg_init(&rng, HASH_SHA256, entropy, 32);

  rsa_pub_init(&k);

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  if (klen < 11 || msg_len > klen - 11)
    goto fail;

  /* EM = 0x00 || 0x02 || PS || 0x00 || M */
  mlen = msg_len;
  plen = klen - mlen - 3;

  em[0] = 0x00;
  em[1] = 0x02;

  drbg_generate(&rng, em + 2, plen);

  for (i = 2; i < 2 + plen; i++) {
    while (em[i] == 0x00)
      drbg_generate(&rng, em + i, 1);
  }

  em[klen - mlen - 1] = 0x00;

  memcpy(em + klen - mlen, msg, msg_len);

  if (!rsa_pub_encrypt(&k, out, em, klen))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  rsa_pub_clear(&k);
  cleanse(&rng, sizeof(rng));
  if (r == 0) cleanse(out, klen);
  return r;
}

int
rsa_decrypt(unsigned char *out,
            size_t *out_len,
            const unsigned char *msg,
            size_t msg_len,
            const unsigned char *key,
            size_t key_len,
            int strict,
            const unsigned char *entropy) {
  /* [RFC8017] Page 29, Section 7.2.2. */
  unsigned char *em = out;
  uint32_t i, zero, two, index, looking;
  uint32_t equals0, validps, valid, offset;
  size_t klen = 0;
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  if (strict && msg_len != klen)
    goto fail;

  if (klen < 11)
    goto fail;

  if (!rsa_priv_decrypt(&k, em, msg, msg_len, entropy))
    goto fail;

  /* EM = 0x00 || 0x02 || PS || 0x00 || M */
  zero = safe_equal(em[0], 0x00);
  two = safe_equal(em[1], 0x02);
  index = 0;
  looking = 1;

  for (i = 2; i < klen; i++) {
    equals0 = safe_equal(em[i], 0x00);
    index = safe_select(index, i, looking & equals0);
    looking = safe_select(looking, 0, equals0);
  }

  validps = safe_lte(2 + 8, index);
  valid = zero & two & (looking ^ 1) & validps;
  offset = safe_select(0, index + 1, valid);

  if (valid == 0)
    goto fail;

  *out_len = klen - offset;
  memmove(out, em + offset, *out_len);

  r = 1;
fail:
  rsa_priv_clear(&k);
  if (r == 0) cleanse(out, klen);
  return r;
}

int
rsa_encrypt_oaep(unsigned char *out,
                 size_t *out_len,
                 int type,
                 const unsigned char *msg,
                 size_t msg_len,
                 const unsigned char *key,
                 size_t key_len,
                 const unsigned char *label,
                 size_t label_len,
                 const unsigned char *entropy) {
  /* [RFC8017] Page 22, Section 7.1.1. */
  unsigned char lhash[HASH_MAX_OUTPUT_SIZE];
  unsigned char *em = out;
  unsigned char *seed, *db;
  size_t hlen = hash_output_size(type);
  size_t klen = 0;
  size_t mlen = msg_len;
  size_t slen, dlen;
  hash_t hash;
  drbg_t rng;
  int r = 0;
  rsa_pub_t k;

  rsa_pub_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  if (klen < 2 * hlen + 2)
    goto fail;

  if (msg_len > klen - 2 * hlen - 2)
    goto fail;

  hash_init(&hash, type);
  hash_update(&hash, label, label_len);
  hash_final(&hash, lhash, hlen);

  /* EM = 0x00 || (seed) || (Hash(L) || PS || 0x01 || M) */
  seed = &em[1];
  slen = hlen;
  db = &em[1 + hlen];
  dlen = klen - (1 + hlen);

  em[0] = 0x00;

  drbg_init(&rng, HASH_SHA256, entropy, 32);
  drbg_generate(&rng, seed, slen);

  memcpy(&db[0], lhash, hlen);
  memset(&db[hlen], 0x00, (dlen - mlen - 1) - hlen);

  db[dlen - mlen - 1] = 0x01;
  memcpy(&db[dlen - mlen], msg, mlen);

  mgf1xor(type, db, dlen, seed, slen);
  mgf1xor(type, seed, slen, db, dlen);

  if (!rsa_pub_encrypt(&k, out, em, klen))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  rsa_pub_clear(&k);
  cleanse(&rng, sizeof(drbg_t));
  cleanse(&hash, sizeof(hash_t));
  if (r == 0) cleanse(out, klen);
  return r;
}

static int
rsa_decrypt_oaep(unsigned char *out,
                 size_t *out_len,
                 int type,
                 const unsigned char *msg,
                 size_t msg_len,
                 const unsigned char *key,
                 size_t key_len,
                 const unsigned char *label,
                 size_t label_len,
                 int strict,
                 const unsigned char *entropy) {
  /* [RFC8017] Page 25, Section 7.1.2. */
  unsigned char *em = out;
  unsigned char *seed, *db, *rest, *lhash;
  size_t i, slen, dlen, rlen;
  size_t hlen = hash_output_size(type);
  size_t klen = 0;
  uint32_t zero, lvalid, looking, index;
  uint32_t invalid, valid, equals0, equals1;
  unsigned char expect[HASH_MAX_OUTPUT_SIZE];
  hash_t hash;
  int r = 0;
  rsa_priv_t k;

  rsa_priv_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  if (strict && msg_len != klen)
    goto fail;

  if (klen < hlen * 2 + 2)
    goto fail;

  if (!rsa_priv_decrypt(&k, em, msg, msg_len, entropy))
    goto fail;

  hash_init(&hash, type);
  hash_update(&hash, label, label_len);
  hash_final(&hash, expect, hlen);

  /* EM = 0x00 || (seed) || (Hash(L) || PS || 0x01 || M) */
  zero = safe_equal(em[0], 0x00);
  seed = &em[1];
  slen = hlen;
  db = &em[hlen + 1];
  dlen = klen - (hlen + 1);

  mgf1xor(type, seed, slen, db, dlen);
  mgf1xor(type, db, dlen, seed, slen);

  lhash = &db[0];
  lvalid = safe_cmp(lhash, expect, hlen);
  rest = &db[hlen];
  rlen = dlen - hlen;

  looking = 1;
  index = 0;
  invalid = 0;

  for (i = 0; i < rlen; i++) {
    equals0 = safe_equal(rest[i], 0x00);
    equals1 = safe_equal(rest[i], 0x01);
    index = safe_select(index, i, looking & equals1);
    looking = safe_select(looking, 0, equals1);
    invalid = safe_select(invalid, 1, looking & (equals0 ^ 1));
  }

  valid = zero & lvalid & (invalid ^ 1) & (looking ^ 1);

  if (valid == 0)
    goto fail;

  *out_len = rlen - (index + 1);
  memmove(out, rest + index + 1, *out_len);

  r = 1;
fail:
  rsa_priv_clear(&k);
  cleanse(&hash, sizeof(hash));
  if (r == 0) cleanse(out, klen);
  return r;
}

int
rsa_sign_pss(unsigned char *out,
             size_t *out_len,
             int type,
             const unsigned char *msg,
             size_t msg_len,
             const unsigned char *key,
             size_t key_len,
             int salt_len,
             const unsigned char *entropy) {
  /* [RFC8017] Page 33, Section 8.1.1. */
  size_t hlen = hash_output_size(type);
  unsigned char *salt = NULL;
  unsigned char *em = out;
  size_t emlen, bits;
  size_t klen = 0;
  drbg_t rng;
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (msg_len != hlen)
    goto fail;

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  bits = mpz_bitlen(k.n);
  klen = (bits + 7) / 8;

  if (salt_len == RSA_SALT_LENGTH_AUTO)
    salt_len = klen - 2 - hlen;
  else if (salt_len == RSA_SALT_LENGTH_HASH)
    salt_len = hlen;

  if (salt_len < 0)
    goto fail;

  if (salt_len > 0) {
    salt = safe_malloc(salt_len);

    if (salt == NULL)
      goto fail;

    drbg_init(&rng, HASH_SHA512, entropy, 32);
    drbg_generate(&rng, salt, salt_len);
  }

  if (!pss_encode(em, &emlen, type, msg, msg_len, bits - 1, salt, salt_len))
    goto fail;

  /* Note that `em` may be one byte less
   * than the modulus size in the case
   * of (bits - 1) mod 8 == 0.
   */
  if (!rsa_priv_decrypt(&k, out, em, emlen, entropy))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  rsa_priv_clear(&k);
  cleanse(&rng, sizeof(rng));
  safe_free(salt, salt_len < 0 ? 0 : salt_len);
  if (r == 0) cleanse(out, klen);
  return r;
}

int
rsa_verify_pss(int type,
               const unsigned char *msg,
               size_t msg_len,
               const unsigned char *sig,
               size_t sig_len,
               const unsigned char *key,
               size_t key_len,
               int salt_len,
               int strict) {
  /* [RFC8017] Page 34, Section 8.1.2. */
  unsigned char *em = NULL;
  size_t hlen = hash_output_size(type);
  size_t klen = 0;
  size_t bits;
  rsa_pub_t k;
  int r = 0;

  rsa_pub_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  bits = mpz_bitlen(k.n);
  klen = (bits + 7) / 8;

  if (msg_len != hlen)
    goto fail;

  if (strict && sig_len != klen)
    goto fail;

  if (salt_len == RSA_SALT_LENGTH_AUTO)
    salt_len = 0; /* Handled in pssVerify. */
  else if (salt_len == RSA_SALT_LENGTH_HASH)
    salt_len = hlen;

  if (salt_len < 0)
    goto fail;

  em = safe_malloc(klen);

  if (em == NULL)
    goto fail;

  if (!rsa_pub_encrypt(&k, em, sig, sig_len))
    goto fail;

  /* Edge case: the encoding crossed a
   * a byte boundary. Our encryption
   * function pads to the modulus size
   * by default, meaning there's one
   * extra zero byte prepended.
   */
  if (((bits - 1) & 7) == 0) {
    if (em[0] != 0x00)
      goto fail;

    em += 1;
  }

  if (!pss_verify(type, msg, msg_len, em, bits - 1, salt_len))
    goto fail;

  r = 1;
fail:
  rsa_pub_clear(&k);
  safe_free(em, klen);
  return r;
}

int
rsa_encrypt_raw(unsigned char *out,
                size_t *out_len,
                const unsigned char *msg,
                size_t msg_len,
                const unsigned char *key,
                size_t key_len,
                int strict) {
  /* [RFC8017] Page 13, Section 5.1.1.
   *           Page 16, Section 5.2.2.
   */
  size_t klen = 0;
  rsa_pub_t k;
  int r = 0;

  rsa_pub_init(&k);

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  /* OpenSSL behavior for public encryption. */
  /* strict flag is new */
  if (strict && msg_len !== klen)
    goto fail;

  if (!rsa_pub_encrypt(&k, out, msg, msg_len))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  rsa_pub_clear(&k);
  return r;
}

int
rsa_decrypt_raw(unsigned char *out,
                size_t *out_len,
                const unsigned char *msg,
                size_t msg_len,
                const unsigned char *key,
                size_t key_len,
                int strict,
                const unsigned char *entropy) {
  /* [RFC8017] Page 13, Section 5.1.2.
   *           Page 15, Section 5.2.1.
   */
  size_t klen = 0;
  rsa_priv_t k;
  int r = 0;

  rsa_priv_init(&k);

  if (!rsa_priv_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_priv_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  /* new */
  if (strict && msg_len !== klen)
    goto fail;

  if (!rsa_priv_decrypt(&k, out, msg, msg_len, entropy))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  rsa_priv_clear(&k);
  return r;
}

int
rsa_veil(unsigned char *out,
         size_t *out_len,
         const unsigned char *msg,
         size_t msg_len,
         size_t bits,
         const unsigned char *key,
         size_t key_len,
         int strict,
         const unsigned char *entropy) {
  size_t mlen = (bits + 7) / 8;
  mpz_t vmax, rmax, c, v, r;
  rsa_pub_t k;
  drbg_t rng;
  int ret = 0;

  mpz_init(vmax);
  mpz_init(rmax);
  mpz_init(c);
  mpz_init(v);
  mpz_init(r);

  rsa_pub_init(&k);

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  /* new */
  if (strict && msg_len != mpz_bytelen(k.n))
    goto fail;

  if (bits < mpz_bitlen(k.n))
    goto fail;

  mpz_import(c, msg_len, 1, 1, 0, 0, msg);

  if (mpz_cmp(c, k.n) >= 0)
    goto fail;

  /* vmax = 1 << bits */
  mpz_set_ui(vmax, 1);
  mpz_mul_2exp(vmax, vmax, bits);

  /* rmax = (vmax - c + n - 1) / n */
  mpz_sub(rmax, vmax, c);
  mpz_add(rmax, rmax, k.n);
  mpz_sub_ui(rmax, rmax, 1);
  mpz_fdiv_q(rmax, rmax, k.n);

  assert(mpz_sgn(rmax) > 0);

  mpz_set(v, vmax);

  drbg_init(&rng, HASH_SHA256, entropy, 32);

  while (mpz_cmp(v, vmax) >= 0) {
    mpz_random_int(r, rmax, &rng);

    /* v = c + r * n */
    mpz_mul(r, r, k.n);
    mpz_add(v, c, r);
  }

  mpz_mod(r, v, k.n);

  assert(mpz_cmp(r, c) == 0);
  assert(mpz_bitlen(v) <= bits);

  mpz_export_pad(out, v, mlen, 1);

  *out_len = mlen;
  ret = 1;
fail:
  cleanse(&rng, sizeof(rng));
  mpz_cleanse(vmax);
  mpz_cleanse(rmax);
  mpz_cleanse(c);
  mpz_cleanse(v);
  mpz_cleanse(r);
  rsa_pub_clear(&k);
  return ret;
}

int
rsa_unveil(unsigned char *out,
           size_t *out_len,
           const unsigned char *msg,
           size_t msg_len,
           size_t bits,
           const unsigned char *key,
           size_t key_len,
           int strict) {
  size_t klen = 0;
  rsa_pub_t k;
  int r = 0;
  mpz_t v;

  mpz_init(v);
  rsa_pub_init(&k);

  if (!rsa_pub_import(&k, key, key_len, 1))
    goto fail;

  if (!rsa_pub_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  /* strict flag new */
  if (strict && msg_len < klen)
    goto fail;

  mpz_import(v, msg_len, 1, 1, 0, 0, msg);

  if (mpz_bitlen(v) > bits)
    goto fail;

  mpz_mod(v, v, k.n);
  mpz_export_pad(out, v, mlen, 1);

  *out_len = klen;
  r = 1;
fail:
  mpz_cleanse(v);
  rsa_pub_clear(&k);
  return r;
}
