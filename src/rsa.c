#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <torsion/hash.h>
#include <torsion/drbg.h>
#include <torsion/rsa.h>

#define RSA_DEFAULT_MOD_BITS 2048
#define RSA_DEFAULT_EXP 65537
#define RSA_MIN_MOD_BITS 512
#define RSA_MAX_MOD_BITS 16384
#define RSA_MIN_MOD_BYTES ((RSA_MIN_MOD_BITS + 7) / 8)
#define RSA_MAX_MOD_BYTES ((RSA_MAX_MOD_BITS + 7) / 8)
#define RSA_MIN_EXP 3ull
#define RSA_MAX_EXP 0x1ffffffffull
#define RSA_MIN_EXP_BITS 2
#define RSA_MAX_EXP_BITS 33
#define RSA_MIN_EXP_BYTES 1
#define RSA_MAX_EXP_BYTES 5

/* Limits:
 * 4096 = 2614
 * 8192 = 5174
 * 16384 = 10294
 */

#define RSA_MAX_PRIV_SIZE (0                   \
  + 4 /* seq */                                \
  + 3 /* version */                            \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* n */          \
  + 2 + 1 + RSA_MAX_EXP_BYTES /* e */          \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* d */          \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* p */  \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* q */  \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* dp */ \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* dq */ \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* qi */         \
)

/* Limits:
 * 4096 = 529
 * 8192 = 1041
 * 16384 = 2065
 */

#define RSA_MAX_PUB_SIZE (0           \
  + 4 /* seq */                       \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* n */ \
  + 2 + 1 + RSA_MAX_EXP_BYTES /* e */ \
)

#define RSA_MAX_SIG_SIZE RSA_MAX_MOD_BYTES

#define RSA_SALT_LENGTH_AUTO 0
#define RSA_SALT_LENGTH_HASH -1

static const unsigned char pss_prefix[8] = {0, 0, 0, 0, 0, 0, 0, 0};

/*
 * Structs
 */

typedef struct _pub_s {
  mpz_t n;
  mpz_t e;
} pub_t;

typedef struct _priv_s {
  mpz_t n;
  mpz_t e;
  mpz_t d;
  mpz_t p;
  mpz_t q;
  mpz_t dp;
  mpz_t dq;
  mpz_t qi;
} priv_t;

/*
 * Digest Info
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
 * Helpers
 */

static void
cleanse(void *ptr, size_t len) {
#if defined(_WIN32)
  /* https://github.com/jedisct1/libsodium/blob/3b26a5c/src/libsodium/sodium/utils.c#L112 */
  SecureZeroMemory(ptr, len);
#elif defined(__GNUC__)
  /* https://github.com/torvalds/linux/blob/37d4e84/include/linux/string.h#L233 */
  /* https://github.com/torvalds/linux/blob/37d4e84/include/linux/compiler-gcc.h#L21 */
  memset(ptr, 0, len);
  __asm__ __volatile__("": :"r"(ptr) :"memory");
#else
  /* http://www.daemonology.net/blog/2014-09-04-how-to-zero-a-buffer.html */
  static void *(*const volatile memset_ptr)(void *, int, size_t) = memset;
  (memset_ptr)(ptr, 0, len);
#endif
}

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
 * MPZ
 */

/* Alias to avoid future collisions with gmp. */
#define mpz_bitlen(n) (mpz_sgn(n) == 0 ? 0 : mpz_sizeinbase(n, 2))
#define mpz_bytelen(n) ((mpz_bitlen(n) + 7) / 8)
#define mpz_export_pad torsion_mpz_export_pad
#define mpz_zerobits torsion_mpz_zerobits
#define mpz_cleanse torsion_mpz_cleanse
#define mpz_random_bits torsion_mpz_random_bits
#define mpz_random_int torsion_mpz_random_int
#define mpz_is_prime_mr torsion_mpz_is_prime_mr
#define mpz_is_prime_lucas torsion_mpz_is_prime_lucas
#define mpz_is_prime torsion_mpz_is_prime
#define mpz_random_prime torsion_mpz_random_prime

static void
mpz_export_pad(unsigned char *out, const mpz_t n, size_t size) {
  size_t len = mpz_bytelen(n);
  size_t pos = size - len;

  assert(len <= size);

  memset(out, 0x00, pos);

  mpz_export(out + pos, NULL, 1, 1, 0, 0, n);
}

static size_t
mpz_zerobits(const mpz_t n) {
  /* Note: mpz_ptr is undocumented. */
  /* https://gmplib.org/list-archives/gmp-discuss/2009-May/003769.html */
  /* https://gmplib.org/list-archives/gmp-devel/2013-February/002775.html */
  int sgn = mpz_sgn(n);
  unsigned long bits;

  if (sgn == 0)
    return 0;

  if (sgn < 0)
    mpz_neg((mpz_ptr)n, n);

  bits = mpz_scan1(n, 0);

  if (sgn < 0)
    mpz_neg((mpz_ptr)n, n);

  return bits;
}

#ifndef TORSION_HAS_GMP
/* `mpz_jacobi` is not implemented in mini-gmp. */
/* https://github.com/golang/go/blob/aadaec5/src/math/big/int.go#L754 */
static int
mpz_jacobi(const mpz_t x, const mpz_t y) {
  mpz_t a, b, c;
  unsigned long s, bmod8;
  int j;

  /* Undefined behavior. */
  /* if y == 0 or y mod 2 == 0 */
  if (mpz_sgn(y) == 0 || mpz_even_p(y))
    return 0;

  mpz_init(a);
  mpz_init(b);
  mpz_init(c);

  /* a = x */
  mpz_set(a, x);

  /* b = y */
  mpz_set(b, y);

  j = 1;

  /* if b < 0 */
  if (mpz_sgn(b) < 0) {
    /* if a < 0 */
    if (mpz_sgn(a) < 0)
      j = -1;

    /* b = -b */
    mpz_neg(b, b);
  }

  for (;;) {
    /* if b == 1 */
    if (mpz_cmp_ui(b, 1) == 0)
      break;

    /* if a == 0 */
    if (mpz_sgn(a) == 0) {
      j = 0;
      break;
    }

    /* a = a mod b */
    mpz_mod(a, a, b);

    /* if a == 0 */
    if (mpz_sgn(a) == 0) {
      j = 0;
      break;
    }

    /* s = a factors of 2 */
    s = mpz_zerobits(a);

    if (s & 1) {
      /* bmod8 = b mod 8 */
      bmod8 = mpz_getlimbn(b, 0) & 7;

      if (bmod8 == 3 || bmod8 == 5)
        j = -j;
    }

    /* c = a >> s */
    mpz_tdiv_q_2exp(c, a, s);

    /* if b mod 4 == 3 and c mod 4 == 3 */
    if ((mpz_getlimbn(b, 0) & 3) == 3 && (mpz_getlimbn(c, 0) & 3) == 3)
      j = -j;

    /* a = b */
    mpz_set(a, b);

    /* b = c */
    mpz_set(b, c);
  }

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);

  return j;
}
#endif

static void
mpz_cleanse(mpz_t n) {
#ifdef TORSION_HAS_GMP
  /* Using the public API. */
  const mp_limb_t *orig = mpz_limbs_read(n);
  size_t size = mpz_size(n);
  mp_limb_t *limbs = mpz_limbs_modify(n, (mp_size_t)size);

  /* Zero the limbs. */
  cleanse(limbs, size * sizeof(mp_limb_t));

  /* Ensure the integer remains in a valid state. */
  mpz_limbs_finish(n, 0);

  /* Sanity checks. */
  assert(limbs == orig);
  assert(mpz_limbs_read(n) == orig);
  assert(mpz_sgn(n) == 0);
#else
  /* Using the internal API. */
  mpz_ptr x = n;
  cleanse(x->_mp_d, x->_mp_alloc * sizeof(mp_limb_t));
  x->_mp_size = 0;
#endif
  mpz_clear(n);
}

static void
mpz_random_bits(mpz_t ret, size_t bits, drbg_t *rng) {
  size_t num_bits = sizeof(mp_limb_t) * CHAR_BIT;
  size_t size = (bits + num_bits - 1) / num_bits;
  size_t low = num_bits - (size * num_bits - bits);
  mp_limb_t *limbs = mpz_limbs_write(ret, size);

  drbg_generate(rng, limbs, size * sizeof(mp_limb_t));

  if (low != num_bits)
    limbs[size - 1] &= (((mp_limb_t)1 << low) - 1);

  mpz_limbs_finish(ret, size);
}

static void
mpz_random_int(mpz_t ret, const mpz_t max, drbg_t *rng) {
  size_t bits;

  if (mpz_sgn(max) <= 0) {
    mpz_set_ui(ret, 0);
    return;
  }

  mpz_set(ret, max);

  bits = mpz_bitlen(ret);

  assert(bits > 0);

  while (mpz_cmp(ret, max) >= 0)
    mpz_random_bits(ret, bits, rng);
}

/* https://github.com/golang/go/blob/aadaec5/src/math/big/prime.go#L81 */
/* https://github.com/indutny/miller-rabin/blob/master/lib/mr.js */
static int
mpz_is_prime_mr(const mpz_t n, unsigned long reps, int force2, drbg_t *rng) {
  int r = 0;
  mpz_t nm1, nm3, q, x, y;
  unsigned long k, i, j;

  /* if n < 7 */
  if (mpz_cmp_ui(n, 7) < 0) {
    /* n == 2 or n == 3 or n == 5 */
    return mpz_cmp_ui(n, 2) == 0
        || mpz_cmp_ui(n, 3) == 0
        || mpz_cmp_ui(n, 5) == 0;
  }

  /* if n mod 2 == 0 */
  if (mpz_even_p(n))
    return 0;

  mpz_init(nm1);
  mpz_init(nm3);
  mpz_init(q);
  mpz_init(x);
  mpz_init(y);

  /* nm1 = n - 1 */
  mpz_sub_ui(nm1, n, 1);

  /* nm3 = nm1 - 2 */
  mpz_sub_ui(nm3, nm1, 2);

  /* k = nm1 factors of 2 */
  k = mpz_zerobits(nm1);

  /* q = nm1 >> k */
  mpz_tdiv_q_2exp(q, nm1, k);

  for (i = 0; i < reps; i++) {
    if (i == reps - 1 && force2) {
      /* x = 2 */
      mpz_set_ui(x, 2);
    } else {
      /* x = random integer in [2,n-1] */
      mpz_random_int(x, nm3, rng);
      mpz_add_ui(x, x, 2);
    }

    /* y = x^q mod n */
    mpz_powm(y, x, q, n);

    /* if y == 1 or y == -1 mod n */
    if (mpz_cmp_ui(y, 1) == 0 || mpz_cmp(y, nm1) == 0)
      continue;

    for (j = 1; j < k; j++) {
      /* y = y^2 mod n */
      mpz_mul(y, y, y);
      mpz_mod(y, y, n);

      /* if y == -1 mod n */
      if (mpz_cmp(y, nm1) == 0)
        goto next;

      /* if y == 1 mod n */
      if (mpz_cmp_ui(y, 1) == 0)
        goto fail;
    }

    goto fail;
next:
    ;
  }

  r = 1;
fail:
  mpz_clear(nm1);
  mpz_clear(nm3);
  mpz_clear(q);
  mpz_clear(x);
  mpz_clear(y);
  return r;
}

/* https://github.com/golang/go/blob/aadaec5/src/math/big/prime.go#L150 */
static int
mpz_is_prime_lucas(const mpz_t n, unsigned long limit) {
  int ret = 0;
  unsigned long p, r;
  mpz_t d, s, nm2, vk, vk1, t1, t2, t3;
  long i, t;
  int j;

  mpz_init(d);
  mpz_init(s);
  mpz_init(nm2);
  mpz_init(vk);
  mpz_init(vk1);
  mpz_init(t1);
  mpz_init(t2);
  mpz_init(t3);

  /* if n <= 1 */
  if (mpz_cmp_ui(n, 1) <= 0)
    goto fail;

  /* if n mod 2 == 0 */
  if (mpz_even_p(n)) {
    /* if n == 2 */
    if (mpz_cmp_ui(n, 2) == 0)
      goto succeed;
    goto fail;
  }

  /* p = 3 */
  p = 3;

  /* d = 1 */
  mpz_set_ui(d, 1);

  for (;;) {
    if (p > 10000) {
      /* Thought to be impossible. */
      goto fail;
    }

    if (limit != 0 && p > limit) {
      /* Enforce a limit to prevent DoS'ing. */
      goto fail;
    }

    /* d = p * p - 4 */
    mpz_set_ui(d, p * p - 4);

    j = mpz_jacobi(d, n);

    /* if d is not square mod n */
    if (j == -1)
      break;

    /* if d == 0 mod n */
    if (j == 0) {
      /* if n == p + 2 */
      if (mpz_cmp_ui(n, p + 2) == 0)
        goto succeed;
      goto fail;
    }

    if (p == 40) {
      /* if floor(n^(1 / 2))^2 == n */
      if (mpz_perfect_square_p(n))
        goto fail;
    }

    p += 1;
  }

  /* s = n + 1 */
  mpz_add_ui(s, n, 1);

  /* r = s factors of 2 */
  r = mpz_zerobits(s);

  /* nm2 = n - 2 */
  mpz_sub_ui(nm2, n, 2);

  /* vk = 2 */
  mpz_set_ui(vk, 2);

  /* vk1 = p */
  mpz_set_ui(vk1, p);

  /* s >>= r */
  mpz_tdiv_q_2exp(s, s, r);

  for (i = (long)mpz_bitlen(s); i >= 0; i--) {
    /* if floor(s / 2^i) mod 2 == 1 */
    if (mpz_tstbit(s, i)) {
      /* vk = (vk * vk1 + n - p) mod n */
      /* vk1 = (vk1^2 + nm2) mod n */
      mpz_mul(t1, vk, vk1);
      mpz_add(t1, t1, n);
      mpz_sub_ui(t1, t1, p);
      mpz_mod(vk, t1, n);
      mpz_mul(t1, vk1, vk1);
      mpz_add(t1, t1, nm2);
      mpz_mod(vk1, t1, n);
    } else {
      /* vk1 = (vk * vk1 + n - p) mod n */
      /* vk = (vk^2 + nm2) mod n */
      mpz_mul(t1, vk, vk1);
      mpz_add(t1, t1, n);
      mpz_sub_ui(t1, t1, p);
      mpz_mod(vk1, t1, n);
      mpz_mul(t1, vk, vk);
      mpz_add(t1, t1, nm2);
      mpz_mod(vk, t1, n);
    }
  }

  /* if vk == 2 or vk == nm2 */
  if (mpz_cmp_ui(vk, 2) == 0 || mpz_cmp(vk, nm2) == 0) {
    /* t3 = abs(vk * p - vk1 * 2) mod n */
    mpz_mul_ui(t1, vk, p);
    mpz_mul_2exp(t2, vk1, 1);

    if (mpz_cmp(t1, t2) < 0)
      mpz_swap(t1, t2);

    mpz_sub(t1, t1, t2);
    mpz_mod(t3, t1, n);

    /* if t3 == 0 */
    if (mpz_sgn(t3) == 0)
      goto succeed;
  }

  for (t = 0; t < (long)r - 1; t++) {
    /* if vk == 0 */
    if (mpz_sgn(vk) == 0)
      goto succeed;

    /* if vk == 2 */
    if (mpz_cmp_ui(vk, 2) == 0)
      goto fail;

    /* vk = (vk^2 - 2) mod n */
    mpz_mul(t1, vk, vk);
    mpz_sub_ui(t1, t1, 2);
    mpz_mod(vk, t1, n);
  }

  goto fail;
succeed:
  ret = 1;
fail:
  mpz_clear(d);
  mpz_clear(s);
  mpz_clear(nm2);
  mpz_clear(vk);
  mpz_clear(vk1);
  mpz_clear(t1);
  mpz_clear(t2);
  mpz_clear(t3);
  return ret;
}

static int
mpz_is_prime(const mpz_t p, drbg_t *rng) {
  if (!mpz_is_prime_mr(p, 20 + 1, 1, rng))
    return 0;

  if (!mpz_is_prime_lucas(p, 0))
    return 0;

  return 1;
}

static void
mpz_random_prime(mpz_t ret, size_t bits, drbg_t *rng) {
  assert(bits > 1);

  do {
    mpz_random_bits(ret, bits, rng);
    mpz_setbit(ret, 0);
    mpz_setbit(ret, bits - 1);
  } while (!mpz_is_prime(ret, rng));
}

/*
 * ASN1
 */

static int
asn1_read_size(size_t *size,
               const unsigned char **data,
               size_t *len,
               int strict) {
  unsigned char ch;

  assert(sizeof(size_t) * CHAR_BIT >= 32);

  if (*len == 0)
    return 0;

  ch = **data;

  *data += 1;
  *len -= 1;

  if ((ch & 0x80) == 0) {
    /* Short form. */
    *size = ch;
  } else {
    size_t bytes = ch & 0x7f;
    size_t i;

    /* Indefinite form. */
    if (strict && bytes == 0)
      return 0;

    /* Long form. */
    *size = 0;

    for (i = 0; i < bytes; i++) {
      if (*len == 0)
        return 0;

      ch = **data;
      *data += 1;
      *len -= 1;

      if (*size >= (1 << 23))
        return 0;

      *size <<= 8;
      *size |= ch;

      if (strict && *size == 0)
        return 0;
    }

    if (strict && *size < 0x80)
      return 0;
  }

  return 1;
}

static int
asn1_read_seq(const unsigned char **data, size_t *len, int strict) {
  size_t size;

  if (*len == 0 || **data != 0x30)
    return 0;

  *data += 1;
  *len -= 1;

  if (!asn1_read_size(&size, data, len, strict))
    return 0;

  if (strict && size != *len)
    return 0;

  return 1;
}

static int
asn1_read_int(mpz_t n, const unsigned char **data, size_t *len, int strict) {
  size_t size;

  if (*len == 0 || **data != 0x02)
    return 0;

  *data += 1;
  *len -= 1;

  if (!asn1_read_size(&size, data, len, strict))
    return 0;

  if (size > *len)
    return 0;

  if (strict) {
    const unsigned char *num = *data;

    /* No reason to have an integer larger than this. */
    if (size == 0 || size > 1 + RSA_MAX_MOD_BYTES)
      return 0;

    /* No negatives. */
    if (num[0] & 0x80)
      return 0;

    /* Allow zero only if it prefixes a high bit. */
    if (size > 1) {
      if (num[0] == 0x00 && (num[1] & 0x80) == 0x00)
        return 0;
    }
  }

  mpz_import(n, size, 1, 1, 0, 0, *data);

  *data += size;
  *len -= size;

  return 1;
}

static size_t
asn1_size_size(size_t size) {
  if (size <= 0x7f) /* [size] */
    return 1;

  if (size <= 0xff) /* 0x81 [size] */
    return 2;

  return 3; /* 0x82 [size-hi] [size-lo] */
}

static size_t
asn1_size_int(const mpz_t n) {
  /* 0x02 [size] [0x00?] [int] */
  size_t bits = mpz_bitlen(n);
  size_t size = (bits + 7) / 8;

  if ((bits & 7) == 0)
    size += mpz_tstbit(n, bits - 1);

  if (bits == 0)
    size = 1;

  return 1 + asn1_size_size(size) + size;
}

static size_t
asn1_write_size(unsigned char *data, size_t pos, size_t size) {
  if (size <= 0x7f)  {
    /* [size] */
    data[pos++] = size;
  } else if (size <= 0xff) {
    /* 0x81 [size] */
    data[pos++] = 0x81;
    data[pos++] = size;
  } else {
    /* 0x82 [size-hi] [size-lo] */
    assert(size <= 0xffff);
    data[pos++] = 0x82;
    data[pos++] = size >> 8;
    data[pos++] = size & 0xff;
  }
  return pos;
}

static size_t
asn1_write_int(unsigned char *data, size_t pos, const mpz_t n) {
  /* 0x02 [size] [0x00?] [int] */
  size_t bits = mpz_bitlen(n);
  size_t size = (bits + 7) / 8;
  size_t pad = 0;

  if ((bits & 7) == 0)
    pad = mpz_tstbit(n, bits - 1);

  if (bits == 0)
    size = 1;

  data[pos++] = 0x02;

  pos += asn1_write_size(data, pos, pad + size);

  if (pad)
    data[pos++] = 0x00;

  if (bits != 0)
    mpz_export(data + pos, NULL, 1, 1, 0, 0, n);
  else
    data[pos] = 0x00;

  pos += size;

  return pos;
}

/*
 * Private Key
 */

static void
priv_init(priv_t *k) {
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
priv_clear(priv_t *k) {
  mpz_cleanse(k->n);
  mpz_cleanse(k->e);
  mpz_cleanse(k->d);
  mpz_cleanse(k->p);
  mpz_cleanse(k->q);
  mpz_cleanse(k->dp);
  mpz_cleanse(k->dq);
  mpz_cleanse(k->qi);
}

static void
priv_import(priv_t *k, const unsigned char *data, size_t len, int strict) {
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
priv_export(unsigned char *out, size_t *out_len, const priv_t *k) {
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
priv_generate(priv_t *k, size_t bits, uint64_t exp, unsigned char *seed) {
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

  drbg_init(&rng, HASH_SHA256, seed, 32);

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
priv_verify(const priv_t *k) {
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
  mpz_mul(t, e, d);
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
priv_decrypt(const priv_t *k,
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
#ifdef TORSION_HAS_GMP
  if (mpz_sgn(k->d) > 0 && mpz_odd_p(k->n))
    mpz_powm_sec(m, m, k->d, k->n);
  else
    mpz_powm(m, m, k->d, k->n);
#else
  mpz_powm(m, m, k->d, k->n);
#endif

  /* m = m' * bi mod n (unblind) */
  mpz_mul(m, m, bi);
  mpz_mod(m, m, k->n);
  mpz_export_pad(out, m, mpz_bytelen(k->n));

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
pub_init(pub_t *k) {
  mpz_init(k->n);
  mpz_init(k->e);
}

static void
pub_clear(pub_t *k) {
  mpz_cleanse(k->n);
  mpz_cleanse(k->e);
}

static void
pub_import(pub_t *k, const unsigned char *data, size_t len, int strict) {
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
pub_export(unsigned char *out, size_t *out_len, const pub_t *k) {
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
pub_verify(const pub_t *k) {
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
pub_encrypt(const pub_t *k,
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
  mpz_export_pad(out, m, mpz_bytelen(k->n));

  r = 1;
fail:
  mpz_cleanse(m);
  return r;
}

/*
 * RSA
 */

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
  priv_t k;

  priv_init(&k);

  if (!get_digest_info(&prefix, &prefix_len, type))
    goto fail;

  if (type == -1)
    hlen = msg_len;

  if (msg_len != hlen)
    goto fail;

  if (!priv_import(&k, key, key_len, 0))
    goto fail;

  if (!priv_verify(&k))
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

  if (!priv_decrypt(&k, out, em, klen, entropy))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  priv_clear(&k);
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
  pub_t k;

  pub_init(&k);

  if (!get_digest_info(&prefix, &prefix_len, type))
    goto fail;

  if (type == -1)
    hlen = msg_len;

  if (msg_len != hlen)
    goto fail;

  if (!pub_import(&k, key, key_len, 0))
    goto fail;

  if (!pub_verify(&k))
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

  if (!pub_encrypt(&k, em, sig, sig_len))
    goto fail;

  /* EM = 0x00 || 0x01 || PS || 0x00 || T */
  ok = 1;

  ok &= safe_equal(em[0], 0x00);
  ok &= safe_equal(em[1], 0x01);

  for (let i = 2; i < klen - tlen - 1; i++)
    ok &= safe_equal(em[i], 0xff);

  ok &= safe_equal(em[klen - tlen - 1], 0x00);
  ok &= safe_cmp(em + klen - tlen, prefix, prefix_len);
  ok &= safe_cmp(em + klen - hlen, msg, msg_len);

  r = (ok == 1);
fail:
  pub_clear(&k);
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
  pub_t k;
  int r = 0;

  drbg_init(&rng, HASH_SHA256, entropy, 32);

  pub_init(&k);

  if (!pub_import(&k, key, key_len, 0))
    goto fail;

  if (!pub_verify(&k))
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

  if (!pub_encrypt(&k, out, em, klen))
    goto fail;

  *out_len = klen;
  r = 1;
fail:
  pub_clear(&k);
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
  unsigned char *em = out;
  uint32_t i, zero, two, index, looking;
  uint32_t equals0, validps, valid, offset;
  size_t klen = 0;
  priv_t k;
  int r = 0;

  priv_init(&k);

  if (!priv_import(&k, key, key_len, 0))
    goto fail;

  if (!priv_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  if (strict && msg_len != klen)
    goto fail;

  if (klen < 11)
    goto fail;

  if (!priv_decrypt(&k, em, msg, msg_len, entropy))
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

  memmove(out, em + offset, klen - offset);
  /*cleanse(em + klen - offset, offset);*/
  *out_len = klen - offset;

  r = 1;
fail:
  priv_clear(&k);
  if (r == 0) cleanse(out, klen);
  return r;
}

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
  pub_t k;

  pub_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (!pub_import(&k, key, key_len, 0))
    goto fail;

  if (!pub_verify(&k))
    goto fail;

  klen = mpz_bytelen(k->n);

  if (klen < 2 * hlen + 2 || msg_len > klen - 2 * hlen - 2)
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

  drbg_init(&rng, entropy, 32);
  drbg_generate(&rng, seed, slen);

  memcpy(&db[0], lhash, hlen);
  memset(&db[hlen], 0x00, (dlen - mlen - 1) - hlen);

  db[dlen - mlen - 1] = 0x01;
  memcpy(&db[dlen - mlen], msg, mlen);

  mgf1xor(type, db, dlen, seed, slen);
  mgf1xor(type, seed, slen, db, dlen);

  if (!pub_encrypt(&k, out, em, klen))
    goto fail;

  *out_len = klen;

  r = 1;
fail:
  pub_clear(&k);
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
  priv_t k;

  priv_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (!priv_import(&k, key, key_len, 0))
    goto fail;

  if (!priv_verify(&k))
    goto fail;

  klen = mpz_bytelen(k.n);

  if (strict && msg_len != klen)
    goto fail;

  if (klen < hlen * 2 + 2)
    goto fail;

  if (!priv_decrypt(&k, em, msg, msg_len, entropy))
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
  priv_clear(&k);
  cleanse(&hash, sizeof(hash));
  if (r == 0) cleanse(out, klen);
  return r;
}

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
           const unsigned char *em,
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
  priv_t k;
  int r = 0;

  priv_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (msg_len != hlen)
    goto fail;

  if (!priv_import(&k, key, key_len, 0))
    goto fail;

  if (!priv_verify(&k))
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
  if (!priv_decrypt(&k, out, em, emlen, entropy))
    goto fail;

  *outlen = klen;
  r = 1;
fail:
  priv_clear(&k);
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
  size_t emlen, bits;
  size_t klen = 0;
  pub_t k;
  int r = 0;

  pub_init(&k);

  if (!hash_has_backend(type))
    goto fail;

  if (!pub_import(&k, key, key_len, 0))
    goto fail;

  if (!pub_verify(&k))
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

  if (!pub_encrypt(&k, em, sig, sig_len))
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
  pub_clear(&k);
  safe_free(em, klen);
  return r;
}
