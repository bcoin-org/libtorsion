/*!
 * mpi.h - multi-precision integers for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 *
 * A from-scratch reimplementation of GMP.
 */

#ifndef _TORSION_MPI_H
#define _TORSION_MPI_H

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include "internal.h"

/*
 * Types
 */

#if defined(TORSION_HAVE_INT128)
typedef uint64_t mp_limb_t;
typedef torsion_uint128_t mp_wide_t;
#  define MP_LIMB_BITS 64
#  define MP_LIMB_BYTES 8
#  define MP_LIMB_C(x) UINT64_C(x)
#  define MP_LIMB_MAX MP_LIMB_C(0xffffffffffffffff)
#else
typedef uint32_t mp_limb_t;
typedef uint64_t mp_wide_t;
#  define MP_LIMB_BITS 32
#  define MP_LIMB_BYTES 4
#  define MP_LIMB_C(x) UINT32_C(x)
#  define MP_LIMB_MAX MP_LIMB_C(0xffffffff)
#endif

#define MP_LOW_BITS (MP_LIMB_BITS / 2)
#define MP_LOW_MASK (MP_LIMB_MAX >> MP_LOW_BITS)
#define MP_SIZE_MAX (INT_MAX / MP_LIMB_BITS)

TORSION_BARRIER(mp_limb_t, mp_limb)

struct mpz_s {
  mp_limb_t *limbs;
  int alloc;
  int size;
};

typedef struct mpz_s mpz_t[1];

typedef void mp_rng_f(void *out, size_t size, void *arg);

/*
 * Definitions
 */

#define MP_SLIDE_WIDTH 4
#define MP_SLIDE_SIZE (1 << (MP_SLIDE_WIDTH - 1))
#define MP_FIXED_WIDTH 4
#define MP_FIXED_SIZE (1 << MP_FIXED_WIDTH)

/*
 * Itches
 */

#define MPN_SQR_ITCH(n) (2 * (n))
#define MPN_MULSHIFT_ITCH(n) (2 * (n))
#define MPN_REDUCE_WEAK_ITCH(n) (n)
#define MPN_BARRETT_ITCH(n, shift) ((shift) + 1 - (n) + 1)
#define MPN_REDUCE_ITCH(n, shift) (1 + (shift) + ((shift) - (n) + 1))
#define MPN_MONT_ITCH(n) (2 * (n) + 1)
#define MPN_MONTMUL_ITCH(n) (2 * (n))
#define MPN_INVERT_ITCH(n) (4 * ((n) + 1))
#define MPN_JACOBI_ITCH(n) (2 * (n))
#define MPN_SLIDE_ITCH(yn, mn) ((yn) > 2 ? (MP_SLIDE_SIZE * (mn)) : 0)
#define MPN_POWM_ITCH(yn, mn) (6 * (mn) + MPN_SLIDE_ITCH(yn, mn))
#define MPN_SEC_POWM_ITCH(n) (5 * (n) + MP_FIXED_SIZE * (n) + 1)

/* Either Barrett or Montgomery precomputation. */
#define MPN_BARRETT_MONT_ITCH(shift) ((shift) + 2)

/*
 * Allocation
 */

mp_limb_t *
mp_alloc_limbs(int size);

mp_limb_t *
mp_realloc_limbs(mp_limb_t *ptr, int size);

void
mp_free_limbs(mp_limb_t *ptr);

/*
 * MPN Interface
 */

/*
 * Initialization
 */

void
mpn_zero(mp_limb_t *xp, int xn);

/*
 * Uninitialization
 */

void
mpn_cleanse(mp_limb_t *xp, int xn);

/*
 * Assignment
 */

void
mpn_set_1(mp_limb_t *xp, int xn, mp_limb_t y);

void
mpn_copy(mp_limb_t *zp, const mp_limb_t *xp, int xn);

/*
 * Comparison
 */

int
mpn_zero_p(const mp_limb_t *xp, int xn);

int
mpn_cmp(const mp_limb_t *xp, const mp_limb_t *yp, int n);

/*
 * Addition
 */

mp_limb_t
mpn_add_1(mp_limb_t *zp, const mp_limb_t *xp, int xn, mp_limb_t y);

mp_limb_t
mpn_add_n(mp_limb_t *zp, const mp_limb_t *xp, const mp_limb_t *yp, int n);

mp_limb_t
mpn_add(mp_limb_t *zp, const mp_limb_t *xp, int xn,
                       const mp_limb_t *yp, int yn);

/*
 * Subtraction
 */

mp_limb_t
mpn_sub_1(mp_limb_t *zp, const mp_limb_t *xp, int xn, mp_limb_t y);

mp_limb_t
mpn_sub_n(mp_limb_t *zp, const mp_limb_t *xp, const mp_limb_t *yp, int n);

mp_limb_t
mpn_sub(mp_limb_t *zp, const mp_limb_t *xp, int xn,
                       const mp_limb_t *yp, int yn);

/*
 * Multiplication
 */

mp_limb_t
mpn_mul_1(mp_limb_t *zp, const mp_limb_t *xp, int xn, mp_limb_t y);

mp_limb_t
mpn_addmul_1(mp_limb_t *zp, const mp_limb_t *xp, int xn, mp_limb_t y);

mp_limb_t
mpn_submul_1(mp_limb_t *zp, const mp_limb_t *xp, int xn, mp_limb_t y);

void
mpn_mul_n(mp_limb_t *zp, const mp_limb_t *xp, const mp_limb_t *yp, int n);

void
mpn_mul(mp_limb_t *zp, const mp_limb_t *xp, int xn,
                       const mp_limb_t *yp, int yn);

void
mpn_sqr(mp_limb_t *zp, const mp_limb_t *xp, int xn, mp_limb_t *scratch);

/*
 * Multiply + Shift
 */

mp_limb_t
mpn_mulshift(mp_limb_t *zp,
             const mp_limb_t *xp,
             const mp_limb_t *yp,
             int n, int bits,
             mp_limb_t *scratch);

/*
 * Weak Reduction
 */

int
mpn_reduce_weak(mp_limb_t *zp,
                const mp_limb_t *xp,
                const mp_limb_t *np,
                int n, mp_limb_t hi,
                mp_limb_t *scratch);

/*
 * Barrett Reduction
 */

void
mpn_barrett(mp_limb_t *mp, const mp_limb_t *np,
            int n, int shift, mp_limb_t *scratch);

void
mpn_reduce(mp_limb_t *zp, const mp_limb_t *xp,
                          const mp_limb_t *mp,
                          const mp_limb_t *np,
                          int n, int shift,
                          mp_limb_t *scratch);

/*
 * Montgomery Multiplication
 */

void
mpn_mont(mp_limb_t *kp, mp_limb_t *rp,
         const mp_limb_t *mp, int n,
         mp_limb_t *scratch);

void
mpn_montmul(mp_limb_t *zp,
            const mp_limb_t *xp,
            const mp_limb_t *yp,
            const mp_limb_t *mp,
            int n, mp_limb_t k,
            mp_limb_t *scratch);

void
mpn_montmul_var(mp_limb_t *zp,
                const mp_limb_t *xp,
                const mp_limb_t *yp,
                const mp_limb_t *mp,
                int n, mp_limb_t k,
                mp_limb_t *scratch);

/*
 * Division
 */

mp_limb_t
mpn_divmod_1(mp_limb_t *qp, const mp_limb_t *np, int nn, mp_limb_t d);

void
mpn_div_1(mp_limb_t *qp, const mp_limb_t *np, int nn, mp_limb_t d);

mp_limb_t
mpn_mod_1(const mp_limb_t *np, int nn, mp_limb_t d);

mp_wide_t
mpn_mod_2(const mp_limb_t *np, int nn, mp_wide_t d);

void
mpn_divmod(mp_limb_t *qp, mp_limb_t *rp,
           const mp_limb_t *np, int nn,
           const mp_limb_t *dp, int dn);

void
mpn_div(mp_limb_t *qp, const mp_limb_t *np, int nn,
                       const mp_limb_t *dp, int dn);

void
mpn_mod(mp_limb_t *rp, const mp_limb_t *np, int nn,
                       const mp_limb_t *dp, int dn);

/*
 * Left Shift
 */

mp_limb_t
mpn_lshift(mp_limb_t *zp, const mp_limb_t *xp, int xn, int bits);

/*
 * Right Shift
 */

mp_limb_t
mpn_rshift(mp_limb_t *zp, const mp_limb_t *xp, int xn, int bits);

/*
 * Bit Manipulation
 */

mp_limb_t
mpn_get_bit(const mp_limb_t *xp, int xn, int pos);

mp_limb_t
mpn_get_bits(const mp_limb_t *xp, int xn, int pos, int width);

/*
 * Number Theoretic Functions
 */

int
mpn_invert(mp_limb_t *zp,
           const mp_limb_t *xp, int xn,
           const mp_limb_t *yp, int yn,
           mp_limb_t *scratch);

int
mpn_invert_n(mp_limb_t *zp,
             const mp_limb_t *xp,
             const mp_limb_t *yp,
             int n,
             mp_limb_t *scratch);

int
mpn_jacobi(const mp_limb_t *xp, int xn,
           const mp_limb_t *yp, int yn,
           mp_limb_t *scratch);

int
mpn_jacobi_n(const mp_limb_t *xp,
             const mp_limb_t *yp,
             int n,
             mp_limb_t *scratch);

void
mpn_powm(mp_limb_t *zp, const mp_limb_t *xp, int xn,
                        const mp_limb_t *yp, int yn,
                        const mp_limb_t *mp, int mn,
                        mp_limb_t *scratch);

void
mpn_sec_powm(mp_limb_t *zp,
             const mp_limb_t *xp, int xn,
             const mp_limb_t *yp, int yn,
             const mp_limb_t *mp, int mn,
             mp_limb_t *scratch);

/*
 * Helpers
 */

int
mpn_strip(const mp_limb_t *xp, int xn);

int
mpn_odd_p(const mp_limb_t *xp, int xn);

int
mpn_even_p(const mp_limb_t *xp, int xn);

int
mpn_ctz(const mp_limb_t *xp, int xn);

int
mpn_bitlen(const mp_limb_t *xp, int xn);

size_t
mpn_bytelen(const mp_limb_t *xp, int xn);

void
mpn_swap(mp_limb_t **xp, int *xn,
         mp_limb_t **yp, int *yn);

/*
 * Constant Time
 */

void
mpn_select(mp_limb_t *zp,
           const mp_limb_t *xp,
           const mp_limb_t *yp,
           int n, int flag);

void
mpn_select_zero(mp_limb_t *zp, const mp_limb_t *xp, int n, int flag);

int
mpn_sec_zero_p(const mp_limb_t *xp, int xn);

int
mpn_sec_equal(const mp_limb_t *xp, const mp_limb_t *yp, int n);

int
mpn_sec_cmp(const mp_limb_t *xp, const mp_limb_t *yp, int n);

int
mpn_sec_lt(const mp_limb_t *xp, const mp_limb_t *yp, int n);

int
mpn_sec_lte(const mp_limb_t *xp, const mp_limb_t *yp, int n);

int
mpn_sec_gt(const mp_limb_t *xp, const mp_limb_t *yp, int n);

int
mpn_sec_gte(const mp_limb_t *xp, const mp_limb_t *yp, int n);

/*
 * Import
 */

void
mpn_import(mp_limb_t *zp, int zn,
           const unsigned char *raw,
           size_t len, int endian);

/*
 * Export
 */

void
mpn_export(unsigned char *raw, size_t len,
           const mp_limb_t *xp, int xn, int endian);

/*
 * MPZ Interface
 */

/*
 * Initialization
 */

void
mpz_init(mpz_t x);

/*
 * Uninitialization
 */

void
mpz_clear(mpz_t x);

void
mpz_cleanse(mpz_t x);

/*
 * Assignment
 */

void
mpz_set(mpz_t z, const mpz_t x);

void
mpz_roset(mpz_t z, const mpz_t x);

void
mpz_set_ui(mpz_t z, mp_limb_t x);

void
mpz_set_u64(mpz_t z, uint64_t x);

/*
 * Conversion
 */

mp_limb_t
mpz_get_ui(const mpz_t x);

uint64_t
mpz_get_u64(const mpz_t x);

/*
 * Comparison
 */

int
mpz_sgn(const mpz_t x);

int
mpz_cmp(const mpz_t x, const mpz_t y);

int
mpz_cmp_ui(const mpz_t x, mp_limb_t y);

/*
 * Unsigned Comparison
 */

int
mpz_cmpabs(const mpz_t x, const mpz_t y);

int
mpz_cmpabs_ui(const mpz_t x, mp_limb_t y);

/*
 * Addition
 */

void
mpz_add(mpz_t z, const mpz_t x, const mpz_t y);

void
mpz_add_ui(mpz_t z, const mpz_t x, mp_limb_t y);

/*
 * Subtraction
 */

void
mpz_sub(mpz_t z, const mpz_t x, const mpz_t y);

void
mpz_sub_ui(mpz_t z, const mpz_t x, mp_limb_t y);

/*
 * Multiplication
 */

void
mpz_mul(mpz_t z, const mpz_t x, const mpz_t y);

void
mpz_mul_ui(mpz_t z, const mpz_t x, mp_limb_t y);

void
mpz_sqr(mpz_t z, const mpz_t x);

/*
 * Truncation Division
 */

void
mpz_quorem(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d);

void
mpz_quo(mpz_t q, const mpz_t n, const mpz_t d);

void
mpz_rem(mpz_t r, const mpz_t n, const mpz_t d);

mp_limb_t
mpz_quo_ui(mpz_t q, const mpz_t n, mp_limb_t d);

mp_limb_t
mpz_rem_ui(const mpz_t n, mp_limb_t d);

/*
 * Euclidean Division
 */

void
mpz_divmod(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d);

void
mpz_div(mpz_t q, const mpz_t n, const mpz_t d);

void
mpz_mod(mpz_t r, const mpz_t n, const mpz_t d);

mp_limb_t
mpz_div_ui(mpz_t q, const mpz_t n, mp_limb_t d);

mp_limb_t
mpz_mod_ui(const mpz_t n, mp_limb_t d);

/*
 * Exact Division
 */

void
mpz_divexact(mpz_t q, const mpz_t n, const mpz_t d);

void
mpz_divexact_ui(mpz_t q, const mpz_t n, mp_limb_t d);

/*
 * Left Shift
 */

void
mpz_lshift(mpz_t z, const mpz_t x, int bits);

/*
 * Right Shift
 */

void
mpz_rshift(mpz_t z, const mpz_t x, int bits);

/*
 * Bit Manipulation
 */

mp_limb_t
mpz_get_bit(const mpz_t x, int pos);

mp_limb_t
mpz_get_bits(const mpz_t x, int pos, int width);

void
mpz_set_bit(mpz_t x, int pos);

void
mpz_clr_bit(mpz_t x, int pos);

/*
 * Negation
 */

void
mpz_abs(mpz_t z, const mpz_t x);

void
mpz_neg(mpz_t z, const mpz_t x);

/*
 * Number Theoretic Functions
 */

void
mpz_gcd(mpz_t z, const mpz_t x, const mpz_t y);

void
mpz_lcm(mpz_t z, const mpz_t x, const mpz_t y);

void
mpz_gcdext(mpz_t g, mpz_t s, mpz_t t, const mpz_t x, const mpz_t y);

int
mpz_invert(mpz_t z, const mpz_t x, const mpz_t y);

int
mpz_jacobi(const mpz_t x, const mpz_t y);

void
mpz_powm(mpz_t z, const mpz_t x, const mpz_t y, const mpz_t m);

void
mpz_powm_sec(mpz_t z, const mpz_t x, const mpz_t y, const mpz_t m);

/*
 * Primality Testing
 */

int
mpz_is_prime_mr(const mpz_t n, int reps, int force2, mp_rng_f *rng, void *arg);

int
mpz_is_prime_lucas(const mpz_t n, mp_limb_t limit);

int
mpz_is_prime(const mpz_t n, int rounds, mp_rng_f *rng, void *arg);

void
mpz_random_prime(mpz_t z, int bits, mp_rng_f *rng, void *arg);

/*
 * Helpers
 */

int
mpz_odd_p(const mpz_t x);

int
mpz_even_p(const mpz_t x);

int
mpz_ctz(const mpz_t x);

int
mpz_bitlen(const mpz_t x);

size_t
mpz_bytelen(const mpz_t x);

void
mpz_swap(mpz_t x, mpz_t y);

/*
 * Import
 */

void
mpz_import(mpz_t z, const unsigned char *raw, size_t size, int endian);

/*
 * Export
 */

void
mpz_export(unsigned char *raw, const mpz_t x, size_t size, int endian);

/*
 * RNG
 */

void
mpz_random_bits(mpz_t z, int bits, mp_rng_f *rng, void *arg);

void
mpz_random_int(mpz_t z, const mpz_t max, mp_rng_f *rng, void *arg);

#endif /* _TORSION_MPI_H */
