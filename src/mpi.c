/*!
 * mpi.c - multi-precision integers for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 *
 * Parts of this software are based on GMP:
 *   https://gmplib.org/
 *   Copyright (C) 1989, 1991 Free Software Foundation, Inc.
 *
 * Our MPI implementation is a fork of mini-gmp, modified to include
 * inline x86-64 assembly from gmp. In addition, all unnecessary
 * functionality has been stripped out. A notable change is that
 * we have switch to fixed/predictable integer widths via stdint.
 * Furthermore, we care less about signed integer arithmetic with
 * regards to bit operations.
 *
 * See the mini-gmp license below.
 *
 * mini-gmp, a minimalistic implementation of a GNU GMP subset.
 *
 * Contributed to the GNU project by Niels Möller
 *
 * Copyright 1991-1997, 1999-2019 Free Software Foundation, Inc.
 *
 * This file is part of the GNU MP Library.
 *
 * The GNU MP Library is free software; you can redistribute it and/or modify
 * it under the terms of either:
 *
 *   * the GNU Lesser General Public License as published by the Free
 *     Software Foundation; either version 3 of the License, or (at your
 *     option) any later version.
 *
 * or
 *
 *   * the GNU General Public License as published by the Free Software
 *     Foundation; either version 2 of the License, or (at your option) any
 *     later version.
 *
 * or both in parallel, as here.
 *
 * The GNU MP Library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received copies of the GNU General Public License and the
 * GNU Lesser General Public License along with the GNU MP Library.  If not,
 * see https://www.gnu.org/licenses/.
 */

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "mpi.h"

/*
 * Prototypes
 */

void
torsion_cleanse(void *, size_t);

/*
 * Types
 */

struct mpi_div_inverse {
  /* Normalization shift count. */
  unsigned shift;
  /* Normalized divisor (d0 unused for mpn_div_qr_1) */
  mp_limb_t d1, d0;
  /* Inverse, for 2/1 or 3/2. */
  mp_limb_t di;
};

enum mpz_div_round_mode { MPI_DIV_FLOOR, MPI_DIV_CEIL, MPI_DIV_TRUNC };

/*
 * Macros
 */

#define MPI_LIMB_MAX ((mp_limb_t)~(mp_limb_t)0)
#define MPI_LIMB_HIGHBIT ((mp_limb_t)1 << (MPI_LIMB_BITS - 1))

#define MPI_HLIMB_BIT ((mp_limb_t)1 << (MPI_LIMB_BITS / 2))
#define MPI_LLIMB_MASK (MPI_HLIMB_BIT - 1)

#define MPI_ABS(x) ((x) >= 0 ? (x) : -(x))
#define MPI_NEG_CAST(T, x) (-((T)((x) + 1) - 1))

#define MPI_MIN(a, b) ((a) < (b) ? (a) : (b))
#define MPI_MAX(a, b) ((a) > (b) ? (a) : (b))

#define MPI_CMP(a, b) (((a) > (b)) - ((a) < (b)))

#define MPI_MPN_OVERLAP_P(xp, xsize, yp, ysize) \
  ((xp) + (xsize) > (yp) && (yp) + (ysize) > (xp))

#define mpi_assert_nocarry(x) do { \
  mp_limb_t __cy = (x);            \
  assert(__cy == 0);               \
} while (0)

#define mpn_invert_limb(x) mpn_invert_3by2((x), 0)

#define MPZ_REALLOC(z, n) \
  ((n) > (z)->_mp_alloc ? mpz_realloc(z, n) : (z)->_mp_d)

/*
 * Low-level Arithmetic Macros
 * See: https://gmplib.org/repo/gmp-6.2/file/tip/longlong.h#l1044
 */

#ifdef TORSION_USE_ASM

#define mpi_clz(count, x) do { \
  uint64_t __cbtmp;            \
  assert((x) != 0);            \
  __asm__ (                    \
    "bsr %1, %0\n"             \
    : "=r" (__cbtmp)           \
    : "rm" ((uint64_t)(x))     \
  );                           \
  (count) = __cbtmp ^ 63;      \
} while (0)

#define mpi_ctz(count, x) do { \
  uint64_t __cbtmp;            \
  assert((x) != 0);            \
  __asm__ (                    \
    "bsf %1, %q0\n"            \
    : "=r" (__cbtmp)           \
    : "rm" ((uint64_t)(x))     \
  );                           \
  (count) = __cbtmp;           \
} while (0)

#define mpi_add_ssaaaa(sh, sl, ah, al, bh, bl)      \
  __asm__ (                                         \
    "addq %5, %q1\n"                                \
    "adcq %3, %q0\n"                                \
    : "=r" (sh), "=&r" (sl)                         \
    : "0" ((uint64_t)(ah)), "rme" ((uint64_t)(bh)), \
      "%1" ((uint64_t)(al)), "rme" ((uint64_t)(bl)) \
  )

#define mpi_sub_ddmmss(sh, sl, ah, al, bh, bl)      \
  __asm__ (                                         \
    "subq %5, %q1\n"                                \
    "sbbq %3, %q0\n"                                \
    : "=r" (sh), "=&r" (sl)                         \
    : "0" ((uint64_t)(ah)), "rme" ((uint64_t)(bh)), \
      "1" ((uint64_t)(al)), "rme" ((uint64_t)(bl))  \
  )

#define mpi_umul_ppmm(w1, w0, u, v) \
  __asm__ (                         \
    "mulq %3\n"                     \
    : "=a" (w0), "=d" (w1)          \
    : "%0" ((uint64_t)(u)),         \
      "rm" ((uint64_t)(v))          \
  )

#else /* TORSION_USE_ASM */

#define mpi_clz(count, x) do {                                      \
  mp_limb_t __clz_x = (x);                                          \
  unsigned __clz_c = 0;                                             \
  for (; (__clz_x & ((mp_limb_t)0xff << (MPI_LIMB_BITS - 8))) == 0; \
         __clz_c += 8) {                                            \
    __clz_x <<= 8;                                                  \
  }                                                                 \
  for (; (__clz_x & MPI_LIMB_HIGHBIT) == 0; __clz_c++)              \
    __clz_x <<= 1;                                                  \
  (count) = __clz_c;                                                \
} while (0)

#define mpi_ctz(count, x) do {           \
  mp_limb_t __ctz_x = (x);               \
  unsigned __ctz_c = 0;                  \
  mpi_clz(__ctz_c, __ctz_x & -__ctz_x);  \
  (count) = MPI_LIMB_BITS - 1 - __ctz_c; \
} while (0)

#define mpi_add_ssaaaa(sh, sl, ah, al, bh, bl) do { \
  mp_limb_t __x;                                    \
  __x = (al) + (bl);                                \
  (sh) = (ah) + (bh) + (__x < (al));                \
  (sl) = __x;                                       \
} while (0)

#define mpi_sub_ddmmss(sh, sl, ah, al, bh, bl) do { \
  mp_limb_t __x;                                    \
  __x = (al) - (bl);                                \
  (sh) = (ah) - (bh) - ((al) < (bl));               \
  (sl) = __x;                                       \
} while (0)

#ifdef MPI_HAS_WIDE
#define mpi_umul_ppmm(w1, w0, u, v) do {   \
  mp_wide_t __ww = (mp_wide_t)(u) * (v);   \
  w0 = (mp_limb_t)__ww;                    \
  w1 = (mp_limb_t)(__ww >> MPI_LIMB_BITS); \
} while (0)
#else
#define mpi_umul_ppmm(w1, w0, u, v) do {                                    \
  mp_limb_t __x0, __x1, __x2, __x3;                                         \
  unsigned __ul, __vl, __uh, __vh;                                          \
  mp_limb_t __u = (u), __v = (v);                                           \
                                                                            \
  __ul = __u & MPI_LLIMB_MASK;                                              \
  __uh = __u >> (MPI_LIMB_BITS / 2);                                        \
  __vl = __v & MPI_LLIMB_MASK;                                              \
  __vh = __v >> (MPI_LIMB_BITS / 2);                                        \
                                                                            \
  __x0 = (mp_limb_t)__ul * __vl;                                            \
  __x1 = (mp_limb_t)__ul * __vh;                                            \
  __x2 = (mp_limb_t)__uh * __vl;                                            \
  __x3 = (mp_limb_t)__uh * __vh;                                            \
                                                                            \
  __x1 += __x0 >> (MPI_LIMB_BITS / 2); /* this can't give carry */          \
  __x1 += __x2;                        /* but this indeed can */            \
                                                                            \
  if (__x1 < __x2)                     /* did we get it? */                 \
    __x3 += MPI_HLIMB_BIT;             /* yes, add it in the proper pos. */ \
                                                                            \
  (w1) = __x3 + (__x1 >> (MPI_LIMB_BITS / 2));                              \
  (w0) = (__x1 << (MPI_LIMB_BITS / 2)) + (__x0 & MPI_LLIMB_MASK);           \
} while (0)
#endif

#endif /* TORSION_USE_ASM */

#define mpi_udiv_qrnnd_preinv(q, r, nh, nl, d, di) do {      \
  mp_limb_t _qh, _ql, _r, _mask;                             \
                                                             \
  mpi_umul_ppmm(_qh, _ql, (nh), (di));                       \
  mpi_add_ssaaaa(_qh, _ql, _qh, _ql, (nh) + 1, (nl));        \
                                                             \
  _r = (nl) - _qh * (d);                                     \
  _mask = -(mp_limb_t)(_r > _ql); /* both > and >= are OK */ \
  _qh += _mask;                                              \
  _r += _mask & (d);                                         \
                                                             \
  if (_r >= (d)) {                                           \
    _r -= (d);                                               \
    _qh++;                                                   \
  }                                                          \
                                                             \
  (r) = _r;                                                  \
  (q) = _qh;                                                 \
} while (0)

#define mpi_udiv_qr_3by2(q, r1, r0, n2, n1, n0, d1, d0, dinv) do {    \
  mp_limb_t _q0, _t1, _t0, _mask;                                     \
  mpi_umul_ppmm((q), _q0, (n2), (dinv));                              \
  mpi_add_ssaaaa((q), _q0, (q), _q0, (n2), (n1));                     \
                                                                      \
  /* Compute the two most significant limbs of n - q'd */             \
  (r1) = (n1) - (d1) * (q);                                           \
  mpi_sub_ddmmss((r1), (r0), (r1), (n0), (d1), (d0));                 \
  mpi_umul_ppmm(_t1, _t0, (d0), (q));                                 \
  mpi_sub_ddmmss((r1), (r0), (r1), (r0), _t1, _t0);                   \
  (q)++;                                                              \
                                                                      \
  /* Conditionally adjust q and the remainders */                     \
  _mask = - (mp_limb_t)((r1) >= _q0);                                 \
  (q) += _mask;                                                       \
  mpi_add_ssaaaa((r1), (r0), (r1), (r0), _mask & (d1), _mask & (d0)); \
                                                                      \
  if ((r1) >= (d1)) {                                                 \
    if ((r1) > (d1) || (r0) >= (d0)) {                                \
      (q)++;                                                          \
      mpi_sub_ddmmss((r1), (r0), (r1), (r0), (d1), (d0));             \
    }                                                                 \
  }                                                                   \
} while (0)

/*
 * Swap Macros
 */

#define MP_LIMB_T_SWAP(x, y) do {        \
  mp_limb_t __mp_limb_t_swap__tmp = (x); \
  (x) = (y);                             \
  (y) = __mp_limb_t_swap__tmp;           \
} while (0)

#define MP_SIZE_T_SWAP(x, y) do {        \
  mp_size_t __mp_size_t_swap__tmp = (x); \
  (x) = (y);                             \
  (y) = __mp_size_t_swap__tmp;           \
} while (0)

#define MP_BITCNT_T_SWAP(x, y) do {          \
  mp_bitcnt_t __mp_bitcnt_t_swap__tmp = (x); \
  (x) = (y);                                 \
  (y) = __mp_bitcnt_t_swap__tmp;             \
} while (0)

#define MP_PTR_SWAP(x, y) do {     \
  mp_ptr __mp_ptr_swap__tmp = (x); \
  (x) = (y);                       \
  (y) = __mp_ptr_swap__tmp;        \
} while (0)

#define MP_SRCPTR_SWAP(x, y) do {        \
  mp_srcptr __mp_srcptr_swap__tmp = (x); \
  (x) = (y);                             \
  (y) = __mp_srcptr_swap__tmp;           \
} while (0)

#define MPN_PTR_SWAP(xp, xs, yp, ys) do { \
  MP_PTR_SWAP(xp, yp);                    \
  MP_SIZE_T_SWAP(xs, ys);                 \
} while (0)

#define MPN_SRCPTR_SWAP(xp, xs, yp, ys) do { \
  MP_SRCPTR_SWAP(xp, yp);                    \
  MP_SIZE_T_SWAP(xs, ys);                    \
} while (0)

#define MPZ_PTR_SWAP(x, y) do {      \
  mpz_ptr __mpz_ptr_swap__tmp = (x); \
  (x) = (y);                         \
  (y) = __mpz_ptr_swap__tmp;         \
} while (0)

#define MPZ_SRCPTR_SWAP(x, y) do {         \
  mpz_srcptr __mpz_srcptr_swap__tmp = (x); \
  (x) = (y);                               \
  (y) = __mpz_srcptr_swap__tmp;            \
} while (0)

/*
 * MPN Arithmetic Macros
 */

#define MPN_MAKE_ODD(bits, rp, rn) do {                        \
  mp_size_t __bits = mpn_ctz(rp, rn);                          \
  mp_size_t __limbs = __bits / MPI_LIMB_BITS;                  \
  mp_size_t __shift = __bits % MPI_LIMB_BITS;                  \
                                                               \
  if (__shift != 0) {                                          \
    mpn_rshift((rp), (rp) + __limbs, (rn) - __limbs, __shift); \
    (rn) -= __limbs;                                           \
    (rn) -= ((rp)[(rn) - 1] == 0);                             \
  } else if (__limbs != 0) {                                   \
    mpn_copyi((rp), (rp) + __limbs, (rn) - __limbs);           \
    (rn) -= __limbs;                                           \
  }                                                            \
                                                               \
  (bits) = __bits;                                             \
} while (0)

#define MPN_ADD(ap, as, bp, bs) do {        \
  mp_size_t __an;                           \
                                            \
  if (((as) ^ (bs)) >= 0)                   \
    __an = mpn_abs_add(ap, ap, as, bp, bs); \
  else                                      \
    __an = mpn_abs_sub(ap, ap, as, bp, bs); \
                                            \
  (as) = (as) >= 0 ? __an : -__an;          \
} while (0)

#define MPN_SUB(ap, as, bp, bs) do {        \
  mp_size_t __an;                           \
                                            \
  if (((as) ^ (bs)) >= 0)                   \
    __an = mpn_abs_sub(ap, ap, as, bp, bs); \
  else                                      \
    __an = mpn_abs_add(ap, ap, as, bp, bs); \
                                            \
  (as) = (as) >= 0 ? __an : -__an;          \
} while (0)

#define MPN_ODD_P(ap, as) \
  (MPI_ABS(as) > 0 && ((ap)[0] & 1) == 1)

#define MPN_RSHIFT(ap, as, bits) do { \
  mp_size_t __an = MPI_ABS(as);       \
                                      \
  if (__an > 0 && (bits) != 0) {      \
    mpn_rshift(ap, ap, __an, bits);   \
    __an -= ((ap)[__an - 1] == 0);    \
  }                                   \
                                      \
  (as) = (as) >= 0 ? __an : -__an;    \
} while (0)

/* itch = (nn - dn + 1) + dn */
#define MPN_MOD(np, nn, dp, dn, sp) do {                      \
  if (mpn_cmp4(np, nn, dp, dn) >= 0) {                        \
    mpn_div_tp(sp, np, nn, dp, dn, (sp) + ((nn) - (dn) + 1)); \
    (nn) = mpn_normalized_size(np, dn);                       \
  }                                                           \
} while (0)

/*
 * Allocation
 */

static void
mpi_die(const char *msg) {
  fprintf(stderr, "%s\n", msg);
  abort();
}

static void *
mpi_xalloc(size_t size) {
  void *p;

  assert(size > 0);

  p = malloc(size);

  if (!p)
    mpi_die("mpi_xalloc: Virtual memory exhausted.");

  return p;
}

static void *
mpi_xrealloc(void *old, size_t size) {
  void *p;

  p = realloc(old, size);

  if (!p)
    mpi_die("mpi_xrealloc: Virtual memory exhausted.");

  return p;
}

static void
mpi_free(void *p) {
  free(p);
}

static mp_ptr
mpi_xalloc_limbs(mp_size_t size) {
  return mpi_xalloc(size * sizeof(mp_limb_t));
}

static mp_ptr
mpi_xrealloc_limbs(mp_ptr old, mp_size_t size) {
  assert(size > 0);
  return mpi_xrealloc(old, size * sizeof(mp_limb_t));
}

/*
 * MPN Interface
 */

/*
 * Initialization
 */

void
mpn_zero(mp_ptr rp, mp_size_t n) {
  while (--n >= 0)
    rp[n] = 0;
}

/*
 * Uninitialization
 */

void
mpn_cleanse(mp_ptr xp, mp_size_t xn) {
  torsion_cleanse(xp, xn * sizeof(mp_limb_t));
}

/*
 * Assignment
 */

void
mpn_copyi(mp_ptr d, mp_srcptr s, mp_size_t n) {
#ifdef TORSION_USE_ASM
  /* From:
   * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/copyi.asm
   *
   * Registers:
   *
   *   %rdi = rp
   *   %rsi = up
   *   %rdx = n
   */
  __asm__ __volatile__(
    "lea -8(%%rdi), %%rdi\n"
    "sub $4, %%rdx\n"
    "jc 2f\n" /* end */

    "1:\n" /* top */
    "mov (%%rsi), %%rax\n"
    "mov 8(%%rsi), %%r9\n"
    "lea 32(%%rdi), %%rdi\n"
    "mov 16(%%rsi), %%r10\n"
    "mov 24(%%rsi), %%r11\n"
    "lea 32(%%rsi), %%rsi\n"
    "mov %%rax, -24(%%rdi)\n"
    "mov %%r9, -16(%%rdi)\n"
    "sub $4, %%rdx\n"
    "mov %%r10, -8(%%rdi)\n"
    "mov %%r11, (%%rdi)\n"
    "jnc 1b\n" /* top */

    "2:\n" /* end */
    "shr %%edx\n"
    "jnc 3f\n"
    "mov (%%rsi), %%rax\n"
    "mov %%rax, 8(%%rdi)\n"
    "lea 8(%%rdi), %%rdi\n"
    "lea 8(%%rsi), %%rsi\n"
    "3:\n"
    "shr %%edx\n"
    "jnc 4f\n"
    "mov (%%rsi), %%rax\n"
    "mov 8(%%rsi), %%r9\n"
    "mov %%rax, 8(%%rdi)\n"
    "mov %%r9, 16(%%rdi)\n"
    "4:\n"
    : "+D" (d), "+S" (s), "+d" (n)
    :
    : "rax", "r9", "r10", "r11",
      "cc", "memory"
  );
#else
  mp_size_t i;
  for (i = 0; i < n; i++)
    d[i] = s[i];
#endif
}

void
mpn_copyd(mp_ptr d, mp_srcptr s, mp_size_t n) {
#ifdef TORSION_USE_ASM
  /* From:
   * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/copyd.asm
   *
   * Registers:
   *
   *   %rdi = rp
   *   %rsi = up
   *   %rdx = n
   */
  __asm__ __volatile__(
    "lea -8(%%rsi,%%rdx,8), %%rsi\n"
    "lea (%%rdi,%%rdx,8), %%rdi\n"
    "sub $4, %%rdx\n"
    "jc 2f\n" /* end */
    "nop\n"

    "1:\n" /* top */
    "mov (%%rsi), %%rax\n"
    "mov -8(%%rsi), %%r9\n"
    "lea -32(%%rdi), %%rdi\n"
    "mov -16(%%rsi), %%r10\n"
    "mov -24(%%rsi), %%r11\n"
    "lea -32(%%rsi), %%rsi\n"
    "mov %%rax, 24(%%rdi)\n"
    "mov %%r9, 16(%%rdi)\n"
    "sub $4, %%rdx\n"
    "mov %%r10, 8(%%rdi)\n"
    "mov %%r11, (%%rdi)\n"
    "jnc 1b\n" /* top */

    "2:\n" /* end */
    "shr %%edx\n"
    "jnc 3f\n"
    "mov (%%rsi), %%rax\n"
    "mov %%rax, -8(%%rdi)\n"
    "lea -8(%%rdi), %%rdi\n"
    "lea -8(%%rsi), %%rsi\n"
    "3:\n"
    "shr %%edx\n"
    "jnc 4f\n"
    "mov (%%rsi), %%rax\n"
    "mov -8(%%rsi), %%r9\n"
    "mov %%rax, -8(%%rdi)\n"
    "mov %%r9, -16(%%rdi)\n"
    "4:\n"
    : "+D" (d), "+S" (s), "+d" (n)
    :
    : "rax", "r9", "r10", "r11",
      "cc", "memory"
  );
#else
  while (--n >= 0)
    d[n] = s[n];
#endif
}

/*
 * Comparison
 */

int
mpn_zero_p(mp_srcptr rp, mp_size_t n) {
  return mpn_normalized_size(rp, n) == 0;
}

int
mpn_cmp(mp_srcptr ap, mp_srcptr bp, mp_size_t n) {
  while (--n >= 0) {
    if (ap[n] != bp[n])
      return ap[n] > bp[n] ? 1 : -1;
  }
  return 0;
}

int
mpn_cmp4(mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn) {
  if (an != bn)
    return an < bn ? -1 : 1;
  else
    return mpn_cmp(ap, bp, an);
}

static int
mpn_cmp_1(mp_srcptr ap, mp_size_t an, mp_limb_t b) {
  mp_limb_t a;

  if (an > 1)
    return 1;

  if (an > 0)
    a = ap[0];
  else
    a = 0;

  return MPI_CMP(a, b);
}

/*
 * Addition/Subtraction Engines (X86-64 ASM)
 */

/* Probably more constant-time than the
 * C code below. This is only half necessary
 * as modern versions of GCC will compile
 * mpn_add_1 as constant-time, which is not
 * the case for mpn_sub_1.
 *
 * Registers:
 *
 *   %rdi = rp
 *   %rsi = ap
 *   %rcx = n
 *   %rdx = b
 */
#define AORS_1(ADDSUB, ADCSBB) \
  __asm__ __volatile__(        \
    "mov (%%rsi), %%rax\n"     \
    ADDSUB " %%rdx, %%rax\n"   \
    "mov %%rax, (%%rdi)\n"     \
    "dec %%rcx\n"              \
    "jz 2f\n"                  \
                               \
    ".align 16\n"              \
    "1:\n"                     \
    "lea 8(%%rsi), %%rsi\n"    \
    "lea 8(%%rdi), %%rdi\n"    \
    "mov (%%rsi), %%rax\n"     \
    ADCSBB " $0, %%rax\n"      \
    "mov %%rax, (%%rdi)\n"     \
    "dec %%rcx\n"              \
    "jnz 1b\n"                 \
                               \
    "2:\n"                     \
    "mov $0, %%rdx\n"          \
    "adc $0, %%rdx\n"          \
    : "+D" (rp), "+S" (ap),    \
      "+c" (n), "+d" (b)       \
    :                          \
    : "rax", "cc", "memory"    \
  );

/* From:
 * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/aors_n.asm
 *
 * Registers:
 *
 *   %rdi = rp (rcx)
 *   %rsi = up (rdx)
 *   %rdx = vp (r8)
 *   %rcx = n (r9)
 */
#define AORS_N(ADCSBB)                      \
  __asm__ __volatile__(                     \
    "mov %%ecx, %%eax\n"                    \
    /* guard for n <= 0 */                  \
    "test %%eax, %%eax\n"                   \
    "jle 7f\n" /* exit */                   \
    /* end guard */                         \
    "shr $2, %%rcx\n"                       \
    "and $3, %%eax\n"                       \
    "jrcxz 1f\n" /* lt4 */                  \
                                            \
    "mov (%%rsi), %%r8\n"                   \
    "mov 8(%%rsi), %%r9\n"                  \
    "dec %%rcx\n"                           \
    "jmp 5f\n" /* mid */                    \
                                            \
    "1:\n" /* lt4 */                        \
    "dec %%eax\n"                           \
    "mov (%%rsi), %%r8\n"                   \
    "jnz 2f\n" /* 2 */                      \
    ADCSBB " (%%rdx), %%r8\n"               \
    "mov %%r8, (%%rdi)\n"                   \
    "adc %%eax, %%eax\n"                    \
    "jmp 7f\n" /* exit */                   \
                                            \
    "2:\n" /* 2 */                          \
    "dec %%eax\n"                           \
    "mov 8(%%rsi), %%r9\n"                  \
    "jnz 3f\n" /* 3 */                      \
    ADCSBB " (%%rdx), %%r8\n"               \
    ADCSBB " 8(%%rdx), %%r9\n"              \
    "mov %%r8, (%%rdi)\n"                   \
    "mov %%r9, 8(%%rdi)\n"                  \
    "adc %%eax, %%eax\n"                    \
    "jmp 7f\n" /* exit */                   \
                                            \
    "3:\n" /* 3 */                          \
    "mov 16(%%rsi), %%r10\n"                \
    ADCSBB " (%%rdx), %%r8\n"               \
    ADCSBB " 8(%%rdx), %%r9\n"              \
    ADCSBB " 16(%%rdx), %%r10\n"            \
    "mov %%r8, (%%rdi)\n"                   \
    "mov %%r9, 8(%%rdi)\n"                  \
    "mov %%r10, 16(%%rdi)\n"                \
    "setc %%al\n"                           \
    "jmp 7f\n" /* exit */                   \
                                            \
    ".align 16\n"                           \
    "4:\n" /* top */                        \
    ADCSBB " (%%rdx), %%r8\n"               \
    ADCSBB " 8(%%rdx), %%r9\n"              \
    ADCSBB " 16(%%rdx), %%r10\n"            \
    ADCSBB " 24(%%rdx), %%r11\n"            \
    "mov %%r8, (%%rdi)\n"                   \
    "lea 32(%%rsi), %%rsi\n"                \
    "mov %%r9, 8(%%rdi)\n"                  \
    "mov %%r10, 16(%%rdi)\n"                \
    "dec %%rcx\n"                           \
    "mov %%r11, 24(%%rdi)\n"                \
    "lea 32(%%rdx), %%rdx\n"                \
    "mov (%%rsi), %%r8\n"                   \
    "mov 8(%%rsi), %%r9\n"                  \
    "lea 32(%%rdi), %%rdi\n"                \
    "5:\n" /* mid */                        \
    "mov 16(%%rsi), %%r10\n"                \
    "mov 24(%%rsi), %%r11\n"                \
    "jnz 4b\n" /* top */                    \
                                            \
    "6:\n" /* end */                        \
    "lea 32(%%rsi), %%rsi\n"                \
    ADCSBB " (%%rdx), %%r8\n"               \
    ADCSBB " 8(%%rdx), %%r9\n"              \
    ADCSBB " 16(%%rdx), %%r10\n"            \
    ADCSBB " 24(%%rdx), %%r11\n"            \
    "lea 32(%%rdx), %%rdx\n"                \
    "mov %%r8, (%%rdi)\n"                   \
    "mov %%r9, 8(%%rdi)\n"                  \
    "mov %%r10, 16(%%rdi)\n"                \
    "mov %%r11, 24(%%rdi)\n"                \
    "lea 32(%%rdi), %%rdi\n"                \
                                            \
    "inc %%eax\n"                           \
    "dec %%eax\n"                           \
    "jnz 1b\n" /* lt4 */                    \
    "adc %%eax, %%eax\n"                    \
    "7:\n" /* exit */                       \
    "movq $0, %q0\n"                        \
    "movb %%al, %q0\n"                      \
    : "=m" (cy),                            \
      "+D" (rp), "+S" (ap),                 \
      "+d" (bp), "+c" (n)                   \
    :                                       \
    : "al", "eax", "ebx",                   \
      "rax", "rbx", "r8", "r9",             \
      "r10", "r11", "cc", "memory"          \
  );                                        \

/*
 * Addition
 */

mp_limb_t
mpn_add_1(mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b) {
#ifdef TORSION_USE_ASM
  AORS_1("add", "adc")
  return b;
#else
  mp_size_t i;

  assert(n > 0);

  i = 0;

  do {
    mp_limb_t r = ap[i] + b;
    /* Carry out */
    b = (r < b);
    rp[i] = r;
  } while (++i < n);

  return b;
#endif
}

mp_limb_t
mpn_add_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n) {
#ifdef TORSION_USE_ASM
  mp_limb_t cy;
  AORS_N("adc")
  return cy;
#else
  mp_size_t i;
  mp_limb_t cy;

  for (i = 0, cy = 0; i < n; i++) {
    mp_limb_t a, b, r;
    a = ap[i]; b = bp[i];
    r = a + cy;
    cy = (r < cy);
    r += b;
    cy += (r < b);
    rp[i] = r;
  }

  return cy;
#endif
}

mp_limb_t
mpn_add(mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn) {
  mp_limb_t cy;

  assert(an >= bn);

  cy = mpn_add_n(rp, ap, bp, bn);

  if (an > bn)
    cy = mpn_add_1(rp + bn, ap + bn, an - bn, cy);

  return cy;
}

/*
 * Subtraction
 */

mp_limb_t
mpn_sub_1(mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b) {
#ifdef TORSION_USE_ASM
  AORS_1("sub", "sbb")
  return b;
#else
  mp_size_t i;

  assert(n > 0);

  i = 0;

  do {
    mp_limb_t a = ap[i];
    /* Carry out */
    mp_limb_t cy = a < b;
    rp[i] = a - b;
    b = cy;
  } while (++i < n);

  return b;
#endif
}

mp_limb_t
mpn_sub_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n) {
#ifdef TORSION_USE_ASM
  mp_limb_t cy;
  AORS_N("sbb")
  return cy;
#else
  mp_size_t i;
  mp_limb_t cy;

  for (i = 0, cy = 0; i < n; i++) {
    mp_limb_t a, b;
    a = ap[i]; b = bp[i];
    b += cy;
    cy = (b < cy);
    cy += (a < b);
    rp[i] = a - b;
  }

  return cy;
#endif
}

mp_limb_t
mpn_sub(mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn) {
  mp_limb_t cy;

  assert(an >= bn);

  cy = mpn_sub_n(rp, ap, bp, bn);

  if (an > bn)
    cy = mpn_sub_1(rp + bn, ap + bn, an - bn, cy);

  return cy;
}

/*
 * High-level Addition/Subtraction
 */

static mp_size_t
mpn_abs_add(mp_ptr rp, mp_srcptr ap, mp_size_t as, mp_srcptr bp, mp_size_t bs) {
  mp_size_t an = MPI_ABS(as);
  mp_size_t bn = MPI_ABS(bs);
  mp_limb_t cy;

  if (an < bn)
    MPN_SRCPTR_SWAP(ap, an, bp, bn);

  cy = mpn_add(rp, ap, an, bp, bn);

  rp[an] = cy;

  return an + cy;
}

static mp_size_t
mpn_abs_sub(mp_ptr rp, mp_srcptr ap, mp_size_t as, mp_srcptr bp, mp_size_t bs) {
  mp_size_t an = MPI_ABS(as);
  mp_size_t bn = MPI_ABS(bs);
  int cmp;

  cmp = mpn_cmp4(ap, an, bp, bn);

  if (cmp > 0) {
    mpi_assert_nocarry(mpn_sub(rp, ap, an, bp, bn));
    return mpn_normalized_size(rp, an);
  } else if (cmp < 0) {
    mpi_assert_nocarry(mpn_sub(rp, bp, bn, ap, an));
    return -mpn_normalized_size(rp, bn);
  } else {
    return 0;
  }
}

/*
 * Multiplication
 */

mp_limb_t
mpn_mul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl) {
#ifdef TORSION_USE_ASM
  /* From:
   * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/mul_1.asm
   *
   * Registers:
   *
   *   %rdi = rp (rcx)
   *   %rsi = up (rdx)
   *   %rdx = n_param (r8)
   *   %rcx = vl (r9)
   *   %r11 = n
   */
  mp_limb_t cy;

  __asm__ __volatile__(
    "xor %%r10, %%r10\n"
    "1:\n" /* common */
    "mov (%%rsi), %%rax\n" /* read first u limb early */
    "mov %%rdx, %%rbx\n" /* move away n from rdx, mul uses it */
    "mul %%rcx\n"
    "mov %%rbx, %%r11\n"

    "add %%r10, %%rax\n"
    "adc $0, %%rdx\n"

    "and $3, %%ebx\n"
    "jz 4f\n" /* b0 */
    "cmp $2, %%ebx\n"
    "jz 6f\n" /* b2 */
    "jg 5f\n" /* b3 */

    "2:\n" /* b1 */
    "dec %%r11\n"
    "jne 3f\n" /* gt1 */
    "mov %%rax, (%%rdi)\n"
    "jmp 12f\n" /* ret */
    "3:\n" /* gt1 */
    "lea 8(%%rsi,%%r11,8), %%rsi\n"
    "lea -8(%%rdi,%%r11,8), %%rdi\n"
    "neg %%r11\n"
    "xor %%r10, %%r10\n"
    "xor %%ebx, %%ebx\n"
    "mov %%rax, %%r9\n"
    "mov (%%rsi,%%r11,8), %%rax\n"
    "mov %%rdx, %%r8\n"
    "jmp 8f\n" /* L1 */

    "4:\n" /* b0 */
    "lea (%%rsi,%%r11,8), %%rsi\n"
    "lea -16(%%rdi,%%r11,8), %%rdi\n"
    "neg %%r11\n"
    "xor %%r10, %%r10\n"
    "mov %%rax, %%r8\n"
    "mov %%rdx, %%rbx\n"
    "jmp 9f\n" /* L0 */

    "5:\n" /* b3 */
    "lea -8(%%rsi,%%r11,8), %%rsi\n"
    "lea -24(%%rdi,%%r11,8), %%rdi\n"
    "neg %%r11\n"
    "mov %%rax, %%rbx\n"
    "mov %%rdx, %%r10\n"
    "jmp 10f\n" /* L3 */

    "6:\n" /* b2 */
    "lea -16(%%rsi,%%r11,8), %%rsi\n"
    "lea -32(%%rdi,%%r11,8), %%rdi\n"
    "neg %%r11\n"
    "xor %%r8, %%r8\n"
    "xor %%ebx, %%ebx\n"
    "mov %%rax, %%r10\n"
    "mov 24(%%rsi,%%r11,8), %%rax\n"
    "mov %%rdx, %%r9\n"
    "jmp 11f\n" /* L2 */

    ".align 16\n"
    "7:\n" /* top */
    "mov %%r10, (%%rdi,%%r11,8)\n"
    "add %%rax, %%r9\n"
    "mov (%%rsi,%%r11,8), %%rax\n"
    "adc %%rdx, %%r8\n"
    "mov $0, %%r10d\n"
    "8:\n" /* L1 */
    "mul %%rcx\n"
    "mov %%r9, 8(%%rdi,%%r11,8)\n"
    "add %%rax, %%r8\n"
    "adc %%rdx, %%rbx\n"
    "9:\n" /* L0 */
    "mov 8(%%rsi,%%r11,8), %%rax\n"
    "mul %%rcx\n"
    "mov %%r8, 16(%%rdi,%%r11,8)\n"
    "add %%rax, %%rbx\n"
    "adc %%rdx, %%r10\n"
    "10:\n" /* L3 */
    "mov 16(%%rsi,%%r11,8), %%rax\n"
    "mul %%rcx\n"
    "mov %%rbx, 24(%%rdi,%%r11,8)\n"
    "mov $0, %%r8d\n" /* zero */
    "mov %%r8, %%rbx\n" /* zero */
    "add %%rax, %%r10\n"
    "mov 24(%%rsi,%%r11,8), %%rax\n"
    "mov %%r8, %%r9\n" /* zero */
    "adc %%rdx, %%r9\n"
    "11:\n" /* L2 */
    "mul %%rcx\n"
    "add $4, %%r11\n"
    "js 7b\n" /* top */

    "mov %%r10, (%%rdi,%%r11,8)\n"
    "add %%rax, %%r9\n"
    "adc %%r8, %%rdx\n"
    "mov %%r9, 8(%%rdi,%%r11,8)\n"
    "add %%r8, %%rdx\n"
    "12:\n" /* ret */
    "movq %%rdx, %q0\n"
    : "=m" (cy), "+D" (rp), "+S" (up), "+d" (n)
    : "c" (vl)
    : "rax", "rbx", "ebx",
      "r8", "r9", "r10", "r11",
      "cc", "memory"
  );

  return cy;
#else
  mp_limb_t ul, cl, hpl, lpl;

  assert(n >= 1);

  cl = 0;

  do {
    ul = *up++;
    mpi_umul_ppmm(hpl, lpl, ul, vl);

    lpl += cl;
    cl = (lpl < cl) + hpl;

    *rp++ = lpl;
  } while (--n != 0);

  return cl;
#endif
}

/* From:
 * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/aorsmul_1.asm
 *
 * Registers:
 *
 *   %rdi = rp (rcx)
 *   %rsi = up (rdx)
 *   %rdx = n_param (r8)
 *   %rcx = vl (r9)
 *   %r11 = n
 */
#define AORSMUL_1(ADDSUB)                                        \
  __asm__ __volatile__(                                          \
    "mov (%%rsi), %%rax\n" /* read first u limb early */         \
    "mov %%rdx, %%rbx\n" /* move away n from rdx, mul uses it */ \
    "mul %%rcx\n"                                                \
    "mov %%rbx, %%r11\n"                                         \
                                                                 \
    "and $3, %%ebx\n"                                            \
    "jz 3f\n" /* b0 */                                           \
    "cmp $2, %%ebx\n"                                            \
    "jz 5f\n" /* b2 */                                           \
    "jg 4f\n" /* b3 */                                           \
                                                                 \
    "1:\n" /* b1 */                                              \
    "dec %%r11\n"                                                \
    "jne 2f\n" /* gt1 */                                         \
    ADDSUB " %%rax, (%%rdi)\n"                                   \
    "jmp 11f\n" /* ret */                                        \
    "2:\n" /* gt1 */                                             \
    "lea 8(%%rsi,%%r11,8), %%rsi\n"                              \
    "lea -8(%%rdi,%%r11,8), %%rdi\n"                             \
    "neg %%r11\n"                                                \
    "xor %%r10, %%r10\n"                                         \
    "xor %%ebx, %%ebx\n"                                         \
    "mov %%rax, %%r9\n"                                          \
    "mov (%%rsi,%%r11,8), %%rax\n"                               \
    "mov %%rdx, %%r8\n"                                          \
    "jmp 7f\n" /* L1 */                                          \
                                                                 \
    "3:\n" /* b0 */                                              \
    "lea (%%rsi,%%r11,8), %%rsi\n"                               \
    "lea -16(%%rdi,%%r11,8), %%rdi\n"                            \
    "neg %%r11\n"                                                \
    "xor %%r10, %%r10\n"                                         \
    "mov %%rax, %%r8\n"                                          \
    "mov %%rdx, %%rbx\n"                                         \
    "jmp 8f\n" /* L0 */                                          \
                                                                 \
    "4:\n" /* b3 */                                              \
    "lea -8(%%rsi,%%r11,8), %%rsi\n"                             \
    "lea -24(%%rdi,%%r11,8), %%rdi\n"                            \
    "neg %%r11\n"                                                \
    "mov %%rax, %%rbx\n"                                         \
    "mov %%rdx, %%r10\n"                                         \
    "jmp 9f\n" /* L3 */                                          \
                                                                 \
    "5:\n" /* b2 */                                              \
    "lea -16(%%rsi,%%r11,8), %%rsi\n"                            \
    "lea -32(%%rdi,%%r11,8), %%rdi\n"                            \
    "neg %%r11\n"                                                \
    "xor %%r8, %%r8\n"                                           \
    "xor %%ebx, %%ebx\n"                                         \
    "mov %%rax, %%r10\n"                                         \
    "mov 24(%%rsi,%%r11,8), %%rax\n"                             \
    "mov %%rdx, %%r9\n"                                          \
    "jmp 10f\n" /* L2 */                                         \
                                                                 \
    ".align 16\n"                                                \
    "6:\n" /* top */                                             \
    ADDSUB " %%r10, (%%rdi,%%r11,8)\n"                           \
    "adc %%rax, %%r9\n"                                          \
    "mov (%%rsi,%%r11,8), %%rax\n"                               \
    "adc %%rdx, %%r8\n"                                          \
    "mov $0, %%r10d\n"                                           \
    "7:\n" /* L1 */                                              \
    "mul %%rcx\n"                                                \
    ADDSUB " %%r9, 8(%%rdi,%%r11,8)\n"                           \
    "adc %%rax, %%r8\n"                                          \
    "adc %%rdx, %%rbx\n"                                         \
    "8:\n" /* L0 */                                              \
    "mov 8(%%rsi,%%r11,8), %%rax\n"                              \
    "mul %%rcx\n"                                                \
    ADDSUB " %%r8, 16(%%rdi,%%r11,8)\n"                          \
    "adc %%rax, %%rbx\n"                                         \
    "adc %%rdx, %%r10\n"                                         \
    "9:\n" /* L3 */                                              \
    "mov 16(%%rsi,%%r11,8), %%rax\n"                             \
    "mul %%rcx\n"                                                \
    ADDSUB " %%rbx, 24(%%rdi,%%r11,8)\n"                         \
    "mov $0, %%r8d\n"                                            \
    "mov %%r8, %%rbx\n" /* zero */                               \
    "adc %%rax, %%r10\n"                                         \
    "mov 24(%%rsi,%%r11,8), %%rax\n"                             \
    "mov %%r8, %%r9\n" /* zero */                                \
    "adc %%rdx, %%r9\n"                                          \
    "10:\n" /* L2 */                                             \
    "mul %%rcx\n"                                                \
    "add $4, %%r11\n"                                            \
    "js 6b\n" /* top */                                          \
                                                                 \
    ADDSUB " %%r10, (%%rdi,%%r11,8)\n"                           \
    "adc %%rax, %%r9\n"                                          \
    "adc %%r8, %%rdx\n"                                          \
    ADDSUB " %%r9, 8(%%rdi,%%r11,8)\n"                           \
    "11:\n" /* ret */                                            \
    "adc $0, %%rdx\n"                                            \
    "movq %%rdx, %q0\n"                                          \
    : "=m" (cy), "+D" (rp), "+S" (up), "+d" (n)                  \
    : "c" (vl)                                                   \
    : "rax", "rbx", "ebx",                                       \
      "r8", "r9", "r10", "r11",                                  \
      "cc", "memory"                                             \
  );                                                             \

mp_limb_t
mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl) {
#ifdef TORSION_USE_ASM
  mp_limb_t cy;
  AORSMUL_1("add")
  return cy;
#else
  mp_limb_t ul, cl, hpl, lpl, rl;

  assert(n >= 1);

  cl = 0;

  do {
    ul = *up++;
    mpi_umul_ppmm(hpl, lpl, ul, vl);

    lpl += cl;
    cl = (lpl < cl) + hpl;

    rl = *rp;
    lpl = rl + lpl;
    cl += lpl < rl;
    *rp++ = lpl;
  } while (--n != 0);

  return cl;
#endif
}

mp_limb_t
mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl) {
#ifdef TORSION_USE_ASM
  mp_limb_t cy;
  AORSMUL_1("sub")
  return cy;
#else
  mp_limb_t ul, cl, hpl, lpl, rl;

  assert(n >= 1);

  cl = 0;

  do {
    ul = *up++;
    mpi_umul_ppmm(hpl, lpl, ul, vl);

    lpl += cl;
    cl = (lpl < cl) + hpl;

    rl = *rp;
    lpl = rl - lpl;
    cl += lpl > rl;
    *rp++ = lpl;
  } while (--n != 0);

  return cl;
#endif
}

mp_limb_t
mpn_mul(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn) {
  assert(un >= vn);
  assert(vn >= 1);
  assert(!MPI_MPN_OVERLAP_P(rp, un + vn, up, un));
  assert(!MPI_MPN_OVERLAP_P(rp, un + vn, vp, vn));

  /* We first multiply by the low order limb.
     This result can be stored, not added, to
     rp. We also avoid a loop for zeroing this
     way. */

  rp[un] = mpn_mul_1(rp, up, un, vp[0]);

  /* Now accumulate the product of up[] and
     the next higher limb from vp[]. */

  while (--vn >= 1) {
    rp += 1, vp += 1;
    rp[un] = mpn_addmul_1(rp, up, un, vp[0]);
  }

  return rp[un];
}

void
mpn_mul_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n) {
  mpn_mul(rp, ap, n, bp, n);
}

#ifdef TORSION_USE_ASM
static void
mpn_sqr_diag_addlsh1(mp_ptr rp, mp_srcptr tp, mp_srcptr up, mp_size_t n) {
  /* From:
   * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/sqr_diag_addlsh1.asm
   *
   * Registers:
   *
   *   %rdi = rp
   *   %rsi = tp
   *   %rdx = up
   *   %rcx = n
   */
  __asm__ __volatile__(
    "dec %%rcx\n"
    "shl %%rcx\n"

    "mov (%%rdx),%%rax\n"

    "lea (%%rdi,%%rcx,8),%%rdi\n"
    "lea (%%rsi,%%rcx,8),%%rsi\n"
    "lea (%%rdx,%%rcx,4),%%r11\n"
    "neg %%rcx\n"

    "mul %%rax\n"
    "mov %%rax,(%%rdi,%%rcx,8)\n"

    "xor %%ebx,%%ebx\n"
    "jmp 2f\n" /* mid */

    ".align 16\n"
    "1:\n" /* top */
    "add %%r10,%%r8\n"
    "adc %%rax,%%r9\n"
    "mov %%r8,-8(%%rdi,%%rcx,8)\n"
    "mov %%r9,(%%rdi,%%rcx,8)\n"
    "2:\n" /* mid */
    "mov 8(%%r11,%%rcx,4),%%rax\n"
    "mov (%%rsi,%%rcx,8),%%r8\n"
    "mov 8(%%rsi,%%rcx,8),%%r9\n"
    "adc %%r8,%%r8\n"
    "adc %%r9,%%r9\n"
    "lea (%%rdx,%%rbx,1),%%r10\n"
    "setc %%bl\n"
    "mul %%rax\n"
    "add $2,%%rcx\n"
    "js 1b\n" /* top */

    "3:\n" /* end */
    "add %%r10,%%r8\n"
    "adc %%rax,%%r9\n"
    "mov %%r8,-8(%%rdi)\n"
    "mov %%r9,(%%rdi)\n"
    "adc %%rbx,%%rdx\n"
    "mov %%rdx,8(%%rdi)\n"
    : "+D" (rp), "+S" (tp),
      "+d" (up), "+c" (n)
    :
    : "rax", "rbx", "ebx", "bl",
      "r8", "r9", "r10", "r11",
      "cc", "memory"
  );
}
#endif

void
mpn_sqr(mp_ptr rp, mp_srcptr up, mp_size_t n) {
#ifdef TORSION_USE_ASM
  /* https://gmplib.org/repo/gmp-6.2/file/tip/mpn/generic/sqr_basecase.c */
  assert(n >= 1);
  assert(!MPI_MPN_OVERLAP_P(rp, 2 * n, up, n));

  if (n == 1) {
    mp_limb_t ul, lpl;
    ul = up[0];
    mpi_umul_ppmm(rp[1], lpl, ul, ul);
    rp[0] = lpl;
  } else {
    mp_size_t i;
    mp_ptr xp;

    rp += 1;
    rp[n - 1] = mpn_mul_1(rp, up + 1, n - 1, up[0]);

    for (i = n - 2; i != 0; i--) {
      up += 1;
      rp += 2;
      rp[i] = mpn_addmul_1(rp, up + 1, i, up[0]);
    }

    xp = rp - 2 * n + 3;
    mpn_sqr_diag_addlsh1(xp, xp + 1, up - n + 2, n);
  }
#else
  mpn_mul_n(rp, up, up, n);
#endif
}

void
mpn_mont(mp_ptr zp,
         mp_srcptr xp,
         mp_srcptr yp,
         mp_srcptr mp,
         mp_limb_t k,
         mp_size_t n) {
  /* Montgomery multiplication. */
  /* Guarantees zn <= mn, but not necessarily z < m. */
  mp_limb_t c2, c3, cx, cy;
  mp_limb_t c1 = 0;
  mp_size_t i;

  mpn_zero(zp, n * 2);

  for (i = 0; i < n; i++) {
    c2 = mpn_addmul_1(zp + i, xp, n, yp[i]);
    c3 = mpn_addmul_1(zp + i, mp, n, zp[i] * k);
    cx = c1 + c2;
    cy = cx + c3;
    zp[n + i] = cy;
    c1 = (cx < c2) | (cy < c3);
  }

  mpn_sub_n(zp, zp + n, mp, n);
  mpn_cnd_select(c1 != 0, zp, zp + n, zp, n);
}

/*
 * Division Engine
 */

static mp_limb_t
mpn_invert_3by2(mp_limb_t u1, mp_limb_t u0) {
  /* The 3/2 inverse is defined as
   *
   *   m = floor((B^3 - 1) / (B u1 + u0)) - B
   */
  mp_limb_t r, m;

  {
    mp_limb_t p, ql;
    unsigned ul, uh, qh;

    /* For notation, let b denote the half-limb base,
       so that B = b^2. Split u1 = b uh + ul. */
    ul = u1 & MPI_LLIMB_MASK;
    uh = u1 >> (MPI_LIMB_BITS / 2);

    /* Approximation of the high half of quotient.
       Differs from the 2/1 inverse of the half limb
       uh, since we have already subtracted u0. */
    qh = (u1 ^ MPI_LIMB_MAX) / uh;

    /* Adjust to get a half-limb 3/2 inverse,
     * i.e., we want
     *
     *   qh' = floor((b^3 - 1) / u) - b
     *       = floor((b^3 - b u - 1) / u
     *       = floor((b (~u) + b-1) / u),
     *
     * and the remainder
     *
     *   r = b (~u) + b-1 - qh (b uh + ul)
     *     = b (~u - qh uh) + b-1 - qh ul
     *
     * Subtraction of qh ul may underflow, which
     * implies adjustments.  But by normalization,
     * 2 u >= B > qh ul, so we need to adjust by
     * at most 2.
     */
    r = ((~u1 - (mp_limb_t)qh * uh) << (MPI_LIMB_BITS / 2)) | MPI_LLIMB_MASK;

    p = (mp_limb_t)qh * ul;

    /* Adjustment steps taken from udiv_qrnnd_c */
    if (r < p) {
      qh--;
      r += u1;

      if (r >= u1) { /* i.e. we didn't get carry when adding to r */
        if (r < p) {
          qh--;
          r += u1;
        }
      }
    }

    r -= p;

    /* Low half of the quotient is
     *
     *   ql = floor((b r + b-1) / u1).
     *
     * This is a 3/2 division (on half-limbs),
     * for which qh is a suitable inverse.
     */
    p = (r >> (MPI_LIMB_BITS / 2)) * qh + r;

    /* Unlike full-limb 3/2, we can add 1 without
       overflow. For this to work, it is essential
       that ql is a full mp_limb_t. */
    ql = (p >> (MPI_LIMB_BITS / 2)) + 1;

    /* By the 3/2 trick, we don't need the high half limb. */
    r = (r << (MPI_LIMB_BITS / 2)) + MPI_LLIMB_MASK - ql * u1;

    if (r >= (MPI_LIMB_MAX & (p << (MPI_LIMB_BITS / 2)))) {
      ql--;
      r += u1;
    }

    m = ((mp_limb_t)qh << (MPI_LIMB_BITS / 2)) + ql;

    if (r >= u1) {
      m++;
      r -= u1;
    }
  }

  /* Now m is the 2/1 inverse of u1. If u0 > 0,
     adjust it to become a 3/2 inverse. */
  if (u0 > 0) {
    mp_limb_t th, tl;

    r = ~r;
    r += u0;

    if (r < u0) {
      m--;

      if (r >= u1) {
        m--;
        r -= u1;
      }

      r -= u1;
    }

    mpi_umul_ppmm(th, tl, u0, m);

    r += th;

    if (r < th) {
      m--;
      m -= ((r > u1) | ((r == u1) & (tl > u0)));
    }
  }

  return m;
}

static void
mpn_div_qr_1_invert(struct mpi_div_inverse *inv, mp_limb_t d) {
  unsigned shift;

  assert(d > 0);

  mpi_clz(shift, d);

  inv->shift = shift;
  inv->d1 = d << shift;
  inv->di = mpn_invert_limb(inv->d1);
}

static void
mpn_div_qr_2_invert(struct mpi_div_inverse *inv,
                    mp_limb_t d1, mp_limb_t d0) {
  unsigned shift;

  assert(d1 > 0);

  mpi_clz(shift, d1);

  inv->shift = shift;

  if (shift > 0) {
    d1 = (d1 << shift) | (d0 >> (MPI_LIMB_BITS - shift));
    d0 <<= shift;
  }

  inv->d1 = d1;
  inv->d0 = d0;
  inv->di = mpn_invert_3by2(d1, d0);
}

static void
mpn_div_qr_invert(struct mpi_div_inverse *inv,
                  mp_srcptr dp, mp_size_t dn) {
  assert(dn > 0);

  if (dn == 1) {
    mpn_div_qr_1_invert(inv, dp[0]);
  } else if (dn == 2) {
    mpn_div_qr_2_invert(inv, dp[1], dp[0]);
  } else {
    unsigned shift;
    mp_limb_t d1, d0;

    d1 = dp[dn - 1];
    d0 = dp[dn - 2];

    assert(d1 > 0);

    mpi_clz(shift, d1);

    inv->shift = shift;

    if (shift > 0) {
      d1 = (d1 << shift) | (d0 >> (MPI_LIMB_BITS - shift));
      d0 = (d0 << shift) | (dp[dn - 3] >> (MPI_LIMB_BITS - shift));
    }

    inv->d1 = d1;
    inv->d0 = d0;
    inv->di = mpn_invert_3by2(d1, d0);
  }
}

static mp_limb_t
mpn_div_qr_1_preinv(mp_ptr qp, mp_srcptr np, mp_size_t nn,
                    const struct mpi_div_inverse *inv) {
  mp_limb_t d, di;
  mp_limb_t r;
  mp_ptr tp = NULL;

  if (inv->shift > 0) {
    /* Shift, reusing qp area if possible. */
    /* In-place shift if qp == np. */
    tp = qp ? qp : mpi_xalloc_limbs(nn);
    r = mpn_lshift(tp, np, nn, inv->shift);
    np = tp;
  } else {
    r = 0;
  }

  d = inv->d1;
  di = inv->di;

  while (--nn >= 0) {
    mp_limb_t q;

    mpi_udiv_qrnnd_preinv(q, r, r, np[nn], d, di);

    if (qp)
      qp[nn] = q;
  }

  if ((inv->shift > 0) && (tp != qp))
    mpi_free(tp);

  return r >> inv->shift;
}

static void
mpn_div_qr_2_preinv(mp_ptr qp, mp_ptr np, mp_size_t nn,
                    const struct mpi_div_inverse *inv) {
  unsigned shift;
  mp_size_t i;
  mp_limb_t d1, d0, di, r1, r0;

  assert(nn >= 2);

  shift = inv->shift;
  d1 = inv->d1;
  d0 = inv->d0;
  di = inv->di;

  if (shift > 0)
    r1 = mpn_lshift(np, np, nn, shift);
  else
    r1 = 0;

  r0 = np[nn - 1];
  i = nn - 2;

  do {
    mp_limb_t n0, q;

    n0 = np[i];
    mpi_udiv_qr_3by2(q, r1, r0, r1, r0, n0, d1, d0, di);

    if (qp)
      qp[i] = q;
  } while (--i >= 0);

  if (shift > 0) {
    assert((r0 & (MPI_LIMB_MAX >> (MPI_LIMB_BITS - shift))) == 0);
    r0 = (r0 >> shift) | (r1 << (MPI_LIMB_BITS - shift));
    r1 >>= shift;
  }

  np[1] = r1;
  np[0] = r0;
}

static void
mpn_div_qr_pi1(mp_ptr qp,
               mp_ptr np, mp_size_t nn, mp_limb_t n1,
               mp_srcptr dp, mp_size_t dn,
               mp_limb_t dinv) {
  mp_size_t i;

  mp_limb_t d1, d0;
  mp_limb_t cy, cy1;
  mp_limb_t q;

  assert(dn > 2);
  assert(nn >= dn);

  d1 = dp[dn - 1];
  d0 = dp[dn - 2];

  assert((d1 & MPI_LIMB_HIGHBIT) != 0);

  /* Iteration variable is the index of the q limb.
   *
   * We divide <n1, np[dn-1+i], np[dn-2+i], np[dn-3+i],..., np[i]>
   * by            <d1,          d0,        dp[dn-3],  ..., dp[0]>
   */
  i = nn - dn;

  do {
    mp_limb_t n0 = np[dn - 1 + i];

    if (n1 == d1 && n0 == d0) {
      q = MPI_LIMB_MAX;
      mpn_submul_1(np + i, dp, dn, q);
      n1 = np[dn - 1 + i]; /* Update n1, last loop's value is now invalid. */
    } else {
      mpi_udiv_qr_3by2(q, n1, n0, n1, n0, np[dn - 2 + i], d1, d0, dinv);

      cy = mpn_submul_1(np + i, dp, dn - 2, q);

      cy1 = n0 < cy;
      n0 = n0 - cy;
      cy = n1 < cy1;
      n1 = n1 - cy1;
      np[dn - 2 + i] = n0;

      if (cy != 0) {
        n1 += d1 + mpn_add_n(np + i, np + i, dp, dn - 1);
        q--;
      }
    }

    if (qp)
      qp[i] = q;
  } while (--i >= 0);

  np[dn - 1] = n1;
}

static void
mpn_div_qr_preinv(mp_ptr qp, mp_ptr np, mp_size_t nn,
                  mp_srcptr dp, mp_size_t dn,
                  const struct mpi_div_inverse *inv) {
  assert(dn > 0);
  assert(nn >= dn);

  if (dn == 1) {
    np[0] = mpn_div_qr_1_preinv(qp, np, nn, inv);
  } else if (dn == 2) {
    mpn_div_qr_2_preinv(qp, np, nn, inv);
  } else {
    mp_limb_t nh;
    unsigned shift;

    assert(inv->d1 == dp[dn - 1]);
    assert(inv->d0 == dp[dn - 2]);
    assert((inv->d1 & MPI_LIMB_HIGHBIT) != 0);

    shift = inv->shift;

    if (shift > 0)
      nh = mpn_lshift(np, np, nn, shift);
    else
      nh = 0;

    mpn_div_qr_pi1(qp, np, nn, nh, dp, dn, inv->di);

    if (shift > 0)
      mpi_assert_nocarry(mpn_rshift(np, np, dn, shift));
  }
}

static void
mpn_div_qr(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn) {
  struct mpi_div_inverse inv;
  mp_ptr tp = NULL;

  assert(dn > 0);
  assert(nn >= dn);

  mpn_div_qr_invert(&inv, dp, dn);

  if (dn > 2 && inv.shift > 0) {
    tp = mpi_xalloc_limbs(dn);
    mpi_assert_nocarry(mpn_lshift(tp, dp, dn, inv.shift));
    dp = tp;
  }

  mpn_div_qr_preinv(qp, np, nn, dp, dn, &inv);

  if (tp)
    mpi_free(tp);
}

static void
mpn_div_tp(mp_ptr qp,
           mp_ptr np,
           mp_size_t nn,
           mp_srcptr dp,
           mp_size_t dn,
           mp_ptr tp) {
  struct mpi_div_inverse inv;

  assert(dn > 0);
  assert(nn >= dn);

  mpn_div_qr_invert(&inv, dp, dn);

  if (dn > 2 && inv.shift > 0) {
    mpi_assert_nocarry(mpn_lshift(tp, dp, dn, inv.shift));
    dp = tp;
  }

  mpn_div_qr_preinv(qp, np, nn, dp, dn, &inv);
}

/*
 * Truncation Division
 */

void
mpn_quorem(mp_ptr qp, mp_ptr rp,
           mp_srcptr np, mp_size_t nn,
           mp_srcptr dp, mp_size_t dn) {
  assert(nn >= dn);
  assert(dn > 0);
  assert(dp[dn - 1] != 0);

  if (rp == np) {
    mpn_div_qr(qp, rp, nn, dp, dn);
  } else {
    mp_ptr tp = mpi_xalloc_limbs(nn);

    mpn_copyi(tp, np, nn);
    mpn_div_qr(qp, tp, nn, dp, dn);

    if (rp)
      mpn_copyi(rp, tp, dn);

    mpi_free(tp);
  }
}

/*
 * Left Shift
 */

mp_limb_t
mpn_lshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt) {
#ifdef TORSION_USE_ASM
  /* From:
   * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/lshift.asm
   *
   *
   * Registers:
   *
   *   %rdi = rp
   *   %rsi = up
   *   %rdx = n
   *   %rcx = cnt
   */
  mp_limb_t cy;

  __asm__ __volatile__(
    "neg %%ecx\n" /* put rsh count in cl */
    "mov -8(%%rsi,%%rdx,8), %%rax\n"
    "shr %%cl, %%rax\n" /* function return value */

    "neg %%ecx\n" /* put lsh count in cl */
    "lea 1(%%rdx), %%r8d\n"
    "and $3, %%r8d\n"
    "je 4f\n" /* (rlx) jump for n = 3, 7, 11, ... */

    "dec %%r8d\n"
    "jne 1f\n" /* 1 */
    /* n = 4, 8, 12, ... */
    "mov -8(%%rsi,%%rdx,8), %%r10\n"
    "shl %%cl, %%r10\n"
    "neg %%ecx\n" /* put rsh count in cl */
    "mov -16(%%rsi,%%rdx,8), %%r8\n"
    "shr %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "mov %%r10, -8(%%rdi,%%rdx,8)\n"
    "dec %%rdx\n"
    "jmp 3f\n" /* rll */

    "1:\n" /* 1 */
    "dec %%r8d\n"
    "je 2f\n" /* (1x) jump for n = 1, 5, 9, 13, ... */
    /* n = 2, 6, 10, 16, ... */
    "mov -8(%%rsi,%%rdx,8), %%r10\n"
    "shl %%cl, %%r10\n"
    "neg %%ecx\n" /* put rsh count in cl */
    "mov -16(%%rsi,%%rdx,8), %%r8\n"
    "shr %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "mov %%r10, -8(%%rdi,%%rdx,8)\n"
    "dec %%rdx\n"
    "neg %%ecx\n" /* put lsh count in cl */
    "2:\n" /* 1x */
    "cmp $1, %%rdx\n"
    "je 7f\n" /* ast */
    "mov -8(%%rsi,%%rdx,8), %%r10\n"
    "shl %%cl, %%r10\n"
    "mov -16(%%rsi,%%rdx,8), %%r11\n"
    "shl %%cl, %%r11\n"
    "neg %%ecx\n" /* put rsh count in cl */
    "mov -16(%%rsi,%%rdx,8), %%r8\n"
    "mov -24(%%rsi,%%rdx,8), %%r9\n"
    "shr %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "shr %%cl, %%r9\n"
    "or %%r9, %%r11\n"
    "mov %%r10, -8(%%rdi,%%rdx,8)\n"
    "mov %%r11, -16(%%rdi,%%rdx,8)\n"
    "sub $2, %%rdx\n"

    "3:\n" /* rll */
    "neg %%ecx\n" /* put lsh count in cl */
    "4:\n" /* rlx */
    "mov -8(%%rsi,%%rdx,8), %%r10\n"
    "shl %%cl, %%r10\n"
    "mov -16(%%rsi,%%rdx,8), %%r11\n"
    "shl %%cl, %%r11\n"

    "sub $4, %%rdx\n" /* 4 */
    "jb 6f\n" /* (end) 2 */
    ".align 16\n"
    "5:\n" /* top */
    /* finish stuff from lsh block */
    "neg %%ecx\n" /* put rsh count in cl */
    "mov 16(%%rsi,%%rdx,8), %%r8\n"
    "mov 8(%%rsi,%%rdx,8), %%r9\n"
    "shr %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "shr %%cl, %%r9\n"
    "or %%r9, %%r11\n"
    "mov %%r10, 24(%%rdi,%%rdx,8)\n"
    "mov %%r11, 16(%%rdi,%%rdx,8)\n"
    /* start two new rsh */
    "mov 0(%%rsi,%%rdx,8), %%r8\n"
    "mov -8(%%rsi,%%rdx,8), %%r9\n"
    "shr %%cl, %%r8\n"
    "shr %%cl, %%r9\n"

    /* finish stuff from rsh block */
    "neg %%ecx\n" /* put lsh count in cl */
    "mov 8(%%rsi,%%rdx,8), %%r10\n"
    "mov 0(%%rsi,%%rdx,8), %%r11\n"
    "shl %%cl, %%r10\n"
    "or %%r10, %%r8\n"
    "shl %%cl, %%r11\n"
    "or %%r11, %%r9\n"
    "mov %%r8, 8(%%rdi,%%rdx,8)\n"
    "mov %%r9, 0(%%rdi,%%rdx,8)\n"
    /* start two new lsh */
    "mov -8(%%rsi,%%rdx,8), %%r10\n"
    "mov -16(%%rsi,%%rdx,8), %%r11\n"
    "shl %%cl, %%r10\n"
    "shl %%cl, %%r11\n"

    "sub $4, %%rdx\n"
    "jae 5b\n" /* (top) 2 */
    "6:\n" /* end */
    "neg %%ecx\n" /* put rsh count in cl */
    "mov 8(%%rsi), %%r8\n"
    "shr %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "mov (%%rsi), %%r9\n"
    "shr %%cl, %%r9\n"
    "or %%r9, %%r11\n"
    "mov %%r10, 16(%%rdi)\n"
    "mov %%r11, 8(%%rdi)\n"

    "neg %%ecx\n" /* put lsh count in cl */
    "7:\n" /* ast */
    "mov (%%rsi), %%r10\n"
    "shl %%cl, %%r10\n"
    "mov %%r10, (%%rdi)\n"
    "movq %%rax, %q0\n"
    : "=m" (cy), "+d" (n), "+c" (cnt)
    : "D" (rp), "S" (up)
    : "rax", "r8", "r9", "r10", "r11",
      "cc", "memory"
  );

  return cy;
#else
  mp_limb_t high_limb, low_limb;
  unsigned int tnc;
  mp_limb_t retval;

  assert(n >= 1);
  assert(cnt >= 1);
  assert(cnt < MPI_LIMB_BITS);

  up += n;
  rp += n;

  tnc = MPI_LIMB_BITS - cnt;
  low_limb = *--up;
  retval = low_limb >> tnc;
  high_limb = (low_limb << cnt);

  while (--n != 0) {
    low_limb = *--up;
    *--rp = high_limb | (low_limb >> tnc);
    high_limb = (low_limb << cnt);
  }

  *--rp = high_limb;

  return retval;
#endif
}

/*
 * Right Shift
 */

mp_limb_t
mpn_rshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt) {
#ifdef TORSION_USE_ASM
  /* From:
   * https://gmplib.org/repo/gmp-6.2/file/tip/mpn/x86_64/rshift.asm
   *
   * Registers:
   *
   *   %rdi = rp
   *   %rsi = up
   *   %rdx = n
   *   %rcx = cnt
   */
  mp_limb_t cy;

  __asm__ __volatile__(
    "neg %%ecx\n" /* put rsh count in cl */
    "mov (%%rsi), %%rax\n"
    "shl %%cl, %%rax\n" /* function return value */
    "neg %%ecx\n" /* put lsh count in cl */

    "lea 1(%%rdx), %%r8d\n"

    "lea -8(%%rsi,%%rdx,8), %%rsi\n"
    "lea -8(%%rdi,%%rdx,8), %%rdi\n"
    "neg %%rdx\n"

    "and $3, %%r8d\n"
    "je 4f\n" /* (rlx) jump for n = 3, 7, 11, ... */

    "dec %%r8d\n"
    "jne 1f\n" /* 1 */
    /* n = 4, 8, 12, ... */
    "mov 8(%%rsi,%%rdx,8), %%r10\n"
    "shr %%cl, %%r10\n"
    "neg %%ecx\n" /* put rsh count in cl */
    "mov 16(%%rsi,%%rdx,8), %%r8\n"
    "shl %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "mov %%r10, 8(%%rdi,%%rdx,8)\n"
    "inc %%rdx\n"
    "jmp 3f\n" /* rll */

    "1:\n" /* 1 */
    "dec %%r8d\n"
    "je 2f\n" /* (1x) jump for n = 1, 5, 9, 13, ... */
    /* n = 2, 6, 10, 16, ... */
    "mov 8(%%rsi,%%rdx,8), %%r10\n"
    "shr %%cl, %%r10\n"
    "neg %%ecx\n" /* put rsh count in cl */
    "mov 16(%%rsi,%%rdx,8), %%r8\n"
    "shl %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "mov %%r10, 8(%%rdi,%%rdx,8)\n"
    "inc %%rdx\n"
    "neg %%ecx\n" /* put lsh count in cl */
    "2:\n" /* 1x */
    "cmp $-1, %%rdx\n"
    "je 7f\n" /* ast */
    "mov 8(%%rsi,%%rdx,8), %%r10\n"
    "shr %%cl, %%r10\n"
    "mov 16(%%rsi,%%rdx,8), %%r11\n"
    "shr %%cl, %%r11\n"
    "neg %%ecx\n" /* put rsh count in cl */
    "mov 16(%%rsi,%%rdx,8), %%r8\n"
    "mov 24(%%rsi,%%rdx,8), %%r9\n"
    "shl %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "shl %%cl, %%r9\n"
    "or %%r9, %%r11\n"
    "mov %%r10, 8(%%rdi,%%rdx,8)\n"
    "mov %%r11, 16(%%rdi,%%rdx,8)\n"
    "add $2, %%rdx\n"

    "3:\n" /* rll */
    "neg %%ecx\n" /* put lsh count in cl */
    "4:\n" /* rlx */
    "mov 8(%%rsi,%%rdx,8), %%r10\n"
    "shr %%cl, %%r10\n"
    "mov 16(%%rsi,%%rdx,8), %%r11\n"
    "shr %%cl, %%r11\n"

    "add $4, %%rdx\n" /* 4 */
    "jb 6f\n" /* (end) 2 */
    ".align 16\n"
    "5:\n" /* top */
    /* finish stuff from lsh block */
    "neg %%ecx\n" /* put rsh count in cl */
    "mov -16(%%rsi,%%rdx,8), %%r8\n"
    "mov -8(%%rsi,%%rdx,8), %%r9\n"
    "shl %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "shl %%cl, %%r9\n"
    "or %%r9, %%r11\n"
    "mov %%r10, -24(%%rdi,%%rdx,8)\n"
    "mov %%r11, -16(%%rdi,%%rdx,8)\n"
    /* start two new rsh */
    "mov (%%rsi,%%rdx,8), %%r8\n"
    "mov 8(%%rsi,%%rdx,8), %%r9\n"
    "shl %%cl, %%r8\n"
    "shl %%cl, %%r9\n"

    /* finish stuff from rsh block */
    "neg %%ecx\n" /* put lsh count in cl */
    "mov -8(%%rsi,%%rdx,8), %%r10\n"
    "mov 0(%%rsi,%%rdx,8), %%r11\n"
    "shr %%cl, %%r10\n"
    "or %%r10, %%r8\n"
    "shr %%cl, %%r11\n"
    "or %%r11, %%r9\n"
    "mov %%r8, -8(%%rdi,%%rdx,8)\n"
    "mov %%r9, 0(%%rdi,%%rdx,8)\n"
    /* start two new lsh */
    "mov 8(%%rsi,%%rdx,8), %%r10\n"
    "mov 16(%%rsi,%%rdx,8), %%r11\n"
    "shr %%cl, %%r10\n"
    "shr %%cl, %%r11\n"

    "add $4, %%rdx\n"
    "jae 5b\n" /* (top) 2 */
    "6:\n" /* end */
    "neg %%ecx\n" /* put rsh count in cl */
    "mov -8(%%rsi), %%r8\n"
    "shl %%cl, %%r8\n"
    "or %%r8, %%r10\n"
    "mov (%%rsi), %%r9\n"
    "shl %%cl, %%r9\n"
    "or %%r9, %%r11\n"
    "mov %%r10, -16(%%rdi)\n"
    "mov %%r11, -8(%%rdi)\n"

    "neg %%ecx\n" /* put lsh count in cl */
    "7:\n" /* ast */
    "mov (%%rsi), %%r10\n"
    "shr %%cl, %%r10\n"
    "mov %%r10, (%%rdi)\n"
    "movq %%rax, %q0\n"
    : "=m" (cy),
      "+D" (rp), "+S" (up),
      "+d" (n), "+c" (cnt)
    :
    : "rax", "r8", "r9", "r10", "r11",
      "cc", "memory"
  );

  return cy;
#else
  mp_limb_t high_limb, low_limb;
  unsigned int tnc;
  mp_limb_t retval;

  assert(n >= 1);
  assert(cnt >= 1);
  assert(cnt < MPI_LIMB_BITS);

  tnc = MPI_LIMB_BITS - cnt;
  high_limb = *up++;
  retval = (high_limb << tnc);
  low_limb = high_limb >> cnt;

  while (--n != 0) {
    high_limb = *up++;
    *rp++ = low_limb | (high_limb << tnc);
    low_limb = high_limb >> cnt;
  }

  *rp = low_limb;

  return retval;
#endif
}

/*
 * Bit Manipulation
 */

mp_limb_t
mpn_get_bit(mp_srcptr xp, mp_size_t xn, mp_bitcnt_t pos) {
  mp_size_t index = pos / MPI_LIMB_BITS;

  if (index >= xn)
    return 0;

  return (xp[index] >> (pos % MPI_LIMB_BITS)) & 1;
}

mp_limb_t
mpn_get_bits(mp_srcptr xp, mp_size_t xn, mp_bitcnt_t pos, mp_bitcnt_t width) {
  mp_size_t index = pos / MPI_LIMB_BITS;
  mp_size_t shift;
  mp_limb_t bits;

  if (index >= xn)
    return 0;

  shift = pos % MPI_LIMB_BITS;
  bits = (xp[index] >> shift) & (((mp_limb_t)1 << width) - 1);

  if (shift + width > MPI_LIMB_BITS && index + 1 < xn) {
    mp_size_t more = shift + width - MPI_LIMB_BITS;
    mp_limb_t next = xp[index + 1] & (((mp_limb_t)1 << more) - 1);

    bits |= next << (MPI_LIMB_BITS - shift);
  }

  return bits;
}

void
mpn_set_bit(mp_ptr xp, mp_size_t xn, mp_bitcnt_t pos) {
  mp_size_t index = pos / MPI_LIMB_BITS;
  mp_size_t shift = pos % MPI_LIMB_BITS;

  assert(index < xn);

  xp[index] |= ((mp_limb_t)1 << shift);
}

void
mpn_clr_bit(mp_ptr xp, mp_size_t xn, mp_bitcnt_t pos) {
  mp_size_t index = pos / MPI_LIMB_BITS;
  mp_size_t shift = pos % MPI_LIMB_BITS;

  assert(index < xn);

  xp[index] &= ~((mp_limb_t)1 << shift);
}

/*
 * Number Theoretic Functions
 */

static mp_limb_t
mpn_gcd_11(mp_limb_t u, mp_limb_t v) {
  unsigned shift;

  assert((u | v) > 0);

  if (u == 0)
    return v;
  else if (v == 0)
    return u;

  mpi_ctz(shift, u | v);

  u >>= shift;
  v >>= shift;

  if ((u & 1) == 0)
    MP_LIMB_T_SWAP(u, v);

  while ((v & 1) == 0)
    v >>= 1;

  while (u != v) {
    if (u > v) {
      u -= v;
      do {
        u >>= 1;
      } while ((u & 1) == 0);
    } else {
      v -= u;
      do {
        v >>= 1;
      } while ((v & 1) == 0);
    }
  }

  return u << shift;
}

int
mpn_invert_n(mp_ptr rp, mp_srcptr xp, mp_srcptr mp, mp_size_t n, mp_ptr scratch) {
  /* Fast binary EGCD to invert over a prime field. */
  mp_ptr ap = &scratch[0 * (n + 2)];
  mp_ptr bp = &scratch[1 * (n + 2)];
  mp_ptr up = &scratch[2 * (n + 2)];
  mp_ptr vp = &scratch[3 * (n + 2)];
  mp_size_t an, bn, us, vs;
  mp_bitcnt_t shift;

  assert(n > 0);
  assert(mp[n - 1] != 0);

  if (mpn_zero_p(xp, n)) {
    mpn_zero(rp, n);
    return 0;
  }

  mpn_zero(ap + n, 2);
  mpn_zero(bp + n, 2);
  mpn_zero(up + 1, n + 1);
  mpn_zero(vp, n + 2);

  mpn_copyi(ap, xp, n);
  mpn_copyi(bp, mp, n);

  up[0] = 1;

  an = mpn_normalized_size(ap, n);
  bn = n;
  us = 1;
  vs = 0;

  while (mpn_cmp_1(ap, an, 1) > 0 && mpn_cmp_1(bp, bn, 1) > 0) {
    MPN_MAKE_ODD(shift, ap, an);

    while (shift--) {
      if (MPN_ODD_P(up, us))
        MPN_ADD(up, us, mp, n);

      MPN_RSHIFT(up, us, 1);
    }

    MPN_MAKE_ODD(shift, bp, bn);

    while (shift--) {
      if (MPN_ODD_P(vp, vs))
        MPN_ADD(vp, vs, mp, n);

      MPN_RSHIFT(vp, vs, 1);
    }

    if (mpn_cmp4(ap, an, bp, bn) >= 0) {
      MPN_SUB(ap, an, bp, bn);
      MPN_SUB(up, us, vp, vs);
    } else {
      MPN_SUB(bp, bn, ap, an);
      MPN_SUB(vp, vs, up, us);
    }
  }

  if (mpn_cmp_1(ap, an, 1) != 0) {
    up = vp;
    us = vs;
  }

  while (us < 0)
    MPN_ADD(up, us, mp, n);

  while (mpn_cmp4(up, us, mp, n) >= 0)
    MPN_SUB(up, us, mp, n);

  assert(us <= n);

  mpn_copyi(rp, up, us);
  mpn_zero(rp + us, n - us);

  return 1;
}

int
mpn_jacobi_n(mp_srcptr xp, mp_srcptr yp, mp_size_t n, mp_ptr scratch) {
  mp_ptr ap = &scratch[0 * n];
  mp_ptr bp = &scratch[1 * n];
  mp_ptr sp = &scratch[2 * n];
  mp_size_t an, bn, bits;
  mp_limb_t bmod8;
  int j = 1;

  assert(n > 0);
  assert(yp[n - 1] != 0);

  if ((yp[0] & 1) == 0)
    mpi_die("mpn_jacobi: Even argument.");

  mpn_copyi(ap, xp, n);
  mpn_copyi(bp, yp, n);

  an = mpn_normalized_size(ap, n);
  bn = n;

  for (;;) {
    if (bn == 1 && bp[0] == 1)
      break;

    if (an == 0) {
      j = 0;
      break;
    }

    MPN_MOD(ap, an, bp, bn, sp);

    if (an == 0) {
      j = 0;
      break;
    }

    MPN_MAKE_ODD(bits, ap, an);

    if (bits & 1) {
      bmod8 = bp[0] & 7;

      if (bmod8 == 3 || bmod8 == 5)
        j = -j;
    }

    if ((ap[0] & 3) == 3 && (bp[0] & 3) == 3)
      j = -j;

    MPN_PTR_SWAP(ap, an, bp, bn);
  }

  return j;
}

void
mpn_powm_sec(mp_ptr zp,
             mp_srcptr xp, mp_size_t xn,
             mp_srcptr yp, mp_size_t yn,
             mp_srcptr mp, mp_size_t mn,
             mp_ptr scratch) {
  /* Scratch Layout:
   *
   *   up = mod_limbs
   *   z1 = 2 * mod_limbs
   *   z2 = 2 * mod_limbs
   *   one = mod_limbs
   *   tmp = mod_limbs
   *   wnds = ((1 << 4) + 1) * mod_limbs
   *   total = 24 * mod_limbs
   *
   * Precomputation:
   *
   *   k = -m^-1 mod 2^limb_width
   *   rr = 2^(2 * mod_limbs) mod m
   *
   * We assume the modulus is not secret.
   */
  mp_ptr up = &scratch[0];
  mp_ptr z1 = &scratch[mn];
  mp_ptr z2 = &scratch[3 * mn];
  mp_ptr one = &scratch[5 * mn];
  mp_ptr tmp = &scratch[6 * mn];
  mp_ptr wnds = &scratch[7 * mn];
  mp_ptr rr = &wnds[3 * mn];
  mp_ptr wnd[1 << 4];
  mp_size_t yb = yn * MPI_LIMB_BITS;
  mp_size_t start = (yb + MPI_WND_WIDTH - 1) / MPI_WND_WIDTH - 1;
  mp_limb_t k, t, b, j, cy;
  mp_size_t i;

  assert(xn >= 0);
  assert(yn >= 0);
  assert(mn > 0);
  assert(xn <= mn);
  assert(mp[mn - 1] != 0);
  assert((mp[0] & 1) != 0);

  mpn_copyi(up, xp, xn);
  mpn_zero(up + xn, mn - xn);

  k = 2 - mp[0];
  t = mp[0] - 1;

  for (i = 1; i < MPI_LIMB_BITS; i <<= 1) {
    t *= t;
    k *= (t + 1);
  }

  k = -k;

  mpn_zero(rr, mn * 2);
  rr[mn * 2] = 1;
  mpn_div_qr(NULL, rr, mn * 2 + 1, mp, mn);

  one[0] = 1;
  mpn_zero(one + 1, mn - 1);

  for (i = 0; i < MPI_WND_SIZE; i++)
    wnd[i] = &wnds[i * mn];

  mpn_mont(wnd[0], one, rr, mp, k, mn);
  mpn_mont(wnd[1], up, rr, mp, k, mn);

  for (i = 2; i < MPI_WND_SIZE; i++)
    mpn_mont(wnd[i], wnd[i - 1], wnd[1], mp, k, mn);

  mpn_copyi(z1, wnd[0], mn);

  for (i = start; i >= 0; i--) {
    b = mpn_get_bits(yp, yn, i * MPI_WND_WIDTH, MPI_WND_WIDTH);

    for (j = 0; j < MPI_WND_SIZE; j++)
      mpn_cnd_select(j == b, tmp, tmp, wnd[j], mn);

    if (i == start) {
      mpn_copyi(z1, tmp, mn);
    } else {
      for (j = 0; j < MPI_WND_WIDTH; j++) {
        if (j & 1)
          mpn_mont(z1, z2, z2, mp, k, mn);
        else
          mpn_mont(z2, z1, z1, mp, k, mn);
      }

#if MPI_WND_WIDTH % 2 == 0
      MP_PTR_SWAP(z1, z2);
#endif

      mpn_mont(z1, z2, tmp, mp, k, mn);
    }
  }

  mpn_mont(z2, z1, one, mp, k, mn);

  cy = mpn_sub_n(z1, z2, mp, mn);
  mpn_cnd_select(cy == 0, zp, z2, z1, mn);
}

/*
 * Helpers
 */

mp_size_t
mpn_normalized_size(mp_srcptr xp, mp_size_t n) {
  while (n > 0 && xp[n - 1] == 0)
    --n;
  return n;
}

mp_bitcnt_t
mpn_bitlen(mp_srcptr xp, mp_size_t xn) {
  mp_size_t bits;

  xn = mpn_normalized_size(xp, xn);

  if (xn == 0)
    return 0;

  mpi_clz(bits, xp[xn - 1]);

  return xn * MPI_LIMB_BITS - bits;
}

mp_bitcnt_t
mpn_ctz(mp_srcptr xp, mp_size_t xn) {
  mp_size_t i, cnt;

  assert(xn >= 0);

  for (i = 0; i < xn; i++) {
    if (xp[i] != 0)
      break;
  }

  if (i == xn)
    return 0;

  mpi_ctz(cnt, xp[i]);

  return cnt + i * MPI_LIMB_BITS;
}

/*
 * Constant Time
 */

void
mpn_cnd_select(mp_limb_t cnd,
               mp_ptr zp,
               mp_srcptr xp,
               mp_srcptr yp,
               mp_size_t n) {
  mp_limb_t cond = (cnd != 0);
  mp_limb_t mask0 = cond - 1;
  mp_limb_t mask1 = ~mask0;
  mp_size_t i;

  for (i = 0; i < n; i++)
    zp[i] = (xp[i] & mask0) | (yp[i] & mask1);
}

void
mpn_cnd_swap(mp_limb_t cnd, mp_ptr ap, mp_ptr bp, mp_size_t n) {
  mp_limb_t mask = -(mp_limb_t)(cnd != 0);
  mp_size_t i;

  for (i = 0; i < n; i++) {
    mp_limb_t a = ap[i];
    mp_limb_t b = bp[i];
    mp_limb_t w = (a ^ b) & mask;

    ap[i] = a ^ w;
    bp[i] = b ^ w;
  }
}

void
mpn_cnd_zero(mp_limb_t cnd, mp_ptr rp, mp_srcptr ap, mp_size_t n) {
  mp_limb_t cond = (cnd != 0);
  mp_limb_t mask = cond - 1;
  mp_size_t i;

  for (i = 0; i < n; i++)
    rp[i] = ap[i] & mask;
}

int
mpn_sec_zero_p(mp_srcptr xp, mp_size_t xn) {
  /* Compute (x == 0) in constant time. */
  mp_limb_t w = 0;

  while (xn--)
    w |= xp[xn];

  w = (w >> 1) | (w & 1);

  return (w - 1) >> (MPI_LIMB_BITS - 1);
}

int
mpn_sec_eq(mp_srcptr xp, mp_srcptr yp, mp_size_t n) {
  /* Compute (x == y) in constant time. */
  mp_limb_t w = 0;

  while (n--)
    w |= xp[n] ^ yp[n];

  w = (w >> 1) | (w & 1);

  return (w - 1) >> (MPI_LIMB_BITS - 1);
}

int
mpn_sec_lt(mp_srcptr xp, mp_srcptr yp, mp_size_t n) {
  /* Compute (x < y) in constant time. */
  size_t shift = MPI_LIMB_BITS - 1;
  mp_size_t i = n * 2;
  mp_limb_t eq = 1;
  mp_limb_t lt = 0;
  mp_limb_t a, b;

  while (i--) {
    a = (xp[i / 2] >> ((i % 2) * (MPI_LIMB_BITS / 2))) & MPI_LLIMB_MASK;
    b = (yp[i / 2] >> ((i % 2) * (MPI_LIMB_BITS / 2))) & MPI_LLIMB_MASK;
    lt = ((eq ^ 1) & lt) | (eq & ((a - b) >> shift));
    eq &= ((a ^ b) - 1) >> shift;
  }

  return lt & (eq ^ 1);
}

int
mpn_sec_lte(mp_srcptr xp, mp_srcptr yp, mp_size_t n) {
  /* Compute (x <= y) in constant time. */
  size_t shift = MPI_LIMB_BITS - 1;
  mp_size_t i = n * 2;
  mp_limb_t eq = 1;
  mp_limb_t lt = 0;
  mp_limb_t a, b;

  while (i--) {
    a = (xp[i / 2] >> ((i % 2) * (MPI_LIMB_BITS / 2))) & MPI_LLIMB_MASK;
    b = (yp[i / 2] >> ((i % 2) * (MPI_LIMB_BITS / 2))) & MPI_LLIMB_MASK;
    lt = ((eq ^ 1) & lt) | (eq & ((a - b) >> shift));
    eq &= ((a ^ b) - 1) >> shift;
  }

  return lt | eq;
}

int
mpn_sec_gt(mp_srcptr xp, mp_srcptr yp, mp_size_t n) {
  /* Compute (x > y) in constant time. */
  return mpn_sec_lte(xp, yp, n) ^ 1;
}

int
mpn_sec_gte(mp_srcptr xp, mp_srcptr yp, mp_size_t n) {
  /* Compute (x >= y) in constant time. */
  return mpn_sec_lt(xp, yp, n) ^ 1;
}

/*
 * Import
 */

static void
mpn_import_be(mp_ptr rp, mp_size_t rn, const unsigned char *xp, size_t xn) {
  unsigned int bits = 0;
  mp_limb_t out = 0;
  size_t xi = xn;

  while (xi > 0 && rn > 0) {
    mp_limb_t in = xp[--xi];

    out |= (in << bits) & MPI_LIMB_MAX;
    bits += 8;

    if (bits >= MPI_LIMB_BITS) {
      *rp++ = out;
      rn--;

      bits -= MPI_LIMB_BITS;
      out = in >> (8 - bits);
    }
  }

  if (rn > 0) {
    *rp++ = out;
    if (--rn > 0)
      mpn_zero(rp, rn);
  }
}

static void
mpn_import_le(mp_ptr rp, mp_size_t rn, const unsigned char *xp, size_t xn) {
  unsigned int bits = 0;
  mp_limb_t out = 0;
  size_t xi = 0;

  while (xi < xn && rn > 0) {
    mp_limb_t in = xp[xi++];

    out |= (in << bits) & MPI_LIMB_MAX;
    bits += 8;

    if (bits >= MPI_LIMB_BITS) {
      *rp++ = out;
      rn--;

      bits -= MPI_LIMB_BITS;
      out = in >> (8 - bits);
    }
  }

  if (rn > 0) {
    *rp++ = out;
    if (--rn > 0)
      mpn_zero(rp, rn);
  }
}

void
mpn_import(mp_ptr rp, mp_size_t rn,
           const unsigned char *xp, size_t xn, int endian) {
  if (endian == 1)
    mpn_import_be(rp, rn, xp, xn);
  else if (endian == -1)
    mpn_import_le(rp, rn, xp, xn);
  else
    mpi_die("mpn_import: invalid endianness.");
}

/*
 * Export
 */

static void
mpn_export_be(unsigned char *rp, size_t rn, mp_srcptr xp, mp_size_t xn) {
  unsigned int bits = 0;
  mp_limb_t in = 0;
  unsigned char old;

  while (xn > 0 && rn > 0) {
    if (bits >= 8) {
      rp[--rn] = in;
      in >>= 8;
      bits -= 8;
    } else {
      old = in;
      in = *xp++;
      xn--;
      rp[--rn] = old | (in << bits);
      in >>= (8 - bits);
      bits += MPI_LIMB_BITS - 8;
    }
  }

  while (rn > 0) {
    rp[--rn] = in;
    in >>= 8;
  }
}

static void
mpn_export_le(unsigned char *rp, size_t rn, mp_srcptr xp, mp_size_t xn) {
  unsigned int bits = 0;
  mp_limb_t in = 0;
  unsigned char old;

  while (xn > 0 && rn > 0) {
    if (bits >= 8) {
      *rp++ = in;
      rn--;
      in >>= 8;
      bits -= 8;
    } else {
      old = in;
      in = *xp++;
      xn--;
      *rp++ = old | (in << bits);
      rn--;
      in >>= (8 - bits);
      bits += MPI_LIMB_BITS - 8;
    }
  }

  while (rn > 0) {
    *rp++ = in;
    rn--;
    in >>= 8;
  }
}

void
mpn_export(unsigned char *rp, size_t rn,
           mp_srcptr xp, mp_size_t xn, int endian) {
  if (endian == 1)
    mpn_export_be(rp, rn, xp, xn);
  else if (endian == -1)
    mpn_export_le(rp, rn, xp, xn);
  else
    mpi_die("mpn_export: invalid endianness.");
}

/*
 * Stringification
 */

size_t
mpn_out_str(FILE *stream, int base, mp_srcptr xp, mp_size_t xn) {
  mp_size_t bytes = 0;
  mp_limb_t ch, hi, lo;
  mp_size_t i;

  assert(base == 16);

  if (xn < 0) {
    fputc('-', stream);
    xn = -xn;
  }

  xn = mpn_normalized_size(xp, xn);

  while (xn--) {
    i = MPI_LIMB_BITS / 8;

    while (i--) {
      ch = (xp[xn] >> (i * 8)) & 0xff;

      if (bytes == 0 && ch == 0)
        continue;

      hi = ch >> 4;
      lo = ch & 0x0f;

      if (hi < 0x0a)
        hi += '0';
      else
        hi += 'a' - 0x0a;

      if (lo < 0x0a)
        lo += '0';
      else
        lo += 'a' - 0x0a;

      fputc(hi, stream);
      fputc(lo, stream);

      bytes += 2;
    }
  }

  return bytes;
}

/*
 * MPZ Interface
 */

/*
 * Initialization
 */

void
mpz_init(mpz_t r) {
  static const mp_limb_t dummy_limb = 0xc1a0;

  r->_mp_alloc = 0;
  r->_mp_size = 0;
  r->_mp_d = (mp_ptr)&dummy_limb;
}

void
mpz_init2(mpz_t r, mp_bitcnt_t bits) {
  /* The utility of this function is a bit limited,
     since many functions assign the result variable
     using mpz_swap. */
  mp_size_t rn;

  bits -= (bits != 0); /* Round down, except if 0. */
  rn = 1 + bits / MPI_LIMB_BITS;

  r->_mp_alloc = rn;
  r->_mp_size = 0;
  r->_mp_d = mpi_xalloc_limbs(rn);
}

void
mpz_init_set(mpz_t r, const mpz_t x) {
  mpz_init(r);
  mpz_set(r, x);
}

void
mpz_init_set_ui(mpz_t r, mp_limb_t x) {
  mpz_init(r);
  mpz_set_ui(r, x);
}

void
mpz_init_set_si(mpz_t r, mp_long_t x) {
  mpz_init(r);
  mpz_set_si(r, x);
}

void
mpz_init_set_u64(mpz_t r, uint64_t x) {
  mpz_init(r);
  mpz_set_u64(r, x);
}

/*
 * Uninitialization
 */

void
mpz_clear(mpz_t r) {
  if (r->_mp_alloc)
    mpi_free(r->_mp_d);
}

void
mpz_cleanse(mpz_t r) {
  if (r->_mp_alloc) {
    mpn_cleanse(r->_mp_d, r->_mp_alloc);
    mpi_free(r->_mp_d);
  }
}

/*
 * Internal
 */

static mp_ptr
mpz_realloc(mpz_t r, mp_size_t size) {
  size = MPI_MAX(size, 1);

  if (r->_mp_alloc)
    r->_mp_d = mpi_xrealloc_limbs(r->_mp_d, size);
  else
    r->_mp_d = mpi_xalloc_limbs(size);

  r->_mp_alloc = size;

  if (MPI_ABS(r->_mp_size) > size)
    r->_mp_size = 0;

  return r->_mp_d;
}

static mpz_srcptr
mpz_roinit_normal_n(mpz_t x, mp_srcptr xp, mp_size_t xs) {
  x->_mp_alloc = 0;
  x->_mp_d = (mp_ptr)xp;
  x->_mp_size = xs;
  return x;
}

/*
 * Assignment
 */

void
mpz_set(mpz_t r, const mpz_t x) {
  /* Allow the NOP r == x. */
  if (r != x) {
    mp_size_t n;
    mp_ptr rp;

    n = MPI_ABS(x->_mp_size);
    rp = MPZ_REALLOC(r, n);

    mpn_copyi(rp, x->_mp_d, n);
    r->_mp_size = x->_mp_size;
  }
}

void
mpz_roset(mpz_t r, const mpz_t x) {
  r->_mp_alloc = 0;
  r->_mp_size = x->_mp_size;
  r->_mp_d = (mp_ptr)x->_mp_d;
}

void
mpz_set_ui(mpz_t r, mp_limb_t x) {
  if (x > 0) {
    r->_mp_size = 1;
    MPZ_REALLOC(r, 1)[0] = x;
  } else {
    r->_mp_size = 0;
  }
}

void
mpz_set_si(mpz_t r, mp_long_t x) {
  if (x >= 0) {
    mpz_set_ui(r, x);
  } else { /* (x < 0) */
    r->_mp_size = -1;
    MPZ_REALLOC(r, 1)[0] = MPI_NEG_CAST(mp_limb_t, x);
  }
}

void
mpz_set_u64(mpz_t r, uint64_t x) {
#if MPI_LIMB_BITS == 32
  if (x > 0) {
    MPZ_REALLOC(r, 2);
    r->_mp_d[0] = x;
    r->_mp_d[1] = x >> 32;
    r->_mp_size = mpn_normalized_size(r->_mp_d, 2);
  } else {
    r->_mp_size = 0;
  }
#else
  mpz_set_ui(r, x);
#endif
}

/*
 * Conversion
 */

mp_limb_t
mpz_get_ui(const mpz_t u) {
  return u->_mp_size == 0 ? 0 : u->_mp_d[0];
}

mp_long_t
mpz_get_si(const mpz_t u) {
  mp_long_t r = mpz_get_ui(u) & (MPI_LIMB_HIGHBIT - 1);
  return u->_mp_size < 0 ? -r : r;
}

uint64_t
mpz_get_u64(const mpz_t u) {
#if MPI_LIMB_BITS == 32
  mp_size_t un = MPI_ABS(u->_mp_size);

  if (un == 0)
    return 0;

  if (un == 1)
    return u->_mp_d[0];

  return ((uint64_t)u->_mp_d[1] << 32) | u->_mp_d[0];
#else
  return mpz_get_ui(u);
#endif
}

/*
 * Conversion Testing
 */

int
mpz_fits_ulong_p(const mpz_t u) {
  return u->_mp_size == 0 || u->_mp_size == 1;
}

int
mpz_fits_slong_p(const mpz_t u) {
  return mpz_cmpabs_ui(u, MPI_LIMB_HIGHBIT - 1) <= 0;
}

int
mpz_fits_u64_p(const mpz_t u) {
#if MPI_LIMB_BITS == 32
  return MPI_ABS(u->_mp_size) <= 2;
#else
  return MPI_ABS(u->_mp_size) <= 1;
#endif
}

/*
 * Comparison
 */

int
mpz_sgn(const mpz_t u) {
  return MPI_CMP(u->_mp_size, 0);
}

int
mpz_cmp(const mpz_t a, const mpz_t b) {
  mp_size_t asize = a->_mp_size;
  mp_size_t bsize = b->_mp_size;

  if (asize != bsize)
    return (asize < bsize) ? -1 : 1;
  else if (asize >= 0)
    return mpn_cmp(a->_mp_d, b->_mp_d, asize);
  else
    return mpn_cmp(b->_mp_d, a->_mp_d, -asize);
}

int
mpz_cmp_ui(const mpz_t u, mp_limb_t v) {
  mp_size_t usize = u->_mp_size;

  if (usize < 0)
    return -1;
  else
    return mpz_cmpabs_ui(u, v);
}

int
mpz_cmp_si(const mpz_t u, mp_long_t v) {
  mp_size_t usize = u->_mp_size;

  if (v >= 0)
    return mpz_cmp_ui(u, v);
  else if (usize >= 0)
    return 1;
  else
    return -mpz_cmpabs_ui(u, MPI_NEG_CAST(mp_limb_t, v));
}

/*
 * Unsigned Comparison
 */

int
mpz_cmpabs(const mpz_t u, const mpz_t v) {
  return mpn_cmp4(u->_mp_d, MPI_ABS(u->_mp_size),
                  v->_mp_d, MPI_ABS(v->_mp_size));
}

int
mpz_cmpabs_ui(const mpz_t u, mp_limb_t v) {
  mp_size_t un = MPI_ABS(u->_mp_size);

  if (un > 1) {
    return 1;
  } else {
    mp_limb_t uu = mpz_get_ui(u);
    return MPI_CMP(uu, v);
  }
}

/*
 * Addition/Subtraction Engine
 */

static mp_size_t
mpz_abs_add(mpz_t r, const mpz_t a, const mpz_t b) {
  mp_size_t an = MPI_ABS(a->_mp_size);
  mp_size_t bn = MPI_ABS(b->_mp_size);
  mp_ptr rp;
  mp_limb_t cy;

  if (an < bn) {
    MPZ_SRCPTR_SWAP(a, b);
    MP_SIZE_T_SWAP(an, bn);
  }

  rp = MPZ_REALLOC(r, an + 1);
  cy = mpn_add(rp, a->_mp_d, an, b->_mp_d, bn);

  rp[an] = cy;

  return an + cy;
}

static mp_size_t
mpz_abs_sub(mpz_t r, const mpz_t a, const mpz_t b) {
  mp_size_t an = MPI_ABS(a->_mp_size);
  mp_size_t bn = MPI_ABS(b->_mp_size);
  int cmp;
  mp_ptr rp;

  cmp = mpn_cmp4(a->_mp_d, an, b->_mp_d, bn);

  if (cmp > 0) {
    rp = MPZ_REALLOC(r, an);
    mpi_assert_nocarry(mpn_sub(rp, a->_mp_d, an, b->_mp_d, bn));
    return mpn_normalized_size(rp, an);
  } else if (cmp < 0) {
    rp = MPZ_REALLOC(r, bn);
    mpi_assert_nocarry(mpn_sub(rp, b->_mp_d, bn, a->_mp_d, an));
    return -mpn_normalized_size(rp, bn);
  } else {
    return 0;
  }
}

/*
 * Addition
 */

void
mpz_add(mpz_t r, const mpz_t a, const mpz_t b) {
  mp_size_t rn;

  if ((a->_mp_size ^ b->_mp_size) >= 0)
    rn = mpz_abs_add(r, a, b);
  else
    rn = mpz_abs_sub(r, a, b);

  r->_mp_size = a->_mp_size >= 0 ? rn : -rn;
}

void
mpz_add_ui(mpz_ptr w, mpz_srcptr u, mp_limb_t vval) {
  mp_srcptr up;
  mp_ptr wp;
  mp_size_t usize, wsize;
  mp_size_t abs_usize;

  usize = u->_mp_size;

  if (usize == 0) {
    MPZ_REALLOC(w, 1)[0] = vval;
    w->_mp_size = (vval != 0);
    return;
  }

  abs_usize = MPI_ABS(usize);

  /* If not space for W (and possible carry), increase space. */
  wp = MPZ_REALLOC(w, abs_usize + 1);

  /* These must be after realloc (U may be the same as W). */
  up = u->_mp_d;

  if (usize >= 0) {
    mp_limb_t cy;
    cy = mpn_add_1(wp, up, abs_usize, vval);
    wp[abs_usize] = cy;
    wsize = (abs_usize + cy);
  } else {
    /* The signs are different. Need exact comparison to
       determine which operand to subtract from which. */
    if (abs_usize == 1 && up[0] < vval) {
      wp[0] = vval - up[0];
      wsize = 1;
    } else {
      mpn_sub_1(wp, up, abs_usize, vval);
      /* Size can decrease with at most one limb. */
      wsize = -(abs_usize - (wp[abs_usize - 1] == 0));
    }
  }

  w->_mp_size = wsize;
}

/*
 * Subtraction
 */

void
mpz_sub(mpz_t r, const mpz_t a, const mpz_t b) {
  mp_size_t rn;

  if ((a->_mp_size ^ b->_mp_size) >= 0)
    rn = mpz_abs_sub(r, a, b);
  else
    rn = mpz_abs_add(r, a, b);

  r->_mp_size = a->_mp_size >= 0 ? rn : -rn;
}

void
mpz_sub_ui(mpz_ptr w, mpz_srcptr u, mp_limb_t vval) {
  mp_srcptr up;
  mp_ptr wp;
  mp_size_t usize, wsize;
  mp_size_t abs_usize;

  usize = u->_mp_size;

  if (usize == 0) {
    MPZ_REALLOC(w, 1)[0] = vval;
    w->_mp_size = -(vval != 0);
    return;
  }

  abs_usize = MPI_ABS(usize);

  /* If not space for W (and possible carry), increase space. */
  wp = MPZ_REALLOC(w, abs_usize + 1);

  /* These must be after realloc (U may be the same as W). */
  up = u->_mp_d;

  if (usize < 0) {
    mp_limb_t cy;
    cy = mpn_add_1(wp, up, abs_usize, vval);
    wp[abs_usize] = cy;
    wsize = -(abs_usize + cy);
  } else {
    /* The signs are different. Need exact comparison to
       determine which operand to subtract from which. */
    if (abs_usize == 1 && up[0] < vval) {
      wp[0] = vval - up[0];
      wsize = -1;
    } else {
      mpn_sub_1(wp, up, abs_usize, vval);
      /* Size can decrease with at most one limb. */
      wsize = (abs_usize - (wp[abs_usize - 1] == 0));
    }
  }

  w->_mp_size = wsize;
}

/*
 * Multiplication
 */

void
mpz_mul(mpz_t r, const mpz_t u, const mpz_t v) {
  int sign;
  mp_size_t un, vn, rn;
  mpz_t t;
  mp_ptr tp;

  un = u->_mp_size;
  vn = v->_mp_size;

  if (un == 0 || vn == 0) {
    r->_mp_size = 0;
    return;
  }

  sign = (un ^ vn) < 0;

  un = MPI_ABS(un);
  vn = MPI_ABS(vn);

  mpz_init2(t, (un + vn) * MPI_LIMB_BITS);

  tp = t->_mp_d;

  if (u == v)
    mpn_sqr(tp, u->_mp_d, un);
  else if (un >= vn)
    mpn_mul(tp, u->_mp_d, un, v->_mp_d, vn);
  else
    mpn_mul(tp, v->_mp_d, vn, u->_mp_d, un);

  rn = un + vn;
  rn -= tp[rn - 1] == 0;

  t->_mp_size = sign ? -rn : rn;

  mpz_swap(r, t);
  mpz_clear(t);
}

void
mpz_mul_ui(mpz_ptr prod, mpz_srcptr mult, mp_limb_t small_mult) {
  mp_size_t size;
  mp_size_t sign_product;
  mp_limb_t cy;
  mp_ptr pp;

  sign_product = mult->_mp_size;

  if (sign_product == 0 || small_mult == 0) {
    prod->_mp_size = 0;
    return;
  }

  size = MPI_ABS(sign_product);

  pp = MPZ_REALLOC(prod, size + 1);
  cy = mpn_mul_1(pp, mult->_mp_d, size, small_mult);
  pp[size] = cy;
  size += cy != 0;

  prod->_mp_size = (sign_product < 0) ? -size : size;
}

void
mpz_mul_si(mpz_t r, const mpz_t u, mp_long_t v) {
  if (v < 0) {
    mpz_mul_ui(r, u, MPI_NEG_CAST(mp_limb_t, v));
    mpz_neg(r, r);
  } else {
    mpz_mul_ui(r, u, v);
  }
}

/*
 * Division Engine
 */

static int
mpz_div_qr(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d,
           enum mpz_div_round_mode mode) {
  /* Allows q or r to be zero. */
  /* Returns 1 iff remainder is non-zero. */
  mp_size_t ns, ds, nn, dn, qs;

  ns = n->_mp_size;
  ds = d->_mp_size;

  if (ds == 0)
    mpi_die("mpz_div_qr: Divide by zero.");

  if (ns == 0) {
    if (q)
      q->_mp_size = 0;

    if (r)
      r->_mp_size = 0;

    return 0;
  }

  nn = MPI_ABS(ns);
  dn = MPI_ABS(ds);

  qs = ds ^ ns;

  if (nn < dn) {
    if (mode == MPI_DIV_CEIL && qs >= 0) {
      /* q = 1, r = n - d */
      if (r)
        mpz_sub(r, n, d);

      if (q)
        mpz_set_ui(q, 1);
    } else if (mode == MPI_DIV_FLOOR && qs < 0) {
      /* q = -1, r = n + d */
      if (r)
        mpz_add(r, n, d);

      if (q)
        mpz_set_si(q, -1);
    } else {
      /* q = 0, r = d */
      if (r)
        mpz_set(r, n);

      if (q)
        q->_mp_size = 0;
    }

    return 1;
  } else {
    mp_ptr np, qp;
    mp_size_t qn, rn;
    mpz_t tq, tr;

    mpz_init_set(tr, n);
    np = tr->_mp_d;

    qn = nn - dn + 1;

    if (q) {
      mpz_init2(tq, qn * MPI_LIMB_BITS);
      qp = tq->_mp_d;
    } else {
      qp = NULL;
    }

    mpn_div_qr(qp, np, nn, d->_mp_d, dn);

    if (qp) {
      qn -= (qp[qn - 1] == 0);

      tq->_mp_size = qs < 0 ? -qn : qn;
    }

    rn = mpn_normalized_size(np, dn);
    tr->_mp_size = ns < 0 ? -rn : rn;

    if (mode == MPI_DIV_FLOOR && qs < 0 && rn != 0) {
      if (q)
        mpz_sub_ui(tq, tq, 1);

      if (r)
        mpz_add(tr, tr, d);
    } else if (mode == MPI_DIV_CEIL && qs >= 0 && rn != 0) {
      if (q)
        mpz_add_ui(tq, tq, 1);

      if (r)
        mpz_sub(tr, tr, d);
    }

    if (q) {
      mpz_swap(tq, q);
      mpz_clear(tq);
    }

    if (r)
      mpz_swap(tr, r);

    mpz_clear(tr);

    return rn != 0;
  }
}

static mp_limb_t
mpz_div_qr_ui(mpz_t q, mpz_t r, const mpz_t n, mp_limb_t d,
              enum mpz_div_round_mode mode) {
  mp_limb_t ret;
  mpz_t rr, dd;

  mpz_init(rr);
  mpz_init_set_ui(dd, d);
  mpz_div_qr(q, rr, n, dd, mode);
  mpz_clear(dd);

  ret = mpz_get_ui(rr);

  if (r)
    mpz_swap(r, rr);

  mpz_clear(rr);

  return ret;
}

/*
 * Truncation Division
 */

void
mpz_quorem(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d) {
  mpz_div_qr(q, r, n, d, MPI_DIV_TRUNC);
}

void
mpz_quo(mpz_t q, const mpz_t n, const mpz_t d) {
  mpz_div_qr(q, NULL, n, d, MPI_DIV_TRUNC);
}

void
mpz_rem(mpz_t r, const mpz_t n, const mpz_t d) {
  mpz_div_qr(NULL, r, n, d, MPI_DIV_TRUNC);
}

mp_limb_t
mpz_quo_ui(mpz_t q, const mpz_t n, mp_limb_t d) {
  return mpz_div_qr_ui(q, NULL, n, d, MPI_DIV_TRUNC);
}

mp_limb_t
mpz_rem_ui(const mpz_t n, mp_limb_t d) {
  return mpz_div_qr_ui(NULL, NULL, n, d, MPI_DIV_TRUNC);
}

/*
 * Euclidean Division
 */

void
mpz_divmod(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d) {
  mpz_div_qr(q, r, n, d, d->_mp_size >= 0 ? MPI_DIV_FLOOR : MPI_DIV_CEIL);
}

void
mpz_div(mpz_t q, const mpz_t n, const mpz_t d) {
  mpz_div_qr(q, NULL, n, d, d->_mp_size >= 0 ? MPI_DIV_FLOOR : MPI_DIV_CEIL);
}

void
mpz_mod(mpz_t r, const mpz_t n, const mpz_t d) {
  if (mpz_cmpabs(n, d) < 0) {
    if (n->_mp_size < 0) {
      mp_size_t nn = MPI_ABS(n->_mp_size);
      mp_size_t dn = MPI_ABS(d->_mp_size);
      mp_ptr rp = MPZ_REALLOC(r, dn);

      mpi_assert_nocarry(mpn_sub(rp, d->_mp_d, dn, n->_mp_d, nn));

      r->_mp_size = mpn_normalized_size(rp, dn);
    } else {
      mpz_set(r, n);
    }
  } else {
    mpz_div_qr(NULL, r, n, d, d->_mp_size >= 0 ? MPI_DIV_FLOOR : MPI_DIV_CEIL);
  }
}

mp_limb_t
mpz_div_ui(mpz_t q, const mpz_t n, mp_limb_t d) {
  return mpz_div_qr_ui(q, NULL, n, d, MPI_DIV_FLOOR);
}

mp_limb_t
mpz_mod_ui(const mpz_t n, mp_limb_t d) {
  return mpz_div_qr_ui(NULL, NULL, n, d, MPI_DIV_FLOOR);
}

/*
 * Exact Division
 */

void
mpz_divexact(mpz_t q, const mpz_t n, const mpz_t d) {
  mpi_assert_nocarry(mpz_div_qr(q, NULL, n, d, MPI_DIV_TRUNC));
}

void
mpz_divexact_ui(mpz_t q, const mpz_t n, mp_limb_t d) {
  mpi_assert_nocarry(mpz_div_qr_ui(q, NULL, n, d, MPI_DIV_TRUNC));
}

/*
 * Left Shift
 */

void
mpz_lshift(mpz_t r, const mpz_t u, mp_bitcnt_t bits) {
  mp_size_t un, rn;
  mp_size_t limbs;
  unsigned shift;
  mp_ptr rp;

  un = MPI_ABS(u->_mp_size);

  if (un == 0) {
    r->_mp_size = 0;
    return;
  }

  limbs = bits / MPI_LIMB_BITS;
  shift = bits % MPI_LIMB_BITS;

  rn = un + limbs + (shift > 0);
  rp = MPZ_REALLOC(r, rn);

  if (shift > 0) {
    mp_limb_t cy = mpn_lshift(rp + limbs, u->_mp_d, un, shift);
    rp[rn - 1] = cy;
    rn -= (cy == 0);
  } else {
    mpn_copyd(rp + limbs, u->_mp_d, un);
  }

  mpn_zero(rp, limbs);

  r->_mp_size = (u->_mp_size < 0) ? -rn : rn;
}

/*
 * Right Shift
 */

void
mpz_rshift(mpz_t r, const mpz_t u, mp_bitcnt_t bits) {
  mp_size_t un, rn, limbs, shift;
  mp_ptr rp;

  un = MPI_ABS(u->_mp_size);

  if (un == 0) {
    r->_mp_size = 0;
    return;
  }

  limbs = bits / MPI_LIMB_BITS;
  shift = bits % MPI_LIMB_BITS;
  rn = un - limbs;

  if (rn <= 0) {
    rn = 0;
  } else {
    rp = MPZ_REALLOC(r, rn);

    if (shift != 0) {
      mpn_rshift(rp, u->_mp_d + limbs, rn, shift);
      rn -= rp[rn - 1] == 0;
    } else {
      mpn_copyi(rp, u->_mp_d + limbs, rn);
    }
  }

  r->_mp_size = (u->_mp_size < 0) ? -rn : rn;
}

/*
 * Bit Manipulation
 */

mp_limb_t
mpz_get_bit(const mpz_t d, mp_bitcnt_t pos) {
  return mpn_get_bit(d->_mp_d, MPI_ABS(d->_mp_size), pos);
}

mp_limb_t
mpz_get_bits(const mpz_t d, mp_bitcnt_t pos, mp_bitcnt_t width) {
  return mpn_get_bits(d->_mp_d, MPI_ABS(d->_mp_size), pos, width);
}

void
mpz_set_bit(mpz_t d, mp_bitcnt_t pos) {
  mp_size_t index = pos / MPI_LIMB_BITS;
  mp_size_t dn = MPI_ABS(d->_mp_size);
  mp_ptr dp = d->_mp_d;

  if (index >= dn) {
    dp = MPZ_REALLOC(d, index + 1);

    while (dn < index + 1)
      dp[dn++] = 0;

    d->_mp_size = (d->_mp_size < 0) ? -dn : dn;
  }

  dp[index] |= (mp_limb_t)1 << (pos % MPI_LIMB_BITS);
}

void
mpz_clr_bit(mpz_t d, mp_bitcnt_t pos) {
  mp_size_t index = pos / MPI_LIMB_BITS;

  if (index < MPI_ABS(d->_mp_size))
    d->_mp_d[index] &= ~((mp_limb_t)1 << (pos % MPI_LIMB_BITS));
}

/*
 * Negation
 */

void
mpz_abs(mpz_t r, const mpz_t u) {
  mpz_set(r, u);
  r->_mp_size = MPI_ABS(r->_mp_size);
}

void
mpz_neg(mpz_t r, const mpz_t u) {
  mpz_set(r, u);
  r->_mp_size = -r->_mp_size;
}

/*
 * Number Theoretic Functions
 */

static mp_bitcnt_t
mpz_make_odd(mpz_t r) {
  mp_bitcnt_t shift;

  assert(r->_mp_size > 0);

  /* Count trailing zeros, equivalent to mpn_scan1,
     because we know that there is a 1. */
  shift = mpn_ctz(r->_mp_d, r->_mp_size);

  mpz_rshift(r, r, shift);

  return shift;
}

void
mpz_gcd(mpz_t g, const mpz_t u, const mpz_t v) {
  mpz_t tu, tv;
  mp_bitcnt_t uz, vz, gz;

  if (u->_mp_size == 0) {
    mpz_abs(g, v);
    return;
  }

  if (v->_mp_size == 0) {
    mpz_abs(g, u);
    return;
  }

  mpz_init(tu);
  mpz_init(tv);

  mpz_abs(tu, u);
  uz = mpz_make_odd(tu);

  mpz_abs(tv, v);
  vz = mpz_make_odd(tv);
  gz = MPI_MIN(uz, vz);

  if (tu->_mp_size < tv->_mp_size)
    mpz_swap(tu, tv);

  mpz_rem(tu, tu, tv);

  if (tu->_mp_size == 0) {
    mpz_swap(g, tv);
  } else {
    for (;;) {
      int c;

      mpz_make_odd(tu);
      c = mpz_cmp(tu, tv);

      if (c == 0) {
        mpz_swap(g, tu);
        break;
      }

      if (c < 0)
        mpz_swap(tu, tv);

      if (tv->_mp_size == 1) {
        mp_limb_t vl = tv->_mp_d[0];
        mp_limb_t ul = mpz_rem_ui(tu, vl);
        mpz_set_ui(g, mpn_gcd_11(ul, vl));
        break;
      }

      mpz_sub(tu, tu, tv);
    }
  }

  mpz_clear(tu);
  mpz_clear(tv);
  mpz_lshift(g, g, gz);
}

void
mpz_lcm(mpz_t r, const mpz_t u, const mpz_t v) {
  mpz_t g;

  if (u->_mp_size == 0 || v->_mp_size == 0) {
    r->_mp_size = 0;
    return;
  }

  mpz_init(g);

  mpz_gcd(g, u, v);
  mpz_divexact(g, u, g);
  mpz_mul(r, g, v);

  mpz_clear(g);
  mpz_abs(r, r);
}

void
mpz_gcdext(mpz_t g, mpz_t s, mpz_t t, const mpz_t u, const mpz_t v) {
  mpz_t tu, tv, s0, s1, t0, t1;
  mp_bitcnt_t uz, vz, gz;
  mp_bitcnt_t power;

  if (u->_mp_size == 0) {
    /* g = 0 u + sgn(v) v */
    mp_long_t sign = mpz_sgn(v);

    mpz_abs(g, v);

    if (s)
      s->_mp_size = 0;

    if (t)
      mpz_set_si(t, sign);

    return;
  }

  if (v->_mp_size == 0) {
    /* g = sgn(u) u + 0 v */
    mp_long_t sign = mpz_sgn(u);

    mpz_abs(g, u);

    if (s)
      mpz_set_si(s, sign);

    if (t)
      t->_mp_size = 0;

    return;
  }

  mpz_init(tu);
  mpz_init(tv);
  mpz_init(s0);
  mpz_init(s1);
  mpz_init(t0);
  mpz_init(t1);

  mpz_abs(tu, u);
  uz = mpz_make_odd(tu);

  mpz_abs(tv, v);
  vz = mpz_make_odd(tv);

  gz = MPI_MIN(uz, vz);

  uz -= gz;
  vz -= gz;

  /* Cofactors corresponding to odd gcd. gz handled later. */
  if (tu->_mp_size < tv->_mp_size) {
    mpz_swap(tu, tv);
    MPZ_SRCPTR_SWAP(u, v);
    MPZ_PTR_SWAP(s, t);
    MP_BITCNT_T_SWAP(uz, vz);
  }

  /* Maintain
   *
   *   u = t0 tu + t1 tv
   *   v = s0 tu + s1 tv
   *
   * where u and v denote the inputs with common factors of
   * two eliminated, and det (s0, t0; s1, t1) = 2^p. Then
   *
   *   2^p tu =  s1 u - t1 v
   *   2^p tv = -s0 u + t0 v
   */

  /* After initial division, tu = q tv + tu', we have
   *
   *   u = 2^uz (tu' + q tv)
   *   v = 2^vz tv
   *
   * or
   *
   *   t0 = 2^uz, t1 = 2^uz q
   *   s0 = 0,    s1 = 2^vz
   */

  mpz_set_bit(t0, uz);
  mpz_quorem(t1, tu, tu, tv);
  mpz_lshift(t1, t1, uz);

  mpz_set_bit(s1, vz);
  power = uz + vz;

  if (tu->_mp_size > 0) {
    mp_bitcnt_t shift;

    shift = mpz_make_odd(tu);

    mpz_lshift(t0, t0, shift);
    mpz_lshift(s0, s0, shift);

    power += shift;

    for (;;) {
      int c;

      c = mpz_cmp(tu, tv);

      if (c == 0)
        break;

      if (c < 0) {
        /* tv = tv' + tu
         *
         * u = t0 tu + t1 (tv' + tu) = (t0 + t1) tu + t1 tv'
         * v = s0 tu + s1 (tv' + tu) = (s0 + s1) tu + s1 tv'
         */

        mpz_sub(tv, tv, tu);
        mpz_add(t0, t0, t1);
        mpz_add(s0, s0, s1);

        shift = mpz_make_odd(tv);

        mpz_lshift(t1, t1, shift);
        mpz_lshift(s1, s1, shift);
      } else {
        mpz_sub(tu, tu, tv);
        mpz_add(t1, t0, t1);
        mpz_add(s1, s0, s1);

        shift = mpz_make_odd(tu);

        mpz_lshift(t0, t0, shift);
        mpz_lshift(s0, s0, shift);
      }

      power += shift;
    }
  }

  /* Now tv = odd part of gcd, and -s0 and t0
     are corresponding cofactors. */

  mpz_lshift(tv, tv, gz);
  mpz_neg(s0, s0);

  /* 2^p g = s0 u + t0 v. Eliminate one factor
     of two at a time. To adjust cofactors, we
     need u / g and v / g */

  mpz_divexact(s1, v, tv);
  mpz_abs(s1, s1);
  mpz_divexact(t1, u, tv);
  mpz_abs(t1, t1);

  while (power-- > 0) {
    /* s0 u + t0 v = (s0 - v/g) u - (t0 + u/g) v */
    if (mpz_odd_p(s0) || mpz_odd_p(t0)) {
      mpz_sub(s0, s0, s1);
      mpz_add(t0, t0, t1);
    }

    assert(mpz_even_p(t0) && mpz_even_p(s0));

    mpz_rshift(s0, s0, 1);
    mpz_rshift(t0, t0, 1);
  }

  /* Arrange so that |s| < |u| / 2g */
  mpz_add(s1, s0, s1);

  if (mpz_cmpabs(s0, s1) > 0) {
    mpz_swap(s0, s1);
    mpz_sub(t0, t0, t1);
  }

  if (u->_mp_size < 0)
    mpz_neg(s0, s0);

  if (v->_mp_size < 0)
    mpz_neg(t0, t0);

  mpz_swap(g, tv);

  if (s)
    mpz_swap(s, s0);

  if (t)
    mpz_swap(t, t0);

  mpz_clear(tu);
  mpz_clear(tv);
  mpz_clear(s0);
  mpz_clear(s1);
  mpz_clear(t0);
  mpz_clear(t1);
}

int
mpz_invert(mpz_t r, const mpz_t u, const mpz_t m) {
  mpz_t g, tr;
  int invertible;

  if (u->_mp_size == 0 || mpz_cmpabs_ui(m, 1) <= 0)
    return 0;

  mpz_init(g);
  mpz_init(tr);

  mpz_gcdext(g, tr, NULL, u, m);

  invertible = (mpz_cmp_ui(g, 1) == 0);

  if (invertible) {
    if (tr->_mp_size < 0) {
      if (m->_mp_size >= 0)
        mpz_add(tr, tr, m);
      else
        mpz_sub(tr, tr, m);
    }
    mpz_swap(r, tr);
  }

  mpz_clear(g);
  mpz_clear(tr);

  return invertible;
}

int
mpz_jacobi(const mpz_t x, const mpz_t y) {
  mp_bitcnt_t bits;
  mp_limb_t bmod8;
  mpz_t a, b;
  int j = 1;

  if (!mpz_odd_p(y))
    mpi_die("mpz_jacobi: Even argument.");

  mpz_init_set(a, x);
  mpz_init_set(b, y);

  if (b->_mp_size < 0) {
    if (a->_mp_size < 0)
      j = -1;

    mpz_neg(b, b);
  }

  for (;;) {
    if (mpz_cmp_ui(b, 1) == 0)
      break;

    if (a->_mp_size == 0) {
      j = 0;
      break;
    }

    mpz_mod(a, a, b);

    if (a->_mp_size == 0) {
      j = 0;
      break;
    }

    bits = mpz_make_odd(a);

    if (bits & 1) {
      bmod8 = b->_mp_d[0] & 7;

      if (bmod8 == 3 || bmod8 == 5)
        j = -j;
    }

    if ((a->_mp_d[0] & 3) == 3 && (b->_mp_d[0] & 3) == 3)
      j = -j;

    mpz_swap(a, b);
  }

  mpz_clear(a);
  mpz_clear(b);

  return j;
}

void
mpz_powm(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t m) {
  mpz_t tr;
  mpz_t base;
  mp_size_t en, mn;
  mp_srcptr mp;
  struct mpi_div_inverse minv;
  unsigned shift;
  mp_ptr tp = NULL;

  en = MPI_ABS(e->_mp_size);
  mn = MPI_ABS(m->_mp_size);

  if (mn == 0)
    mpi_die("mpz_powm: Zero modulo.");

  if (en == 0) {
    mpz_set_ui(r, 1);
    return;
  }

  mp = m->_mp_d;
  mpn_div_qr_invert(&minv, mp, mn);

  shift = minv.shift;

  if (shift > 0) {
    /* To avoid shifts, we do all our reductions, except
       the final one, using a *normalized* m. */
    minv.shift = 0;

    tp = mpi_xalloc_limbs(mn);
    mpi_assert_nocarry(mpn_lshift(tp, mp, mn, shift));
    mp = tp;
  }

  mpz_init(base);

  if (e->_mp_size < 0) {
    if (!mpz_invert(base, b, m))
      mpi_die("mpz_powm: Negative exponent and non-invertible base.");
  } else {
    mp_size_t bn;

    mpz_abs(base, b);

    bn = base->_mp_size;

    if (bn >= mn) {
      mpn_div_qr_preinv(NULL, base->_mp_d, base->_mp_size, mp, mn, &minv);
      bn = mn;
    }

    /* We have reduced the absolute value. Now take
       care of the sign. Note that we get zero represented
       non-canonically as m. */
    if (b->_mp_size < 0) {
      mp_ptr bp = MPZ_REALLOC(base, mn);
      mpi_assert_nocarry(mpn_sub(bp, mp, mn, bp, bn));
      bn = mn;
    }

    base->_mp_size = mpn_normalized_size(base->_mp_d, bn);
  }

  mpz_init_set_ui(tr, 1);

  while (--en >= 0) {
    mp_limb_t w = e->_mp_d[en];
    mp_limb_t bit;

    bit = MPI_LIMB_HIGHBIT;

    do {
      mpz_mul(tr, tr, tr);

      if (w & bit)
        mpz_mul(tr, tr, base);

      if (tr->_mp_size > mn) {
        mpn_div_qr_preinv(NULL, tr->_mp_d, tr->_mp_size, mp, mn, &minv);
        tr->_mp_size = mpn_normalized_size(tr->_mp_d, mn);
      }

      bit >>= 1;
    } while (bit > 0);
  }

  /* Final reduction */
  if (tr->_mp_size >= mn) {
    minv.shift = shift;
    mpn_div_qr_preinv(NULL, tr->_mp_d, tr->_mp_size, mp, mn, &minv);
    tr->_mp_size = mpn_normalized_size(tr->_mp_d, mn);
  }

  if (tp)
    mpi_free(tp);

  mpz_swap(r, tr);
  mpz_clear(tr);
  mpz_clear(base);
}

void
mpz_powm_sec(mpz_ptr r, mpz_srcptr b, mpz_srcptr e, mpz_srcptr m) {
  mp_ptr rp, scratch;
  mp_size_t mn, itch;
  mpz_t t;

  if (e->_mp_size < 0)
    mpi_die("mpz_powm_sec: Negative exponent.");

  if (e->_mp_size == 0)
    mpi_die("mpz_powm_sec: Zero exponent.");

  if (m->_mp_size == 0)
    mpi_die("mpz_powm_sec: Zero modulus.");

  if ((m->_mp_d[0] & 1) == 0)
    mpi_die("mpz_powm_sec: Even modulus.");

  mn = MPI_ABS(m->_mp_size);

  mpz_init(t);

  if (mpz_sgn(b) < 0 || mpz_cmpabs(b, m) >= 0)
    mpz_mod(t, b, m);
  else
    mpz_set(t, b);

  itch = MPN_POWM_SEC_ITCH(mn);
  scratch = mpi_xalloc_limbs(itch);

  rp = MPZ_REALLOC(r, mn);

  mpn_powm_sec(rp, t->_mp_d, t->_mp_size,
                   e->_mp_d, e->_mp_size,
                   m->_mp_d, mn,
                   scratch);

  r->_mp_size = mpn_normalized_size(rp, mn);

  mpi_free(scratch);
  mpz_clear(t);
}

/*
 * Primality Testing
 */

int
mpz_is_prime_mr(const mpz_t n, unsigned long reps,
                int force2, mp_rng_t rng, void *arg) {
  mpz_t nm1, nm3, q, x, y;
  unsigned long k, i, j;
  int ret = 0;

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
  k = mpz_ctz(nm1);

  /* q = nm1 >> k */
  mpz_rshift(q, nm1, k);

  for (i = 0; i < reps; i++) {
    if (i == reps - 1 && force2) {
      /* x = 2 */
      mpz_set_ui(x, 2);
    } else {
      /* x = random integer in [2,n-1] */
      mpz_random_int(x, nm3, rng, arg);
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

  ret = 1;
fail:
  mpz_clear(nm1);
  mpz_clear(nm3);
  mpz_clear(q);
  mpz_clear(x);
  mpz_clear(y);
  return ret;
}

int
mpz_is_prime_lucas(const mpz_t n, unsigned long limit) {
  mpz_t d, s, nm2, vk, vk1, t1, t2, t3;
  unsigned long i, p, r, t;
  int ret = 0;
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
      mpz_set_bit(t2, mpz_bitlen(n) / 2 + 1);

      do {
        mpz_swap(t1, t2);
        mpz_quo(t2, n, t1);
        mpz_add(t2, t2, t1);
        mpz_rshift(t2, t2, 1);
      } while (mpz_cmpabs(t2, t1) < 0);

      mpz_mul(t2, t1, t1);

      if (mpz_cmp(t2, n) == 0)
        goto fail;
    }

    p += 1;
  }

  /* s = n + 1 */
  mpz_add_ui(s, n, 1);

  /* r = s factors of 2 */
  r = mpz_ctz(s);

  /* nm2 = n - 2 */
  mpz_sub_ui(nm2, n, 2);

  /* vk = 2 */
  mpz_set_ui(vk, 2);

  /* vk1 = p */
  mpz_set_ui(vk1, p);

  /* s >>= r */
  mpz_rshift(s, s, r);

  for (i = mpz_bitlen(s) + 1; i-- > 0;) {
    /* if floor(s / 2^i) mod 2 == 1 */
    if (mpz_get_bit(s, i)) {
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
    mpz_lshift(t2, vk1, 1);

    if (mpz_cmp(t1, t2) < 0)
      mpz_swap(t1, t2);

    mpz_sub(t1, t1, t2);
    mpz_mod(t3, t1, n);

    /* if t3 == 0 */
    if (mpz_sgn(t3) == 0)
      goto succeed;
  }

  for (t = 1; t < r; t++) {
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

int
mpz_is_prime(const mpz_t p, unsigned long rounds, mp_rng_t rng, void *arg) {
  static const mp_limb_t primes_a =
    3ul * 5ul * 7ul * 11ul * 13ul * 17ul * 19ul * 23ul * 37ul;
  static const mp_limb_t primes_b = 29ul * 31ul * 41ul * 43ul * 47ul * 53ul;
  mp_limb_t ra, rb;

  if (mpz_sgn(p) <= 0)
    return 0;

  if (mpz_cmp_ui(p, 64) < 0) {
    static const uint64_t prime_mask = 0ull
      | 1ull <<  2 | 1ull <<  3 | 1ull <<  5 | 1ull << 7
      | 1ull << 11 | 1ull << 13 | 1ull << 17 | 1ull << 19
      | 1ull << 23 | 1ull << 29 | 1ull << 31 | 1ull << 37
      | 1ull << 41 | 1ull << 43 | 1ull << 47 | 1ull << 53
      | 1ull << 59 | 1ull << 61;

    return (prime_mask >> mpz_get_ui(p)) & 1;
  }

  if (mpz_even_p(p))
    return 0;

  ra = mpz_rem_ui(p, primes_a);
  rb = mpz_rem_ui(p, primes_b);

  if (ra % 3 == 0
      || ra % 5 == 0
      || ra % 7 == 0
      || ra % 11 == 0
      || ra % 13 == 0
      || ra % 17 == 0
      || ra % 19 == 0
      || ra % 23 == 0
      || ra % 37 == 0
      || rb % 29 == 0
      || rb % 31 == 0
      || rb % 41 == 0
      || rb % 43 == 0
      || rb % 47 == 0
      || rb % 53 == 0) {
    return 0;
  }

  if (!mpz_is_prime_mr(p, rounds + 1, 1, rng, arg))
    return 0;

  if (!mpz_is_prime_lucas(p, 0))
    return 0;

  return 1;
}

void
mpz_random_prime(mpz_t ret, mp_bitcnt_t bits, mp_rng_t rng, void *arg) {
  static const uint64_t primes[15] =
    { 3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 };
  static const uint64_t product = 16294579238595022365ull;
  uint64_t mod, delta, m, p;
  mpz_t prod, tmp;
  size_t i;

  assert(bits > 1);

  mpz_init_set_u64(prod, product);
  mpz_init(tmp);

  for (;;) {
    mpz_random_bits(ret, bits, rng, arg);

    mpz_set_bit(ret, bits - 1);
    mpz_set_bit(ret, bits - 2);
    mpz_set_bit(ret, 0);

    mpz_mod(tmp, ret, prod);

    mod = mpz_get_u64(tmp);

    for (delta = 0; delta < (1ull << 20); delta += 2) {
      m = mod + delta;

      for (i = 0; i < sizeof(primes) / sizeof(primes[0]); i++) {
        p = primes[i];

        if ((m % p) == 0 && (bits > 6 || m != p))
          goto next;
      }

      mpz_add_ui(ret, ret, (mp_limb_t)delta);

      break;
next:
      ;
    }

    if (mpz_bitlen(ret) != bits)
      continue;

    if (!mpz_is_prime(ret, 20, rng, arg))
      continue;

    break;
  }

  mpz_cleanse(prod);
  mpz_cleanse(tmp);
}

/*
 * Helpers
 */

int
mpz_odd_p(const mpz_t z) {
  if (z->_mp_size == 0)
    return 0;

  return z->_mp_d[0] & 1;
}

int
mpz_even_p(const mpz_t z) {
  return !mpz_odd_p(z);
}

mp_bitcnt_t
mpz_bitlen(const mpz_t u) {
  return mpn_bitlen(u->_mp_d, MPI_ABS(u->_mp_size));
}

mp_bitcnt_t
mpz_ctz(const mpz_t u) {
  return mpn_ctz(u->_mp_d, MPI_ABS(u->_mp_size));
}

size_t
mpz_bytelen(const mpz_t u) {
  return (mpz_bitlen(u) + 7) / 8;
}

void
mpz_swap(mpz_t u, mpz_t v) {
  MP_SIZE_T_SWAP(u->_mp_size, v->_mp_size);
  MP_SIZE_T_SWAP(u->_mp_alloc, v->_mp_alloc);
  MP_PTR_SWAP(u->_mp_d, v->_mp_d);
}

/*
 * Limb Manipulation
 */

mp_size_t
mpz_size(const mpz_t u) {
  return MPI_ABS(u->_mp_size);
}

mp_limb_t
mpz_getlimbn(const mpz_t u, mp_size_t n) {
  if (n >= 0 && n < MPI_ABS(u->_mp_size))
    return u->_mp_d[n];
  else
    return 0;
}

mp_srcptr
mpz_limbs_read(mpz_srcptr x) {
  return x->_mp_d;
}

mp_ptr
mpz_limbs_modify(mpz_t x, mp_size_t n) {
  assert(n > 0);
  return MPZ_REALLOC(x, n);
}

mp_ptr
mpz_limbs_write(mpz_t x, mp_size_t n) {
  return mpz_limbs_modify(x, n);
}

void
mpz_limbs_finish(mpz_t x, mp_size_t xs) {
  mp_size_t xn;
  xn = mpn_normalized_size(x->_mp_d, MPI_ABS(xs));
  x->_mp_size = xs < 0 ? -xn : xn;
}

mpz_srcptr
mpz_roinit_n(mpz_t x, mp_srcptr xp, mp_size_t xs) {
  mpz_roinit_normal_n(x, xp, xs);
  mpz_limbs_finish(x, xs);
  return x;
}

/*
 * Import
 */

void
mpz_import(mpz_t r, const unsigned char *u, size_t size, int endian) {
  mp_size_t rn;
  mp_ptr rp;

  if (size == 0) {
    r->_mp_size = 0;
    return;
  }

  rn = (size + sizeof(mp_limb_t) - 1) / sizeof(mp_limb_t);
  rp = MPZ_REALLOC(r, rn);

  mpn_import(rp, rn, u, size, endian);

  r->_mp_size = mpn_normalized_size(rp, rn);
}

/*
 * Export
 */

void
mpz_export(unsigned char *r, const mpz_t u, size_t size, int endian) {
  mpn_export(r, size, u->_mp_d, MPI_ABS(u->_mp_size), endian);
}

/*
 * Stringification
 */

size_t
mpz_out_str(FILE *stream, int base, const mpz_t x) {
  return mpn_out_str(stream, base, x->_mp_d, x->_mp_size);
}

/*
 * RNG
 */

void
mpz_random_bits(mpz_t r, mp_bitcnt_t bits, mp_rng_t rng, void *arg) {
  mp_size_t size = (bits + MPI_LIMB_BITS - 1) / MPI_LIMB_BITS;
  mp_size_t low = bits % MPI_LIMB_BITS;
  mp_ptr rp = MPZ_REALLOC(r, size);

  (*rng)(rp, size * sizeof(mp_limb_t), arg);

  if (low != 0)
    rp[size - 1] &= ((mp_limb_t)1 << low) - 1;

  r->_mp_size = mpn_normalized_size(rp, size);

  assert(mpz_bitlen(r) <= bits);
}

void
mpz_random_int(mpz_t ret, const mpz_t max, mp_rng_t rng, void *arg) {
  mp_bitcnt_t bits = mpz_bitlen(max);

  mpz_set(ret, max);

  if (bits > 0) {
    while (mpz_cmpabs(ret, max) >= 0)
      mpz_random_bits(ret, bits, rng, arg);

    if (mpz_sgn(max) < 0)
      mpz_neg(ret, ret);
  }
}
