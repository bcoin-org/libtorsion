/*!
 * internal.h - internal utils for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef TORSION_INTERNAL_H
#define TORSION_INTERNAL_H

/*
 * Clang Compat
 */

#if defined(__has_builtin)
#  define TORSION_HAS_BUILTIN __has_builtin
#else
#  define TORSION_HAS_BUILTIN(x) 0
#endif

/*
 * GNUC Compat
 */

#if defined(__GNUC__) && defined(__GNUC_MINOR__)
#  define TORSION_GNUC_PREREQ(maj, min) \
    ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
#  define TORSION_GNUC_PREREQ(maj, min) 0
#endif

/*
 * Builtins
 */

#undef LIKELY
#undef UNLIKELY
#undef UNPREDICTABLE

#if TORSION_GNUC_PREREQ(3, 0) || TORSION_HAS_BUILTIN(__builtin_expect)
#  define LIKELY(x) __builtin_expect(x, 1)
#  define UNLIKELY(x) __builtin_expect(x, 0)
#else
#  define LIKELY(x) (x)
#  define UNLIKELY(x) (x)
#endif

#if TORSION_HAS_BUILTIN(__builtin_unpredictable)
#  define UNPREDICTABLE __builtin_unpredictable
#else
#  define UNPREDICTABLE(x) (x)
#endif

/*
 * Sanity Checks
 */

#undef CHECK_ALWAYS
#undef CHECK_NEVER
#undef CHECK

#define CHECK_ALWAYS(expr) do { \
  if (UNLIKELY(!(expr)))        \
    torsion__abort();           \
} while (0)

#define CHECK_NEVER(expr) do { \
  (void)(expr);                \
} while (0)

#if !defined(TORSION_COVERAGE)
#  define CHECK CHECK_ALWAYS
#else
#  define CHECK CHECK_NEVER
#endif

/*
 * Assertions
 */

#undef ASSERT_ALWAYS
#undef ASSERT_NEVER
#undef ASSERT

#define ASSERT_ALWAYS(expr) do {                     \
  if (UNLIKELY(!(expr)))                             \
    torsion__assert_fail(__FILE__, __LINE__, #expr); \
} while (0)

#define ASSERT_NEVER(expr) do { \
  (void)(expr);                 \
} while (0)

#if defined(TORSION_DEBUG) && !defined(TORSION_COVERAGE)
#  define ASSERT ASSERT_ALWAYS
#else
#  define ASSERT ASSERT_NEVER
#endif

/*
 * Static Assertions
 */

#undef STATIC_ASSERT

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L \
                              && !defined(__cplusplus)
#  undef _Static_assert
#  define STATIC_ASSERT(expr) _Static_assert(expr, "check failed")
#elif defined(__cplusplus) && (__cplusplus + 0L) >= 201703L
#  define STATIC_ASSERT(expr) static_assert(expr)
#elif defined(__cplusplus) && (__cplusplus + 0L) >= 201103L
#  define STATIC_ASSERT(expr) static_assert(expr, "check failed")
#elif TORSION_GNUC_PREREQ(2, 7)
#  define STATIC_ASSERT_2(x, y) \
     typedef char torsion__assert_ ## y[(x) ? 1 : -1] __attribute__((unused))
#  define STATIC_ASSERT_1(x, y) STATIC_ASSERT_2(x, y)
#  define STATIC_ASSERT(expr) STATIC_ASSERT_1(expr, __LINE__)
#else
#  define STATIC_ASSERT(expr) struct torsion__assert_empty
#endif

/*
 * Keywords/Attributes
 */

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L \
                              && !defined(__cplusplus)
#  define TORSION_INLINE inline
#elif defined(__cplusplus) && (__cplusplus + 0L) >= 199711L
#  define TORSION_INLINE inline
#elif TORSION_GNUC_PREREQ(2, 7)
#  define TORSION_INLINE __inline__
#elif defined(_MSC_VER) && _MSC_VER >= 900
#  define TORSION_INLINE __inline
#elif (defined(__SUNPRO_C) && __SUNPRO_C >= 0x560) \
   || (defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x560)
#  define TORSION_INLINE inline
#else
#  define TORSION_INLINE
#endif

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L \
                              && !defined(__cplusplus)
#  define TORSION_RESTRICT restrict
#elif TORSION_GNUC_PREREQ(3, 0)
#  define TORSION_RESTRICT __restrict__
#elif defined(_MSC_VER) && _MSC_VER >= 1400
#  define TORSION_RESTRICT __restrict
#elif defined(__SUNPRO_C) && __SUNPRO_C >= 0x530
#  define TORSION_RESTRICT _Restrict
#else
#  define TORSION_RESTRICT
#endif

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L \
                              && !defined(__cplusplus)
#  define TORSION_NORETURN _Noreturn
#elif defined(__cplusplus) && (__cplusplus + 0L) >= 201103L
#  undef noreturn
#  define TORSION_NORETURN [[noreturn]]
#elif TORSION_GNUC_PREREQ(2, 7)
#  undef noreturn
#  define TORSION_NORETURN __attribute__((noreturn))
#elif defined(_MSC_VER) && _MSC_VER >= 1200
#  undef noreturn
#  define TORSION_NORETURN __declspec(noreturn)
#elif (defined(__SUNPRO_C) && __SUNPRO_C >= 0x590) \
   || (defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x590)
#  undef noreturn
#  define TORSION_NORETURN __attribute__((noreturn))
#else
#  define TORSION_NORETURN
#endif

#if defined(__cplusplus) && (__cplusplus + 0L) >= 201703L
#  define TORSION_UNUSED [[maybe_unused]]
#elif TORSION_GNUC_PREREQ(2, 7)
#  define TORSION_UNUSED __attribute__((unused))
#else
#  define TORSION_UNUSED
#endif

#if defined(__GNUC__)
#  define TORSION_EXTENSION __extension__
#else
#  define TORSION_EXTENSION
#endif

/*
 * Endianness
 */

/* Any decent compiler should be able to optimize this out. */
static const unsigned long torsion__endian_check TORSION_UNUSED = 1;

#define TORSION_BIGENDIAN \
  (*((const unsigned char *)&torsion__endian_check) == 0)

/*
 * Configuration
 */

#ifndef TORSION_HAVE_CONFIG
/* TORSION_HAVE_CONFIG signals that the config
 * will be passed in via the commandline (-D).
 * Otherwise, auto configuration is useful if
 * you're using an awful build system like gyp.
 *
 * Start by clearing everything...
 */
#undef TORSION_HAVE_ASM
#undef TORSION_HAVE_ASM_X86
#undef TORSION_HAVE_ASM_X64
#undef TORSION_HAVE_INT128
#undef TORSION_HAVE_TLS
#undef TORSION_TLS

/* Detect inline ASM support for x86-64.
 *
 * GCC inline assembly has been documented as
 * far back as 2.95[1]. It appears in the GCC
 * codebase as early as 2.0. However, early
 * implementations may not have the features
 * we require, so to be practical, we require
 * GNUC version 4.0.
 *
 * [1] https://gcc.gnu.org/onlinedocs/gcc-2.95.3/gcc_4.html#SEC93
 */
#if TORSION_GNUC_PREREQ(4, 0)
#  define TORSION_HAVE_ASM
#  if defined(__amd64__) || defined(__x86_64__)
#    define TORSION_HAVE_ASM_X64
#  elif defined(__i386__)
#    define TORSION_HAVE_ASM_X86
#  endif
#endif

/* Detect __int128 support.
 *
 * Support (verified on godbolt):
 *
 *   x86-64:
 *     gcc 4.6.4 (gnuc 4.6.4)
 *     clang 3.1 (gnuc 4.2.1) (__SIZEOF_INT128__ defined in 3.3.0)
 *     icc <=13.0.1 (gnuc 4.7.0) (__SIZEOF_INT128__ defined in 16.0.3)
 *
 *   arm64:
 *     gcc <=5.4.0 (gnuc 5.4.9)
 *     clang <=9.0 (gnuc 4.2.1)
 *
 *   mips64:
 *     gcc <=5.4.0 (gnuc 5.4.9)
 *
 *   power64:
 *     gcc <=6.3.0 (gnuc 6.3.0)
 *     clang <=12.0.0 (gnuc 4.2.1)
 *     at <=12.0.0 (gnuc 8.2.1)
 *
 *   riscv64:
 *     gcc <=8.2.0 (gnuc 8.2.0)
 *     clang <=12.0.0 (gnuc 4.2.1)
 *
 *   wasm32/wasm64:
 *     clang <=7.0 (gnuc 4.2.1)
 *
 * See: https://stackoverflow.com/a/54815033
 */
#if defined(__GNUC__) && defined(__SIZEOF_INT128__)  \
                      && defined(__SIZEOF_POINTER__)
#  if __SIZEOF_POINTER__ >= 8
#    define TORSION_HAVE_INT128
#  endif
#endif

/* Basically a stripped down version of our old file[1].
 * It only includes the compilers we for sure know work.
 *
 * [1] https://github.com/bcoin-org/libtorsion/blob/2fe6cd3/src/tls.h
 */
#if defined(__clang__) || defined(__llvm__)
#  ifdef __has_extension
#    if __has_extension(c_thread_local)
#      if defined(_MSC_VER) || defined(__BORLANDC__)
#        define TORSION_TLS __declspec(thread)
#      elif defined(__ANDROID__)
#        if defined(__clang_major__) && __clang_major__ >= 5
#          define TORSION_TLS __thread
#        endif
#      else
#        define TORSION_TLS __thread
#      endif
#    endif
#  endif
#elif defined(__INTEL_COMPILER)
#  if defined(_WIN32) && __INTEL_COMPILER >= 1000
#    define TORSION_TLS __declspec(thread)
#  elif defined(__linux__) && __INTEL_COMPILER >= 810
#    if TORSION_GNUC_PREREQ(3, 3)
#      define TORSION_TLS __thread
#    endif
#  elif defined(__APPLE__) && __INTEL_COMPILER >= 1500
#    define TORSION_TLS __thread
#  endif
#elif defined(__GNUC__) && !defined(__CC_ARM) \
                        && !defined(__PCC__)  \
                        && !defined(__NWCC__)
#  if TORSION_GNUC_PREREQ(4, 3)
#    define TORSION_TLS __thread
#  elif TORSION_GNUC_PREREQ(3, 3)
#    if defined(__ELF__) && (defined(__i386__) || defined(__x86_64__))
#      define TORSION_TLS __thread
#    endif
#  endif
#elif (defined(_MSC_VER) && _MSC_VER >= 1200)          \
   || (defined(__WATCOMC__) && __WATCOMC__ >= 1200)    \
   || (defined(__BORLANDC__) && __BORLANDC__ >= 0x610) \
   || (defined(__DMC__) && __DMC__ >= 0x822)
#  define TORSION_TLS __declspec(thread)
#elif (defined(__SUNPRO_C) && __SUNPRO_C >= 0x590)     \
   || (defined(__SUNPRO_CC) && __SUNPRO_CC >= 0x590)   \
   || (defined(__HP_cc) && __HP_cc >= 53600)           \
   || (defined(__HP_aCC) && __HP_aCC >= 53600)         \
   || (defined(__CC_ARM) && __ARMCC_VERSION >= 510000) \
   || (defined(__PCC__) && __PCC__ >= 1)               \
   || (defined(__NWCC__))
#  define TORSION_TLS __thread
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L \
                                && !defined(__cplusplus)
#  ifndef __STDC_NO_THREADS__
#    define TORSION_TLS _Thread_local
#  endif
#elif defined(__cplusplus) && (__cplusplus + 0L) >= 201103L
#  define TORSION_TLS thread_local
#endif

#if defined(TORSION_TLS)
#  define TORSION_HAVE_TLS
#else
#  define TORSION_TLS
#endif

/* Allow some overrides. */
#ifdef TORSION_NO_ASM
#  undef TORSION_HAVE_ASM
#  undef TORSION_HAVE_ASM_X86
#  undef TORSION_HAVE_ASM_X64
#endif

#ifdef TORSION_NO_INT128
#  undef TORSION_HAVE_INT128
#endif

#ifdef TORSION_NO_TLS
#  undef TORSION_HAVE_TLS
#  undef TORSION_TLS
#  define TORSION_TLS
#endif

#ifdef TORSION_FORCE_32BIT
#  undef TORSION_HAVE_ASM_X64
#  undef TORSION_HAVE_INT128
#endif

#endif /* !TORSION_HAVE_CONFIG */

/*
 * Types
 */

#ifdef TORSION_HAVE_INT128
TORSION_EXTENSION typedef unsigned __int128 torsion_uint128_t;
TORSION_EXTENSION typedef signed __int128 torsion_int128_t;
#endif

/*
 * Value Barrier
 */

#if defined(TORSION_HAVE_ASM)
#define TORSION_BARRIER(type, prefix) \
static TORSION_INLINE type            \
prefix ## _barrier(type x) {          \
  __asm__ ("" : "+r" (x));            \
  return x;                           \
}
#else
#define TORSION_BARRIER(type, prefix) \
static TORSION_INLINE type            \
prefix ## _barrier(type x) {          \
  return x;                           \
}
#endif

/*
 * Sanity Checks
 */

#if (-1 & 3) != 3
#  error "Two's complement is required."
#endif

/*
 * Macros
 */

#define ENTROPY_SIZE 32
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

/*
 * Helpers
 */

#define torsion_abort torsion__abort

TORSION_NORETURN void
torsion__assert_fail(const char *file, int line, const char *expr);

TORSION_NORETURN void
torsion__abort(void);

/*
 * Character Transcoding
 */

extern const int torsion__ascii[256];
extern const int torsion__native[256];

/* We could check the character set in preprocessor, but the
 * standard has some very strange wording around character
 * constants in preprocessor. Specifically, the standard says,
 *
 *   "Whether the numeric value for these character constants
 *    matches the value obtained when an identical character
 *    constant occurs in an expression (other than within a
 *    #if or #elif directive) is implementation-defined."[1]
 *
 * I suppose this can be taken to mean that the preprocessor
 * may use the source character set instead of the execution
 * character set. Vague wording like this has often been the
 * justification for compiler developers to do wacky stuff,
 * so we instead check the character set at "runtime". Every
 * compiler should treat this as a constant expression and
 * optimize it out.
 *
 * [1] ANSI/ISO 9899-1990, Page 87, Section 6.8.1 ("Conditional Inclusion")
 */
#define torsion_a (' ' == 32 && '0' == 48 && 'A' == 65 && 'a' == 97)
#define torsion_ascii(c) (torsion_a ? ((c) & 0xff) : torsion__ascii[(c) & 0xff])
#define torsion_native(c) (torsion_a ? (c) : torsion__native[(c) & 0xff])

#endif /* TORSION_INTERNAL_H */
