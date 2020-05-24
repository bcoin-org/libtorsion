AC_PREREQ(2.59)
AC_INIT([libtorsion], [0.0.0], [https://github.com/bcoin-org/libtorsion/issues])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([include/torsion/ecc.h])
m4_include([m4/attributes.m4])
m4_include([m4/as_case.m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CANONICAL_HOST
AC_ENABLE_SHARED
AC_ENABLE_STATIC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C89

if test x"$ac_cv_prog_cc_c89" = x"no"; then
  AC_MSG_ERROR([CC does not support ANSI C89])
fi

CC_CHECK_CFLAGS_APPEND([-std=c89])
CC_CHECK_CFLAGS_APPEND([-pedantic])
CC_CHECK_CFLAGS_APPEND([-Wall])
CC_CHECK_CFLAGS_APPEND([-Wextra])
CC_CHECK_CFLAGS_APPEND([-Wshadow])
CC_CHECK_CFLAGS_APPEND([-Wcast-align])
CC_CHECK_CFLAGS_APPEND([-Wno-implicit-fallthrough])
CC_CHECK_CFLAGS_APPEND([-Wno-declaration-after-statement])
CC_CHECK_CFLAGS_APPEND([-Wno-overlength-strings])
CC_CHECK_CFLAGS_APPEND([-Wno-long-long])

AC_MSG_CHECKING(for 64 bit support)
AS_CASE([$host_cpu],
  [aarch64*|ia64|mips64*|ppc64*|sparc64|x86_64], [has_64bit=yes],
  [*], [has_64bit=no])
AC_MSG_RESULT([$has_64bit])

AC_MSG_CHECKING(for x86_64 inline asm support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdint.h>]],[[
#ifdef __EMSCRIPTEN__
#  error "no asm"
#endif
#if !defined(__amd64__) && !defined(__x86_64__)
#  error "no asm"
#endif
#if defined(__CYGWIN__) || defined(__MINGW32__)
#  error "no asm"
#endif
  uint64_t u = 2, v = 3, w0, w1;
  __asm__ __volatile__(
    "mulq %q3\n"
    : "=a" (w0), "=d" (w1)
    : "%0" (u), "rm" (v)
  );
]])],[has_asm=yes],[has_asm=no])
AC_MSG_RESULT([$has_asm])

AC_C_BIGENDIAN()

AC_MSG_CHECKING(for int128 support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#ifdef __EMSCRIPTEN__
#  error "no int128"
#endif
#ifdef __GNUC__
#  define TORSION_EXTENSION __extension__
#else
#  define TORSION_EXTENSION
#endif
TORSION_EXTENSION typedef signed __int128 int128_t;
TORSION_EXTENSION typedef unsigned __int128 uint128_t;
void compute_int128(uint128_t *r) {
  *r = ((uint128_t)1) << 64;
  *r *= 256;
  *r >>= 65;
}
]],[[
  uint128_t r;
  compute_int128(&r);
]])],[has_int128=yes],[has_int128=no])
AC_MSG_RESULT([$has_int128])

if test x"$has_int128" = x"yes" -a x"$has_64bit" != x"yes"; then
  has_int128=no
fi

AC_MSG_CHECKING(for TI mode)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#ifdef __EMSCRIPTEN__
#  error "no mode(TI)"
#endif
typedef unsigned uint128_t __attribute__((mode(TI)));
typedef signed int128_t __attribute__((mode(TI)));
void compute_modeti(uint128_t *r) {
  *r = ((uint128_t)1) << 64;
  *r *= 256;
  *r >>= 65;
}
]],[[
  uint128_t r;
  compute_modeti(&r);
]])],[has_modeti=yes],[has_modeti=no])
AC_MSG_RESULT([$has_modeti])

if test x"$has_modeti" = x"yes" -a x"$has_64bit" != x"yes"; then
  has_modeti=no
fi

AC_MSG_CHECKING(for sse2 support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdint.h>]],[[
#ifdef __EMSCRIPTEN__
#  error "no asm"
#endif
#if !defined(__amd64__) && !defined(__x86_64__)
#  error "no asm"
#endif
#if defined(__CYGWIN__) || defined(__MINGW32__)
#  error "no asm"
#endif
  uint32_t stream[8], state[8];
  __asm__ __volatile__(
    "movups (%%rsi), %%xmm0\n"
    "movups 16(%%rsi), %%xmm1\n"
    "pxor %%xmm1, %%xmm0\n"
    "movups %%xmm0,(%%rdi)\n"
    "movups %%xmm1,16(%%rdi)\n"
    :
    : "D" (stream), "S" (state)
    : "xmm0", "xmm1", "cc", "memory"
  );
]])],[has_sse2=yes],[has_sse2=no])
AC_MSG_RESULT([$has_sse2])

AC_CHECK_HEADER([valgrind/memcheck.h], [has_valgrind=yes], [has_valgrind=no], [])

AC_ARG_WITH([64bit], [AS_HELP_STRING([--with-64bit=yes|no|auto],
[enable to use 64 bit backends [default=auto]])],[use_64bit=$withval], [use_64bit=auto])

AC_ARG_WITH([asm], [AS_HELP_STRING([--with-asm=yes|no|auto],
[enable to use x86_64 asm [default=auto]])],[use_asm=$withval], [use_asm=auto])

AC_ARG_ENABLE(assert,
  AS_HELP_STRING([--enable-assert],[enable assertions [default=yes]]),
  [enable_assert=$enableval],
  [enable_assert=yes])

AC_ARG_ENABLE(crt,
  AS_HELP_STRING([--enable-crt],[enable chinese remainder theorem for RSA [default=no]]),
  [enable_crt=$enableval],
  [enable_crt=no])

AC_ARG_ENABLE(ctime,
  AS_HELP_STRING([--enable-ctime],[enable valgrind constant-time test [default=no]]),
  [enable_ctime=$enableval],
  [enable_ctime=no])

AC_ARG_ENABLE(libsecp256k1,
  AS_HELP_STRING([--enable-libsecp256k1],[enable libsecp256k1 field element backend [default=no]]),
  [enable_libsecp256k1=$enableval],
  [enable_libsecp256k1=no])

AC_ARG_ENABLE(rng,
  AS_HELP_STRING([--enable-rng],[enable rng [default=no]]),
  [enable_rng=$enableval],
  [enable_rng=no])

AC_ARG_ENABLE(tests,
  AS_HELP_STRING([--enable-tests],[enable tests [default=no]]),
  [enable_tests=$enableval],
  [enable_tests=no])

AC_ARG_ENABLE(verify,
  AS_HELP_STRING([--enable-verify],[enable scalar bounds checks [default=no]]),
  [enable_verify=$enableval],
  [enable_verify=no])

if test x"$use_64bit" = x"auto"; then
  if test x"$EMSCRIPTEN" != x""; then
    use_64bit=no
  else
    use_64bit=$has_64bit
  fi
fi

if test x"$use_asm" = x"auto"; then
  if test x"$EMSCRIPTEN" != x""; then
    use_asm=no
  else
    use_asm=$has_asm
  fi
fi

if test x"$enable_ctime" = x"yes"; then
  enable_assert=no
  enable_verify=no
  enable_tests=yes
fi

AC_DEFINE(TORSION_HAVE_CONFIG, 1, [Define to pass in config])

if test x"$use_64bit" = x"yes"; then
  if test x"$has_64bit" != x"yes"; then
    AC_MSG_ERROR([64 bit support not available])
  fi

  AC_DEFINE(TORSION_HAVE_64BIT, 1, [Define to enable 64 bit backends])
fi

if test x"$use_asm" = x"yes"; then
  if test x"$has_asm" != x"yes"; then
    AC_MSG_ERROR([x86_64 asm not available])
  fi

  if test x"$host_cpu" != x"x86_64"; then
    AC_MSG_ERROR([asm only supported for x86_64])
  fi

  AC_DEFINE(TORSION_HAVE_ASM_X64, 1, [Define to enable x86_64 asm])
fi

if test x"$enable_assert" != x"yes"; then
  AC_DEFINE(TORSION_NO_ASSERT, 1, [Define to disable assertions])
fi

if test x"$ac_cv_c_bigendian" = x"yes"; then
  AC_DEFINE(TORSION_BIGENDIAN, 1, [Define if machine is big endian])
fi

if test x"$enable_crt" = x"yes"; then
  AC_DEFINE(TORSION_USE_CRT, 1, [Define to enable RSA CRT])
fi

if test x"$enable_ctime" = x"yes"; then
  if test x"$has_valgrind" != x"yes"; then
    AC_MSG_ERROR([valgrind must be installed for constant time test])
  fi
  AC_DEFINE(TORSION_VALGRIND, 1, [Define to make the ECC code valgrind-aware])
fi

if test x"$has_int128" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_INT128, 1, [Define to enable __int128 backends])
fi

if test x"$enable_libsecp256k1" = x"yes"; then
  AC_DEFINE(TORSION_USE_LIBSECP256K1, 1, [Define to enable the libsecp256k1 field element backend])
fi

if test x"$has_modeti" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_MODETI, 1, [Define to enable TI mode])
fi

if test x"$enable_rng" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_RNG, 1, [Define to notify tests of RNG support])
fi

if test x"$has_sse2" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_ASM_SSE2, 1, [Define to enable SSE2 instructions])
fi

if test x"$enable_tests" = x"yes"; then
  AC_DEFINE(TORSION_TEST, 1, [Define to enable internal tests])
fi

if test x"$enable_verify" = x"yes"; then
  if test x"$enable_tests" != x"yes"; then
    AC_MSG_ERROR([scalar bounds verification for testing only])
  fi
  AC_DEFINE(TORSION_VERIFY, 1, [Define to enable scalar bounds checks])
fi

AM_CONDITIONAL([ENABLE_RNG], [test x"$enable_rng" = x"yes"])
AM_CONDITIONAL([ENABLE_TESTS], [test x"$enable_tests" = x"yes"])
AM_CONDITIONAL([ENABLE_CTIME], [test x"$enable_ctime" = x"yes"])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
LT_INIT
AC_CONFIG_FILES([Makefile libtorsion.pc])
AC_OUTPUT

echo ""
echo "Build Options:"
echo ""
echo "  64bit        = $use_64bit"
echo "  asm          = $use_asm"
echo "  assert       = $enable_assert"
echo "  big endian   = $ac_cv_c_bigendian"
echo "  crt          = $enable_crt"
echo "  ctime        = $enable_ctime"
echo "  int128       = $has_int128"
echo "  libsecp256k1 = $enable_libsecp256k1"
echo "  modeti       = $has_modeti"
echo "  rng          = $enable_rng"
echo "  sse2         = $has_sse2"
echo "  tests        = $enable_tests"
echo "  verify       = $enable_verify"
echo ""
echo "  CC           = $CC"
echo "  CFLAGS       = $CFLAGS"
echo "  LDFLAGS      = $LDFLAGS"
echo ""
