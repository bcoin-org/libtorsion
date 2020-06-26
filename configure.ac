AC_PREREQ(2.59)
AC_INIT([libtorsion], [0.0.0], [https://github.com/bcoin-org/libtorsion/issues])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([include/torsion/ecc.h])
m4_include([m4/attributes.m4])
m4_include([m4/as_case.m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
: ${CFLAGS=""}
AC_CANONICAL_HOST
AC_ENABLE_SHARED
AC_ENABLE_STATIC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C89

if test x"$ac_cv_prog_cc_c89" = x"no"; then
  AC_MSG_ERROR([CC does not support ANSI C89])
fi

AC_ARG_ENABLE(64bit,
  AS_HELP_STRING([--enable-64bit],[use 64-bit backends if available [default=yes]]),
  [enable_64bit=$enableval],
  [enable_64bit=yes])

AC_ARG_ENABLE(asm,
  AS_HELP_STRING([--enable-asm],[use inline x86-64 assembly if available [default=yes]]),
  [enable_asm=$enableval],
  [enable_asm=yes])

AC_ARG_ENABLE(int128,
  AS_HELP_STRING([--enable-int128],[use __int128 if available [default=yes]]),
  [enable_int128=$enableval],
  [enable_int128=yes])

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
  AS_HELP_STRING([--enable-rng],[enable rng [default=yes]]),
  [enable_rng=$enableval],
  [enable_rng=yes])

AC_ARG_ENABLE(rand,
  AS_HELP_STRING([--enable-rand],[enable global rng [default=yes]]),
  [enable_rand=$enableval],
  [enable_rand=yes])

AC_ARG_ENABLE(tests,
  AS_HELP_STRING([--enable-tests],[enable tests [default=yes]]),
  [enable_tests=$enableval],
  [enable_tests=yes])

AC_ARG_ENABLE(verify,
  AS_HELP_STRING([--enable-verify],[enable scalar bounds checks [default=no]]),
  [enable_verify=$enableval],
  [enable_verify=no])

if test x"$enable_ctime" = x"yes"; then
  enable_assert=no
  enable_verify=no
fi

if test x"$EMSCRIPTEN" != x""; then
  CC_CHECK_CFLAGS_APPEND([-std=gnu11])
else
  CC_CHECK_CFLAGS_APPEND([-std=c89])
  CC_CHECK_CFLAGS_APPEND([-pedantic])
  CC_CHECK_CFLAGS_APPEND([-fvisibility=hidden])
fi

CC_CHECK_CFLAGS_APPEND([-Wall])
CC_CHECK_CFLAGS_APPEND([-Wextra])
CC_CHECK_CFLAGS_APPEND([-Wcast-align])
CC_CHECK_CFLAGS_APPEND([-Wno-declaration-after-statement])
CC_CHECK_CFLAGS_APPEND([-Wno-implicit-fallthrough])
CC_CHECK_CFLAGS_APPEND([-Wno-long-long])
CC_CHECK_CFLAGS_APPEND([-Wno-overlength-strings])
CC_CHECK_CFLAGS_APPEND([-Wshadow])
CC_CHECK_CFLAGS_APPEND([-O3])

AC_MSG_CHECKING(for 64 bit support)
AS_CASE([$host_cpu],
  [aarch64*|ia64|mips64*|ppc64*|sparc64|x86_64], [has_64bit=yes],
  [*], [has_64bit=no])
AC_MSG_RESULT([$has_64bit])

AC_MSG_CHECKING(for x86_64 inline asm support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdint.h>]],[[
#if !defined(__amd64__) && !defined(__x86_64__)
#  error "not an x64 platform"
#endif
#if defined(__CYGWIN__) || defined(__MINGW32__)
#  error "no asm on unix lookalikes"
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
]])],[has_asm=yes],[has_asm=no])
AC_MSG_RESULT([$has_asm])

AC_MSG_CHECKING(for int128 support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
typedef signed __int128 int128_t;
typedef unsigned __int128 uint128_t;
]],[[
uint128_t r;
r = ((uint128_t)1) << 64;
r *= 113;
r >>= 65;
]])],[has_int128=yes],[has_int128=no])
AC_MSG_RESULT([$has_int128])

AC_CHECK_HEADER([valgrind/memcheck.h], [has_valgrind=yes], [has_valgrind=no], [])

AC_DEFINE(TORSION_HAVE_CONFIG)

if test x"$EMSCRIPTEN" != x""; then
  has_64bit=no
  has_asm=no
  has_int128=no
  enable_assert=no
fi

if test x"$enable_64bit" = x"yes" -a x"$has_64bit" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_64BIT)
  if test x"$enable_asm" = x"yes" -a x"$has_asm" = x"yes"; then
    AC_DEFINE(TORSION_HAVE_ASM_X64)
  fi
  if test x"$enable_int128" = x"yes" -a x"$has_int128" = x"yes"; then
    AC_DEFINE(TORSION_HAVE_INT128)
  fi
fi

if test x"$enable_assert" != x"yes"; then
  AC_DEFINE(TORSION_NO_ASSERT)
fi

if test x"$enable_crt" = x"yes"; then
  AC_DEFINE(TORSION_USE_CRT)
fi

if test x"$enable_ctime" = x"yes"; then
  if test x"$has_valgrind" != x"yes"; then
    AC_MSG_ERROR([valgrind must be installed for constant time test])
  fi
  AC_DEFINE(TORSION_VALGRIND)
fi

if test x"$enable_libsecp256k1" = x"yes"; then
  AC_DEFINE(TORSION_USE_LIBSECP256K1)
fi

if test x"$enable_rng" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_RNG)
  if test x"$enable_rand" = x"yes"; then
    AC_DEFINE(TORSION_HAVE_RAND)
  fi
else
  enable_rand=no
fi

if test x"$enable_tests" = x"yes"; then
  if test x"$EMSCRIPTEN" = x""; then
    CC_CHECK_CFLAGS_APPEND([-g])
  fi
  AC_DEFINE(TORSION_TEST)
fi

if test x"$enable_verify" = x"yes"; then
  if test x"$enable_tests" != x"yes"; then
    AC_MSG_ERROR([scalar bounds verification for testing only])
  fi
  AC_DEFINE(TORSION_VERIFY)
fi

AM_CONDITIONAL([ENABLE_RNG], [test x"$enable_rng" = x"yes"])
AM_CONDITIONAL([ENABLE_RAND], [test x"$enable_rand" = x"yes"])
AM_CONDITIONAL([ENABLE_TESTS], [test x"$enable_tests" = x"yes"])
AM_CONDITIONAL([ENABLE_CTIME], [test x"$enable_ctime" = x"yes"])
AM_CONDITIONAL([EMSCRIPTEN], [test x"$EMSCRIPTEN" != x""])

AM_CONDITIONAL([DARWIN], [AS_CASE([$host_os],[darwin*], [true], [false])])
AM_CONDITIONAL([LINUX], [AS_CASE([$host_os],[linux*], [true], [false])])

m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
LT_INIT
AC_CONFIG_FILES([Makefile libtorsion.pc])
AC_OUTPUT

echo ""
echo "Build Options:"
echo ""
echo "  64bit        = $enable_64bit (support=$has_64bit)"
echo "  asm          = $enable_asm (support=$has_asm)"
echo "  assert       = $enable_assert"
echo "  crt          = $enable_crt"
echo "  ctime        = $enable_ctime"
echo "  emscripten   = $EMSCRIPTEN"
echo "  int128       = $enable_int128 (support=$has_int128)"
echo "  libsecp256k1 = $enable_libsecp256k1"
echo "  rng          = $enable_rng"
echo "  rand         = $enable_rand"
echo "  tests        = $enable_tests"
echo "  verify       = $enable_verify"
echo ""
echo "  CC           = $CC"
echo "  CFLAGS       = $CFLAGS"
echo "  LDFLAGS      = $LDFLAGS"
echo ""
