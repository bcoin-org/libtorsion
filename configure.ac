dnl configure.ac - autoconf script for libtorsion
dnl Copyright (c) 2020, Christopher Jeffrey (MIT License).
dnl https://github.com/bcoin-org/libtorsion

AC_PREREQ(2.59)
AC_INIT([libtorsion], [0.0.0], [https://github.com/bcoin-org/libtorsion/issues])
LIB_VERSION=0:0:0
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([include/torsion/ecc.h])
m4_include([m4/attributes.m4])
m4_include([m4/as_case.m4])
m4_include([m4/ax_pthread.m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
: ${CFLAGS=""}
AC_CANONICAL_HOST
AC_ENABLE_SHARED
AC_ENABLE_STATIC
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C89
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
LT_INIT

if test x"$ac_cv_prog_cc_c89" = x"no"; then
  AC_MSG_ERROR([CC does not support ANSI C89])
fi

if test x"$EMSCRIPTEN" != x"" -o x"$WASI" != x""; then
  WASM=yes
else
  WASM=no
fi

AC_ARG_ENABLE(asm,
  AS_HELP_STRING([--enable-asm], [use inline x86-64 assembly if available [default=yes]]),
  [enable_asm=$enableval],
  [enable_asm=yes])

AC_ARG_ENABLE(coverage,
  AS_HELP_STRING([--enable-coverage], [enable coverage [default=no]]),
  [enable_coverage=$enableval],
  [enable_coverage=no])

AC_ARG_ENABLE(crt,
  AS_HELP_STRING([--enable-crt], [enable chinese remainder theorem for RSA [default=no]]),
  [enable_crt=$enableval],
  [enable_crt=no])

AC_ARG_ENABLE(debug,
  AS_HELP_STRING([--enable-debug], [enable tests [default=no]]),
  [enable_debug=$enableval],
  [enable_debug=no])

AC_ARG_ENABLE(int128,
  AS_HELP_STRING([--enable-int128], [use __int128 if available [default=yes]]),
  [enable_int128=$enableval],
  [enable_int128=yes])

AC_ARG_ENABLE(libsecp256k1,
  AS_HELP_STRING([--enable-libsecp256k1], [use libsecp256k1 field element backend [default=no]]),
  [enable_libsecp256k1=$enableval],
  [enable_libsecp256k1=no])

AC_ARG_ENABLE(pthread,
  AS_HELP_STRING([--enable-pthread], [use pthread as a fallback for tls [default=yes]]),
  [enable_pthread=$enableval],
  [enable_pthread=yes])

AC_ARG_ENABLE(rng,
  AS_HELP_STRING([--enable-rng], [enable rng [default=yes]]),
  [enable_rng=$enableval],
  [enable_rng=yes])

AC_ARG_ENABLE(tls,
  AS_HELP_STRING([--enable-tls], [enable tls [default=yes]]),
  [enable_tls=$enableval],
  [enable_tls=yes])

AC_ARG_ENABLE(verify,
  AS_HELP_STRING([--enable-verify], [enable scalar bounds checks [default=no]]),
  [enable_verify=$enableval],
  [enable_verify=no])

if test x"$EMSCRIPTEN" != x""; then
  CC_CHECK_CFLAGS_APPEND([-std=gnu11])
else
  CC_CHECK_CFLAGS_APPEND([-std=c89])
  CC_CHECK_CFLAGS_APPEND([-pedantic])
fi

CC_CHECK_CFLAGS_APPEND([-Wall])
CC_CHECK_CFLAGS_APPEND([-Wextra])
CC_CHECK_CFLAGS_APPEND([-Wcast-align])
CC_CHECK_CFLAGS_APPEND([-Wno-implicit-fallthrough])
CC_CHECK_CFLAGS_APPEND([-Wno-long-long])
CC_CHECK_CFLAGS_APPEND([-Wno-overlength-strings])
CC_CHECK_CFLAGS_APPEND([-Wshadow])
CC_CHECK_CFLAGS_APPEND([-O3])

if test x"$WASM" != x"yes"; then
  CC_CHECK_CFLAGS_APPEND([-fvisibility=hidden])

  AS_CASE([$host_os], [darwin*],
    [CC_CHECK_CFLAGS_APPEND([-mmacosx-version-min=10.7])], [])

  AS_CASE([$host_os], [linux*],
    [CFLAGS="$CFLAGS -D_POSIX_C_SOURCE=200112"], [])
fi

AC_MSG_CHECKING(for inline asm support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], [[
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  error "cannot use inline asm with wasm"
#endif
unsigned char ptr[32];
int x = 0;
__asm__ ("" : "+r" (x) ::);
__asm__ __volatile__("" :: "r" (ptr) : "memory");
]])], [has_asm=yes], [has_asm=no])
AC_MSG_RESULT([$has_asm])

AC_MSG_CHECKING(for x86_64 inline asm support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <stdint.h>]], [[
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  error "cannot use inline asm with wasm"
#endif
#if !defined(__amd64__) && !defined(__x86_64__)
#  error "not an x64 platform"
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
]])], [has_asm_x64=yes], [has_asm_x64=no])
AC_MSG_RESULT([$has_asm_x64])

AC_MSG_CHECKING(for int128 support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  error "cannot use __int128 with wasm"
#endif
typedef char check_64bit_t[sizeof(void *) == 8 ? 1 : -1];
typedef signed __int128 int128_t;
typedef unsigned __int128 uint128_t;
]], [[
uint128_t r = 1;
r <<= 64;
r *= 113;
r >>= 65;
]])], [has_int128=yes], [has_int128=no])
AC_MSG_RESULT([$has_int128])

AC_MSG_CHECKING(for tls support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  define TLS
#elif defined(_WIN32) && !defined(__MINGW32__)
#  define TLS __declspec(thread)
#else
#  define TLS __thread
#endif
static TLS int value;
]], [[
value = 1;
]])], [has_tls=yes], [has_tls=no])
AC_MSG_RESULT([$has_tls])

if test x"$WASM" != x"yes"; then
  AX_PTHREAD([has_threads=yes], [has_threads=no])
else
  has_threads=no
fi

AC_MSG_CHECKING(for fork support)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  error "cannot use fork with wasm"
#endif
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
]], [[
unsigned char data[32];
int pfds[2];
int status;
pid_t pid;
pipe(pfds);
pid = fork();
close(pfds[1]);
read(pfds[0], data, 32);
close(pfds[0]);
waitpid(pid, &status, 0);
WIFEXITED(status);
WEXITSTATUS(status);
]])], [has_fork=yes], [has_fork=no])
AC_MSG_RESULT([$has_fork])

AC_CHECK_HEADER([valgrind/memcheck.h], [has_valgrind=yes], [has_valgrind=no], [])

AC_DEFINE(TORSION_HAVE_CONFIG)

if test x"$enable_asm" = x"yes" -a x"$has_asm" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_ASM)
fi

if test x"$enable_asm" = x"yes" -a x"$has_asm_x64" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_ASM_X64)
fi

if test x"$enable_coverage" = x"yes"; then
  if test x"$WASM" = x"yes"; then
    AC_MSG_ERROR([coverage not supported with wasm])
  fi
  CFLAGS="$CFLAGS -O0 --coverage"
  LDFLAGS="$LDFLAGS --coverage"
  AC_DEFINE(TORSION_COVERAGE)
fi

if test x"$enable_crt" = x"yes"; then
  AC_DEFINE(TORSION_USE_CRT)
fi

if test x"$enable_debug" = x"yes"; then
  CC_CHECK_CFLAGS_APPEND([-g])
  AC_DEFINE(TORSION_DEBUG)
fi

if test x"$has_fork" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_FORK)
fi

if test x"$enable_int128" = x"yes" -a x"$has_int128" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_INT128)
fi

if test x"$enable_libsecp256k1" = x"yes"; then
  AC_DEFINE(TORSION_USE_LIBSECP256K1)
fi

if test x"$enable_pthread" = x"yes" -a x"$has_threads" = x"yes"; then
  if test x"$enable_tls" != x"yes" -o x"$has_tls" != x"yes"; then
    AC_DEFINE(TORSION_HAVE_PTHREAD)
    AC_SUBST(LIBS, [$PTHREAD_LIBS])
    pthread_fallback=yes
  fi
fi

if test x"$enable_rng" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_RNG)
fi

if test x"$has_threads" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_THREADS)
fi

if test x"$enable_tls" = x"yes" -a x"$has_tls" = x"yes"; then
  AC_DEFINE(TORSION_HAVE_TLS)
fi

if test x"$enable_verify" = x"yes"; then
  AC_DEFINE(TORSION_VERIFY)
fi

AS_VAR_SET(PKG_VERSION, [AC_PACKAGE_VERSION])

AC_SUBST(PKG_VERSION)
AC_SUBST(LIB_VERSION)

AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)

AC_SUBST(WASI)

AM_CONDITIONAL([EMSCRIPTEN], [test x"$EMSCRIPTEN" != x""])
AM_CONDITIONAL([ENABLE_CTGRIND], [test x"$has_valgrind" = x"yes"])
AM_CONDITIONAL([ENABLE_PTHREAD], [test x"$pthread_fallback" = x"yes"])
AM_CONDITIONAL([ENABLE_RNG], [test x"$enable_rng" = x"yes"])
AM_CONDITIONAL([ENABLE_THREADS], [test x"$has_threads" = x"yes"])
AM_CONDITIONAL([NATIVE], [test x"$WASM" != x"yes"])
AM_CONDITIONAL([WASI], [test x"$WASI" != x""])

AC_CONFIG_FILES([Makefile libtorsion.pc])

AC_OUTPUT

echo ""
echo "Build Options:"
echo ""
echo "  pkg version  = $PKG_VERSION"
echo "  lib version  = $LIB_VERSION"
echo ""
echo "  asm          = $enable_asm (support=$has_asm)"
echo "  asm_x64      = $enable_asm (support=$has_asm_x64)"
echo "  coverage     = $enable_coverage"
echo "  crt          = $enable_crt"
echo "  debug        = $enable_debug"
echo "  emscripten   = $EMSCRIPTEN"
echo "  fork         = $has_fork"
echo "  int128       = $enable_int128 (support=$has_int128)"
echo "  libsecp256k1 = $enable_libsecp256k1"
echo "  rng          = $enable_rng"
echo "  threads      = $has_threads"
echo "  tls          = $enable_tls (support=$has_tls)"
echo "  valgrind     = $has_valgrind"
echo "  verify       = $enable_verify"
echo "  wasi         = $WASI"
echo ""
echo "  CC           = $CC"
echo "  CFLAGS       = $CFLAGS"
echo "  LDFLAGS      = $LDFLAGS"
echo ""
