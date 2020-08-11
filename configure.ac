dnl configure.ac - autoconf script for libtorsion
dnl Copyright (c) 2020, Christopher Jeffrey (MIT License).
dnl https://github.com/bcoin-org/libtorsion

dnl
dnl Initialization
dnl

AC_PREREQ(2.64)
AC_INIT([libtorsion], [0.0.0], [https://github.com/bcoin-org/libtorsion/issues])
LIB_VERSION=0:0:0
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([include/torsion/ecc.h])
AM_INIT_AUTOMAKE([foreign subdir-objects no-dist])
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

dnl
dnl Sanity Checks
dnl

AS_IF([test x"$ac_cv_prog_cc_c89" = x"no"], [
  AC_MSG_ERROR([CC does not support ANSI C89])
])

AC_CHECK_DEFINE([__EMSCRIPTEN__], [
  AS_IF([test x"$EMSCRIPTEN" = x""], [
    AC_MSG_ERROR([./configure must be run with emconfigure])
  ])
], [
  AS_IF([test x"$EMSCRIPTEN" != x""], [
    AC_MSG_WARN([EMSCRIPTEN env var set but environment is not emscripten])
    AS_UNSET(EMSCRIPTEN)
  ])
])

AC_CHECK_DEFINE([__wasi__], [
  AS_IF([test x"$WASI" = x""], [
    AC_MSG_ERROR([./configure must be run with ./scripts/wasi-configure])
  ])
], [
  AS_IF([test x"$WASI" != x""], [
    AC_MSG_WARN([WASI env var set but environment is not wasi])
    AS_UNSET(WASI)
  ])
])

AS_IF([test x"$EMSCRIPTEN" != x"" -o x"$WASI" != x""], [WASM=yes], [WASM=no])

dnl
dnl Options
dnl

AC_ARG_ENABLE(
  asm,
  AS_HELP_STRING([--enable-asm],
                 [use inline x86-64 assembly if available [default=yes]]),
  [enable_asm=$enableval],
  [enable_asm=yes]
)

AC_ARG_ENABLE(
  coverage,
  AS_HELP_STRING([--enable-coverage],
                 [enable coverage [default=no]]),
  [enable_coverage=$enableval],
  [enable_coverage=no]
)

AC_ARG_ENABLE(
  crt,
  AS_HELP_STRING([--enable-crt],
                 [enable chinese remainder theorem for RSA [default=no]]),
  [enable_crt=$enableval],
  [enable_crt=no]
)

AC_ARG_ENABLE(
  debug,
  AS_HELP_STRING([--enable-debug],
                 [enable debug build [default=no]]),
  [enable_debug=$enableval],
  [enable_debug=no]
)

AC_ARG_ENABLE(
  int128,
  AS_HELP_STRING([--enable-int128],
                 [use __int128 if available [default=yes]]),
  [enable_int128=$enableval],
  [enable_int128=yes]
)

AC_ARG_ENABLE(
  libsecp256k1,
  AS_HELP_STRING([--enable-libsecp256k1],
                 [use libsecp256k1 field element backend [default=no]]),
  [enable_libsecp256k1=$enableval],
  [enable_libsecp256k1=no]
)

AC_ARG_ENABLE(
  pthread,
  AS_HELP_STRING([--enable-pthread],
                 [use pthread as a fallback for tls [default=yes]]),
  [enable_pthread=$enableval],
  [enable_pthread=yes]
)

AC_ARG_ENABLE(
  rng,
  AS_HELP_STRING([--enable-rng],
                 [enable rng [default=yes]]),
  [enable_rng=$enableval],
  [enable_rng=yes]
)

AC_ARG_ENABLE(
  tests,
  AS_HELP_STRING([--enable-tests],
                 [enable tests [default=yes]]),
  [enable_tests=$enableval],
  [enable_tests=yes]
)

AC_ARG_ENABLE(
  tls,
  AS_HELP_STRING([--enable-tls],
                 [enable tls [default=yes]]),
  [enable_tls=$enableval],
  [enable_tls=yes]
)

AC_ARG_ENABLE(
  verify,
  AS_HELP_STRING([--enable-verify],
                 [enable scalar bounds checks [default=no]]),
  [enable_verify=$enableval],
  [enable_verify=no]
)

dnl
dnl Flags
dnl

AX_CHECK_COMPILE_FLAG([-std=c89], [CFLAGS="$CFLAGS -std=c89"])
AX_CHECK_COMPILE_FLAG([-pedantic], [CFLAGS="$CFLAGS -pedantic"])
AX_CHECK_COMPILE_FLAG([-Wall], [CFLAGS="$CFLAGS -Wall"])
AX_CHECK_COMPILE_FLAG([-Wextra], [CFLAGS="$CFLAGS -Wextra"])
AX_CHECK_COMPILE_FLAG([-Wcast-align], [CFLAGS="$CFLAGS -Wcast-align"])
AX_CHECK_COMPILE_FLAG([-Wno-implicit-fallthrough],
                      [CFLAGS="$CFLAGS -Wno-implicit-fallthrough"])
AX_CHECK_COMPILE_FLAG([-Wno-long-long], [CFLAGS="$CFLAGS -Wno-long-long"])
AX_CHECK_COMPILE_FLAG([-Wno-overlength-strings],
                      [CFLAGS="$CFLAGS -Wno-overlength-strings"])
AX_CHECK_COMPILE_FLAG([-Wshadow], [CFLAGS="$CFLAGS -Wshadow"])
AX_CHECK_COMPILE_FLAG([-O3], [CFLAGS="$CFLAGS -O3"])

AS_IF([test x"$WASM" != x"yes"], [
  AX_CHECK_COMPILE_FLAG([-fvisibility=hidden],
                        [CFLAGS="$CFLAGS -fvisibility=hidden"])

  AS_CASE([$host_os], [darwin*], [
    AX_CHECK_COMPILE_FLAG([-mmacosx-version-min=10.7],
                          [CFLAGS="$CFLAGS -mmacosx-version-min=10.7"])
  ])

  AS_CASE([$host_os], [linux*], [CFLAGS="$CFLAGS -D_POSIX_C_SOURCE=200112"])
])

dnl
dnl Feature Testing
dnl

AC_MSG_CHECKING(for inline asm support)
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[]], [[
    unsigned char ptr[32];
    int x = 0;
    __asm__ ("" : "+r" (x) ::);
    __asm__ __volatile__("" :: "r" (ptr) : "memory");
  ]])
], [
  has_asm=yes
], [
  has_asm=no
])
AC_MSG_RESULT([$has_asm])

AC_MSG_CHECKING(for x86_64 inline asm support)
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[
#if !defined(__amd64__) && !defined(__x86_64__)
#  error "not an x64 platform"
#endif
#include <stdint.h>
  ]], [[
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
  ]])
], [
  has_asm_x64=yes
], [
  has_asm_x64=no
])
AC_MSG_RESULT([$has_asm_x64])

AC_MSG_CHECKING(for int128 support)
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[
#ifndef __wasm__
    typedef char check_64bit_t[sizeof(void *) == 8 ? 1 : -1];
#endif
    typedef signed __int128 int128_t;
    typedef unsigned __int128 uint128_t;
  ]], [[
    uint128_t r = 1;
    r <<= 64;
    r *= 113;
    r >>= 65;
    r &= 1;
  ]])
], [
  has_int128=yes
], [
  has_int128=no
])
AC_MSG_RESULT([$has_int128])

AC_MSG_CHECKING(for tls support)
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[
#if defined(_WIN32) && !defined(__MINGW32__)
#  define TLS __declspec(thread)
#else
#  define TLS __thread
#endif
    static TLS int value;
  ]], [[
    value = 1;
  ]])
], [
  has_tls=yes
], [
  has_tls=no
])
AC_MSG_RESULT([$has_tls])

has_threads=no

AS_IF([test x"$WASM" != x"yes"], [
  AX_PTHREAD([has_threads=yes])
])

AC_MSG_CHECKING(for fork support)
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[
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
  ]])
], [
  has_fork=yes
], [
  has_fork=no
])
AC_MSG_RESULT([$has_fork])

has_valgrind=no
has_zlib=no

AS_IF([test x"$WASM" != x"yes"], [
  AC_CHECK_HEADER([valgrind/memcheck.h], [has_valgrind=yes])

  AS_IF([test x"$enable_rng" = x"yes"], [
    AC_CHECK_HEADER([zlib.h], [
      AC_CHECK_LIB(z, compress2, [has_zlib=yes])
    ])
  ])
])

dnl
dnl Defines
dnl

AC_DEFINE(TORSION_HAVE_CONFIG)

AS_IF([test x"$enable_asm" = x"yes" -a x"$has_asm" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_ASM)
])

AS_IF([test x"$enable_asm" = x"yes" -a x"$has_asm_x64" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_ASM_X64)
])

AS_IF([test x"$enable_coverage" = x"yes"], [
  AS_IF([test x"$WASM" = x"yes"], [
    AC_MSG_ERROR([coverage not supported with wasm])
  ])
  CFLAGS="$CFLAGS -O0 --coverage"
  AC_DEFINE(TORSION_COVERAGE)
])

AS_IF([test x"$enable_crt" = x"yes"], [
  AC_DEFINE(TORSION_USE_CRT)
])

AS_IF([test x"$enable_debug" = x"yes"], [
  AX_CHECK_COMPILE_FLAG([-g], [CFLAGS="$CFLAGS -g"])
  AC_DEFINE(TORSION_DEBUG)
])

AS_IF([test x"$has_fork" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_FORK)
])

AS_IF([test x"$enable_int128" = x"yes" -a x"$has_int128" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_INT128)
])

AS_IF([test x"$enable_libsecp256k1" = x"yes"], [
  AC_DEFINE(TORSION_USE_LIBSECP256K1)
])

pthread_fallback=no

AS_IF([test x"$enable_pthread" = x"yes" -a x"$has_threads" = x"yes"], [
  AS_IF([test x"$enable_tls" != x"yes" -o x"$has_tls" != x"yes"], [
    AC_DEFINE(TORSION_HAVE_PTHREAD)
    AC_SUBST(LIBS, [$PTHREAD_LIBS])
    pthread_fallback=yes
  ])
])

AS_IF([test x"$enable_rng" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_RNG)
])

AS_IF([test x"$has_threads" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_THREADS)
])

AS_IF([test x"$enable_tls" = x"yes" -a x"$has_tls" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_TLS)
])

AS_IF([test x"$enable_verify" = x"yes"], [
  AC_DEFINE(TORSION_VERIFY)
])

AS_IF([test x"$has_zlib" = x"yes"], [
  AC_DEFINE(TORSION_HAVE_ZLIB)
])

dnl
dnl Configuration
dnl

AC_SUBST(LIB_VERSION)

AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)

AC_SUBST(WASI)

AM_CONDITIONAL([EMSCRIPTEN], [test x"$EMSCRIPTEN" != x""])
AM_CONDITIONAL([ENABLE_CTGRIND], [test x"$has_valgrind" = x"yes"])
AM_CONDITIONAL([ENABLE_PTHREAD], [test x"$pthread_fallback" = x"yes"])
AM_CONDITIONAL([ENABLE_RNG], [test x"$enable_rng" = x"yes"])
AM_CONDITIONAL([ENABLE_TESTS], [test x"$enable_tests" = x"yes"])
AM_CONDITIONAL([ENABLE_THREADS], [test x"$has_threads" = x"yes"])
AM_CONDITIONAL([ENABLE_ZLIB], [test x"$has_zlib" = x"yes"])
AM_CONDITIONAL([NATIVE], [test x"$WASM" != x"yes"])
AM_CONDITIONAL([WASI], [test x"$WASI" != x""])

AC_CONFIG_FILES([Makefile libtorsion.pc])

dnl
dnl Output
dnl

AC_OUTPUT

AC_MSG_NOTICE([Build Options:

  asm          = $enable_asm (support=$has_asm)
  asm_x64      = $enable_asm (support=$has_asm_x64)
  coverage     = $enable_coverage
  crt          = $enable_crt
  debug        = $enable_debug
  emscripten   = $EMSCRIPTEN
  fork         = $has_fork
  int128       = $enable_int128 (support=$has_int128)
  libsecp256k1 = $enable_libsecp256k1
  rng          = $enable_rng
  threads      = $has_threads
  tests        = $enable_tests
  tls          = $enable_tls (support=$has_tls)
  valgrind     = $has_valgrind
  verify       = $enable_verify
  wasi         = $WASI
  zlib         = $has_zlib

  CC           = $CC
  CPPFLAGS     = $CPPFLAGS
  CFLAGS       = $CFLAGS
  LDFLAGS      = $LDFLAGS
])
