dnl configure.ac - autoconf script for libtorsion
dnl Copyright (c) 2020, Christopher Jeffrey (MIT License).
dnl https://github.com/bcoin-org/libtorsion

m4_define(TORSION_PKG_VERSION, 0.0.0)
m4_define(TORSION_ABI_VERSION, 0:0:0)

#
# Initialization
#

AC_PREREQ([2.64])
AC_INIT([libtorsion], [m4_substr(TORSION_PKG_VERSION, 0)],
        [https://github.com/bcoin-org/libtorsion/issues])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([include/torsion/ecc.h])
AC_CANONICAL_HOST
AC_LANG_ASSERT([C])
AX_CHECK_CC_E([__PCC__], [cross_compiling=yes])
AC_PROG_CC

AM_INIT_AUTOMAKE([1.11 foreign subdir-objects tar-v7 dist-zip])
AM_PROG_CC_C_O
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AM_SILENT_RULES([yes])

LT_PREREQ([2.2.8])
LT_INIT([shared static])

#
# Sanity Checks
#

AC_CHECK_DEFINE([__wasi__], [
  AS_IF([test x"$WASI" = x''], [
    AC_MSG_ERROR([./configure must be run with ./scripts/wasi-configure])
  ])
], [
  AS_IF([test x"$WASI" != x''], [
    AC_MSG_WARN([WASI env var set but environment is not wasi])
    AS_UNSET([WASI])
  ])
])

AC_CHECK_DEFINE([__EMSCRIPTEN__], [
  AS_IF([test x"$EMSCRIPTEN" = x''], [
    AC_MSG_ERROR([./configure must be run with emconfigure])
  ])
], [
  AS_IF([test x"$EMSCRIPTEN" != x''], [
    AC_MSG_WARN([EMSCRIPTEN env var set but environment is not emscripten])
    AS_UNSET([EMSCRIPTEN])
  ])
])

#
# Options
#

AC_ARG_ENABLE(
  asm,
  AS_HELP_STRING([--enable-asm],
                 [use inline assembly if available [default=yes]]),
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
  mpi,
  AS_HELP_STRING([--enable-mpi],
                 [export mpi functions [default=no]]),
  [enable_mpi=$enableval],
  [enable_mpi=no]
)

AC_ARG_ENABLE(
  pthread,
  AS_HELP_STRING([--enable-pthread],
                 [use pthread if present in libc [default=yes]]),
  [enable_pthread=$enableval],
  [enable_pthread=yes]
)

AC_ARG_ENABLE(
  rng,
  AS_HELP_STRING([--enable-rng],
                 [enable RNG [default=yes]]),
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
                 [use thread-local storage if available [default=yes]]),
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

#
# Global Flags
#

AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[]], [[
#   if !defined(__GNUC__) || __GNUC__ < 4
      choke me
#   endif
  ]])
], [
  AX_CHECK_COMPILE_FLAG([-fvisibility=hidden],
                        [CFLAGS="$CFLAGS -fvisibility=hidden"])
])

AS_CASE([$host_os], [darwin*], [
  AX_CHECK_COMPILE_FLAG([-mmacosx-version-min=10.7],
                        [CFLAGS="$CFLAGS -mmacosx-version-min=10.7"])
])

AS_CASE([$host_os], [mingw*], [
  LIBS="-lkernel32 -ladvapi32 $LIBS"
])

#
# Compiler Fixes
#

# Encourage the user to build with xlc_r.
AC_COMPILE_IFELSE([
  AC_LANG_PROGRAM([[]], [[
#   if defined(__xlC__) && defined(_AIX) && !defined(_THREAD_SAFE)
      choke me
#   endif
  ]])
], [], [
  AC_MSG_WARN([Please use a thread-safe compiler invocation.])
  AC_MSG_WARN([For example, './configure CC=xlc_r'.])
])

# Modern glibc (stdlib.h in particular) breaks pcc.
AC_CHECK_DEFINE([__PCC__], [
  AX_CHECK_DEFINE([limits.h], [__GLIBC__], [
    AC_CHECK_TYPE([__float128], [], [
      AC_DEFINE_UNQUOTED([__float128], [long double])
      AC_DEFINE_UNQUOTED([__TC__], [__SC__])
    ])
  ])
])

# NWCC has issues with -fPIC and shared libraries.
AC_CHECK_DEFINE([__NWCC__], [
  pic_mode=no
  enable_shared=no
])

# chibicc doesn't support .so versioning.
AC_CHECK_DEFINE([__chibicc__], [
  enable_shared=no
])

# Hack to get OpenWatcom's owcc working.
AC_CHECK_DEFINE([__WATCOMC__], [
  grep -v '#define PACKAGE_STRING' confdefs.h > confdefs.h.tmp
  mv -f confdefs.h.tmp confdefs.h
])

# OpenWatcom-Linux can't create shared libraries.
AC_CHECK_DEFINE([__WATCOMC__], [
  AC_CHECK_DEFINE([__LINUX__], [
    enable_shared=no
  ])
])

# dietlibc is for static linking only.
AC_CHECK_DEFINE([__dietlibc__], [
  enable_shared=no
])

#
# Feature Testing
#

has_asm=no
has_clock_gettime=no
has_emutls=no
has_fork=no
has_gettimeofday=no
has_int128=no
has_pthread=no
has_time=no
has_tls=no
has_memcheck=no
has_zlib=no

AS_IF([test x"$enable_asm" = x'yes'], [
  AC_MSG_CHECKING(for inline asm support)
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
      int main(void) {
        unsigned long z = 953;
        unsigned long x = 109;
        unsigned long y = 577;
        unsigned long c;
        __asm__ __volatile__ (
#     if defined(__amd64__) || defined(__amd64)
          "movq \@S|@0, %q1\n"
#     elif defined(__x86_64__) || defined(__x86_64)
          "movq \@S|@0, %q1\n"
#     elif defined(__i386__) || defined(__i386) || defined(i386)
          "movl \@S|@0, %k1\n"
#     elif defined(__arm__) || defined(__aarch64__)
          "mov %1, #0\n"
#     elif defined(__powerpc__) || defined(__powerpc64__) || defined(__PPC__)
          "li %1, 0\n"
#     elif defined(__riscv)
          "lui %1, 0\n"
#     else
          ""
#     endif
          : "+r" (z), "=&r" (c)
#     if defined(__wasm__) || defined(__EMSCRIPTEN__)
          :
#     elif defined(__TINYC__)
          : "rm" (x), "rm" (y)
#     else
          : "%rm" (x), "rm" (y)
#     endif
          : "cc", "memory"
        );
        return z & 0x7f;
      }
    ]])
  ], [
    has_asm=yes
  ])
  AC_MSG_RESULT([$has_asm])
])

AS_IF([test x"$enable_tests" = x'yes'], [
  AS_CASE([$host_os], [darwin*], [], [
    AC_MSG_CHECKING(for clock_gettime support)
    AC_LINK_IFELSE([
      AC_LANG_SOURCE([[
#       include <time.h>
        int main(void) {
          struct timespec ts;
          (void)clock_gettime(CLOCK_REALTIME, &ts);
          (void)clock_gettime(CLOCK_MONOTONIC, &ts);
          return ts.tv_sec & ts.tv_nsec & 0x7f;
        }
      ]])
    ], [
      has_clock_gettime=yes
    ])
    AC_MSG_RESULT([$has_clock_gettime])
  ])
])

AS_IF([test x"$enable_tests" = x'yes' -a x"$EMSCRIPTEN" = x''], [
  AC_MSG_CHECKING(for fork support)
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
#     include <sys/types.h>
#     include <sys/wait.h>
#     include <unistd.h>
      int main(void) {
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
        return status;
      }
    ]])
  ], [
    has_fork=yes
  ])
  AC_MSG_RESULT([$has_fork])
])

AS_IF([test x"$enable_tests" = x'yes'], [
  AC_MSG_CHECKING(for gettimeofday support)
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
#     include <stddef.h>
#     include <sys/time.h>
      int main(void) {
        struct timeval tv;
        (void)gettimeofday(&tv, NULL);
        return tv.tv_sec & tv.tv_usec & 0x7f;
      }
    ]])
  ], [
    has_gettimeofday=yes
  ])
  AC_MSG_RESULT([$has_gettimeofday])
])

AS_IF([test x"$enable_int128" = x'yes'], [
  AC_MSG_CHECKING(for int128 support)
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
      typedef signed __int128 xint128_t;
      typedef unsigned __int128 xuint128_t;
      typedef char check_voidptr_t[sizeof(void *) >= 8 ? 1 : -1];
      typedef char check_int128_t[sizeof(xint128_t) == 16 ? 1 : -1];
      typedef char check_uint128_t[sizeof(xuint128_t) == 16 ? 1 : -1];
      int main(int argc, char **argv) {
        xint128_t c = argv[0][0];
        xuint128_t r = argc + c;
        while (argc--) r *= r;
        return r >> 121;
      }
    ]])
  ], [
    has_int128=yes
  ])
  AC_MSG_RESULT([$has_int128])
])

AS_IF([test x"$enable_pthread" = x'yes' -a x"$EMSCRIPTEN" = x''], [
  AC_MSG_CHECKING(for pthread support)
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
#     include <pthread.h>
      static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
      int main(void) {
        (void)pthread_mutex_lock(&lock);
        (void)pthread_mutex_unlock(&lock);
        return 0;
      }
    ]])
  ], [
    has_pthread=yes
  ])
  AC_MSG_RESULT([$has_pthread])
])

AS_IF([test x"$enable_tests" = x'yes'], [
  AC_MSG_CHECKING(for time support)
  AC_LINK_IFELSE([
    AC_LANG_SOURCE([[
#     include <stddef.h>
#     include <time.h>
      int main(void) {
        return time(NULL) & 0x7f;
      }
    ]])
  ], [
    has_time=yes
  ])
  AC_MSG_RESULT([$has_time])
])

AS_IF([test x"$enable_tls" = x'yes'], [
  AX_TLS([has_tls=yes], [has_tls=no])
  AS_IF([test x"$ax_cv_tls_emulated$has_pthread" = x'yesyes'], [has_tls=no])
])

AS_IF([test x"$enable_tests" = x'yes'], [
  AC_CHECK_HEADER([valgrind/memcheck.h], [has_memcheck=yes])

  AS_IF([test x"$enable_rng" = x'yes'], [
    AC_CHECK_HEADER([zlib.h], [
      AC_CHECK_LIB([z], [compress2], [has_zlib=yes])
    ])
  ])
])

#
# Flags
#

AX_CHECK_COMPILE_FLAG([-Werror=unknown-warning-option],
                      [CFLAGS="$CFLAGS -Werror=unknown-warning-option"])
AX_CHECK_COMPILE_FLAG([-pedantic], [CFLAGS="$CFLAGS -pedantic"])
AX_CHECK_COMPILE_FLAG([-Wall], [CFLAGS="$CFLAGS -Wall"])
AX_CHECK_COMPILE_FLAG([-Wextra], [CFLAGS="$CFLAGS -Wextra"])
AX_CHECK_COMPILE_FLAG([-Wcast-align], [CFLAGS="$CFLAGS -Wcast-align"])
AX_CHECK_COMPILE_FLAG([-Wcast-align=strict],
                      [CFLAGS="$CFLAGS -Wcast-align=strict"])
AX_CHECK_COMPILE_FLAG([-Wconditional-uninitialized],
                      [CFLAGS="$CFLAGS -Wconditional-uninitialized"])
AX_CHECK_COMPILE_FLAG([-Wmissing-prototypes],
                      [CFLAGS="$CFLAGS -Wmissing-prototypes"])
AX_CHECK_COMPILE_FLAG([-Wno-implicit-fallthrough],
                      [CFLAGS="$CFLAGS -Wno-implicit-fallthrough"])
AX_CHECK_COMPILE_FLAG([-Wno-long-long], [CFLAGS="$CFLAGS -Wno-long-long"])
AX_CHECK_COMPILE_FLAG([-Wno-overlength-strings],
                      [CFLAGS="$CFLAGS -Wno-overlength-strings"])
AX_CHECK_COMPILE_FLAG([-Wshadow], [CFLAGS="$CFLAGS -Wshadow"])
AX_CHECK_COMPILE_FLAG([-Wstrict-prototypes],
                      [CFLAGS="$CFLAGS -Wstrict-prototypes"])
AX_CHECK_COMPILE_FLAG([-Wundef], [CFLAGS="$CFLAGS -Wundef"])

AC_CHECK_DEFINE([__WATCOMC__], [
  AX_CHECK_COMPILE_FLAG([-Wno-n201], [CFLAGS="$CFLAGS -Wno-n201"])
  AX_CHECK_COMPILE_FLAG([-Wno-n202], [CFLAGS="$CFLAGS -Wno-n202"])
])

AS_IF([test x"$enable_coverage" = x'yes'], [
  CFLAGS="$CFLAGS -O0 --coverage"
])

AS_IF([test x"$enable_debug" = x'yes'], [
  AX_CHECK_COMPILE_FLAG([-g], [CFLAGS="$CFLAGS -g"])
])

AS_IF([test x"$has_tls" = x'yes'], [
  CFLAGS="$CFLAGS $ax_cv_tls_cflags"
])

AS_CASE([$host_os], [mingw*], [
  # Ensure we are redistributable on windows.
  LDFLAGS="$LDFLAGS -static-libgcc"
])

#
# Defines
#

AC_DEFINE([TORSION_HAVE_CONFIG])

AS_IF([test x"$has_asm" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_ASM])
])

AS_IF([test x"$has_clock_gettime" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_CLOCK_GETTIME])
])

AS_IF([test x"$enable_coverage" = x'yes'], [
  AC_DEFINE([TORSION_COVERAGE])
])

AS_IF([test x"$enable_debug" = x'yes'], [
  AC_DEFINE([TORSION_DEBUG])
])

AS_IF([test x"$has_fork" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_FORK])
])

AS_IF([test x"$has_gettimeofday" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_GETTIMEOFDAY])
])

AS_IF([test x"$has_int128" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_INT128])
])

AS_IF([test x"$enable_mpi" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_MPI])
])

AS_IF([test x"$has_pthread" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_PTHREAD])
])

AS_IF([test x"$enable_rng" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_RNG])
])

AS_IF([test x"$has_time" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_TIME])
])

AS_IF([test x"$has_tls" = x'yes'], [
  AC_DEFINE_UNQUOTED([TORSION_TLS], [$ax_cv_tls_keyword])
])

AS_IF([test x"$enable_verify" = x'yes'], [
  AC_DEFINE([TORSION_VERIFY])
])

AS_IF([test x"$has_zlib" = x'yes'], [
  AC_DEFINE([TORSION_HAVE_ZLIB])
])

#
# Feature Test Macros
#

AS_CASE([$host_os], [solaris*], [
  AC_DEFINE([_TS_ERRNO])
])

AS_CASE([$host_os], [aix*], [
  AC_DEFINE([_THREAD_SAFE_ERRNO])
])

#
# Configuration
#

AC_SUBST([ABI_VERSION], [TORSION_ABI_VERSION])

AM_CONDITIONAL([EMSCRIPTEN], [test x"$EMSCRIPTEN" != x''])
AM_CONDITIONAL([ENABLE_CTGRIND], [test x"$has_memcheck" = x'yes'])
AM_CONDITIONAL([ENABLE_RNG], [test x"$enable_rng" = x'yes'])
AM_CONDITIONAL([ENABLE_SHARED], [test x"$enable_shared" = x'yes'])
AM_CONDITIONAL([ENABLE_TESTS], [test x"$enable_tests" = x'yes'])
AM_CONDITIONAL([ENABLE_ZLIB], [test x"$has_zlib" = x'yes'])
AM_CONDITIONAL([MINGW], [AS_CASE([$host_os], [mingw*], [true], [false])])
AM_CONDITIONAL([NATIVE], [test x"$WASI$EMSCRIPTEN" = x''])
AM_CONDITIONAL([WASI], [test x"$WASI" != x''])

AC_CONFIG_FILES([Makefile libtorsion.pc])

#
# Output
#

AC_OUTPUT

AC_MSG_NOTICE([Build Options:

  asm           = $has_asm
  clock_gettime = $has_clock_gettime
  coverage      = $enable_coverage
  debug         = $enable_debug
  emscripten    = $EMSCRIPTEN
  fork          = $has_fork
  gettimeofday  = $has_gettimeofday
  int128        = $has_int128
  memcheck      = $has_memcheck
  mpi           = $enable_mpi
  pthread       = $has_pthread
  rng           = $enable_rng
  tests         = $enable_tests
  time          = $has_time
  tls           = $has_tls
  verify        = $enable_verify
  wasi          = $WASI
  zlib          = $has_zlib

  PREFIX        = $prefix
  HOST          = $host
  CC            = $CC
  CPPFLAGS      = $CPPFLAGS
  CFLAGS        = $CFLAGS
  LDFLAGS       = $LDFLAGS
  DEFS          = $DEFS
  LIBS          = $LIBS
])
