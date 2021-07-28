# ax_tls.m4 - thread-local storage test for autoconf
# Copyright (c) 2021, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion
#
# SYNOPSIS
#
#   AX_TLS([action-if-found], [action-if-not-found])
#
# DESCRIPTION
#
#   TODO

AC_DEFUN([AX_TLS_RUN_IFELSE], [
  AS_IF([test x"$cross_compiling" != x"no"],
        [AC_LINK_IFELSE([$1], [m4_default([$2],[[:]])],
                              [m4_default([$3],[[:]])])],
        [AC_RUN_IFELSE([$1], [m4_default([$2],[[:]])],
                             [m4_default([$3],[[:]])],
                             [:])])
])

AC_DEFUN([AX_TLS_CHECK_CFLAG], [
  ax_tlscf_has_flag=no
  ax_tlscf_save_CFLAGS="$CFLAGS"
  CFLAGS="$CFLAGS $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [ax_tlscf_has_flag=yes])
  CFLAGS="$ax_tlscf_save_CFLAGS"
  AS_IF([test x"$ax_tlscf_has_flag" = x"yes"],
        [m4_default([$2],[[:]])],
        [m4_default([$3],[[:]])])
])

AC_DEFUN([AX_TLS], [
  AC_MSG_CHECKING([for thread-local storage support])

  ax_tls_check=no

  AC_CACHE_VAL([ax_cv_tls_keyword], [ax_tls_check=yes])
  AC_CACHE_VAL([ax_cv_tls_cflags], [ax_tls_check=yes])

  AS_IF([test x"$ax_tls_check" = x"yes"], [
    ax_tls_save_CFLAGS="$CFLAGS"
    ax_cv_tls_keyword=none
    ax_cv_tls_cflags=none
    ax_tls_cflag=

    # Try to deoptimize.
    AX_TLS_CHECK_CFLAG([-O0], [CFLAGS="$CFLAGS -O0"])

    # XL requires a special flag. Don't ask me why.
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([], [[
#       if !defined(__xlC__) && !defined(__ibmxl__)
#         error "not an ibm compiler"
#       endif
      ]])
    ], [
      AX_TLS_CHECK_CFLAG([-qtls], [
        ax_tls_cflag="-qtls"
        CFLAGS="$CFLAGS -qtls"
      ])
    ])

    # Various TLS keywords. We prepend or append
    # _Thread_local depending on the C standard.
    # The last keyword, __declspec(__thread), is
    # not widely known, but there is evidence
    # that Compaq C for Tru64 UNIX supported it
    # at one point.
    ax_tls_keywords="__thread __declspec(thread) __declspec(__thread)"

    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([], [[
#       ifndef __cplusplus
#         error "not c++"
#       endif
      ]])
    ], [
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [[
#         if !defined(__cplusplus) || (__cplusplus + 0L) < 201103L
#           error "not c++11"
#         endif
        ]])
      ], [
        ax_tls_keywords="thread_local $ax_tls_keywords"
      ], [
        ax_tls_keywords="$ax_tls_keywords thread_local"
      ])
    ], [
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [[
#         if !defined(__STDC_VERSION__) || (__STDC_VERSION__ + 0L) < 201112L
#           error "not c11"
#         endif
        ]])
      ], [
        ax_tls_keywords="_Thread_local $ax_tls_keywords"
      ], [
        ax_tls_keywords="$ax_tls_keywords _Thread_local"
      ])
    ])

    for ax_tls_keyword in $ax_tls_keywords; do
      AX_TLS_RUN_IFELSE([
        AC_LANG_SOURCE([[
          static $ax_tls_keyword int x;
          int main(void) {
            x = 1;
            return !x;
          }
        ]])
      ], [
        ax_cv_tls_keyword="$ax_tls_keyword"
        ax_cv_tls_cflags="$ax_tls_cflag"
        break
      ])
    done

    CFLAGS="$ax_tls_save_CFLAGS"
  ])

  AS_IF([test x"$ax_cv_tls_keyword" != x"none"], [
    TLS_KEYWORD="$ax_cv_tls_keyword"
    TLS_CFLAGS="$ax_cv_tls_cflags"
  ], [
    TLS_KEYWORD=""
    TLS_CFLAGS=""
  ])

  AC_SUBST([TLS_KEYWORD])
  AC_SUBST([TLS_CFLAGS])

  AC_MSG_RESULT([$ax_cv_tls_keyword])

  AS_IF([test x"$ax_cv_tls_keyword" != x"none"],
        [m4_default([$1],[[:]])],
        [m4_default([$2],[[:]])])
])
