# ax_check_cpp.m4 - preprocessor checker for autoconf
# Copyright (c) 2021, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion
#
# SYNOPSIS
#
#   AX_CHECK_CPP([action-if-found], [action-if-not-found])
#
# DESCRIPTION
#
#   Like AC_CHECK_DEFINE but only invokes the preprocessor.
#   Useful if AC_PROG_CC hasn't been expanded yet.

AC_DEFUN([AX_CHECK_CPP], [
  AC_CACHE_VAL([ax_cv_check_cpp_cc], [
    AS_IF([test x"$CC" != x],
          [ax_cv_check_cpp_cc="$CC"],
          [ax_cv_check_cpp_cc=cc])
  ])

  AC_CACHE_VAL([ax_cv_check_cpp_cc_e], [
    ax_cv_check_cpp_cc_e=no
    AS_IF([echo foobar1337 | $ax_cv_check_cpp_cc -E - 2> /dev/null | fgrep foobar1337 > /dev/null 2>& 1], [
      AS_IF([! echo __LINE__ | $ax_cv_check_cpp_cc -E - 2> /dev/null | fgrep __LINE__ > /dev/null 2>& 1], [
        ax_cv_check_cpp_cc_e=yes
      ])
    ])
  ])

  ax_check_cpp_has=no

  AS_IF([test x"$ax_cv_check_cpp_cc_e" = x"yes"], [
    AS_IF([! echo $1 | $ax_cv_check_cpp_cc -E - 2> /dev/null | fgrep $1 > /dev/null 2>& 1], [
      ax_check_cpp_has=yes
    ])
  ])

  AS_IF([test x"$ax_check_cpp_has" = x"yes"],
        [m4_default([$2],[[:]])],
        [m4_default([$3],[[:]])])
])
