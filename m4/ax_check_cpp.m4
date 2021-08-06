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
  AC_CACHE_VAL([ax_cv_check_cpp_cc_e], [
    ax_cv_check_cpp_cc_e=no
    AS_IF([echo foobar1337 | ${CC-cc} -E - 2> /dev/null | grep foobar1337 > /dev/null], [
      AS_IF([! echo __LINE__ | ${CC-cc} -E - 2> /dev/null | grep __LINE__ > /dev/null], [
        ax_cv_check_cpp_cc_e=yes
      ])
    ])
  ])

  ax_check_cpp_has=no

  AS_IF([test x"$ax_cv_check_cpp_cc_e" = x"yes"], [
    AS_IF([! echo $1 | ${CC-cc} -E - 2> /dev/null | grep $1 > /dev/null], [
      ax_check_cpp_has=yes
    ])
  ])

  AS_IF([test x"$ax_check_cpp_has" = x"yes"],
        [m4_default([$2],[[:]])],
        [m4_default([$3],[[:]])])
])
