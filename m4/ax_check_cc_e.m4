# ax_check_cc_e.m4 - preprocessor checker for autoconf
# Copyright (c) 2021, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion
#
# SYNOPSIS
#
#   AX_CHECK_CC_E([action-if-found], [action-if-not-found], [action-if-fail])
#
# DESCRIPTION
#
#   Like AC_CHECK_DEFINE but only invokes the preprocessor.
#   Useful if AC_PROG_CC hasn't been expanded yet.

AC_DEFUN([AX_CHECK_CC_E], [
  AC_CACHE_CHECK([for ${CC-cc} -E -], [ax_cv_cce_works], [
    ax_cv_cce_works=no
    AS_IF([echo helloworld | ${CC-cc} -E - 2> /dev/null | grep helloworld > /dev/null], [
      AS_IF([! echo __LINE__ | ${CC-cc} -E - 2> /dev/null | grep __LINE__ > /dev/null],
            [ax_cv_cce_works=yes])
    ])
  ])
  AS_IF([test x"$ax_cv_cce_works" = x'yes'], [
    AC_CACHE_CHECK([for $1 defined], [ac_cv_defined_$1], [
      AS_IF([! echo $1 | ${CC-cc} -E - 2> /dev/null | grep $1 > /dev/null],
            [ac_cv_defined_$1=yes], [ac_cv_defined_$1=no])
    ])
    AS_IF([test x"$ac_cv_defined_$1" = x'yes'], [$2], [$3])
  ], [$4])
])
