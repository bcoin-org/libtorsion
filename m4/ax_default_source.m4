# ===========================================================================
# https://github.com/bcoin-org/libtorsion/blob/master/m4/ax_default_source.m4
# ===========================================================================
#
# SYNOPSIS
#
#   AX_DEFAULT_SOURCE
#
# DESCRIPTION
#
#   This macro is a better version of AC_USE_SYSTEM_EXTENSIONS (with slightly
#   different semantics). Rather than attempting to "enable extensions" it
#   tries to define macros which expose the _default_ source for the system
#   regardless of the language standard used. This involes seeking out C
#   standard libraries which determine the set of POSIX features by language
#   standard. In effect, we divorce the set of POSIX features from the set of
#   language features (this forces more BSD-like behavior).
#
# SEE ALSO
#
#   https://github.com/bcoin-org/libtorsion/blob/62694bf/src/entropy/posix.h
#   https://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Posix-Variants.html
#   https://github.com/autotools-mirror/autoconf/blob/6d38e9f/lib/autoconf/specific.m4#L357
#
# LICENSE
#
#   Copyright (c) 2020, Christopher Jeffrey (MIT License).
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

AC_DEFUN([AX_DEFAULT_SOURCE], [
  AC_REQUIRE([AC_CANONICAL_HOST])

  AS_CASE([$host_os], [linux*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [gnu*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [kfreebsd*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [knetbsd*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [solaris*], [
    AC_DEFINE([__EXTENSIONS__])
    AC_DEFINE([_XOPEN_SOURCE], [500])
  ])

  AS_CASE([$host_os], [cygwin*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [msys*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [hpux*], [
    AC_DEFINE([_HPUX_SOURCE])
    AC_DEFINE([_XOPEN_SOURCE], [500])
  ])

  AS_CASE([$host_os], [nonstopux*], [
    AC_DEFINE([_TANDEM_SOURCE])
  ])

  AS_CASE([$host_os], [aix*], [
    AC_DEFINE([_ALL_SOURCE])
  ])

  AS_CASE([$host_os], [openedition*], [
    AC_DEFINE([_ALL_SOURCE])
    AC_DEFINE([_UNIX03_SOURCE])
    AC_DEFINE([_UNIX03_THREADS])
  ])

  AS_CASE([$host_os], [qnx* | nto-qnx*], [
    AC_DEFINE([_QNX_SOURCE])
  ])

  AS_CASE([$host_os], [haiku*], [
    AC_DEFINE([_BSD_SOURCE])
  ])

  AS_CASE([$host_os], [wasi*], [
    AC_DEFINE([_GNU_SOURCE])
  ])

  AS_CASE([$host_os], [emscripten*], [
    AC_DEFINE([_GNU_SOURCE])
  ])
])dnl AX_DEFAULT_SOURCE
