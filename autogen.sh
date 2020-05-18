#!/bin/sh

cd `dirname "$0"`

set -ex

if test x"`uname`" = x"Darwin" && type glibtoolize > /dev/null 2>& 1; then
  glibtoolize --copy
else
  libtoolize --copy
fi

aclocal -I m4
autoconf
automake --add-missing --copy
