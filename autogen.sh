#!/bin/sh

cd `dirname "$0"`

set -ex

if test x"`uname`" = x"Darwin" && type glibtoolize > /dev/null 2>& 1; then
  glibtoolize --copy --force
else
  libtoolize --copy --force
fi

aclocal --force -I m4
autoconf --force
automake --add-missing --copy --force-missing
