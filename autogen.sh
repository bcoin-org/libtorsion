#!/bin/sh

# autogen.sh - autotools generation script
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

set -ex

cd `dirname "$0"`

if test x"`uname`" = x"Darwin" && type glibtoolize > /dev/null 2>& 1; then
  glibtoolize --copy --force
else
  libtoolize --copy --force
fi

aclocal --force -I m4
autoconf --force
automake --add-missing --copy --force-missing
