#!/bin/sh

# autogen.sh - autotools generation script
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

set -ex

cd `dirname "$0"`

if type glibtoolize > /dev/null 2>& 1; then
  glibtoolize --copy --force
else
  libtoolize --copy --force
fi

aclocal --force -I m4
autoconf --force
automake --add-missing --copy --force-missing

# Hack to get dietlibc working with autotools.
sed -e 's/| uclibc\*)$/| uclibc* | dietlibc*)/' \
  < build-aux/config.sub > build-aux/config.sub.tmp
mv -f build-aux/config.sub.tmp build-aux/config.sub
chmod 0755 build-aux/config.sub
