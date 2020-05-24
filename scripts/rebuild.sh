#!/bin/sh

set -ex

make clean || echo -n ''
./autogen.sh
./configure --enable-tests --enable-rng CFLAGS="-g -O3"
make
