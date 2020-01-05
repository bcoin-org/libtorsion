#!/bin/sh

set -ex

def='-DBCRYPTO_EC_64BIT'

if test x"$1" = x'32'; then
  def='-DBCRYPTO_EC_32BIT'
fi

if test x"$2" = x'mini'; then
  gcc -g \
    -Wall \
    -Wextra \
    -Wno-unused-function \
    -Wno-unused-const-variable \
    -Wno-unused-parameter \
    -std=c89 \
    -O3 \
    "$def" \
    -lcrypto \
    -o test \
    mini-gmp.c \
    test.c
else
  gcc -g \
    -Wall \
    -Wextra \
    -Wno-unused-function \
    -Wno-unused-const-variable \
    -Wno-unused-parameter \
    -std=c89 \
    -O3 \
    "$def" \
    -DBCRYPTO_HAS_GMP \
    -lcrypto \
    -lgmp \
    -o test \
    test.c
fi
