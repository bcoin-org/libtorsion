#!/bin/sh

set -ex

def='-DTORSION_64BIT'

if test x"$1" = x'32'; then
  def='-DTORSION_32BIT'
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
    -I./include \
    -o test \
    src/drbg.c \
    src/hash.c \
    src/kdf.c \
    src/mini-gmp.c \
    src/poly1305.c \
    src/test.c
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
    -DTORSION_HAS_GMP \
    -I./include \
    -lgmp \
    -o test \
    src/drbg.c \
    src/hash.c \
    src/kdf.c \
    src/poly1305.c \
    src/test.c
fi
