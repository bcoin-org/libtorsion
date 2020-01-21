#!/bin/sh

set -ex

def='-DTORSION_USE_64BIT'

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
    -Wno-implicit-fallthrough \
    -std=c89 \
    -O3 \
    "$def" \
    -I./include \
    -o test \
    src/aead.c \
    src/chacha20.c \
    src/drbg.c \
    src/dsa.c \
    src/hash.c \
    src/kdf.c \
    src/mini-gmp.c \
    src/poly1305.c \
    src/rsa.c \
    src/salsa20.c \
    src/siphash.c \
    src/test-internal.c \
    src/util.c
else
  gcc -g \
    -Wall \
    -Wextra \
    -Wno-unused-function \
    -Wno-unused-const-variable \
    -Wno-unused-parameter \
    -Wno-implicit-fallthrough \
    -std=c89 \
    -O3 \
    "$def" \
    -DTORSION_USE_GMP \
    -I./include \
    -lgmp \
    -o test \
    src/aead.c \
    src/chacha20.c \
    src/drbg.c \
    src/dsa.c \
    src/hash.c \
    src/kdf.c \
    src/poly1305.c \
    src/rsa.c \
    src/salsa20.c \
    src/siphash.c \
    src/test-internal.c \
    src/util.c
fi
