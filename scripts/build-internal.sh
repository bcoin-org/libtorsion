#!/bin/sh

set -ex

def='-DTORSION_USE_64BIT'

if test x"$1" = x'32'; then
  def='-DTORSION_32BIT'
fi

if test x"$2" = x'mini'; then
  gcc -g \
    -std=c89 \
    -pedantic \
    -Wall \
    -Wextra \
    -Wshadow \
    -Wno-unused-function \
    -Wno-implicit-fallthrough \
    -Wno-unused-parameter \
    -Wno-sign-compare \
    -Wno-declaration-after-statement \
    -Wno-long-long \
    -Wno-overlength-strings \
    -O3 \
    "$def" \
    -I./include \
    -o test-internal \
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
    -std=c89 \
    -pedantic \
    -Wall \
    -Wextra \
    -Wshadow \
    -Wno-unused-function \
    -Wno-implicit-fallthrough \
    -Wno-declaration-after-statement \
    -Wno-long-long \
    -Wno-overlength-strings \
    -O3 \
    "$def" \
    -DTORSION_USE_GMP \
    -I./include \
    -lgmp \
    -o test-internal \
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
