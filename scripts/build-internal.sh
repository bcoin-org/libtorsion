#!/bin/sh

set -ex

def1='-DTORSION_USE_64BIT'
def2='-DTORSION_USE_ASM'

if test x"$1" = x'32'; then
  def1='-DTORSION_USE_32BIT'
  def2='-DTORSION_NO_ASM'
fi

if test x"$1" = x'64'; then
  def2='-DTORSION_NO_ASM'
fi

gcc -g \
  -std=c89 \
  -pedantic \
  -Wall \
  -Wextra \
  -Wshadow \
  -Wno-implicit-fallthrough \
  -Wno-declaration-after-statement \
  -Wno-long-long \
  -Wno-overlength-strings \
  -O3 \
  "$def1" \
  "$def2" \
  -I./include \
  -o test/test-ecc-internal \
  src/aead.c \
  src/asn1.c \
  src/chacha20.c \
  src/drbg.c \
  src/dsa.c \
  src/hash.c \
  src/internal.c \
  src/kdf.c \
  src/mpi.c \
  src/poly1305.c \
  src/rsa.c \
  src/salsa20.c \
  src/siphash.c \
  test/test-ecc-internal.c \
  src/util.c
