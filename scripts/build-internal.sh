#!/bin/sh

set -ex

def1='-DTORSION_USE_64BIT'
def2='-DTORSION_USE_ASM'

if test x"$1" = x'32'; then
  def1='-DTORSION_USE_32BIT'
  def2='-DTORSION_NO_ASM'
fi

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
  "$def1" \
  "$def2" \
  -I./include \
  -o test-internal \
  src/aead.c \
  src/chacha20.c \
  src/drbg.c \
  src/dsa.c \
  src/hash.c \
  src/kdf.c \
  src/mpi.c \
  src/poly1305.c \
  src/prime.c \
  src/rsa.c \
  src/salsa20.c \
  src/siphash.c \
  src/test-internal.c \
  src/util.c
