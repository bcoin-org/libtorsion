#!/bin/sh

set -ex
gcc -g \
  -Wall \
  -Wextra \
  -Wno-unused-function \
  -Wno-unused-const-variable \
  -Wno-unused-parameter \
  -std=c89 \
  -O3 \
  -DBCRYPTO_EC_64BIT \
  -DBCRYPTO_HAS_GMP \
  -lcrypto \
  -lgmp \
  -o test \
  test.c
