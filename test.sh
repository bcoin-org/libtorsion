#!/bin/sh

def='-DFOOBAR'

if test x"$1" = x'bench'; then
  def='-DBENCH';
fi

set -ex
gcc -g \
  -Wall \
  -Wextra \
  -Wno-unused-function \
  -Wno-unused-const-variable \
  -Wno-unused-parameter \
  -std=c89 \
  -O3 \
  "$def" \
  -lcrypto -lgmp \
  -o test \
  test.c

./test
