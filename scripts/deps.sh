#!/bin/sh

# deps.sh - dependency generator for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

O=o
C=-c

if test x"$1" = x"win32"; then
  O=obj
  C=/c
fi

get_deps() {
  gcc -I ./include -DTORSION_HAVE_RNG -MF - -MM "$1" \
    | cut -d' ' -f3- -z                              \
    | tr -d '\\'                                     \
    | tr '\n' ' '                                    \
    | sed -e 's/ \+/ /g'                             \
    | tr ' ' '\n'                                    \
    | sort -u
}

format_rule() {
  rule="$(basename $1 | sed "s/c$/$O/"): $1 "
  space=$(echo "$rule" | sed 's/./ /g')
  get_deps $1 | while read line; do
    if test x"$line" = x; then
      continue
    fi
    echo "$rule$line                       \\"
    rule="$space"
  done
  printf '\t'
  echo '$(CC) '$C' $(LIB_CFLAGS) '$1
  echo ''
}

echo '#'
echo '# Library Objects'
echo '#'
echo ''

format_rule src/aead.c
format_rule src/asn1.c
format_rule src/cipher.c
format_rule src/ecc.c
format_rule src/encoding.c
format_rule src/entropy/hw.c
format_rule src/entropy/sys.c
format_rule src/drbg.c
format_rule src/dsa.c
format_rule src/hash.c
format_rule src/ies.c
format_rule src/internal.c
format_rule src/kdf.c
format_rule src/mac.c
format_rule src/mpi.c
format_rule src/rand.c
format_rule src/rsa.c
format_rule src/stream.c
format_rule src/util.c

echo '#'
echo '# Test Objects'
echo '#'
echo ''

format_rule test/bench.c
format_rule test/ctgrind.c
format_rule test/hrtime.c
format_rule test/test.c
format_rule test/utils.c
