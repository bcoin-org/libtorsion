#!/bin/sh

set -ex

prefix="$1"
solinas="$prefix/src/ExtractionOCaml/unsaturated_solinas"
montgomery="$prefix/src/ExtractionOCaml/word_by_word_montgomery"

if test -z "$prefix" -o ! -e "$solinas"; then
  echo 'Please specify a valid fiat location.'
  exit 1
fi

# 2^192 - 2^64 - 1
"$solinas" --static 'p192' '9' '2^192 - 2^64 - 1' '32' > ./src/fields/p192_32.h
"$solinas" --static 'p192' '4' '2^192 - 2^64 - 1' '64' > ./src/fields/p192_64.h

# 2^224 - 2^96 + 1
"$montgomery" --static 'p224' '2^224 - 2^96 + 1' '32' > ./src/fields/p224_32.h
"$montgomery" --static 'p224' '2^224 - 2^96 + 1' '64' > ./src/fields/p224_64.h

# 2^256 - 2^224 + 2^192 + 2^96 - 1
"$montgomery" --static 'p256' '2^256 - 2^224 + 2^192 + 2^96 - 1' '32' > ./src/fields/p256_32.h
"$montgomery" --static 'p256' '2^256 - 2^224 + 2^192 + 2^96 - 1' '64' > ./src/fields/p256_64.h

# 2^384 - 2^128 - 2^96 + 2^32 - 1
"$montgomery" --static 'p384' '2^384 - 2^128 - 2^96 + 2^32 - 1' '32' > ./src/fields/p384_32.h
"$montgomery" --static 'p384' '2^384 - 2^128 - 2^96 + 2^32 - 1' '64' > ./src/fields/p384_64.h

# 2^521 - 1
"$solinas" --static 'p521' '19' '2^521 - 1' '32' > ./src/fields/p521_32.h
"$solinas" --static 'p521' '9' '2^521 - 1' '64' > ./src/fields/p521_64.h

# 2^256 - 2^32 - 977
"$montgomery" --static 'secp256k1' '2^256 - 2^32 - 977' '32' > ./src/fields/secp256k1_32.h
"$montgomery" --static 'secp256k1' '2^256 - 2^32 - 977' '64' > ./src/fields/secp256k1_64.h

# 2^256 - 2^32 - 977 (faster for GCC, slower for clang)
# "$solinas" --static 'secp256k1' '12' '2^256 - 2^32 - 977' '32' > ./src/fields/secp256k1_32.h
# "$solinas" --static 'secp256k1' '8' '2^256 - 2^32 - 977' '64' > ./src/fields/secp256k1_64.h

# 2^255 - 19
"$solinas" --static 'p25519' '10' '2^255 - 19' '32' \
  carry_mul \
  carry_square \
  carry_scmul121666 \
  carry \
  add \
  sub \
  opp \
  selectznz \
  to_bytes \
  from_bytes > ./src/fields/p25519_32.h

"$solinas" --static 'p25519' '5' '2^255 - 19' '64' \
  carry_mul \
  carry_square \
  carry_scmul121666 \
  carry \
  add \
  sub \
  opp \
  selectznz \
  to_bytes \
  from_bytes > ./src/fields/p25519_64.h

# 2^448 - 2^224 - 1
"$solinas" --static 'p448' '18' '2^448 - 2^224 - 1' '32' > ./src/fields/p448_32.h
"$solinas" --static 'p448' '8' '2^448 - 2^224 - 1' '64' > ./src/fields/p448_64.h

# 2^251 - 9
"$solinas" --static 'p251' '10' '2^251 - 9' '32' > ./src/fields/p251_32.h
"$solinas" --static 'p251' '5' '2^251 - 9' '64' > ./src/fields/p251_64.h

./scripts/preprocess.js
