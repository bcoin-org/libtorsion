#!/bin/bash

sources=(
  'src/aead.c'
  'src/asn1.c'
  'src/chacha20.c'
  'src/cipher.c'
  'src/ecc.c'
  'src/encoding.c'
  'src/drbg.c'
  'src/dsa.c'
  'src/hash.c'
  'src/internal.c'
  'src/kdf.c'
  'src/mpi.c'
  'src/poly1305.c'
  'src/rsa.c'
  'src/salsa20.c'
  'src/secretbox.c'
  'src/siphash.c'
  'src/util.c'
  'src/entropy/env.c'
  'src/entropy/hw.c'
  'src/entropy/sys.c'
  'src/rand.c'
)

tests=(
  'test/test.c'
  'test/hrtime.c'
)

cflags=(
  '-std=gnu11'
  '-Wall'
  '-Wextra'
  '-Wcast-align'
  '-Wno-implicit-fallthrough'
  '-Wshadow'
  '-I./include'
  '-O3'
)

emflags=(
  '-s' 'WASM=1'
  '-s' 'SINGLE_FILE=1'
  '-s' 'ALLOW_MEMORY_GROWTH=1'
  '-s' 'INITIAL_MEMORY=16777216'
  '-s' 'MAXIMUM_MEMORY=2147483648'
  '-s' 'TOTAL_STACK=5242880'
  '-s' 'EMIT_EMSCRIPTEN_METADATA=1'
  '-s' 'ENVIRONMENT=node'
)

emexports=(
  '-s' 'EXPORTED_FUNCTIONS=["_malloc","_free"]'
  '-s' 'EXTRA_EXPORTED_RUNTIME_METHODS=["stackAlloc","stackSave","stackRestore"]'
)

set -ex

test -f include/torsion/ecc.h

emcc "${sources[@]}" \
  -o torsion.js "${emflags[@]}" "${emexports[@]}" "${cflags[@]}" \
  -DTORSION_BUILD -DTORSION_NO_ASSERT

emcc "${sources[@]}" "${tests[@]}" \
  -o tests.js "${emflags[@]}" "${cflags[@]}" \
  -DTORSION_TEST -DTORSION_HAVE_RNG
