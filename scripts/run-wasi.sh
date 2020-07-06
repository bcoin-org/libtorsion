#!/bin/sh

# run-wasi.sh - wasi runner for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

DIR=`dirname "$0"`

exec node                               \
  --no-warnings                         \
  --experimental-wasi-unstable-preview1 \
  --experimental-wasm-bigint            \
  $DIR/run-wasi.js "$@"
