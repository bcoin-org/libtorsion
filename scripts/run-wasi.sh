#!/bin/sh

DIR=`dirname "$0"`

exec node                               \
  --no-warnings                         \
  --experimental-wasi-unstable-preview1 \
  --experimental-wasm-bigint            \
  $DIR/run-wasi.js "$@"
