#!/bin/sh

cd `dirname "$0"`

exec node \
  --no-warnings \
  --experimental-wasi-unstable-preview1 \
  --experimental-wasm-bigint \
  ./wasi-test.js "$@"
