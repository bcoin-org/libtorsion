/*!
 * run-wasi.js - wasi-runner for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

'use strict';

const fs = require('fs');
const {WASI} = require('wasi');
const code = fs.readFileSync(process.argv[2])
const args = process.argv.slice(3);

args.unshift('wasi');

const wasi = new WASI({
  args,
  env: process.env
});

const info = {
  env: {
    emscripten_notify_memory_growth: (memoryIndex) => {}
  },
  wasi_snapshot_preview1: wasi.wasiImport
};

const module_ = new WebAssembly.Module(code);
const instance = new WebAssembly.Instance(module_, info);

wasi.start(instance);
