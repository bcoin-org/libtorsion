'use strict';

const fs = require('fs');
const {WASI} = require('wasi');
const code = fs.readFileSync(process.argv[2])

// Work around a WASI bug in node.js.
let args = [''];

if (process.argv.length > 3)
  args = process.argv.slice(3);

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
