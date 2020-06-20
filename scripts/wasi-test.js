'use strict';

const path = require('path');
const fs = require('fs/promises');
const {WASI} = require('wasi');
const TESTS_WASM = path.resolve(__dirname, '..', 'tests.wasm');

(async () => {
  const wasi = new WASI({
    args: process.argv.slice(2),
    env: process.env
  });

  const info = {
    env: {
      emscripten_notify_memory_growth: (memoryIndex) => {}
    },
    wasi_snapshot_preview1: wasi.wasiImport
  };

  const code = await fs.readFile(TESTS_WASM)
  const wasm = await WebAssembly.compile(code);
  const instance = await WebAssembly.instantiate(wasm, info);

  wasi.start(instance);
})().catch((err) => {
  console.error(err.stack);
  process.exit(1);
});
