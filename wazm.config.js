'use strict';

// Build file for wazm.
// Usage: $ wazm wazm.config.js
//        $ wazm-run tests.wasm
// See: https://github.com/chjj/wazm

const sources = [
  'src/aead.c',
  'src/asn1.c',
  'src/chacha20.c',
  'src/cipher.c',
  'src/ecc.c',
  'src/encoding.c',
  'src/drbg.c',
  'src/dsa.c',
  'src/hash.c',
  'src/internal.c',
  'src/kdf.c',
  'src/mpi.c',
  'src/poly1305.c',
  'src/rsa.c',
  'src/salsa20.c',
  'src/secretbox.c',
  'src/siphash.c',
  'src/util.c',
  'src/entropy/env.c',
  'src/entropy/hw.c',
  'src/entropy/sys.c',
  'src/rand.c'
];

const flags = [
  '-std=gnu11',
  '-Wall',
  '-Wextra',
  '-Wcast-align',
  '-Wno-implicit-fallthrough',
  '-Wshadow',
  '-O3'
];

const MB = 2 ** 20;
const GB = 2 ** 30;

module.exports = [
  {
    root: __dirname,
    target: 'torsion',
    base64: true,
    wat: true,
    sources,
    flags,
    includes: [
      './include'
    ],
    defines: [
      'TORSION_BUILD',
      'TORSION_NO_ASSERT'
    ],
    memory: {
      initial: 16 * MB,
      stack: 5 * MB,
      max: 2 * GB,
      grow: true
    }
  },
  {
    root: __dirname,
    target: 'tests',
    executable: true,
    base64: true,
    wat: true,
    sources: [
      ...sources,
      'test/test.c',
      'test/hrtime.c'
    ],
    flags,
    includes: [
      './include'
    ],
    defines: [
      'TORSION_TEST',
      'TORSION_HAVE_RNG'
    ],
    memory: {
      initial: 16 * MB,
      stack: 5 * MB,
      max: 2 * GB,
      grow: true
    }
  }
];
