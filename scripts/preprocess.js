#!/usr/bin/env node

'use strict';

const assert = require('assert');
const fs = require('fs');
const path = require('path');

const TO_BYTES =
  /static void fiat_[^_]+_to_bytes\(uint8_t out1\[(\d+)\][^}]+}/

const FROM_BYTES =
  /static void fiat_[^_]+_from_bytes\(uint(?:32|64)_t out1\[\d+\], const uint8_t arg1\[(\d+)\][^}]+}/

const MASK =
  /(uint(?:32|64)_t x1 = \(\(uint(?:32|64)_t\)(?:\(fiat_[^_]+_uint1\))?\(arg1\[\d+\])(\) << \d+\);)/

const ALL = false;

// Examples:
// uint32_t x1 = ((uint32_t)(fiat_p521_uint1)(arg1[65]) << 26);
// uint64_t x1 = ((uint64_t)(fiat_p521_uint1)(arg1[65]) << 56);
// uint32_t x1 = ((uint32_t)(arg1[31]) << 18);
// uint64_t x1 = ((uint64_t)(arg1[31]) << 44);
// uint32_t x1 = ((uint32_t)(arg1[31]) << 22);
// uint64_t x1 = ((uint64_t)(arg1[31]) << 47);
const masks = {
  'p521_32.h': '0x01',
  'p521_64.h': '0x01',
  'p25519_32.h': '0x7f',
  'p25519_64.h': '0x7f',
  'p251_32.h': '0x07',
  'p251_64.h': '0x07'
};

function preprocess(file) {
  const name = path.basename(file);

  console.log('Rewriting %s.', name);

  let text = fs.readFileSync(file, 'utf8').trim() + '\n';

  // Fix p224 64 bit.
  if (ALL && name === 'p224_64.h') {
    const size1 = TO_BYTES.exec(text)[1] >>> 0;
    const size2 = FROM_BYTES.exec(text)[1] >>> 0;

    if (size1 === 32) {
      text = text.replace('uint8_t out1[32]', 'uint8_t out1[28]');
      text = text.replace(/  out1\[28\] = 0x0;[^}]+/, '');
    }

    if (size2 === 32)
      text = text.replace('uint8_t arg1[32]', 'uint8_t arg1[28]');
  }

  // Mask top byte properly.
  if (ALL && masks[name]) {
    const m = `UINT8_C(${masks[name]})`;

    text = text.replace(MASK, '$1 & ' + m + '$2');
  }

  // Big endian to_bytes.
  if (ALL && !text.includes('_to_bytes_be')) {
    const m = TO_BYTES.exec(text);

    let str = m[0].replace('_to_bytes', '_to_bytes_be');
    let size = m[1] >>> 0;

    str = str.replace(/out1\[\d+\]/g, () => {
      return `out1[${size--}]`;
    });

    assert(size === -1);

    text += '\n';
    text += str;
    text += '\n';
  }

  // Big endian from_bytes.
  if (ALL && !text.includes('_from_bytes_be')) {
    const m = FROM_BYTES.exec(text);

    let str = m[0].replace('_from_bytes', '_from_bytes_be');
    let size = m[1] >>> 0;
    let index = 0;

    str = str.replace(/arg1\[(\d+)\]/g, (s, n) => {
      if ((n >>> 0) === size)
        return s;

      return `arg1[${index++}]`;
    });

    assert(index === size);

    text += '\n';
    text += str;
    text += '\n';
  }

  // Get GCC to stop complaining about __int128.
  if (!text.includes('torsion_uint128_t')
      && !text.includes('torsion_int128_t')) {
    text = text.replace(/unsigned __int128/g, 'torsion_uint128_t');
    text = text.replace(/signed __int128/g, 'torsion_int128_t');
  }

  fs.writeFileSync(file, text, 'utf8');
}

function main() {
  const prefix = path.resolve(__dirname, '..', 'src', 'fields');
  const list = fs.readdirSync(prefix);

  for (const name of list) {
    if (name.includes('_') && name.endsWith('.h'))
      preprocess(path.join(prefix, name));
  }
}

main();
