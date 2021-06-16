#!/usr/bin/env node

/*!
 * prngd.js - simple prngd implementation for testing
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

'use strict';

const crypto = require('crypto');
const fs = require('fs');
const net = require('net');
const server = net.createServer();
const name = process.argv[2] || '/tmp/entropy';

function log(msg) {
  console.log('[%s] %s', new Date().toISOString(), msg);
}

function destroy(socket, err) {
  try {
    socket.destroy();
  } catch (e) {
    ;
  }
  console.error(err.stack);
}

function handle(socket, data) {
  if (data.length === 0)
    return;

  const cmd = data[0];

  switch (cmd) {
    case 0: {
      if (data.length !== 1)
        throw new Error('invalid length');

      log('sending current entropy');

      socket.write(Buffer.from([0, 0, 0xff, 0xff]));

      break;
    }

    case 1:
    case 2: {
      if (data.length !== 2)
        throw new Error('invalid length');

      let len = data[1];
      let ptr = 0;

      const buf = crypto.randomBytes(len);

      log(`sending ${len} bytes`);

      socket.write(Buffer.from([len]));

      while (len > 0) {
        const h = Math.floor((len + 1) / 2);
        const n = Math.floor(Math.random() * (h + 1));

        if (n === 0)
          continue;

        socket.write(buf.slice(ptr, ptr + n));

        ptr += n;
        len -= n;
      }

      break;
    }

    case 3: {
      if (data.length < 4)
        throw new Error('invalid length');

      const bits = (data[1] << 8) | data[2];
      const size = data[3];

      if (data.length !== size + 4)
        throw new Error('invalid length');

      log(`received ${size} bytes (${bits} bits)`);

      break;
    }

    default: {
      throw new Error('invalid cmd');
    }
  }
}

server.on('error', (err) => {
  console.error(err.stack);
  process.exit(1);
});

server.on('connection', (socket) => {
  socket.on('error', (err) => {
    destroy(socket, err);
  });

  socket.on('data', (data) => {
    try {
      handle(socket, data);
    } catch (e) {
      destroy(socket, e);
    }
  });
});

let listening = false;

fs.unlink(name, () => {
  server.listen(name, () => {
    log(`listening on ${name}`);
    listening = true;
  });
});

const signalHandler = (code) => {
  return () => {
    if (listening) {
      log('shutting down');
      server.close(() => {
        log(`removing ${name}`);
        fs.unlink(name, () => {
          process.exit(code);
        });
      });
    } else {
      process.exit(code);
    }
  };
};

process.on('SIGHUP', signalHandler(0x81));
process.on('SIGINT', signalHandler(0x82));
process.on('SIGTERM', signalHandler(0x8f));
