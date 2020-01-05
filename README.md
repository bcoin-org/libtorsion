# ec

Experimental EC library in C (currently lacking a name).

Intended to be the new EC backend for [bcrypto].

## Design Approach

- Generic, without sacrificing speed or security.
- Use field element backends generated by [fiat-crypto] (which have proofs of
  correctness, formal verification, etc).
- Use a generic barrett reduction for scalar arithmetic (leverage GMP's `mpn_`
  functions for this).
- _Simple_ constant-time multiplication algorithms to avoid any timing leakage
  (ladders only -- for short weierstrass, we use the co-z montgomery ladder).
- If that's not enough, use a random blinding factor as well.
- Constant time, stack-based, and so on.
- Dependency-less (aside from a vendored mini-gmp).

## Current Features

- Curve Support:
    - p192
    - p224
    - p256
    - p384
    - p521
    - secp256k1
    - curve25519
    - curve448
    - edwards25519
    - edwards448

- Schemes:
    - ECDSA
    - ECDH
    - EdDSA

## Todo

- Optimize functions to reduce stack usage.
- Batch verification.
- bip-schnorr support.
- Elligators.

[bcrypto]: https://github.com/bcoin-org/bcrypto
[fiat-crypto]: https://github.com/mit-plv/fiat-crypto