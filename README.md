# libtorsion

![cmake](https://github.com/bcoin-org/libtorsion/workflows/cmake/badge.svg)

Crypto library in C, primarily for EC. Still semi-experimental.

Used as the backend for [bcrypto].

## Design Approach

- Generic, without sacrificing speed or security.
- Use field element backends generated by [fiat-crypto] (which have proofs of
  correctness, formal verification, etc).
- Use a generic barrett reduction for scalar arithmetic (we use a
  reimplementation of GMP's [mpn_] functions for this).
- Keep constant-time multiplication algorithms as simple as possible.
- Constant time, stack-based, and so on.
- Expose a simple opaque API which accepts all arguments in wire format
  (this eases implementation for language bindings, wasm, and other ffis).

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
    - curve1174

- Schemes:
    - ECDSA
    - [BIP-Schnorr][schnorr]
    - [BIP340][bip340]
    - ECDH
    - EdDSA
    - [Ristretto][ristretto]

## Example

``` c
#include <assert.h>
#include <stddef.h>
#include <torsion/ecc.h>
#include <torsion/hash.h>
#include <torsion/rand.h>

int main(void) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  const char str[] = "hello world";
  unsigned char priv[32];
  unsigned char entropy[32];
  unsigned char msg[32];
  unsigned char pub[33];
  unsigned char sig[64];
  sha256_t hash;

  /* Hash our message with sha256. */
  sha256_init(&hash);
  sha256_update(&hash, str, sizeof(str) - 1);
  sha256_final(&hash, msg);

  /* Grab some entropy from the OS. */
  assert(torsion_getentropy(entropy, 32));

  /* Generate a key pair using entropy. */
  ecdsa_privkey_generate(ec, priv, entropy);
  assert(ecdsa_pubkey_create(ec, pub, NULL, priv, 1));

  /* Sign in constant time. */
  assert(ecdsa_sign(ec, sig, NULL, msg, sizeof(msg), priv));

  /* Verify our signature. */
  assert(ecdsa_verify(ec, msg, sizeof(msg), sig, pub, sizeof(pub)));

  /* Cleanup. */
  wei_curve_destroy(ec);

  return 0;
}
```

Compile with:

``` bash
$ cc -o example -I/path/to/libtorsion/include example.c /path/to/libtorsion.a
```

## Building

libtorsion supports a CMake build (recommended) and a fallback build for
systems which do not have CMake available. The "fallback build" consists of a
POSIX-conforming Makefile on unix-like OSes and an NMake file for Windows.

### Unix

#### CMake

The CMake build is fairly straightforward and offers the following options:

- `TORSION_ENABLE_ASM=ON` - Use inline assembly if available
- `TORSION_ENABLE_COVERAGE=OFF` - Enable coverage
- `TORSION_ENABLE_DEBUG=ON` - Enable debug build (forces -g or /Zi)
- `TORSION_ENABLE_INT128=ON` - Use `__int128` if available
- `TORSION_ENABLE_PIC=ON` - Enable PIC
- `TORSION_ENABLE_PTHREAD=ON` - Use pthread if present in libc
- `TORSION_ENABLE_RNG=ON` - Enable RNG
- `TORSION_ENABLE_SHARED=ON` - Build shared library
- `TORSION_ENABLE_TESTS=ON` - Build tests
- `TORSION_ENABLE_TLS=ON` - Use thread-local storage if available
- `TORSION_ENABLE_VERIFY=OFF` - Enable scalar bounds checks

To build:

``` sh
$ cmake . -DCMAKE_BUILD_TYPE=Release
$ make
```

#### Make

Our Makefile strives for perfect POSIX conformance. This means flags like
`-fPIC` aren't assumed as position independent code isn't specified by POSIX.

Likewise, warnings and optimization flags are not assumed as POSIX only
specifies `-O`.

Example for Linux/\*BSD/Darwin:

``` sh
$ make -f Makefile.unix CFLAGS='-fPIC -O3'
```

Example for AIX:

``` sh
$ make -f Makefile.unix CC=xlc_r CFLAGS='-qpic -qoptimize=3 -qtls -DHAVE_QTLS'
```

Example for Solaris:

``` sh
$ make -f Makefile.unix CFLAGS='-KPIC -xO3 -D_REENTRANT'
```

Example for HP-UX:

``` sh
$ make -f Makefile.unix CFLAGS='+Z +O3'
```

Linking a shared library is up to you at that point. The makefile will generate
"object libraries" to make things easy. For example:

``` sh
$ gcc -o libtorsion.so -shared -fPIC @libtorsion.obj
$ gcc -o torsion_test @torsion_test.obj libtorsion.so
```

On Darwin:

``` sh
$ gcc -o libtorsion.dylib -dynamiclib @libtorsion.obj
$ gcc -o torsion_test @torsion_test.obj libtorsion.dylib
```

Or, if you have libtool installed, you can let it do the heavy lifting:

``` sh
$ make -f Makefile.unix all libtorsion.la CFLAGS=-fPIC
```

### Windows

Builds on windows will produce both a static and shared library. To deal with
the naming collision, the static import library is called `libtorsion.lib`
whereas the shared import library is named `torsion.lib` (this follows
Microsoft's new naming convention).

#### CMake

See the unix cmake build documentation above for a list of available options.

MSVC:

``` sh
$ cmake .
$ cmake --build . --config Release
```

Clang (better performance):

``` sh
$ cmake . -G 'NMake Makefiles' -DCMAKE_C_COMPILER=clang-cl \
                               -DCMAKE_BUILD_TYPE=Release
$ nmake
```

#### NMake

The NMake build assumes MSVC (cl.exe), but also works with Windows Clang
(clang-cl).

``` sh
$ nmake /F Makefile.nmake
```

``` sh
$ nmake /F Makefile.nmake CC=clang-cl
```

### MinGW

Another way to build for Windows is by cross-compiling with MinGW (useful for
those averse to Windows or who do not have a Windows machine to build on).

#### CMake

``` sh
$ ./scripts/mingw-cmake cmake . -DCMAKE_BUILD_TYPE=Release
$ make
```

#### Make

``` sh
$ ./scripts/mingw-make make -f Makefile.unix CFLAGS=-O3 mingw
```

### WASI

#### CMake

``` sh
$ ./scripts/wasi-cmake cmake . -DCMAKE_BUILD_TYPE=Release
$ make
```

#### Make

``` sh
$ ./scripts/wasi-make make -f Makefile.unix CFLAGS=-O3 wasm
```

### Emscripten

#### CMake

``` sh
$ emcmake cmake . -DCMAKE_BUILD_TYPE=Release
$ make
```

#### Make

``` sh
$ emmake make -f Makefile.unix CFLAGS=-O3 js
```

### CMake Subprojects

When CMakeLists.txt is invoked as a CMake subproject, it will expose a single
static library target named `torsion`, skipping all other targets (tests,
benchmarks, etc).

Example:

``` cmake
add_subdirectory(deps/libtorsion)
target_link_libraries(my_project PRIVATE torsion)
```

Note that `torsion` will be an object library if you are doing a WASI or
Emscripten build.

### Alternate Build Systems

libtorsion was written generically enough to support any build system you might
want to use, with only a few key configuration options necessary.

#### Configuration

To incorporate libtorsion into a different build system you must be aware of
various configuration options that torsion uses internally.

Any of the following _may_ be passed as defines to the preprocessor:

- `TORSION_HAVE_CONFIG` - Disables preprocessor-based autoconfiguration¹.
- `TORSION_HAVE_ASM` - GNU-flavored inline assembly is available².
- `TORSION_HAVE_INT128` - The `__int128` type is available².
- `TORSION_HAVE_PTHREAD` - The pthread API is available².
- `TORSION_TLS=[tls-keyword]` - Thread-local storage is available².
- `TORSION_COVERAGE` - Coverage is enabled via `gcov`. Disable assertions.
- `TORSION_DEBUG` - Enable assertions.
- `TORSION_VERIFY` - Enable extra debugging checks.
- `TORSION_HAVE_CLOCK_GETTIME` - `clock_gettime(3)` is available²³.
- `TORSION_HAVE_FORK` - `fork(2)` and `waitpid(2)` are available²³.
- `TORSION_HAVE_GETTIMEOFDAY` - `gettimeofday(3)` is available²³.
- `TORSION_HAVE_RNG` - libtorsion was compiled with RNG support³.
- `TORSION_HAVE_TIME` - `time(2)` is available²³.
- `TORSION_HAVE_ZLIB` - `<zlib.h>` is available and we are linked with `-lz`³.

Footnotes:

1. By default torsion will attempt to autoconfigure some options by using the
   C preprocessor to detect features.
2. `TORSION_HAVE_CONFIG` must be defined for this option to have any effect.
3. This option affects tests only.

---

If you are using a build system where detecting features is difficult, you can
have libtorsion attempt to auto-detect features using the C preprocessor:
simply do not pass `-DTORSION_HAVE_CONFIG`.

#### Notes

libtorsion is written to be compatible with C89/C90 (with one caveat: the
system must provide an stdint.h header). This means that passing `-std=c89`
should be fine unless you are building with RNG support. Various C standard
libraries base their default features on the language standard and will not
expose certain POSIX APIs (or system-specific APIs) if a strict C version is
passed.

This also means you should not pass any feature test macros that may downgrade
the features apart from what is necessary for the RNG to function. For example,
passing `_POSIX_SOURCE` or `_POSIX_C_SOURCE` will probably break the RNG on
Linux (assuming glibc is used).

Furthermore, the RNG accesses `errno`, which means a thread-local `errno` is
desirable. By defalt, certain OSes do not expose a thread-safe `errno`. This
includes Solaris and AIX. Please ensure that your build produces thread-safe
code (the flags required for this differ from system to system).

The libtorsion codebase is also valid C++ and can be compiled inside your C++
project without any modifications.

## Contribution and License Agreement

If you contribute code to this project, you are implicitly allowing your code
to be distributed under the MIT license. You are also implicitly verifying that
all code is your original work. `</legalese>`

## License

- Copyright (c) 2020, Christopher Jeffrey (MIT License).

See LICENSE for more info.

[bcrypto]: https://github.com/bcoin-org/bcrypto
[fiat-crypto]: https://github.com/mit-plv/fiat-crypto
[mpn_]: https://gmplib.org/manual/Low_002dlevel-Functions.html
[schnorr]: https://github.com/sipa/bips/blob/d194620/bip-schnorr.mediawiki
[bip340]: https://github.com/bitcoin/bips/blob/master/bip-0340.mediawiki
[ristretto]: https://ristretto.group/
