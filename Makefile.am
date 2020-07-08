# Makefile.am - automake file for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(top_srcdir)/include

CLEANFILES =
EXTRA_DIST = LICENSE README.md

TESTS =
noinst_PROGRAMS =

include_HEADERS = include/torsion/aead.h      \
                  include/torsion/cipher.h    \
                  include/torsion/common.h    \
                  include/torsion/drbg.h      \
                  include/torsion/dsa.h       \
                  include/torsion/ecc.h       \
                  include/torsion/encoding.h  \
                  include/torsion/hash.h      \
                  include/torsion/ies.h       \
                  include/torsion/kdf.h       \
                  include/torsion/mac.h       \
                  include/torsion/rsa.h       \
                  include/torsion/stream.h    \
                  include/torsion/util.h

noinst_HEADERS = src/asn1.h                         \
                 src/bio.h                          \
                 src/fields/libsecp256k1_32.h       \
                 src/fields/libsecp256k1_64.h       \
                 src/fields/p192_32.h               \
                 src/fields/p192_64.h               \
                 src/fields/p192.h                  \
                 src/fields/p224_32.h               \
                 src/fields/p224_64.h               \
                 src/fields/p224.h                  \
                 src/fields/p251_32.h               \
                 src/fields/p251_64.h               \
                 src/fields/p251.h                  \
                 src/fields/p25519_32.h             \
                 src/fields/p25519_64.h             \
                 src/fields/p25519.h                \
                 src/fields/p256_32.h               \
                 src/fields/p256_64.h               \
                 src/fields/p256.h                  \
                 src/fields/p384_32.h               \
                 src/fields/p384_64.h               \
                 src/fields/p384.h                  \
                 src/fields/p448_32.h               \
                 src/fields/p448_64.h               \
                 src/fields/p448.h                  \
                 src/fields/p521_32.h               \
                 src/fields/p521_64.h               \
                 src/fields/p521.h                  \
                 src/fields/scalar.h                \
                 src/fields/secp256k1_32.h          \
                 src/fields/secp256k1_64.h          \
                 src/fields/secp256k1.h             \
                 src/internal.h                     \
                 src/mpi.h                          \
                 src/subgroups.h                    \
                 src/tls.h                          \
                 test/data/chacha20-vectors.h       \
                 test/data/chachapoly-vectors.h     \
                 test/data/cipher-aead-vectors.h    \
                 test/data/cipher-mode-vectors.h    \
                 test/data/cipher-vectors.h         \
                 test/data/ctr-drbg-vectors.h       \
                 test/data/dsa-vectors.h            \
                 test/data/eb2k-vectors.h           \
                 test/data/ecdsa-vectors.h          \
                 test/data/eddsa-vectors.h          \
                 test/data/hash-drbg-vectors.h      \
                 test/data/hash-vectors.h           \
                 test/data/hkdf-vectors.h           \
                 test/data/hmac-drbg-vectors.h      \
                 test/data/hmac-vectors.h           \
                 test/data/pbkdf2-vectors.h         \
                 test/data/poly1305-vectors.h       \
                 test/data/rsa-vectors.h            \
                 test/data/schnorr-legacy-vectors.h \
                 test/data/schnorr-vectors.h        \
                 test/ecc_internal.h                \
                 test/hrtime.h                      \
                 test/thread.h                      \
                 test/utils.h

if ENABLE_RNG
include_HEADERS += include/torsion/rand.h
noinst_HEADERS += src/entropy/entropy.h
endif

torsion_sources = src/aead.c      \
                  src/asn1.c      \
                  src/cipher.c    \
                  src/ecc.c       \
                  src/encoding.c  \
                  src/drbg.c      \
                  src/dsa.c       \
                  src/hash.c      \
                  src/ies.c       \
                  src/internal.c  \
                  src/kdf.c       \
                  src/mac.c       \
                  src/mpi.c       \
                  src/rsa.c       \
                  src/stream.c    \
                  src/util.c

if ENABLE_RNG
torsion_sources += src/entropy/env.c
torsion_sources += src/entropy/hw.c
torsion_sources += src/entropy/sys.c
torsion_sources += src/rand.c
endif

bench_sources = test/bench.c test/hrtime.c
test_sources = test/test.c

if EMSCRIPTEN
LOG_COMPILER = ./scripts/run-wasi.sh

torsion_emflags = -s WASM=1                       \
                  -s STANDALONE_WASM=1            \
                  -s WASM_BIGINT=1                \
                  -s ALLOW_MEMORY_GROWTH=1        \
                  -s INITIAL_MEMORY=16777216      \
                  -s MAXIMUM_MEMORY=2147483648    \
                  -s TOTAL_STACK=5242880          \
                  -s ERROR_ON_UNDEFINED_SYMBOLS=0 \
                  -s WARN_ON_UNDEFINED_SYMBOLS=0

torsion_emexports = -s EXPORTED_FUNCTIONS=@etc/exports.json

noinst_PROGRAMS += torsion.wasm
torsion_wasm_SOURCES = $(torsion_sources)
torsion_wasm_CFLAGS = -D__wasi__ -DTORSION_BUILD
torsion_wasm_LDFLAGS = $(torsion_emflags) $(torsion_emexports)

noinst_PROGRAMS += torsion_bench.wasm
torsion_bench_wasm_SOURCES = $(torsion_sources) $(bench_sources)
torsion_bench_wasm_CFLAGS = -D__wasi__
torsion_bench_wasm_LDFLAGS = $(torsion_emflags)

TESTS += torsion_test.wasm
noinst_PROGRAMS += torsion_test.wasm
torsion_test_wasm_SOURCES = $(torsion_sources) $(test_sources)
torsion_test_wasm_CFLAGS = -g -D__wasi__ -DTORSION_DEBUG
torsion_test_wasm_LDFLAGS = $(torsion_emflags)
else !EMSCRIPTEN
lib_LTLIBRARIES = libtorsion.la
libtorsion_la_SOURCES = $(torsion_sources)
libtorsion_la_CFLAGS = -DTORSION_BUILD
libtorsion_la_LDFLAGS = -no-undefined -version-info @TORSION_LIB_VERSION@
libtorsion_la_LIBADD =

if ENABLE_PTHREAD
libtorsion_la_CFLAGS += @PTHREAD_CFLAGS@
libtorsion_la_LIBADD += @PTHREAD_LIBS@
endif

noinst_PROGRAMS += torsion_bench
torsion_bench_SOURCES = $(bench_sources)
torsion_bench_LDFLAGS = -no-install
torsion_bench_LDADD = libtorsion.la

TESTS += torsion_test
noinst_PROGRAMS += torsion_test
torsion_test_SOURCES = $(test_sources)
torsion_test_CFLAGS =
torsion_test_LDFLAGS = -no-install
torsion_test_LDADD = libtorsion.la

if ENABLE_THREADS
torsion_test_SOURCES += test/thread.c
torsion_test_CFLAGS += @PTHREAD_CFLAGS@
torsion_test_LDADD += @PTHREAD_LIBS@
endif

TESTS += torsion_test_static
noinst_PROGRAMS += torsion_test_static
torsion_test_static_SOURCES := $(torsion_test_SOURCES)
torsion_test_static_CFLAGS := $(torsion_test_CFLAGS)
torsion_test_static_LDFLAGS = -static
torsion_test_static_LDADD := $(torsion_test_LDADD)

if ENABLE_COVERAGE
libtorsion_la_CFLAGS += -O0 --coverage
libtorsion_la_LIBADD += -lgcov
torsion_test_static_LDADD += -lgcov
endif

if ENABLE_CTIME
noinst_PROGRAMS += torsion_ctime
torsion_ctime_SOURCES = test/ctime.c
torsion_ctime_LDFLAGS = -no-install
torsion_ctime_LDADD = libtorsion.la
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libtorsion.pc
endif !EMSCRIPTEN
