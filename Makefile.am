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
                 test/data/chacha20_vectors.h       \
                 test/data/chachapoly_vectors.h     \
                 test/data/cipher_aead_vectors.h    \
                 test/data/cipher_mode_vectors.h    \
                 test/data/cipher_vectors.h         \
                 test/data/ctr_drbg_vectors.h       \
                 test/data/dsa_vectors.h            \
                 test/data/eb2k_vectors.h           \
                 test/data/ecdsa_vectors.h          \
                 test/data/eddsa_vectors.h          \
                 test/data/hash_drbg_vectors.h      \
                 test/data/hash_vectors.h           \
                 test/data/hkdf_vectors.h           \
                 test/data/hmac_drbg_vectors.h      \
                 test/data/hmac_vectors.h           \
                 test/data/pbkdf2_vectors.h         \
                 test/data/poly1305_vectors.h       \
                 test/data/rsa_vectors.h            \
                 test/data/schnorr_legacy_vectors.h \
                 test/data/schnorr_vectors.h        \
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

if WASM
LOG_COMPILER = ./scripts/wasi-run

torsion_objects := $(subst .c,.o,$(torsion_sources))

if EMSCRIPTEN
torsion_ldflags = -s WASM=1                       \
                  -s STANDALONE_WASM=1            \
                  -s WASM_BIGINT=1                \
                  -s ALLOW_MEMORY_GROWTH=1        \
                  -s INITIAL_MEMORY=16777216      \
                  -s MAXIMUM_MEMORY=2147483648    \
                  -s TOTAL_STACK=5242880          \
                  -s ERROR_ON_UNDEFINED_SYMBOLS=0 \
                  -s WARN_ON_UNDEFINED_SYMBOLS=0

torsion_sources_lib = etc/empty.c
torsion_ldflags_lib = -s EXPORTED_FUNCTIONS=@etc/exports.json

CFLAGS += -D__wasi__ -DTORSION_BUILD
else
torsion_ldflags = -Wl,--allow-undefined         \
                  -Wl,--initial-memory=16777216 \
                  -Wl,--max-memory=2147483648   \
                  -Wl,-z -Wl,stack-size=5242880

torsion_sources_lib = etc/crt1-reactor.c etc/main.c
torsion_ldflags_lib = -Wl,--export-all
endif

noinst_PROGRAMS += torsion.wasm
torsion_wasm_SOURCES = $(torsion_sources_lib)
torsion_wasm_LDADD = $(torsion_objects)
torsion_wasm_LDFLAGS = $(torsion_ldflags) $(torsion_ldflags_lib)

noinst_PROGRAMS += torsion_bench.wasm
torsion_bench_wasm_SOURCES = $(bench_sources)
torsion_bench_wasm_LDADD = $(torsion_objects)
torsion_bench_wasm_LDFLAGS = $(torsion_ldflags)

TESTS += torsion_test.wasm
noinst_PROGRAMS += torsion_test.wasm
torsion_test_wasm_SOURCES = $(test_sources)
torsion_test_wasm_LDADD = $(torsion_objects)
torsion_test_wasm_LDFLAGS = $(torsion_ldflags)
else !WASM
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

if ENABLE_CTGRIND
noinst_PROGRAMS += torsion_ctgrind
torsion_ctgrind_SOURCES = test/ctgrind.c
torsion_ctgrind_LDFLAGS = -no-install
torsion_ctgrind_LDADD = libtorsion.la
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libtorsion.pc
endif !WASM
