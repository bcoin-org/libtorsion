ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(top_srcdir)/include

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
                 test/data/aead-vectors.h           \
                 test/data/chacha20-vectors.h       \
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
                 test/os.h                          \
                 test/testutil.h

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

CLEANFILES =
EXTRA_DIST = LICENSE README.md

TESTS =
noinst_PROGRAMS =

if EMSCRIPTEN
torsion_emflags = -s WASM=1                       \
                  -s STANDALONE_WASM=1            \
                  -s WASM_BIGINT=1                \
                  -s ALLOW_MEMORY_GROWTH=1        \
                  -s INITIAL_MEMORY=16777216      \
                  -s MAXIMUM_MEMORY=2147483648    \
                  -s TOTAL_STACK=5242880          \
                  -s ERROR_ON_UNDEFINED_SYMBOLS=0

torsion_emexports = -s EXPORTED_FUNCTIONS=@etc/exports.json

noinst_PROGRAMS += torsion.wasm
torsion_wasm_CFLAGS = -DTORSION_BUILD -D__wasi__
torsion_wasm_LDFLAGS = $(torsion_emflags) $(torsion_emexports)
torsion_wasm_SOURCES = $(torsion_sources)

if ENABLE_TEST
noinst_PROGRAMS += tests.wasm
tests_wasm_CFLAGS = -D__wasi__
tests_wasm_LDFLAGS = $(torsion_emflags)
tests_wasm_SOURCES = test/test.c test/os.c $(torsion_sources)
endif
else !EMSCRIPTEN
lib_LTLIBRARIES = libtorsion.la
libtorsion_la_CFLAGS = -DTORSION_BUILD
libtorsion_la_LDFLAGS = -no-undefined -version-info 0:0:0
libtorsion_la_SOURCES = $(torsion_sources)

if ENABLE_COVERAGE
libtorsion_la_CFLAGS += --coverage
libtorsion_la_LDFLAGS += --coverage
endif

if ENABLE_TEST
TESTS += tests
noinst_PROGRAMS += tests
tests_CFLAGS =
tests_LDFLAGS = -static
tests_SOURCES = test/test.c test/os.c
tests_LDADD = libtorsion.la
if ENABLE_RNG
tests_CFLAGS += @PTHREAD_CFLAGS@
tests_LDADD += @PTHREAD_LIBS@
endif
endif

if ENABLE_CTIME
noinst_PROGRAMS += ctime-test
ctime_test_SOURCES = test/ctime-test.c
ctime_test_LDADD = libtorsion.la
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libtorsion.pc
endif !EMSCRIPTEN
