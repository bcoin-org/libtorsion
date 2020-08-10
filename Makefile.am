# Makefile.am - automake file for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

#
# Initialization
#

ACLOCAL_AMFLAGS = -I m4

noinst_PROGRAMS =
lib_LTLIBRARIES =
TESTS =

licensedir = $(datadir)/licenses/$(PACKAGE_NAME)
pkgconfigdir = $(libdir)/pkgconfig

#
# Headers
#

include_HEADERS = include/torsion/aead.h     \
                  include/torsion/cipher.h   \
                  include/torsion/common.h   \
                  include/torsion/drbg.h     \
                  include/torsion/dsa.h      \
                  include/torsion/ecc.h      \
                  include/torsion/encoding.h \
                  include/torsion/hash.h     \
                  include/torsion/ies.h      \
                  include/torsion/kdf.h      \
                  include/torsion/mac.h      \
                  include/torsion/rsa.h      \
                  include/torsion/stream.h   \
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

#
# Sources & Flags
#

torsion_sources = src/aead.c     \
                  src/asn1.c     \
                  src/cipher.c   \
                  src/ecc.c      \
                  src/encoding.c \
                  src/drbg.c     \
                  src/dsa.c      \
                  src/hash.c     \
                  src/ies.c      \
                  src/internal.c \
                  src/kdf.c      \
                  src/mac.c      \
                  src/mpi.c      \
                  src/rsa.c      \
                  src/stream.c   \
                  src/util.c

rng_sources = src/entropy/env.c \
              src/entropy/hw.c  \
              src/entropy/sys.c \
              src/rand.c

if ENABLE_RNG
torsion_sources += $(rng_sources)
endif

torsion_cflags = -I$(top_srcdir)/include

bench_sources = test/bench.c test/hrtime.c
test_sources = test/test.c

#
# Emscripten
#

if EMSCRIPTEN
INITIAL_MEMORY ?= 16777216
MAX_MEMORY ?= 2147483648
STACK_SIZE ?= 5242880
ENVIRONMENT ?= node
LOG_COMPILER = node

torsion_ldflags = -s SINGLE_FILE=1                    \
                  -s ASSERTIONS=0                     \
                  -s NODEJS_CATCH_EXIT=0              \
                  -s NODEJS_CATCH_REJECTION=0         \
                  -s ALLOW_MEMORY_GROWTH=1            \
                  -s INITIAL_MEMORY=$(INITIAL_MEMORY) \
                  -s MAXIMUM_MEMORY=$(MAX_MEMORY)     \
                  -s TOTAL_STACK=$(STACK_SIZE)        \
                  -s ENVIRONMENT=$(ENVIRONMENT)

noinst_PROGRAMS += torsion.js
torsion_js_SOURCES = $(torsion_sources)
torsion_js_CFLAGS = $(torsion_cflags) -DTORSION_BUILD
torsion_js_LDFLAGS = $(torsion_ldflags) -s EXPORTED_FUNCTIONS=@etc/exports.json

if ENABLE_TESTS
noinst_PROGRAMS += torsion_bench.js
torsion_bench_js_SOURCES = $(bench_sources)
torsion_bench_js_CFLAGS = $(torsion_cflags)
torsion_bench_js_LDADD = $(torsion_js_OBJECTS)
torsion_bench_js_LDFLAGS = $(torsion_ldflags)

noinst_PROGRAMS += torsion_test.js
torsion_test_js_SOURCES = $(test_sources)
torsion_test_js_CFLAGS = $(torsion_cflags)
torsion_test_js_LDADD = $(torsion_js_OBJECTS)
torsion_test_js_LDFLAGS = $(torsion_ldflags)
TESTS += torsion_test.js

bench: torsion_bench.js
	@node $^

test: torsion_test.js
	@node $^
endif ENABLE_TESTS
endif EMSCRIPTEN

#
# WASI
#

if WASI
INITIAL_MEMORY ?= 16777216
MAX_MEMORY ?= 2147483648
STACK_SIZE ?= 5242880
LOG_COMPILER = $(top_srcdir)/scripts/wasi-run
WASI_SYSROOT = $(WASI)/share/wasi-sysroot

torsion_ldflags = -Wl,--allow-undefined                  \
                  -Wl,--initial-memory=$(INITIAL_MEMORY) \
                  -Wl,--max-memory=$(MAX_MEMORY)         \
                  -Wl,-z -Wl,stack-size=$(STACK_SIZE)    \
                  -Wl,--stack-first

# Note: We have to manually choose crt1-reactor.o
# until reactor support is widely available via:
#
#   -mexec-model=reactor
#
# See: https://github.com/llvm/llvm-project/commit/d496437
#      https://reviews.llvm.org/D62922
noinst_PROGRAMS += torsion.wasm
torsion_wasm_SOURCES = $(torsion_sources)
torsion_wasm_CFLAGS = $(torsion_cflags) -DTORSION_BUILD
torsion_wasm_LDADD = $(WASI_SYSROOT)/lib/$(host_cpu)-wasi/crt1-reactor.o
torsion_wasm_LDFLAGS = $(torsion_ldflags)      \
                       -nostartfiles           \
                       -Wl,--entry=_initialize \
                       -Wl,--export-dynamic    \
                       -Wl,--export=malloc     \
                       -Wl,--export=free

if ENABLE_TESTS
noinst_PROGRAMS += torsion_bench.wasm
torsion_bench_wasm_SOURCES = $(bench_sources)
torsion_bench_wasm_CFLAGS = $(torsion_cflags)
torsion_bench_wasm_LDADD = $(torsion_wasm_OBJECTS)
torsion_bench_wasm_LDFLAGS = $(torsion_ldflags)

noinst_PROGRAMS += torsion_test.wasm
torsion_test_wasm_SOURCES = $(test_sources)
torsion_test_wasm_CFLAGS = $(torsion_cflags)
torsion_test_wasm_LDADD = $(torsion_wasm_OBJECTS)
torsion_test_wasm_LDFLAGS = $(torsion_ldflags)
TESTS += torsion_test.wasm

bench: torsion_bench.wasm
	@./scripts/wasi-run $^

test: torsion_test.wasm
	@./scripts/wasi-run $^
endif ENABLE_TESTS
endif WASI

#
# Native
#

if NATIVE
lib_LTLIBRARIES += libtorsion.la
libtorsion_la_SOURCES = $(torsion_sources)
libtorsion_la_CFLAGS = $(torsion_cflags) -DTORSION_BUILD
libtorsion_la_LDFLAGS = -no-undefined -version-info @LIB_VERSION@
libtorsion_la_LIBADD =

if ENABLE_PTHREAD
libtorsion_la_CFLAGS += @PTHREAD_CFLAGS@
libtorsion_la_LIBADD += @PTHREAD_LIBS@
endif

dist_doc_DATA = README.md
dist_license_DATA = LICENSE
pkgconfig_DATA = libtorsion.pc

# Add a `package` target, similar to CMake/CPack.
torsion_distdir = $(distdir)-bin
torsion_abs_distdir = $(abs_top_builddir)/$(torsion_distdir)
torsion_tarname = $(PACKAGE)-$(VERSION)-bin.tar.gz
torsion_tar_czf = $(am__tar) | eval GZIP= gzip $(GZIP_ENV) -c

package: all
	rm -rf $(torsion_tarname) $(torsion_distdir)
	$(MAKE) $(AM_MAKEFLAGS) install 'DESTDIR=$(torsion_abs_distdir)' prefix=/
	tardir=$(torsion_distdir) && $(torsion_tar_czf) > $(torsion_tarname)
	rm -rf $(torsion_distdir)

if ENABLE_TESTS
noinst_PROGRAMS += torsion_bench
torsion_bench_SOURCES = $(bench_sources)
torsion_bench_CFLAGS = $(torsion_cflags)
torsion_bench_LDFLAGS = -no-install
torsion_bench_LDADD = libtorsion.la

noinst_PROGRAMS += torsion_test
torsion_test_SOURCES = $(test_sources)
torsion_test_CFLAGS = $(torsion_cflags)
torsion_test_LDFLAGS = -no-install
torsion_test_LDADD = libtorsion.la
TESTS += torsion_test

if ENABLE_THREADS
torsion_test_SOURCES += test/thread.c
torsion_test_CFLAGS += @PTHREAD_CFLAGS@
torsion_test_LDADD += @PTHREAD_LIBS@
endif

if ENABLE_ZLIB
torsion_test_LDADD += -lz
endif

noinst_PROGRAMS += torsion_test_static
torsion_test_static_SOURCES = $(torsion_test_SOURCES)
torsion_test_static_CFLAGS = $(torsion_test_CFLAGS)
torsion_test_static_LDFLAGS = -static
torsion_test_static_LDADD = $(torsion_test_LDADD)
TESTS += torsion_test_static

if ENABLE_CTGRIND
noinst_PROGRAMS += torsion_ctgrind
torsion_ctgrind_SOURCES = test/ctgrind.c
torsion_ctgrind_CFLAGS = $(torsion_cflags)
torsion_ctgrind_LDFLAGS = -no-install
torsion_ctgrind_LDADD = libtorsion.la
endif

bench: torsion_bench
	@./$^

test: torsion_test
	@./$^

test_static: torsion_test_static
	@./$^

if ENABLE_CTGRIND
ctgrind: torsion_ctgrind
	@./libtool --mode=execute valgrind ./$^
endif
endif ENABLE_TESTS
endif NATIVE
