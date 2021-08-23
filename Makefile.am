# Makefile.am - automake file for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

#
# Initialization
#

ACLOCAL_AMFLAGS = -I m4

includedir = $(prefix)/include/torsion
licensedir = $(datadir)/licenses/$(PACKAGE_NAME)
pkgconfigdir = $(libdir)/pkgconfig

dist_doc_DATA = README.md
dist_license_DATA = LICENSE
pkgconfig_DATA = libtorsion.pc

EXTRA_DIST = cmake/                 \
             etc/                   \
             scripts/               \
             autogen.sh             \
             CHANGELOG.md           \
             CMakeLists.txt         \
             libtorsion-cmake.pc.in \
             Makefile.dmc           \
             Makefile.fiat          \
             Makefile.msvc          \
             Makefile.unix

#
# Helpers
#

VALGRIND = ./libtool --mode=execute valgrind --leak-check=full  \
                                             --error-exitcode=1

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
                  include/torsion/mpi.h      \
                  include/torsion/rand.h     \
                  include/torsion/rsa.h      \
                  include/torsion/stream.h   \
                  include/torsion/util.h

noinst_HEADERS = src/asn1.h                      \
                 src/bf.h                        \
                 src/bio.h                       \
                 src/entropy/entropy.h           \
                 src/fields/p192_32.h            \
                 src/fields/p192_64.h            \
                 src/fields/p192.h               \
                 src/fields/p224_32.h            \
                 src/fields/p224_64.h            \
                 src/fields/p224.h               \
                 src/fields/p251_32.h            \
                 src/fields/p251_64.h            \
                 src/fields/p251.h               \
                 src/fields/p25519_32.h          \
                 src/fields/p25519_64.h          \
                 src/fields/p25519.h             \
                 src/fields/p256_32.h            \
                 src/fields/p256_64.h            \
                 src/fields/p256.h               \
                 src/fields/p384_32.h            \
                 src/fields/p384_64.h            \
                 src/fields/p384.h               \
                 src/fields/p448_32.h            \
                 src/fields/p448_64.h            \
                 src/fields/p448.h               \
                 src/fields/p521_32.h            \
                 src/fields/p521_64.h            \
                 src/fields/p521.h               \
                 src/fields/scalar.h             \
                 src/fields/secp256k1_32.h       \
                 src/fields/secp256k1_64.h       \
                 src/fields/secp256k1.h          \
                 src/internal.h                  \
                 src/subgroups.h                 \
                 test/data/chacha20_vectors.h    \
                 test/data/chachapoly_vectors.h  \
                 test/data/cipher_aead_vectors.h \
                 test/data/cipher_mode_vectors.h \
                 test/data/cipher_vectors.h      \
                 test/data/ctr_drbg_vectors.h    \
                 test/data/dsa_vectors.h         \
                 test/data/eb2k_vectors.h        \
                 test/data/ecdsa_vectors.h       \
                 test/data/eddsa_vectors.h       \
                 test/data/hash_drbg_vectors.h   \
                 test/data/hash_vectors.h        \
                 test/data/hkdf_vectors.h        \
                 test/data/hmac_drbg_vectors.h   \
                 test/data/hmac_vectors.h        \
                 test/data/jacobi_vectors.h      \
                 test/data/mpz_vectors.h         \
                 test/data/pbkdf2_vectors.h      \
                 test/data/poly1305_vectors.h    \
                 test/data/prime_vectors.h       \
                 test/data/rsa_vectors.h         \
                 test/data/bipschnorr_vectors.h  \
                 test/data/bip340_vectors.h      \
                 test/ecc_internal.h             \
                 test/mpi_internal.h             \
                 test/utils.h

#
# Sources & Flags
#

if ENABLE_RNG
rng_sources = src/entropy/hw.c  \
              src/entropy/sys.c \
              src/rand.c
endif

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
                  src/util.c     \
                  $(rng_sources)

torsion_cflags = -I$(top_srcdir)/include

if ENABLE_ZLIB
torsion_libs = -lz
endif

bench_sources = test/bench.c test/hrtime.c test/utils.c
test_sources = test/test.c test/utils.c
ctgrind_sources = test/ctgrind.c test/utils.c

#
# Native
#

if NATIVE
if MINGW
torsion_ldflags = -avoid-version -Wl,--output-def,.libs/libtorsion.def
else
torsion_ldflags = -version-info @ABI_VERSION@
endif

libtorsion_la_NAME = libtorsion.la
libtorsion_la_SOURCES = $(torsion_sources)
libtorsion_la_CFLAGS = $(torsion_cflags) -DTORSION_EXPORT
libtorsion_la_LDFLAGS = -no-undefined $(torsion_ldflags)

if ENABLE_TESTS
torsion_bench_NAME = torsion_bench
torsion_bench_SOURCES = $(bench_sources)
torsion_bench_CFLAGS = $(torsion_cflags)
torsion_bench_LDADD = libtorsion.la

torsion_test_NAME = torsion_test
torsion_test_SOURCES = $(test_sources)
torsion_test_CFLAGS = $(torsion_cflags)
torsion_test_LDADD = libtorsion.la $(torsion_libs)

torsion_test_static_NAME = torsion_test_static
torsion_test_static_SOURCES = $(torsion_test_SOURCES)
torsion_test_static_CFLAGS = $(torsion_test_CFLAGS)
torsion_test_static_LDFLAGS = -static
torsion_test_static_LDADD = $(torsion_test_LDADD)

if ENABLE_CTGRIND
torsion_ctgrind_NAME = torsion_ctgrind
torsion_ctgrind_SOURCES = $(ctgrind_sources)
torsion_ctgrind_CFLAGS = $(torsion_cflags)
torsion_ctgrind_LDADD = libtorsion.la
endif
endif ENABLE_TESTS

if MINGW
LOG_COMPILER = wine

libtorsion.lib: libtorsion.la
	$(AM_V_CCLD)
	@cp -f .libs/libtorsion.a $@

if ENABLE_SHARED
torsion.def: libtorsion.la
	$(AM_V_CCLD)
	@cp -f .libs/libtorsion.def $@

torsion.dll: libtorsion.la
	$(AM_V_CCLD)
	@cp -f .libs/libtorsion.dll $@

torsion.lib: torsion.def torsion.dll
	$(AM_V_CCLD)
	@$(DLLTOOL) -l $@ -d torsion.def -D torsion.dll
endif

if ENABLE_TESTS
bench: torsion_bench.exe
	@wine torsion_bench.exe
endif

if ENABLE_SHARED
all-local: libtorsion.lib torsion.lib
CLEANFILES = libtorsion.lib torsion.def torsion.dll torsion.lib
else
all-local: libtorsion.lib
CLEANFILES = libtorsion.lib
endif

else !MINGW

if ENABLE_TESTS
bench: torsion_bench
	@./torsion_bench

if ENABLE_CTGRIND
ctgrind: torsion_ctgrind
	@$(VALGRIND) ./torsion_ctgrind
endif

valgrind: torsion_test
	@$(VALGRIND) ./torsion_test
endif ENABLE_TESTS

endif !MINGW

lib_LTLIBRARIES = $(libtorsion_la_NAME)

noinst_PROGRAMS = $(torsion_bench_NAME)       \
                  $(torsion_test_NAME)        \
                  $(torsion_test_static_NAME) \
                  $(torsion_ctgrind_NAME)

TESTS = $(torsion_test_NAME) $(torsion_test_static_NAME)

endif NATIVE

#
# Static
#

if !NATIVE
libtorsion_a_NAME = libtorsion.a
libtorsion_a_SOURCES = $(torsion_sources)
libtorsion_a_CFLAGS = $(torsion_cflags)

lib_LIBRARIES = $(libtorsion_a_NAME)
endif !NATIVE

#
# WASI
#

if WASI
INITIAL_MEMORY = 16777216
MAX_MEMORY = 2147483648
STACK_SIZE = 5242880
LOG_COMPILER = $(top_srcdir)/scripts/wasi-run

torsion_ldflags = -Wl,--allow-undefined                  \
                  -Wl,--initial-memory=$(INITIAL_MEMORY) \
                  -Wl,--max-memory=$(MAX_MEMORY)         \
                  -Wl,-z,stack-size=$(STACK_SIZE)        \
                  -Wl,--stack-first

torsion_wasm_NAME = torsion.wasm
torsion_wasm_SOURCES = $(torsion_sources)
torsion_wasm_CFLAGS = $(torsion_cflags) -DTORSION_EXPORT
torsion_wasm_LDFLAGS = $(torsion_ldflags)   \
                       -mexec-model=reactor \
                       -Wl,--export-dynamic \
                       -Wl,--export=malloc  \
                       -Wl,--export=free

if ENABLE_TESTS
torsion_bench_wasm_NAME = torsion_bench.wasm
torsion_bench_wasm_SOURCES = $(bench_sources)
torsion_bench_wasm_CFLAGS = $(torsion_cflags)
torsion_bench_wasm_LDADD = libtorsion.a
torsion_bench_wasm_LDFLAGS = $(torsion_ldflags)

torsion_test_wasm_NAME = torsion_test.wasm
torsion_test_wasm_SOURCES = $(test_sources)
torsion_test_wasm_CFLAGS = $(torsion_cflags)
torsion_test_wasm_LDADD = libtorsion.a
torsion_test_wasm_LDFLAGS = $(torsion_ldflags)

bench: torsion_bench.wasm
	@./scripts/wasi-run torsion_bench.wasm
endif ENABLE_TESTS

bin_PROGRAMS = $(torsion_wasm_NAME)

noinst_PROGRAMS = $(torsion_bench_wasm_NAME) \
                  $(torsion_test_wasm_NAME)

TESTS = $(torsion_test_wasm_NAME)

endif WASI

#
# Emscripten
#

if EMSCRIPTEN
INITIAL_MEMORY = 16777216
MAX_MEMORY = 2147483648
STACK_SIZE = 5242880
ENVIRONMENT = node
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

torsion_js_NAME = torsion.js
torsion_js_SOURCES = $(torsion_sources)
torsion_js_CFLAGS = $(torsion_cflags) -DTORSION_EXPORT
torsion_js_LDFLAGS = $(torsion_ldflags) -s EXPORTED_FUNCTIONS=@etc/exports.json

if ENABLE_TESTS
torsion_bench_js_NAME = torsion_bench.js
torsion_bench_js_SOURCES = $(bench_sources)
torsion_bench_js_CFLAGS = $(torsion_cflags)
torsion_bench_js_LDADD = libtorsion.a
torsion_bench_js_LDFLAGS = $(torsion_ldflags)

torsion_test_js_NAME = torsion_test.js
torsion_test_js_SOURCES = $(test_sources)
torsion_test_js_CFLAGS = $(torsion_cflags)
torsion_test_js_LDADD = libtorsion.a
torsion_test_js_LDFLAGS = $(torsion_ldflags)

bench: torsion_bench.js
	@node torsion_bench.js
endif ENABLE_TESTS

bin_PROGRAMS = $(torsion_js_NAME)

noinst_PROGRAMS = $(torsion_bench_js_NAME) \
                  $(torsion_test_js_NAME)

TESTS = $(torsion_test_js_NAME)

endif EMSCRIPTEN
