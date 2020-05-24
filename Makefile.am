ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -I$(top_srcdir)/include

include_HEADERS = include/torsion/aead.h \
                  include/torsion/chacha20.h \
                  include/torsion/cipher.h \
                  include/torsion/drbg.h \
                  include/torsion/dsa.h \
                  include/torsion/ecc.h \
                  include/torsion/encoding.h \
                  include/torsion/hash.h \
                  include/torsion/kdf.h \
                  include/torsion/poly1305.h \
                  include/torsion/rsa.h \
                  include/torsion/salsa20.h \
                  include/torsion/secretbox.h \
                  include/torsion/siphash.h \
                  include/torsion/util.h

noinst_HEADERS = src/asn1.h \
                 src/bio.h \
                 src/fields/libsecp256k1_32.h \
                 src/fields/libsecp256k1_64.h \
                 src/fields/p192_32.h \
                 src/fields/p192_64.h \
                 src/fields/p192.h \
                 src/fields/p224_32.h \
                 src/fields/p224_64.h \
                 src/fields/p224.h \
                 src/fields/p251_32.h \
                 src/fields/p251_64.h \
                 src/fields/p251.h \
                 src/fields/p25519_32.h \
                 src/fields/p25519_64.h \
                 src/fields/p25519.h \
                 src/fields/p256_32.h \
                 src/fields/p256_64.h \
                 src/fields/p256.h \
                 src/fields/p384_32.h \
                 src/fields/p384_64.h \
                 src/fields/p384.h \
                 src/fields/p448_32.h \
                 src/fields/p448_64.h \
                 src/fields/p448.h \
                 src/fields/p521_32.h \
                 src/fields/p521_64.h \
                 src/fields/p521.h \
                 src/fields/scalar.h \
                 src/fields/secp256k1_32.h \
                 src/fields/secp256k1_64.h \
                 src/fields/secp256k1.h \
                 src/internal.h \
                 src/mpi.h \
                 src/subgroups.h

if ENABLE_RNG
include_HEADERS += include/torsion/rand.h
noinst_HEADERS += src/entropy/entropy.h
endif

CLEANFILES =

lib_LTLIBRARIES = libtorsion.la
libtorsion_la_CFLAGS = @CFLAGS@
libtorsion_la_LDFLAGS = -no-undefined -version-info 0:0:0
libtorsion_la_SOURCES = src/aead.c \
                        src/asn1.c \
                        src/chacha20.c \
                        src/cipher.c \
                        src/ecc.c \
                        src/encoding.c \
                        src/drbg.c \
                        src/dsa.c \
                        src/hash.c \
                        src/internal.c \
                        src/kdf.c \
                        src/mpi.c \
                        src/poly1305.c \
                        src/rsa.c \
                        src/salsa20.c \
                        src/secretbox.c \
                        src/siphash.c \
                        src/util.c

if ENABLE_RNG
libtorsion_la_SOURCES += src/entropy/env.c
libtorsion_la_SOURCES += src/entropy/hw.c
libtorsion_la_SOURCES += src/entropy/sys.c
libtorsion_la_SOURCES += src/rand.c
endif

EXTRA_DIST = LICENSE \
             README.md

if ENABLE_TESTS
TESTS = tests
noinst_PROGRAMS = tests
tests_CFLAGS = @CFLAGS@
tests_LDFLAGS = -static
tests_SOURCES = test/test.c
tests_LDADD = libtorsion.la
if ENABLE_CTIME
noinst_PROGRAMS += ctime-test
ctime_test_SOURCES = test/ctime-test.c
ctime_test_LDADD = libtorsion.la
endif
else
TESTS =
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libtorsion.pc
