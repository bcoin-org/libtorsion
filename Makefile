CFLAGS := \
	-g \
	-std=c89 \
	-pedantic \
	-Wall \
	-Wextra \
	-Wshadow \
	-Wno-unused-function \
	-Wno-implicit-fallthrough \
	-Wno-declaration-after-statement \
	-Wno-long-long \
	-Wno-overlength-strings \
	-I./include \
	-fPIC \
	-O3 \
	-DTORSION_TEST \
	-DTORSION_USE_64BIT \
	-DTORSION_USE_ASM

SOURCES := \
	src/aead.c \
	src/chacha20.c \
	src/ecc.c \
	src/drbg.c \
	src/dsa.c \
	src/hash.c \
	src/internal.c \
	src/kdf.c \
	src/mpi.c \
	src/poly1305.c \
	src/rsa.c \
	src/salsa20.c \
	src/siphash.c \
	src/util.c

LDFLAGS :=

OUTPUT = libtorsion.so

CLEANEXTS = o so

all: $(OUTPUT)

$(OUTPUT): $(subst .c,.o,$(SOURCES))
	$(CC) -shared -fPIC $(LDFLAGS) -o $@ $^

test: $(OUTPUT) src/test.c
	$(CC) $(CFLAGS) -L./ $(LDFLAGS) -ltorsion -o test src/test.c

clean:
	for ext in $(CLEANEXTS); do rm -f src/*.$$ext; done
	for ext in $(CLEANEXTS); do rm -f *.$$ext; done
	rm -f test-internal
	rm -f test

.PHONY: all clean
