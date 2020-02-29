CFLAGS := \
	-g \
	-pedantic \
	-Wall \
	-Wextra \
	-Wshadow \
	-Wno-unused-function \
	-Wno-unused-const-variable \
	-Wno-implicit-fallthrough \
	-Wno-declaration-after-statement \
	-Wno-long-long \
	-Wno-overlength-strings \
	-std=c89 \
	-I./include \
	-fPIC \
	-O3 \
	-DTORSION_TEST \
	-DTORSION_USE_64BIT \
	-DTORSION_USE_ASM \
	-DTORSION_USE_GMP

SOURCES := \
	src/aead.c \
	src/chacha20.c \
	src/ecc.c \
	src/drbg.c \
	src/dsa.c \
	src/hash.c \
	src/kdf.c \
	src/poly1305.c \
	src/rsa.c \
	src/salsa20.c \
	src/siphash.c \
	src/util.c

LDFLAGS := -lgmp

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
