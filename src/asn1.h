#ifndef _TORSION_ASN1_H
#define _TORSION_ASN1_H

#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include "mpz.h"

/*
 * ASN1
 */

static int
asn1_read_size(size_t *size,
               const unsigned char **data,
               size_t *len,
               int strict) {
  unsigned char ch;

  assert(sizeof(size_t) * CHAR_BIT >= 32);

  if (*len == 0)
    return 0;

  ch = **data;

  *data += 1;
  *len -= 1;

  if ((ch & 0x80) == 0) {
    /* Short form. */
    *size = ch;
  } else {
    size_t bytes = ch & 0x7f;
    size_t i;

    /* Indefinite form. */
    if (strict && bytes == 0)
      return 0;

    /* Long form. */
    *size = 0;

    for (i = 0; i < bytes; i++) {
      if (*len == 0)
        return 0;

      ch = **data;
      *data += 1;
      *len -= 1;

      if (*size >= (1 << 23))
        return 0;

      *size <<= 8;
      *size |= ch;

      if (strict && *size == 0)
        return 0;
    }

    if (strict && *size < 0x80)
      return 0;
  }

  return 1;
}

static int
asn1_read_seq(const unsigned char **data, size_t *len, int strict) {
  size_t size;

  if (*len == 0 || **data != 0x30)
    return 0;

  *data += 1;
  *len -= 1;

  if (!asn1_read_size(&size, data, len, strict))
    return 0;

  if (strict && size != *len)
    return 0;

  return 1;
}

static int
asn1_read_int(mpz_t n, const unsigned char **data, size_t *len, int strict) {
  size_t size;

  if (*len == 0 || **data != 0x02)
    return 0;

  *data += 1;
  *len -= 1;

  if (!asn1_read_size(&size, data, len, strict))
    return 0;

  if (size > *len)
    return 0;

  if (strict) {
    const unsigned char *num = *data;

    /* No reason to have an integer larger than this. */
    if (size == 0 || size > 1 + 2048)
      return 0;

    /* No negatives. */
    if (num[0] & 0x80)
      return 0;

    /* Allow zero only if it prefixes a high bit. */
    if (size > 1) {
      if (num[0] == 0x00 && (num[1] & 0x80) == 0x00)
        return 0;
    }
  }

  mpz_import(n, size, 1, 1, 0, 0, *data);

  *data += size;
  *len -= size;

  return 1;
}

static size_t
asn1_size_size(size_t size) {
  if (size <= 0x7f) /* [size] */
    return 1;

  if (size <= 0xff) /* 0x81 [size] */
    return 2;

  return 3; /* 0x82 [size-hi] [size-lo] */
}

static size_t
asn1_size_int(const mpz_t n) {
  /* 0x02 [size] [0x00?] [int] */
  size_t bits = mpz_bitlen(n);
  size_t size = (bits + 7) / 8;

  if ((bits & 7) == 0)
    size += (mpz_tstbit(n, bits - 1) != 0);

  if (bits == 0)
    size = 1;

  return 1 + asn1_size_size(size) + size;
}

static size_t
asn1_write_size(unsigned char *data, size_t pos, size_t size) {
  if (size <= 0x7f)  {
    /* [size] */
    data[pos++] = size;
  } else if (size <= 0xff) {
    /* 0x81 [size] */
    data[pos++] = 0x81;
    data[pos++] = size;
  } else {
    /* 0x82 [size-hi] [size-lo] */
    assert(size <= 0xffff);
    data[pos++] = 0x82;
    data[pos++] = size >> 8;
    data[pos++] = size & 0xff;
  }
  return pos;
}

static size_t
asn1_write_int(unsigned char *data, size_t pos, const mpz_t n) {
  /* 0x02 [size] [0x00?] [int] */
  size_t bits = mpz_bitlen(n);
  size_t size = (bits + 7) / 8;
  size_t pad = 0;

  if ((bits & 7) == 0)
    pad = (mpz_tstbit(n, bits - 1) != 0);

  if (bits == 0)
    size = 1;

  data[pos++] = 0x02;

  pos = asn1_write_size(data, pos, pad + size);

  if (pad)
    data[pos++] = 0x00;

  if (bits != 0)
    mpz_export(data + pos, NULL, 1, 1, 0, 0, n);
  else
    data[pos] = 0x00;

  pos += size;

  return pos;
}

#endif /* _TORSION_ASN1_H */
