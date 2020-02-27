/*!
 * asn1.h - asn1 for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

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
               size_t *len, int strict) {
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

      if (*size >= (1ul << 24))
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
asn1_read_int(unsigned char *out, size_t out_len,
              const unsigned char **data, size_t *len, int strict) {
  size_t size;

  if (*len == 0 || **data != 0x02)
    return 0;

  *data += 1;
  *len -= 1;

  if (!asn1_read_size(&size, data, len, strict))
    return 0;

  /* Out of bounds. */
  if (size > *len)
    return 0;

  if (strict) {
    const unsigned char *num = *data;

    /* Zero-length integer. */
    if (size == 0)
      return 0;

    /* No negatives. */
    if (num[0] & 0x80)
      return 0;

    /* Allow zero only if it prefixes a high bit. */
    if (size > 1 && num[0] == 0x00) {
      if ((num[1] & 0x80) == 0x00)
        return 0;
    }
  }

  /* Eat leading zeroes. */
  while (size > 1 && **data == 0x00) {
    *data += 1;
    *len -= 1;
    size -= 1;
  }

  /* Invalid size. */
  if (size > out_len)
    return 0;

  memset(out, 0x00, out_len - size);
  memcpy(out + out_len - size, *data, size);

  *data += size;
  *len -= size;

  return 1;
}

static int
asn1_read_mpz(mpz_t n, const unsigned char **data, size_t *len, int strict) {
  size_t size;

  if (*len == 0 || **data != 0x02)
    return 0;

  *data += 1;
  *len -= 1;

  if (!asn1_read_size(&size, data, len, 1))
    return 0;

  /* Out of bounds. */
  if (size > *len)
    return 0;

  if (strict) {
    const unsigned char *num = *data;

    /* Zero-length integer. */
    if (size == 0)
      return 0;

    /* No negatives. */
    if (num[0] & 0x80)
      return 0;

    /* Allow zero only if it prefixes a high bit. */
    if (size > 1 && num[0] == 0x00) {
      if ((num[1] & 0x80) == 0x00)
        return 0;
    }
  }

  /* Eat leading zeroes. */
  while (size > 1 && **data == 0x00) {
    *data += 1;
    *len -= 1;
    size -= 1;
  }

  /* Invalid size. */
  if (size > 2048)
    return 0;

  mpz_import(n, size, 1, 1, 0, 0, *data);

  *data += size;
  *len -= size;

  return 1;
}

static int
asn1_read_version(const unsigned char **data, size_t *len,
                  unsigned char version, int strict) {
  int ret = 0;
  mpz_t n;

  if (strict) {
    const unsigned char *num = *data;

    if (*len < 3)
      return 0;

    if (num[0] != 0x02 || num[1] != 0x01 || num[2] != version)
      return 0;

    *data += 3;
    *len -= 3;

    return 1;
  }

  mpz_init(n);

  if (asn1_read_mpz(n, data, len, 0))
    ret = (mpz_cmp_ui(n, version) == 0);

  mpz_clear(n);

  return ret;
}

static int
asn1_read_dumb(mpz_t n, const unsigned char **data, size_t *len) {
  const unsigned char *buf = *data;
  size_t size;

  if (*len < 2)
    return 0;

  size = ((size_t)buf[0] << 8) | (size_t)buf[1];

  *data += 2;
  *len -= 2;

  if (size > *len)
    return 0;

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
asn1_size_int(const unsigned char *num, size_t len) {
  /* 0x02 [size] [0x00?] [int] */
  while (len > 1 && num[0] == 0x00) {
    len--;
    num++;
  }

  if (len == 0)
    return 3;

  len += num[0] >> 7;

  return 1 + asn1_size_size(len) + len;
}

static size_t
asn1_size_mpz(const mpz_t n) {
  /* 0x02 [size] [0x00?] [int] */
  size_t bits = mpz_bitlen(n);
  size_t size = (bits + 7) / 8;

  if (bits > 0 && (bits & 7) == 0)
    size += mpz_tstbit(n, bits - 1);

  if (bits == 0)
    size = 1;

  return 1 + asn1_size_size(size) + size;
}

static size_t
asn1_size_version(unsigned char version) {
  (void)version;
  return 3;
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
asn1_write_seq(unsigned char *data, size_t pos, size_t size) {
  data[pos++] = 0x30;
  return asn1_write_size(data, pos, size);
}

static size_t
asn1_write_int(unsigned char *data, size_t pos,
               const unsigned char *num, size_t len) {
  size_t pad = 0;

  /* 0x02 [size] [0x00?] [int] */
  while (len > 1 && num[0] == 0x00) {
    len--;
    num++;
  }

  if (len == 0) {
    data[pos++] = 0x02;
    data[pos++] = 0x01;
    data[pos++] = 0x00;
    return pos;
  }

  pad = num[0] >> 7;

  data[pos++] = 0x02;

  pos = asn1_write_size(data, pos, pad + len);

  if (pad)
    data[pos++] = 0x00;

  memcpy(data + pos, num, len);

  pos += len;

  return pos;
}

static size_t
asn1_write_mpz(unsigned char *data, size_t pos, const mpz_t n) {
  /* 0x02 [size] [0x00?] [int] */
  size_t bits = mpz_bitlen(n);
  size_t size = (bits + 7) / 8;
  size_t pad = 0;

  if (bits > 0 && (bits & 7) == 0)
    pad = mpz_tstbit(n, bits - 1);

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

static size_t
asn1_write_version(unsigned char *data, size_t pos, unsigned char version) {
  data[pos++] = 0x02;
  data[pos++] = 0x01;
  data[pos++] = version;
  return pos;
}

static size_t
asn1_write_dumb(unsigned char *data, size_t pos, const mpz_t n) {
  size_t size = mpz_bytelen(n);

  data[pos++] = size >> 8;
  data[pos++] = size;

  mpz_export(data + pos, NULL, 1, 1, 0, 0, n);

  pos += size;

  return pos;
}

#endif /* _TORSION_ASN1_H */
