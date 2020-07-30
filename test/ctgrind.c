/*!
 * ctgrind.c - constant-time valgrind test
 *
 * Idea taken from the libsecp256k1 tests:
 *   https://github.com/bitcoin-core/secp256k1/commit/3d23022
 *
 * Code rewritten from:
 *   https://github.com/bitcoin-core/secp256k1/blob/3d23022/src/valgrind_ctime_test.c
 *   Copyright (c) 2020 Gregory Maxwell (MIT License)
 *
 * Update (2020-06-04): It seems the original idea for this came from
 * Adam Langley. Curiously, the libsecp256k1 tests do not ascribe credit to
 * him. See: https://www.imperialviolet.org/2010/04/01/ctgrind.html
 *
 * Explanation:
 *
 *   Valgrind can tell us whether jumps depend on uninitialized data.
 *   If we make our private key uninitialized, valgrind can give us
 *   a warning if our code is leaking any secret data through variable
 *   time operations. A successful run of this program should have no
 *   error/warning output.
 *
 *   We have only one part of our ECC code which is variable time:
 *   the ECDSA signing loop. We have code present to explicitly
 *   exclude this.
 */

#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <valgrind/memcheck.h>
#include <torsion/ecc.h>
#include <torsion/util.h>
#include "utils.h"

static void
redefine(void *ptr, size_t size) {
  VALGRIND_MAKE_MEM_DEFINED(ptr, size);
}

static void
test_ecdsa(int curve_name) {
  wei_curve_t *ec = wei_curve_create(curve_name);
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char pub[33];
  unsigned char secret[33];
  size_t pub_len;
  size_t secret_len;
  size_t i;
  int ret;

  printf("Testing ECDSA (%d)...\n", curve_name);

  for (i = 0; i < 32; i++)
    priv[i] = i + 65;

  for (i = 0; i < 32; i++)
    msg[i] = i + 1;

  /* Public key generation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = ecdsa_pubkey_create(ec, pub, &pub_len, priv, 1);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Public key validation */
  VALGRIND_MAKE_MEM_DEFINED(pub, sizeof(pub));
  ASSERT(ecdsa_pubkey_verify(ec, pub, 33));

  /* Signing */
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = ecdsa_sign_internal(ec, sig, NULL, msg, 32, priv, redefine);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);
  VALGRIND_MAKE_MEM_DEFINED(msg, sizeof(msg));
  VALGRIND_MAKE_MEM_DEFINED(sig, sizeof(sig));
  ASSERT(ecdsa_verify(ec, msg, 32, sig, pub, 33));

  /* Schnorr Legacy */
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = schnorr_legacy_sign(ec, sig, msg, 32, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);
  VALGRIND_MAKE_MEM_DEFINED(msg, sizeof(msg));
  VALGRIND_MAKE_MEM_DEFINED(sig, sizeof(sig));
  ASSERT(schnorr_legacy_verify(ec, msg, 32, sig, pub, 33));

  /* ECDH */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = ecdsa_derive(ec, secret, &secret_len, pub, 33, priv, 1);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Secret validation */
  VALGRIND_MAKE_MEM_DEFINED(secret, sizeof(secret));
  ASSERT(ecdsa_pubkey_verify(ec, secret, 33));

  /* Validation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = ecdsa_privkey_verify(ec, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Negation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = ecdsa_privkey_negate(ec, priv, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Addition */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  ret = ecdsa_privkey_tweak_add(ec, priv, priv, msg);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Multiplication */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  ret = ecdsa_privkey_tweak_mul(ec, priv, priv, msg);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  wei_curve_destroy(ec);
}

static void
test_schnorr(void) {
  wei_curve_t *ec = wei_curve_create(WEI_CURVE_SECP256K1);
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char aux[32];
  unsigned char sig[64];
  unsigned char pub[32];
  unsigned char secret[32];
  size_t i;
  int ret;

  printf("Testing Schnorr...\n");

  for (i = 0; i < 32; i++)
    priv[i] = i + 65;

  for (i = 0; i < 32; i++) {
    msg[i] = i + 1;
    aux[i] = i + 2;
  }

  /* Public key generation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = schnorr_pubkey_create(ec, pub, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Public key validation */
  VALGRIND_MAKE_MEM_UNDEFINED(pub, 32);
  ret = schnorr_pubkey_verify(ec, pub);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  VALGRIND_MAKE_MEM_DEFINED(pub, sizeof(pub));
  ASSERT(ret);

  /* Signing */
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(aux, 32);
  ret = schnorr_sign(ec, sig, msg, 32, priv, aux);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);
  VALGRIND_MAKE_MEM_DEFINED(msg, sizeof(msg));
  VALGRIND_MAKE_MEM_DEFINED(sig, sizeof(sig));
  ASSERT(schnorr_verify(ec, msg, 32, sig, pub));

  /* ECDH */
  VALGRIND_MAKE_MEM_UNDEFINED(pub, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = schnorr_derive(ec, secret, pub, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Secret validation */
  VALGRIND_MAKE_MEM_UNDEFINED(secret, 32);
  ret = schnorr_pubkey_verify(ec, secret);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Validation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = schnorr_privkey_verify(ec, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Addition */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  ret = schnorr_privkey_tweak_add(ec, priv, priv, msg);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Multiplication */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  ret = schnorr_privkey_tweak_mul(ec, priv, priv, msg);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  wei_curve_destroy(ec);
}

static void
test_ecdh(void) {
  mont_curve_t *ec = mont_curve_create(MONT_CURVE_X25519);
  unsigned char priv[32];
  unsigned char pub[32];
  unsigned char secret[32];
  size_t i;
  int ret;

  printf("Testing ECDH...\n");

  for (i = 0; i < 32; i++)
    priv[i] = i + 65;

  /* Public key generation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ecdh_pubkey_create(ec, pub, priv);

  /* Public key validation */
  VALGRIND_MAKE_MEM_UNDEFINED(pub, 32);
  ret = ecdh_pubkey_verify(ec, pub);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* ECDH */
  VALGRIND_MAKE_MEM_UNDEFINED(pub, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = ecdh_derive(ec, secret, pub, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Secret validation */
  VALGRIND_MAKE_MEM_UNDEFINED(secret, 32);
  ret = ecdh_pubkey_verify(ec, secret);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  mont_curve_destroy(ec);
}

static void
test_eddsa(void) {
  edwards_curve_t *ec = edwards_curve_create(EDWARDS_CURVE_ED25519);
  unsigned char priv[32];
  unsigned char msg[32];
  unsigned char sig[64];
  unsigned char pub[32];
  unsigned char secret[32];
  size_t i;
  int ret;

  printf("Testing EdDSA...\n");

  for (i = 0; i < 32; i++)
    priv[i] = i + 65;

  for (i = 0; i < 32; i++)
    msg[i] = i + 1;

  /* Public key generation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  eddsa_pubkey_create(ec, pub, priv);

  /* Public key validation */
  VALGRIND_MAKE_MEM_UNDEFINED(pub, 32);
  ret = eddsa_pubkey_verify(ec, pub);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  VALGRIND_MAKE_MEM_DEFINED(pub, sizeof(pub));
  ASSERT(ret);

  /* Signing */
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  eddsa_sign(ec, sig, msg, 32, priv, 0, NULL, 0);
  VALGRIND_MAKE_MEM_DEFINED(msg, sizeof(msg));
  VALGRIND_MAKE_MEM_DEFINED(sig, sizeof(sig));
  ASSERT(eddsa_verify(ec, msg, 32, sig, pub, 0, NULL, 0));

  /* ECDH */
  VALGRIND_MAKE_MEM_UNDEFINED(pub, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  ret = eddsa_derive(ec, secret, pub, priv);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Secret validation */
  VALGRIND_MAKE_MEM_UNDEFINED(secret, 32);
  ret = eddsa_pubkey_verify(ec, secret);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* Negation */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  eddsa_scalar_negate(ec, priv, priv);

  /* Addition */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  eddsa_scalar_tweak_add(ec, priv, priv, msg);

  /* Multiplication */
  VALGRIND_MAKE_MEM_UNDEFINED(priv, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(msg, 32);
  eddsa_scalar_tweak_mul(ec, priv, priv, msg);

  edwards_curve_destroy(ec);
}

static void
test_util(void) {
  unsigned char x[32];
  unsigned char y[32];
  unsigned char z[32];
  size_t i;
  int ret;

  printf("Testing Util...\n");

  for (i = 0; i < 32; i++) {
    x[i] = i + 1;
    y[i] = i + 1;
    z[i] = i + 2;
  }

  /* x == x */
  VALGRIND_MAKE_MEM_UNDEFINED(x, 32);
  ret = torsion_memequal(x, x, 32);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* x == y */
  VALGRIND_MAKE_MEM_UNDEFINED(x, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(y, 32);
  ret = torsion_memequal(x, y, 32);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(ret);

  /* x != z */
  VALGRIND_MAKE_MEM_UNDEFINED(x, 32);
  VALGRIND_MAKE_MEM_UNDEFINED(z, 32);
  ret = torsion_memequal(x, z, 32);
  VALGRIND_MAKE_MEM_DEFINED(&ret, sizeof(ret));
  ASSERT(!ret);
}

int
main(void) {
  if (!RUNNING_ON_VALGRIND) {
    fprintf(stderr, "Test must be run with valgrind.\n");
    return 1;
  }

  test_ecdsa(WEI_CURVE_P256);
  test_ecdsa(WEI_CURVE_SECP256K1);
  test_schnorr();
  test_ecdh();
  test_eddsa();
  test_util();

  return 0;
}