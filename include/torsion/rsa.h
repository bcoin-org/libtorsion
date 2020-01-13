#ifndef _TORSION_RSA_H
#define _TORSION_RSA_H

#ifdef __cplusplus
extern "C" {
#endif

#define RSA_DEFAULT_MOD_BITS 2048
#define RSA_DEFAULT_EXP 65537
#define RSA_MIN_MOD_BITS 512
#define RSA_MAX_MOD_BITS 16384
#define RSA_MIN_MOD_BYTES ((RSA_MIN_MOD_BITS + 7) / 8)
#define RSA_MAX_MOD_BYTES ((RSA_MAX_MOD_BITS + 7) / 8)
#define RSA_MIN_EXP 3ull
#define RSA_MAX_EXP 0x1ffffffffull
#define RSA_MIN_EXP_BITS 2
#define RSA_MAX_EXP_BITS 33
#define RSA_MIN_EXP_BYTES 1
#define RSA_MAX_EXP_BYTES 5

/* Limits:
 * 4096 = 2614
 * 8192 = 5174
 * 16384 = 10294
 */

#define RSA_MAX_PRIV_SIZE (0                   \
  + 4 /* seq */                                \
  + 3 /* version */                            \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* n */          \
  + 2 + 1 + RSA_MAX_EXP_BYTES /* e */          \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* d */          \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* p */  \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* q */  \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* dp */ \
  + 4 + 1 + RSA_MAX_MOD_BYTES / 2 + 1 /* dq */ \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* qi */         \
)

/* Limits:
 * 4096 = 529
 * 8192 = 1041
 * 16384 = 2065
 */

#define RSA_MAX_PUB_SIZE (0           \
  + 4 /* seq */                       \
  + 4 + 1 + RSA_MAX_MOD_BYTES /* n */ \
  + 2 + 1 + RSA_MAX_EXP_BYTES /* e */ \
)

#define RSA_MAX_SIG_SIZE RSA_MAX_MOD_BYTES

#define RSA_SALT_LENGTH_AUTO 0
#define RSA_SALT_LENGTH_HASH -1

#ifdef __cplusplus
}
#endif

#endif /* _TORSION_RSA_H */
