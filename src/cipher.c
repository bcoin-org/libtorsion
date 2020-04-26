/*!
 * cipher.c - ciphers for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 *
 * Parts of this software are based on gnutls/nettle:
 *   Copyright (c) 1998-2019, Niels Möller and Contributors
 *   https://github.com/gnutls/nettle
 *
 * Parts of this software are based on openssl/openssl:
 *   Based on code entered into the public domain by Vincent Rijmen.
 *   https://github.com/openssl/openssl/blob/master/crypto/aes/aes_core.c
 *
 * Parts of this software are based on joyent/node-bcrypt-pbkdf:
 *   Copyright (c) 2016, Joyent Inc
 *   https://github.com/joyent/node-bcrypt-pbkdf
 *
 * Parts of this software are based on aead/camellia:
 *   Copyright (c) 2016, Andreas Auernhammer (MIT License).
 *   https://github.com/aead/camellia
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <torsion/cipher.h>
#include "bio.h"

/*
 * AES
 *
 * Resources:
 *   https://en.wikipedia.org/wiki/Advanced_Encryption_Standard
 *   http://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.197.pdf
 *   https://github.com/openssl/openssl/blob/master/crypto/aes/aes_core.c
 */

static const uint32_t aes_TE0[256] = {
  0xc66363a5, 0xf87c7c84, 0xee777799, 0xf67b7b8d,
  0xfff2f20d, 0xd66b6bbd, 0xde6f6fb1, 0x91c5c554,
  0x60303050, 0x02010103, 0xce6767a9, 0x562b2b7d,
  0xe7fefe19, 0xb5d7d762, 0x4dababe6, 0xec76769a,
  0x8fcaca45, 0x1f82829d, 0x89c9c940, 0xfa7d7d87,
  0xeffafa15, 0xb25959eb, 0x8e4747c9, 0xfbf0f00b,
  0x41adadec, 0xb3d4d467, 0x5fa2a2fd, 0x45afafea,
  0x239c9cbf, 0x53a4a4f7, 0xe4727296, 0x9bc0c05b,
  0x75b7b7c2, 0xe1fdfd1c, 0x3d9393ae, 0x4c26266a,
  0x6c36365a, 0x7e3f3f41, 0xf5f7f702, 0x83cccc4f,
  0x6834345c, 0x51a5a5f4, 0xd1e5e534, 0xf9f1f108,
  0xe2717193, 0xabd8d873, 0x62313153, 0x2a15153f,
  0x0804040c, 0x95c7c752, 0x46232365, 0x9dc3c35e,
  0x30181828, 0x379696a1, 0x0a05050f, 0x2f9a9ab5,
  0x0e070709, 0x24121236, 0x1b80809b, 0xdfe2e23d,
  0xcdebeb26, 0x4e272769, 0x7fb2b2cd, 0xea75759f,
  0x1209091b, 0x1d83839e, 0x582c2c74, 0x341a1a2e,
  0x361b1b2d, 0xdc6e6eb2, 0xb45a5aee, 0x5ba0a0fb,
  0xa45252f6, 0x763b3b4d, 0xb7d6d661, 0x7db3b3ce,
  0x5229297b, 0xdde3e33e, 0x5e2f2f71, 0x13848497,
  0xa65353f5, 0xb9d1d168, 0x00000000, 0xc1eded2c,
  0x40202060, 0xe3fcfc1f, 0x79b1b1c8, 0xb65b5bed,
  0xd46a6abe, 0x8dcbcb46, 0x67bebed9, 0x7239394b,
  0x944a4ade, 0x984c4cd4, 0xb05858e8, 0x85cfcf4a,
  0xbbd0d06b, 0xc5efef2a, 0x4faaaae5, 0xedfbfb16,
  0x864343c5, 0x9a4d4dd7, 0x66333355, 0x11858594,
  0x8a4545cf, 0xe9f9f910, 0x04020206, 0xfe7f7f81,
  0xa05050f0, 0x783c3c44, 0x259f9fba, 0x4ba8a8e3,
  0xa25151f3, 0x5da3a3fe, 0x804040c0, 0x058f8f8a,
  0x3f9292ad, 0x219d9dbc, 0x70383848, 0xf1f5f504,
  0x63bcbcdf, 0x77b6b6c1, 0xafdada75, 0x42212163,
  0x20101030, 0xe5ffff1a, 0xfdf3f30e, 0xbfd2d26d,
  0x81cdcd4c, 0x180c0c14, 0x26131335, 0xc3ecec2f,
  0xbe5f5fe1, 0x359797a2, 0x884444cc, 0x2e171739,
  0x93c4c457, 0x55a7a7f2, 0xfc7e7e82, 0x7a3d3d47,
  0xc86464ac, 0xba5d5de7, 0x3219192b, 0xe6737395,
  0xc06060a0, 0x19818198, 0x9e4f4fd1, 0xa3dcdc7f,
  0x44222266, 0x542a2a7e, 0x3b9090ab, 0x0b888883,
  0x8c4646ca, 0xc7eeee29, 0x6bb8b8d3, 0x2814143c,
  0xa7dede79, 0xbc5e5ee2, 0x160b0b1d, 0xaddbdb76,
  0xdbe0e03b, 0x64323256, 0x743a3a4e, 0x140a0a1e,
  0x924949db, 0x0c06060a, 0x4824246c, 0xb85c5ce4,
  0x9fc2c25d, 0xbdd3d36e, 0x43acacef, 0xc46262a6,
  0x399191a8, 0x319595a4, 0xd3e4e437, 0xf279798b,
  0xd5e7e732, 0x8bc8c843, 0x6e373759, 0xda6d6db7,
  0x018d8d8c, 0xb1d5d564, 0x9c4e4ed2, 0x49a9a9e0,
  0xd86c6cb4, 0xac5656fa, 0xf3f4f407, 0xcfeaea25,
  0xca6565af, 0xf47a7a8e, 0x47aeaee9, 0x10080818,
  0x6fbabad5, 0xf0787888, 0x4a25256f, 0x5c2e2e72,
  0x381c1c24, 0x57a6a6f1, 0x73b4b4c7, 0x97c6c651,
  0xcbe8e823, 0xa1dddd7c, 0xe874749c, 0x3e1f1f21,
  0x964b4bdd, 0x61bdbddc, 0x0d8b8b86, 0x0f8a8a85,
  0xe0707090, 0x7c3e3e42, 0x71b5b5c4, 0xcc6666aa,
  0x904848d8, 0x06030305, 0xf7f6f601, 0x1c0e0e12,
  0xc26161a3, 0x6a35355f, 0xae5757f9, 0x69b9b9d0,
  0x17868691, 0x99c1c158, 0x3a1d1d27, 0x279e9eb9,
  0xd9e1e138, 0xebf8f813, 0x2b9898b3, 0x22111133,
  0xd26969bb, 0xa9d9d970, 0x078e8e89, 0x339494a7,
  0x2d9b9bb6, 0x3c1e1e22, 0x15878792, 0xc9e9e920,
  0x87cece49, 0xaa5555ff, 0x50282878, 0xa5dfdf7a,
  0x038c8c8f, 0x59a1a1f8, 0x09898980, 0x1a0d0d17,
  0x65bfbfda, 0xd7e6e631, 0x844242c6, 0xd06868b8,
  0x824141c3, 0x299999b0, 0x5a2d2d77, 0x1e0f0f11,
  0x7bb0b0cb, 0xa85454fc, 0x6dbbbbd6, 0x2c16163a
};

static const uint32_t aes_TE1[256] = {
  0xa5c66363, 0x84f87c7c, 0x99ee7777, 0x8df67b7b,
  0x0dfff2f2, 0xbdd66b6b, 0xb1de6f6f, 0x5491c5c5,
  0x50603030, 0x03020101, 0xa9ce6767, 0x7d562b2b,
  0x19e7fefe, 0x62b5d7d7, 0xe64dabab, 0x9aec7676,
  0x458fcaca, 0x9d1f8282, 0x4089c9c9, 0x87fa7d7d,
  0x15effafa, 0xebb25959, 0xc98e4747, 0x0bfbf0f0,
  0xec41adad, 0x67b3d4d4, 0xfd5fa2a2, 0xea45afaf,
  0xbf239c9c, 0xf753a4a4, 0x96e47272, 0x5b9bc0c0,
  0xc275b7b7, 0x1ce1fdfd, 0xae3d9393, 0x6a4c2626,
  0x5a6c3636, 0x417e3f3f, 0x02f5f7f7, 0x4f83cccc,
  0x5c683434, 0xf451a5a5, 0x34d1e5e5, 0x08f9f1f1,
  0x93e27171, 0x73abd8d8, 0x53623131, 0x3f2a1515,
  0x0c080404, 0x5295c7c7, 0x65462323, 0x5e9dc3c3,
  0x28301818, 0xa1379696, 0x0f0a0505, 0xb52f9a9a,
  0x090e0707, 0x36241212, 0x9b1b8080, 0x3ddfe2e2,
  0x26cdebeb, 0x694e2727, 0xcd7fb2b2, 0x9fea7575,
  0x1b120909, 0x9e1d8383, 0x74582c2c, 0x2e341a1a,
  0x2d361b1b, 0xb2dc6e6e, 0xeeb45a5a, 0xfb5ba0a0,
  0xf6a45252, 0x4d763b3b, 0x61b7d6d6, 0xce7db3b3,
  0x7b522929, 0x3edde3e3, 0x715e2f2f, 0x97138484,
  0xf5a65353, 0x68b9d1d1, 0x00000000, 0x2cc1eded,
  0x60402020, 0x1fe3fcfc, 0xc879b1b1, 0xedb65b5b,
  0xbed46a6a, 0x468dcbcb, 0xd967bebe, 0x4b723939,
  0xde944a4a, 0xd4984c4c, 0xe8b05858, 0x4a85cfcf,
  0x6bbbd0d0, 0x2ac5efef, 0xe54faaaa, 0x16edfbfb,
  0xc5864343, 0xd79a4d4d, 0x55663333, 0x94118585,
  0xcf8a4545, 0x10e9f9f9, 0x06040202, 0x81fe7f7f,
  0xf0a05050, 0x44783c3c, 0xba259f9f, 0xe34ba8a8,
  0xf3a25151, 0xfe5da3a3, 0xc0804040, 0x8a058f8f,
  0xad3f9292, 0xbc219d9d, 0x48703838, 0x04f1f5f5,
  0xdf63bcbc, 0xc177b6b6, 0x75afdada, 0x63422121,
  0x30201010, 0x1ae5ffff, 0x0efdf3f3, 0x6dbfd2d2,
  0x4c81cdcd, 0x14180c0c, 0x35261313, 0x2fc3ecec,
  0xe1be5f5f, 0xa2359797, 0xcc884444, 0x392e1717,
  0x5793c4c4, 0xf255a7a7, 0x82fc7e7e, 0x477a3d3d,
  0xacc86464, 0xe7ba5d5d, 0x2b321919, 0x95e67373,
  0xa0c06060, 0x98198181, 0xd19e4f4f, 0x7fa3dcdc,
  0x66442222, 0x7e542a2a, 0xab3b9090, 0x830b8888,
  0xca8c4646, 0x29c7eeee, 0xd36bb8b8, 0x3c281414,
  0x79a7dede, 0xe2bc5e5e, 0x1d160b0b, 0x76addbdb,
  0x3bdbe0e0, 0x56643232, 0x4e743a3a, 0x1e140a0a,
  0xdb924949, 0x0a0c0606, 0x6c482424, 0xe4b85c5c,
  0x5d9fc2c2, 0x6ebdd3d3, 0xef43acac, 0xa6c46262,
  0xa8399191, 0xa4319595, 0x37d3e4e4, 0x8bf27979,
  0x32d5e7e7, 0x438bc8c8, 0x596e3737, 0xb7da6d6d,
  0x8c018d8d, 0x64b1d5d5, 0xd29c4e4e, 0xe049a9a9,
  0xb4d86c6c, 0xfaac5656, 0x07f3f4f4, 0x25cfeaea,
  0xafca6565, 0x8ef47a7a, 0xe947aeae, 0x18100808,
  0xd56fbaba, 0x88f07878, 0x6f4a2525, 0x725c2e2e,
  0x24381c1c, 0xf157a6a6, 0xc773b4b4, 0x5197c6c6,
  0x23cbe8e8, 0x7ca1dddd, 0x9ce87474, 0x213e1f1f,
  0xdd964b4b, 0xdc61bdbd, 0x860d8b8b, 0x850f8a8a,
  0x90e07070, 0x427c3e3e, 0xc471b5b5, 0xaacc6666,
  0xd8904848, 0x05060303, 0x01f7f6f6, 0x121c0e0e,
  0xa3c26161, 0x5f6a3535, 0xf9ae5757, 0xd069b9b9,
  0x91178686, 0x5899c1c1, 0x273a1d1d, 0xb9279e9e,
  0x38d9e1e1, 0x13ebf8f8, 0xb32b9898, 0x33221111,
  0xbbd26969, 0x70a9d9d9, 0x89078e8e, 0xa7339494,
  0xb62d9b9b, 0x223c1e1e, 0x92158787, 0x20c9e9e9,
  0x4987cece, 0xffaa5555, 0x78502828, 0x7aa5dfdf,
  0x8f038c8c, 0xf859a1a1, 0x80098989, 0x171a0d0d,
  0xda65bfbf, 0x31d7e6e6, 0xc6844242, 0xb8d06868,
  0xc3824141, 0xb0299999, 0x775a2d2d, 0x111e0f0f,
  0xcb7bb0b0, 0xfca85454, 0xd66dbbbb, 0x3a2c1616
};

static const uint32_t aes_TE2[256] = {
  0x63a5c663, 0x7c84f87c, 0x7799ee77, 0x7b8df67b,
  0xf20dfff2, 0x6bbdd66b, 0x6fb1de6f, 0xc55491c5,
  0x30506030, 0x01030201, 0x67a9ce67, 0x2b7d562b,
  0xfe19e7fe, 0xd762b5d7, 0xabe64dab, 0x769aec76,
  0xca458fca, 0x829d1f82, 0xc94089c9, 0x7d87fa7d,
  0xfa15effa, 0x59ebb259, 0x47c98e47, 0xf00bfbf0,
  0xadec41ad, 0xd467b3d4, 0xa2fd5fa2, 0xafea45af,
  0x9cbf239c, 0xa4f753a4, 0x7296e472, 0xc05b9bc0,
  0xb7c275b7, 0xfd1ce1fd, 0x93ae3d93, 0x266a4c26,
  0x365a6c36, 0x3f417e3f, 0xf702f5f7, 0xcc4f83cc,
  0x345c6834, 0xa5f451a5, 0xe534d1e5, 0xf108f9f1,
  0x7193e271, 0xd873abd8, 0x31536231, 0x153f2a15,
  0x040c0804, 0xc75295c7, 0x23654623, 0xc35e9dc3,
  0x18283018, 0x96a13796, 0x050f0a05, 0x9ab52f9a,
  0x07090e07, 0x12362412, 0x809b1b80, 0xe23ddfe2,
  0xeb26cdeb, 0x27694e27, 0xb2cd7fb2, 0x759fea75,
  0x091b1209, 0x839e1d83, 0x2c74582c, 0x1a2e341a,
  0x1b2d361b, 0x6eb2dc6e, 0x5aeeb45a, 0xa0fb5ba0,
  0x52f6a452, 0x3b4d763b, 0xd661b7d6, 0xb3ce7db3,
  0x297b5229, 0xe33edde3, 0x2f715e2f, 0x84971384,
  0x53f5a653, 0xd168b9d1, 0x00000000, 0xed2cc1ed,
  0x20604020, 0xfc1fe3fc, 0xb1c879b1, 0x5bedb65b,
  0x6abed46a, 0xcb468dcb, 0xbed967be, 0x394b7239,
  0x4ade944a, 0x4cd4984c, 0x58e8b058, 0xcf4a85cf,
  0xd06bbbd0, 0xef2ac5ef, 0xaae54faa, 0xfb16edfb,
  0x43c58643, 0x4dd79a4d, 0x33556633, 0x85941185,
  0x45cf8a45, 0xf910e9f9, 0x02060402, 0x7f81fe7f,
  0x50f0a050, 0x3c44783c, 0x9fba259f, 0xa8e34ba8,
  0x51f3a251, 0xa3fe5da3, 0x40c08040, 0x8f8a058f,
  0x92ad3f92, 0x9dbc219d, 0x38487038, 0xf504f1f5,
  0xbcdf63bc, 0xb6c177b6, 0xda75afda, 0x21634221,
  0x10302010, 0xff1ae5ff, 0xf30efdf3, 0xd26dbfd2,
  0xcd4c81cd, 0x0c14180c, 0x13352613, 0xec2fc3ec,
  0x5fe1be5f, 0x97a23597, 0x44cc8844, 0x17392e17,
  0xc45793c4, 0xa7f255a7, 0x7e82fc7e, 0x3d477a3d,
  0x64acc864, 0x5de7ba5d, 0x192b3219, 0x7395e673,
  0x60a0c060, 0x81981981, 0x4fd19e4f, 0xdc7fa3dc,
  0x22664422, 0x2a7e542a, 0x90ab3b90, 0x88830b88,
  0x46ca8c46, 0xee29c7ee, 0xb8d36bb8, 0x143c2814,
  0xde79a7de, 0x5ee2bc5e, 0x0b1d160b, 0xdb76addb,
  0xe03bdbe0, 0x32566432, 0x3a4e743a, 0x0a1e140a,
  0x49db9249, 0x060a0c06, 0x246c4824, 0x5ce4b85c,
  0xc25d9fc2, 0xd36ebdd3, 0xacef43ac, 0x62a6c462,
  0x91a83991, 0x95a43195, 0xe437d3e4, 0x798bf279,
  0xe732d5e7, 0xc8438bc8, 0x37596e37, 0x6db7da6d,
  0x8d8c018d, 0xd564b1d5, 0x4ed29c4e, 0xa9e049a9,
  0x6cb4d86c, 0x56faac56, 0xf407f3f4, 0xea25cfea,
  0x65afca65, 0x7a8ef47a, 0xaee947ae, 0x08181008,
  0xbad56fba, 0x7888f078, 0x256f4a25, 0x2e725c2e,
  0x1c24381c, 0xa6f157a6, 0xb4c773b4, 0xc65197c6,
  0xe823cbe8, 0xdd7ca1dd, 0x749ce874, 0x1f213e1f,
  0x4bdd964b, 0xbddc61bd, 0x8b860d8b, 0x8a850f8a,
  0x7090e070, 0x3e427c3e, 0xb5c471b5, 0x66aacc66,
  0x48d89048, 0x03050603, 0xf601f7f6, 0x0e121c0e,
  0x61a3c261, 0x355f6a35, 0x57f9ae57, 0xb9d069b9,
  0x86911786, 0xc15899c1, 0x1d273a1d, 0x9eb9279e,
  0xe138d9e1, 0xf813ebf8, 0x98b32b98, 0x11332211,
  0x69bbd269, 0xd970a9d9, 0x8e89078e, 0x94a73394,
  0x9bb62d9b, 0x1e223c1e, 0x87921587, 0xe920c9e9,
  0xce4987ce, 0x55ffaa55, 0x28785028, 0xdf7aa5df,
  0x8c8f038c, 0xa1f859a1, 0x89800989, 0x0d171a0d,
  0xbfda65bf, 0xe631d7e6, 0x42c68442, 0x68b8d068,
  0x41c38241, 0x99b02999, 0x2d775a2d, 0x0f111e0f,
  0xb0cb7bb0, 0x54fca854, 0xbbd66dbb, 0x163a2c16
};

static const uint32_t aes_TE3[256] = {
  0x6363a5c6, 0x7c7c84f8, 0x777799ee, 0x7b7b8df6,
  0xf2f20dff, 0x6b6bbdd6, 0x6f6fb1de, 0xc5c55491,
  0x30305060, 0x01010302, 0x6767a9ce, 0x2b2b7d56,
  0xfefe19e7, 0xd7d762b5, 0xababe64d, 0x76769aec,
  0xcaca458f, 0x82829d1f, 0xc9c94089, 0x7d7d87fa,
  0xfafa15ef, 0x5959ebb2, 0x4747c98e, 0xf0f00bfb,
  0xadadec41, 0xd4d467b3, 0xa2a2fd5f, 0xafafea45,
  0x9c9cbf23, 0xa4a4f753, 0x727296e4, 0xc0c05b9b,
  0xb7b7c275, 0xfdfd1ce1, 0x9393ae3d, 0x26266a4c,
  0x36365a6c, 0x3f3f417e, 0xf7f702f5, 0xcccc4f83,
  0x34345c68, 0xa5a5f451, 0xe5e534d1, 0xf1f108f9,
  0x717193e2, 0xd8d873ab, 0x31315362, 0x15153f2a,
  0x04040c08, 0xc7c75295, 0x23236546, 0xc3c35e9d,
  0x18182830, 0x9696a137, 0x05050f0a, 0x9a9ab52f,
  0x0707090e, 0x12123624, 0x80809b1b, 0xe2e23ddf,
  0xebeb26cd, 0x2727694e, 0xb2b2cd7f, 0x75759fea,
  0x09091b12, 0x83839e1d, 0x2c2c7458, 0x1a1a2e34,
  0x1b1b2d36, 0x6e6eb2dc, 0x5a5aeeb4, 0xa0a0fb5b,
  0x5252f6a4, 0x3b3b4d76, 0xd6d661b7, 0xb3b3ce7d,
  0x29297b52, 0xe3e33edd, 0x2f2f715e, 0x84849713,
  0x5353f5a6, 0xd1d168b9, 0x00000000, 0xeded2cc1,
  0x20206040, 0xfcfc1fe3, 0xb1b1c879, 0x5b5bedb6,
  0x6a6abed4, 0xcbcb468d, 0xbebed967, 0x39394b72,
  0x4a4ade94, 0x4c4cd498, 0x5858e8b0, 0xcfcf4a85,
  0xd0d06bbb, 0xefef2ac5, 0xaaaae54f, 0xfbfb16ed,
  0x4343c586, 0x4d4dd79a, 0x33335566, 0x85859411,
  0x4545cf8a, 0xf9f910e9, 0x02020604, 0x7f7f81fe,
  0x5050f0a0, 0x3c3c4478, 0x9f9fba25, 0xa8a8e34b,
  0x5151f3a2, 0xa3a3fe5d, 0x4040c080, 0x8f8f8a05,
  0x9292ad3f, 0x9d9dbc21, 0x38384870, 0xf5f504f1,
  0xbcbcdf63, 0xb6b6c177, 0xdada75af, 0x21216342,
  0x10103020, 0xffff1ae5, 0xf3f30efd, 0xd2d26dbf,
  0xcdcd4c81, 0x0c0c1418, 0x13133526, 0xecec2fc3,
  0x5f5fe1be, 0x9797a235, 0x4444cc88, 0x1717392e,
  0xc4c45793, 0xa7a7f255, 0x7e7e82fc, 0x3d3d477a,
  0x6464acc8, 0x5d5de7ba, 0x19192b32, 0x737395e6,
  0x6060a0c0, 0x81819819, 0x4f4fd19e, 0xdcdc7fa3,
  0x22226644, 0x2a2a7e54, 0x9090ab3b, 0x8888830b,
  0x4646ca8c, 0xeeee29c7, 0xb8b8d36b, 0x14143c28,
  0xdede79a7, 0x5e5ee2bc, 0x0b0b1d16, 0xdbdb76ad,
  0xe0e03bdb, 0x32325664, 0x3a3a4e74, 0x0a0a1e14,
  0x4949db92, 0x06060a0c, 0x24246c48, 0x5c5ce4b8,
  0xc2c25d9f, 0xd3d36ebd, 0xacacef43, 0x6262a6c4,
  0x9191a839, 0x9595a431, 0xe4e437d3, 0x79798bf2,
  0xe7e732d5, 0xc8c8438b, 0x3737596e, 0x6d6db7da,
  0x8d8d8c01, 0xd5d564b1, 0x4e4ed29c, 0xa9a9e049,
  0x6c6cb4d8, 0x5656faac, 0xf4f407f3, 0xeaea25cf,
  0x6565afca, 0x7a7a8ef4, 0xaeaee947, 0x08081810,
  0xbabad56f, 0x787888f0, 0x25256f4a, 0x2e2e725c,
  0x1c1c2438, 0xa6a6f157, 0xb4b4c773, 0xc6c65197,
  0xe8e823cb, 0xdddd7ca1, 0x74749ce8, 0x1f1f213e,
  0x4b4bdd96, 0xbdbddc61, 0x8b8b860d, 0x8a8a850f,
  0x707090e0, 0x3e3e427c, 0xb5b5c471, 0x6666aacc,
  0x4848d890, 0x03030506, 0xf6f601f7, 0x0e0e121c,
  0x6161a3c2, 0x35355f6a, 0x5757f9ae, 0xb9b9d069,
  0x86869117, 0xc1c15899, 0x1d1d273a, 0x9e9eb927,
  0xe1e138d9, 0xf8f813eb, 0x9898b32b, 0x11113322,
  0x6969bbd2, 0xd9d970a9, 0x8e8e8907, 0x9494a733,
  0x9b9bb62d, 0x1e1e223c, 0x87879215, 0xe9e920c9,
  0xcece4987, 0x5555ffaa, 0x28287850, 0xdfdf7aa5,
  0x8c8c8f03, 0xa1a1f859, 0x89898009, 0x0d0d171a,
  0xbfbfda65, 0xe6e631d7, 0x4242c684, 0x6868b8d0,
  0x4141c382, 0x9999b029, 0x2d2d775a, 0x0f0f111e,
  0xb0b0cb7b, 0x5454fca8, 0xbbbbd66d, 0x16163a2c
};

static const uint32_t aes_TD0[256] = {
  0x51f4a750, 0x7e416553, 0x1a17a4c3, 0x3a275e96,
  0x3bab6bcb, 0x1f9d45f1, 0xacfa58ab, 0x4be30393,
  0x2030fa55, 0xad766df6, 0x88cc7691, 0xf5024c25,
  0x4fe5d7fc, 0xc52acbd7, 0x26354480, 0xb562a38f,
  0xdeb15a49, 0x25ba1b67, 0x45ea0e98, 0x5dfec0e1,
  0xc32f7502, 0x814cf012, 0x8d4697a3, 0x6bd3f9c6,
  0x038f5fe7, 0x15929c95, 0xbf6d7aeb, 0x955259da,
  0xd4be832d, 0x587421d3, 0x49e06929, 0x8ec9c844,
  0x75c2896a, 0xf48e7978, 0x99583e6b, 0x27b971dd,
  0xbee14fb6, 0xf088ad17, 0xc920ac66, 0x7dce3ab4,
  0x63df4a18, 0xe51a3182, 0x97513360, 0x62537f45,
  0xb16477e0, 0xbb6bae84, 0xfe81a01c, 0xf9082b94,
  0x70486858, 0x8f45fd19, 0x94de6c87, 0x527bf8b7,
  0xab73d323, 0x724b02e2, 0xe31f8f57, 0x6655ab2a,
  0xb2eb2807, 0x2fb5c203, 0x86c57b9a, 0xd33708a5,
  0x302887f2, 0x23bfa5b2, 0x02036aba, 0xed16825c,
  0x8acf1c2b, 0xa779b492, 0xf307f2f0, 0x4e69e2a1,
  0x65daf4cd, 0x0605bed5, 0xd134621f, 0xc4a6fe8a,
  0x342e539d, 0xa2f355a0, 0x058ae132, 0xa4f6eb75,
  0x0b83ec39, 0x4060efaa, 0x5e719f06, 0xbd6e1051,
  0x3e218af9, 0x96dd063d, 0xdd3e05ae, 0x4de6bd46,
  0x91548db5, 0x71c45d05, 0x0406d46f, 0x605015ff,
  0x1998fb24, 0xd6bde997, 0x894043cc, 0x67d99e77,
  0xb0e842bd, 0x07898b88, 0xe7195b38, 0x79c8eedb,
  0xa17c0a47, 0x7c420fe9, 0xf8841ec9, 0x00000000,
  0x09808683, 0x322bed48, 0x1e1170ac, 0x6c5a724e,
  0xfd0efffb, 0x0f853856, 0x3daed51e, 0x362d3927,
  0x0a0fd964, 0x685ca621, 0x9b5b54d1, 0x24362e3a,
  0x0c0a67b1, 0x9357e70f, 0xb4ee96d2, 0x1b9b919e,
  0x80c0c54f, 0x61dc20a2, 0x5a774b69, 0x1c121a16,
  0xe293ba0a, 0xc0a02ae5, 0x3c22e043, 0x121b171d,
  0x0e090d0b, 0xf28bc7ad, 0x2db6a8b9, 0x141ea9c8,
  0x57f11985, 0xaf75074c, 0xee99ddbb, 0xa37f60fd,
  0xf701269f, 0x5c72f5bc, 0x44663bc5, 0x5bfb7e34,
  0x8b432976, 0xcb23c6dc, 0xb6edfc68, 0xb8e4f163,
  0xd731dcca, 0x42638510, 0x13972240, 0x84c61120,
  0x854a247d, 0xd2bb3df8, 0xaef93211, 0xc729a16d,
  0x1d9e2f4b, 0xdcb230f3, 0x0d8652ec, 0x77c1e3d0,
  0x2bb3166c, 0xa970b999, 0x119448fa, 0x47e96422,
  0xa8fc8cc4, 0xa0f03f1a, 0x567d2cd8, 0x223390ef,
  0x87494ec7, 0xd938d1c1, 0x8ccaa2fe, 0x98d40b36,
  0xa6f581cf, 0xa57ade28, 0xdab78e26, 0x3fadbfa4,
  0x2c3a9de4, 0x5078920d, 0x6a5fcc9b, 0x547e4662,
  0xf68d13c2, 0x90d8b8e8, 0x2e39f75e, 0x82c3aff5,
  0x9f5d80be, 0x69d0937c, 0x6fd52da9, 0xcf2512b3,
  0xc8ac993b, 0x10187da7, 0xe89c636e, 0xdb3bbb7b,
  0xcd267809, 0x6e5918f4, 0xec9ab701, 0x834f9aa8,
  0xe6956e65, 0xaaffe67e, 0x21bccf08, 0xef15e8e6,
  0xbae79bd9, 0x4a6f36ce, 0xea9f09d4, 0x29b07cd6,
  0x31a4b2af, 0x2a3f2331, 0xc6a59430, 0x35a266c0,
  0x744ebc37, 0xfc82caa6, 0xe090d0b0, 0x33a7d815,
  0xf104984a, 0x41ecdaf7, 0x7fcd500e, 0x1791f62f,
  0x764dd68d, 0x43efb04d, 0xccaa4d54, 0xe49604df,
  0x9ed1b5e3, 0x4c6a881b, 0xc12c1fb8, 0x4665517f,
  0x9d5eea04, 0x018c355d, 0xfa877473, 0xfb0b412e,
  0xb3671d5a, 0x92dbd252, 0xe9105633, 0x6dd64713,
  0x9ad7618c, 0x37a10c7a, 0x59f8148e, 0xeb133c89,
  0xcea927ee, 0xb761c935, 0xe11ce5ed, 0x7a47b13c,
  0x9cd2df59, 0x55f2733f, 0x1814ce79, 0x73c737bf,
  0x53f7cdea, 0x5ffdaa5b, 0xdf3d6f14, 0x7844db86,
  0xcaaff381, 0xb968c43e, 0x3824342c, 0xc2a3405f,
  0x161dc372, 0xbce2250c, 0x283c498b, 0xff0d9541,
  0x39a80171, 0x080cb3de, 0xd8b4e49c, 0x6456c190,
  0x7bcb8461, 0xd532b670, 0x486c5c74, 0xd0b85742
};

static const uint32_t aes_TD1[256] = {
  0x5051f4a7, 0x537e4165, 0xc31a17a4, 0x963a275e,
  0xcb3bab6b, 0xf11f9d45, 0xabacfa58, 0x934be303,
  0x552030fa, 0xf6ad766d, 0x9188cc76, 0x25f5024c,
  0xfc4fe5d7, 0xd7c52acb, 0x80263544, 0x8fb562a3,
  0x49deb15a, 0x6725ba1b, 0x9845ea0e, 0xe15dfec0,
  0x02c32f75, 0x12814cf0, 0xa38d4697, 0xc66bd3f9,
  0xe7038f5f, 0x9515929c, 0xebbf6d7a, 0xda955259,
  0x2dd4be83, 0xd3587421, 0x2949e069, 0x448ec9c8,
  0x6a75c289, 0x78f48e79, 0x6b99583e, 0xdd27b971,
  0xb6bee14f, 0x17f088ad, 0x66c920ac, 0xb47dce3a,
  0x1863df4a, 0x82e51a31, 0x60975133, 0x4562537f,
  0xe0b16477, 0x84bb6bae, 0x1cfe81a0, 0x94f9082b,
  0x58704868, 0x198f45fd, 0x8794de6c, 0xb7527bf8,
  0x23ab73d3, 0xe2724b02, 0x57e31f8f, 0x2a6655ab,
  0x07b2eb28, 0x032fb5c2, 0x9a86c57b, 0xa5d33708,
  0xf2302887, 0xb223bfa5, 0xba02036a, 0x5ced1682,
  0x2b8acf1c, 0x92a779b4, 0xf0f307f2, 0xa14e69e2,
  0xcd65daf4, 0xd50605be, 0x1fd13462, 0x8ac4a6fe,
  0x9d342e53, 0xa0a2f355, 0x32058ae1, 0x75a4f6eb,
  0x390b83ec, 0xaa4060ef, 0x065e719f, 0x51bd6e10,
  0xf93e218a, 0x3d96dd06, 0xaedd3e05, 0x464de6bd,
  0xb591548d, 0x0571c45d, 0x6f0406d4, 0xff605015,
  0x241998fb, 0x97d6bde9, 0xcc894043, 0x7767d99e,
  0xbdb0e842, 0x8807898b, 0x38e7195b, 0xdb79c8ee,
  0x47a17c0a, 0xe97c420f, 0xc9f8841e, 0x00000000,
  0x83098086, 0x48322bed, 0xac1e1170, 0x4e6c5a72,
  0xfbfd0eff, 0x560f8538, 0x1e3daed5, 0x27362d39,
  0x640a0fd9, 0x21685ca6, 0xd19b5b54, 0x3a24362e,
  0xb10c0a67, 0x0f9357e7, 0xd2b4ee96, 0x9e1b9b91,
  0x4f80c0c5, 0xa261dc20, 0x695a774b, 0x161c121a,
  0x0ae293ba, 0xe5c0a02a, 0x433c22e0, 0x1d121b17,
  0x0b0e090d, 0xadf28bc7, 0xb92db6a8, 0xc8141ea9,
  0x8557f119, 0x4caf7507, 0xbbee99dd, 0xfda37f60,
  0x9ff70126, 0xbc5c72f5, 0xc544663b, 0x345bfb7e,
  0x768b4329, 0xdccb23c6, 0x68b6edfc, 0x63b8e4f1,
  0xcad731dc, 0x10426385, 0x40139722, 0x2084c611,
  0x7d854a24, 0xf8d2bb3d, 0x11aef932, 0x6dc729a1,
  0x4b1d9e2f, 0xf3dcb230, 0xec0d8652, 0xd077c1e3,
  0x6c2bb316, 0x99a970b9, 0xfa119448, 0x2247e964,
  0xc4a8fc8c, 0x1aa0f03f, 0xd8567d2c, 0xef223390,
  0xc787494e, 0xc1d938d1, 0xfe8ccaa2, 0x3698d40b,
  0xcfa6f581, 0x28a57ade, 0x26dab78e, 0xa43fadbf,
  0xe42c3a9d, 0x0d507892, 0x9b6a5fcc, 0x62547e46,
  0xc2f68d13, 0xe890d8b8, 0x5e2e39f7, 0xf582c3af,
  0xbe9f5d80, 0x7c69d093, 0xa96fd52d, 0xb3cf2512,
  0x3bc8ac99, 0xa710187d, 0x6ee89c63, 0x7bdb3bbb,
  0x09cd2678, 0xf46e5918, 0x01ec9ab7, 0xa8834f9a,
  0x65e6956e, 0x7eaaffe6, 0x0821bccf, 0xe6ef15e8,
  0xd9bae79b, 0xce4a6f36, 0xd4ea9f09, 0xd629b07c,
  0xaf31a4b2, 0x312a3f23, 0x30c6a594, 0xc035a266,
  0x37744ebc, 0xa6fc82ca, 0xb0e090d0, 0x1533a7d8,
  0x4af10498, 0xf741ecda, 0x0e7fcd50, 0x2f1791f6,
  0x8d764dd6, 0x4d43efb0, 0x54ccaa4d, 0xdfe49604,
  0xe39ed1b5, 0x1b4c6a88, 0xb8c12c1f, 0x7f466551,
  0x049d5eea, 0x5d018c35, 0x73fa8774, 0x2efb0b41,
  0x5ab3671d, 0x5292dbd2, 0x33e91056, 0x136dd647,
  0x8c9ad761, 0x7a37a10c, 0x8e59f814, 0x89eb133c,
  0xeecea927, 0x35b761c9, 0xede11ce5, 0x3c7a47b1,
  0x599cd2df, 0x3f55f273, 0x791814ce, 0xbf73c737,
  0xea53f7cd, 0x5b5ffdaa, 0x14df3d6f, 0x867844db,
  0x81caaff3, 0x3eb968c4, 0x2c382434, 0x5fc2a340,
  0x72161dc3, 0x0cbce225, 0x8b283c49, 0x41ff0d95,
  0x7139a801, 0xde080cb3, 0x9cd8b4e4, 0x906456c1,
  0x617bcb84, 0x70d532b6, 0x74486c5c, 0x42d0b857
};

static const uint32_t aes_TD2[256] = {
  0xa75051f4, 0x65537e41, 0xa4c31a17, 0x5e963a27,
  0x6bcb3bab, 0x45f11f9d, 0x58abacfa, 0x03934be3,
  0xfa552030, 0x6df6ad76, 0x769188cc, 0x4c25f502,
  0xd7fc4fe5, 0xcbd7c52a, 0x44802635, 0xa38fb562,
  0x5a49deb1, 0x1b6725ba, 0x0e9845ea, 0xc0e15dfe,
  0x7502c32f, 0xf012814c, 0x97a38d46, 0xf9c66bd3,
  0x5fe7038f, 0x9c951592, 0x7aebbf6d, 0x59da9552,
  0x832dd4be, 0x21d35874, 0x692949e0, 0xc8448ec9,
  0x896a75c2, 0x7978f48e, 0x3e6b9958, 0x71dd27b9,
  0x4fb6bee1, 0xad17f088, 0xac66c920, 0x3ab47dce,
  0x4a1863df, 0x3182e51a, 0x33609751, 0x7f456253,
  0x77e0b164, 0xae84bb6b, 0xa01cfe81, 0x2b94f908,
  0x68587048, 0xfd198f45, 0x6c8794de, 0xf8b7527b,
  0xd323ab73, 0x02e2724b, 0x8f57e31f, 0xab2a6655,
  0x2807b2eb, 0xc2032fb5, 0x7b9a86c5, 0x08a5d337,
  0x87f23028, 0xa5b223bf, 0x6aba0203, 0x825ced16,
  0x1c2b8acf, 0xb492a779, 0xf2f0f307, 0xe2a14e69,
  0xf4cd65da, 0xbed50605, 0x621fd134, 0xfe8ac4a6,
  0x539d342e, 0x55a0a2f3, 0xe132058a, 0xeb75a4f6,
  0xec390b83, 0xefaa4060, 0x9f065e71, 0x1051bd6e,
  0x8af93e21, 0x063d96dd, 0x05aedd3e, 0xbd464de6,
  0x8db59154, 0x5d0571c4, 0xd46f0406, 0x15ff6050,
  0xfb241998, 0xe997d6bd, 0x43cc8940, 0x9e7767d9,
  0x42bdb0e8, 0x8b880789, 0x5b38e719, 0xeedb79c8,
  0x0a47a17c, 0x0fe97c42, 0x1ec9f884, 0x00000000,
  0x86830980, 0xed48322b, 0x70ac1e11, 0x724e6c5a,
  0xfffbfd0e, 0x38560f85, 0xd51e3dae, 0x3927362d,
  0xd9640a0f, 0xa621685c, 0x54d19b5b, 0x2e3a2436,
  0x67b10c0a, 0xe70f9357, 0x96d2b4ee, 0x919e1b9b,
  0xc54f80c0, 0x20a261dc, 0x4b695a77, 0x1a161c12,
  0xba0ae293, 0x2ae5c0a0, 0xe0433c22, 0x171d121b,
  0x0d0b0e09, 0xc7adf28b, 0xa8b92db6, 0xa9c8141e,
  0x198557f1, 0x074caf75, 0xddbbee99, 0x60fda37f,
  0x269ff701, 0xf5bc5c72, 0x3bc54466, 0x7e345bfb,
  0x29768b43, 0xc6dccb23, 0xfc68b6ed, 0xf163b8e4,
  0xdccad731, 0x85104263, 0x22401397, 0x112084c6,
  0x247d854a, 0x3df8d2bb, 0x3211aef9, 0xa16dc729,
  0x2f4b1d9e, 0x30f3dcb2, 0x52ec0d86, 0xe3d077c1,
  0x166c2bb3, 0xb999a970, 0x48fa1194, 0x642247e9,
  0x8cc4a8fc, 0x3f1aa0f0, 0x2cd8567d, 0x90ef2233,
  0x4ec78749, 0xd1c1d938, 0xa2fe8cca, 0x0b3698d4,
  0x81cfa6f5, 0xde28a57a, 0x8e26dab7, 0xbfa43fad,
  0x9de42c3a, 0x920d5078, 0xcc9b6a5f, 0x4662547e,
  0x13c2f68d, 0xb8e890d8, 0xf75e2e39, 0xaff582c3,
  0x80be9f5d, 0x937c69d0, 0x2da96fd5, 0x12b3cf25,
  0x993bc8ac, 0x7da71018, 0x636ee89c, 0xbb7bdb3b,
  0x7809cd26, 0x18f46e59, 0xb701ec9a, 0x9aa8834f,
  0x6e65e695, 0xe67eaaff, 0xcf0821bc, 0xe8e6ef15,
  0x9bd9bae7, 0x36ce4a6f, 0x09d4ea9f, 0x7cd629b0,
  0xb2af31a4, 0x23312a3f, 0x9430c6a5, 0x66c035a2,
  0xbc37744e, 0xcaa6fc82, 0xd0b0e090, 0xd81533a7,
  0x984af104, 0xdaf741ec, 0x500e7fcd, 0xf62f1791,
  0xd68d764d, 0xb04d43ef, 0x4d54ccaa, 0x04dfe496,
  0xb5e39ed1, 0x881b4c6a, 0x1fb8c12c, 0x517f4665,
  0xea049d5e, 0x355d018c, 0x7473fa87, 0x412efb0b,
  0x1d5ab367, 0xd25292db, 0x5633e910, 0x47136dd6,
  0x618c9ad7, 0x0c7a37a1, 0x148e59f8, 0x3c89eb13,
  0x27eecea9, 0xc935b761, 0xe5ede11c, 0xb13c7a47,
  0xdf599cd2, 0x733f55f2, 0xce791814, 0x37bf73c7,
  0xcdea53f7, 0xaa5b5ffd, 0x6f14df3d, 0xdb867844,
  0xf381caaf, 0xc43eb968, 0x342c3824, 0x405fc2a3,
  0xc372161d, 0x250cbce2, 0x498b283c, 0x9541ff0d,
  0x017139a8, 0xb3de080c, 0xe49cd8b4, 0xc1906456,
  0x84617bcb, 0xb670d532, 0x5c74486c, 0x5742d0b8
};

static const uint32_t aes_TD3[256] = {
  0xf4a75051, 0x4165537e, 0x17a4c31a, 0x275e963a,
  0xab6bcb3b, 0x9d45f11f, 0xfa58abac, 0xe303934b,
  0x30fa5520, 0x766df6ad, 0xcc769188, 0x024c25f5,
  0xe5d7fc4f, 0x2acbd7c5, 0x35448026, 0x62a38fb5,
  0xb15a49de, 0xba1b6725, 0xea0e9845, 0xfec0e15d,
  0x2f7502c3, 0x4cf01281, 0x4697a38d, 0xd3f9c66b,
  0x8f5fe703, 0x929c9515, 0x6d7aebbf, 0x5259da95,
  0xbe832dd4, 0x7421d358, 0xe0692949, 0xc9c8448e,
  0xc2896a75, 0x8e7978f4, 0x583e6b99, 0xb971dd27,
  0xe14fb6be, 0x88ad17f0, 0x20ac66c9, 0xce3ab47d,
  0xdf4a1863, 0x1a3182e5, 0x51336097, 0x537f4562,
  0x6477e0b1, 0x6bae84bb, 0x81a01cfe, 0x082b94f9,
  0x48685870, 0x45fd198f, 0xde6c8794, 0x7bf8b752,
  0x73d323ab, 0x4b02e272, 0x1f8f57e3, 0x55ab2a66,
  0xeb2807b2, 0xb5c2032f, 0xc57b9a86, 0x3708a5d3,
  0x2887f230, 0xbfa5b223, 0x036aba02, 0x16825ced,
  0xcf1c2b8a, 0x79b492a7, 0x07f2f0f3, 0x69e2a14e,
  0xdaf4cd65, 0x05bed506, 0x34621fd1, 0xa6fe8ac4,
  0x2e539d34, 0xf355a0a2, 0x8ae13205, 0xf6eb75a4,
  0x83ec390b, 0x60efaa40, 0x719f065e, 0x6e1051bd,
  0x218af93e, 0xdd063d96, 0x3e05aedd, 0xe6bd464d,
  0x548db591, 0xc45d0571, 0x06d46f04, 0x5015ff60,
  0x98fb2419, 0xbde997d6, 0x4043cc89, 0xd99e7767,
  0xe842bdb0, 0x898b8807, 0x195b38e7, 0xc8eedb79,
  0x7c0a47a1, 0x420fe97c, 0x841ec9f8, 0x00000000,
  0x80868309, 0x2bed4832, 0x1170ac1e, 0x5a724e6c,
  0x0efffbfd, 0x8538560f, 0xaed51e3d, 0x2d392736,
  0x0fd9640a, 0x5ca62168, 0x5b54d19b, 0x362e3a24,
  0x0a67b10c, 0x57e70f93, 0xee96d2b4, 0x9b919e1b,
  0xc0c54f80, 0xdc20a261, 0x774b695a, 0x121a161c,
  0x93ba0ae2, 0xa02ae5c0, 0x22e0433c, 0x1b171d12,
  0x090d0b0e, 0x8bc7adf2, 0xb6a8b92d, 0x1ea9c814,
  0xf1198557, 0x75074caf, 0x99ddbbee, 0x7f60fda3,
  0x01269ff7, 0x72f5bc5c, 0x663bc544, 0xfb7e345b,
  0x4329768b, 0x23c6dccb, 0xedfc68b6, 0xe4f163b8,
  0x31dccad7, 0x63851042, 0x97224013, 0xc6112084,
  0x4a247d85, 0xbb3df8d2, 0xf93211ae, 0x29a16dc7,
  0x9e2f4b1d, 0xb230f3dc, 0x8652ec0d, 0xc1e3d077,
  0xb3166c2b, 0x70b999a9, 0x9448fa11, 0xe9642247,
  0xfc8cc4a8, 0xf03f1aa0, 0x7d2cd856, 0x3390ef22,
  0x494ec787, 0x38d1c1d9, 0xcaa2fe8c, 0xd40b3698,
  0xf581cfa6, 0x7ade28a5, 0xb78e26da, 0xadbfa43f,
  0x3a9de42c, 0x78920d50, 0x5fcc9b6a, 0x7e466254,
  0x8d13c2f6, 0xd8b8e890, 0x39f75e2e, 0xc3aff582,
  0x5d80be9f, 0xd0937c69, 0xd52da96f, 0x2512b3cf,
  0xac993bc8, 0x187da710, 0x9c636ee8, 0x3bbb7bdb,
  0x267809cd, 0x5918f46e, 0x9ab701ec, 0x4f9aa883,
  0x956e65e6, 0xffe67eaa, 0xbccf0821, 0x15e8e6ef,
  0xe79bd9ba, 0x6f36ce4a, 0x9f09d4ea, 0xb07cd629,
  0xa4b2af31, 0x3f23312a, 0xa59430c6, 0xa266c035,
  0x4ebc3774, 0x82caa6fc, 0x90d0b0e0, 0xa7d81533,
  0x04984af1, 0xecdaf741, 0xcd500e7f, 0x91f62f17,
  0x4dd68d76, 0xefb04d43, 0xaa4d54cc, 0x9604dfe4,
  0xd1b5e39e, 0x6a881b4c, 0x2c1fb8c1, 0x65517f46,
  0x5eea049d, 0x8c355d01, 0x877473fa, 0x0b412efb,
  0x671d5ab3, 0xdbd25292, 0x105633e9, 0xd647136d,
  0xd7618c9a, 0xa10c7a37, 0xf8148e59, 0x133c89eb,
  0xa927eece, 0x61c935b7, 0x1ce5ede1, 0x47b13c7a,
  0xd2df599c, 0xf2733f55, 0x14ce7918, 0xc737bf73,
  0xf7cdea53, 0xfdaa5b5f, 0x3d6f14df, 0x44db8678,
  0xaff381ca, 0x68c43eb9, 0x24342c38, 0xa3405fc2,
  0x1dc37216, 0xe2250cbc, 0x3c498b28, 0x0d9541ff,
  0xa8017139, 0x0cb3de08, 0xb4e49cd8, 0x56c19064,
  0xcb84617b, 0x32b670d5, 0x6c5c7448, 0xb85742d0
};

static const uint8_t aes_TD4[256] = {
  0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38,
  0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
  0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87,
  0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
  0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d,
  0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
  0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2,
  0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
  0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16,
  0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
  0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda,
  0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
  0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a,
  0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
  0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02,
  0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
  0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea,
  0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
  0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85,
  0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
  0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89,
  0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
  0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20,
  0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
  0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31,
  0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
  0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d,
  0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
  0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0,
  0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
  0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26,
  0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d
};

static const uint32_t aes_RCON[10] = {
  0x01000000, 0x02000000, 0x04000000, 0x08000000,
  0x10000000, 0x20000000, 0x40000000, 0x80000000,
  0x1b000000, 0x36000000
};

#define TE0 aes_TE0
#define TE1 aes_TE1
#define TE2 aes_TE2
#define TE3 aes_TE3
#define TD0 aes_TD0
#define TD1 aes_TD1
#define TD2 aes_TD2
#define TD3 aes_TD3
#define TD4 aes_TD4
#define RCON aes_RCON
#define K (ctx->key)

void
aes_init_encrypt(aes_t *ctx, unsigned int bits, const unsigned char *key) {
  size_t p = 0;
  size_t i = 0;
  uint32_t tmp;

  if (bits == 128)
    ctx->rounds = 10;
  else if (bits == 192)
    ctx->rounds = 12;
  else if (bits == 256)
    ctx->rounds = 14;
  else
    ASSERT(0);

  memset(K, 0, sizeof(K));

  K[0] = read32be(key + 0);
  K[1] = read32be(key + 4);
  K[2] = read32be(key + 8);
  K[3] = read32be(key + 12);

  if (bits == 128) {
    for (;;) {
      tmp = K[p + 3];

      K[p + 4] = K[p]
        ^ (TE2[(tmp >> 16) & 0xff] & 0xff000000)
        ^ (TE3[(tmp >>  8) & 0xff] & 0x00ff0000)
        ^ (TE0[(tmp >>  0) & 0xff] & 0x0000ff00)
        ^ (TE1[(tmp >> 24) & 0xff] & 0x000000ff)
        ^ RCON[i];

      K[p + 5] = K[p + 1] ^ K[p + 4];
      K[p + 6] = K[p + 2] ^ K[p + 5];
      K[p + 7] = K[p + 3] ^ K[p + 6];

      i += 1;

      if (i == 10)
        break;

      p += 4;
    }

    return;
  }

  K[p + 4] = read32be(key + 16);
  K[p + 5] = read32be(key + 20);

  if (bits == 192) {
    for (;;) {
      tmp = K[p + 5];

      K[p + 6] = K[p]
        ^ (TE2[(tmp >> 16) & 0xff] & 0xff000000)
        ^ (TE3[(tmp >>  8) & 0xff] & 0x00ff0000)
        ^ (TE0[(tmp >>  0) & 0xff] & 0x0000ff00)
        ^ (TE1[(tmp >> 24) & 0xff] & 0x000000ff)
        ^ RCON[i];

      K[p + 7] = K[p + 1] ^ K[p + 6];
      K[p + 8] = K[p + 2] ^ K[p + 7];
      K[p + 9] = K[p + 3] ^ K[p + 8];

      i += 1;

      if (i == 8)
        break;

      K[p + 10] = K[p + 4] ^ K[p +  9];
      K[p + 11] = K[p + 5] ^ K[p + 10];
      p += 6;
    }

    return;
  }

  K[p + 6] = read32be(key + 24);
  K[p + 7] = read32be(key + 28);

  if (bits == 256) {
    for (;;) {
      tmp = K[p + 7];

      K[p + 8] = K[p]
        ^ (TE2[(tmp >> 16) & 0xff] & 0xff000000)
        ^ (TE3[(tmp >>  8) & 0xff] & 0x00ff0000)
        ^ (TE0[(tmp >>  0) & 0xff] & 0x0000ff00)
        ^ (TE1[(tmp >> 24) & 0xff] & 0x000000ff)
        ^ RCON[i];

      K[p +  9] = K[p + 1] ^ K[p +  8];
      K[p + 10] = K[p + 2] ^ K[p +  9];
      K[p + 11] = K[p + 3] ^ K[p + 10];

      i += 1;

      if (i == 7)
        break;

      tmp = K[p + 11];

      K[p + 12] = K[p + 4]
        ^ (TE2[(tmp >> 24) & 0xff] & 0xff000000)
        ^ (TE3[(tmp >> 16) & 0xff] & 0x00ff0000)
        ^ (TE0[(tmp >>  8) & 0xff] & 0x0000ff00)
        ^ (TE1[(tmp >>  0) & 0xff] & 0x000000ff);

      K[p + 13] = K[p +  5] ^ K[p + 12];
      K[p + 14] = K[p +  6] ^ K[p + 13];
      K[p + 15] = K[p +  7] ^ K[p + 14];

      p += 8;
    }

    return;
  }

  ASSERT(0);
}

void
aes_init_decrypt(aes_t *ctx, unsigned int bits, const unsigned char *key) {
  size_t i, j, p;
  uint32_t tmp;

  aes_init_encrypt(ctx, bits, key);

  /* Invert the order of the round keys. */
  for (i = 0, j = 4 * ctx->rounds; i < j; i += 4, j -= 4) {
    tmp = K[i + 0];
    K[i + 0] = K[j + 0];
    K[j + 0] = tmp;

    tmp = K[i + 1];
    K[i + 1] = K[j + 1];
    K[j + 1] = tmp;

    tmp = K[i + 2];
    K[i + 2] = K[j + 2];
    K[j + 2] = tmp;

    tmp = K[i + 3];
    K[i + 3] = K[j + 3];
    K[j + 3] = tmp;
  }

  p = 0;

  /* Apply the inverse MixColumn transform to
     all round keys but the first and the last. */
  for (i = 1; i < ctx->rounds; i++) {
    p += 4;

    K[p + 0] = TD0[TE1[(K[p + 0] >> 24) & 0xff] & 0xff]
             ^ TD1[TE1[(K[p + 0] >> 16) & 0xff] & 0xff]
             ^ TD2[TE1[(K[p + 0] >>  8) & 0xff] & 0xff]
             ^ TD3[TE1[(K[p + 0] >>  0) & 0xff] & 0xff];

    K[p + 1] = TD0[TE1[(K[p + 1] >> 24) & 0xff] & 0xff]
             ^ TD1[TE1[(K[p + 1] >> 16) & 0xff] & 0xff]
             ^ TD2[TE1[(K[p + 1] >>  8) & 0xff] & 0xff]
             ^ TD3[TE1[(K[p + 1] >>  0) & 0xff] & 0xff];

    K[p + 2] = TD0[TE1[(K[p + 2] >> 24) & 0xff] & 0xff]
             ^ TD1[TE1[(K[p + 2] >> 16) & 0xff] & 0xff]
             ^ TD2[TE1[(K[p + 2] >>  8) & 0xff] & 0xff]
             ^ TD3[TE1[(K[p + 2] >>  0) & 0xff] & 0xff];

    K[p + 3] = TD0[TE1[(K[p + 3] >> 24) & 0xff] & 0xff]
             ^ TD1[TE1[(K[p + 3] >> 16) & 0xff] & 0xff]
             ^ TD2[TE1[(K[p + 3] >>  8) & 0xff] & 0xff]
             ^ TD3[TE1[(K[p + 3] >>  0) & 0xff] & 0xff];
  }
}

void
aes_encrypt(const aes_t *ctx, unsigned char *dst, const unsigned char *src) {
  uint32_t s0, s1, s2, s3;
  uint32_t t0, t1, t2, t3;
  size_t r, p;

  /* Map byte array block to cipher
     state and add initial round key. */
  s0 = read32be(src +  0) ^ K[0];
  s1 = read32be(src +  4) ^ K[1];
  s2 = read32be(src +  8) ^ K[2];
  s3 = read32be(src + 12) ^ K[3];

  /* Nr - 1 full rounds */
  r = ctx->rounds >> 1;
  p = 0;

  for (;;) {
    t0 = TE0[(s0 >> 24) & 0xff]
       ^ TE1[(s1 >> 16) & 0xff]
       ^ TE2[(s2 >>  8) & 0xff]
       ^ TE3[(s3 >>  0) & 0xff]
       ^ K[p + 4];

    t1 = TE0[(s1 >> 24) & 0xff]
       ^ TE1[(s2 >> 16) & 0xff]
       ^ TE2[(s3 >>  8) & 0xff]
       ^ TE3[(s0 >>  0) & 0xff]
       ^ K[p + 5];

    t2 = TE0[(s2 >> 24) & 0xff]
       ^ TE1[(s3 >> 16) & 0xff]
       ^ TE2[(s0 >>  8) & 0xff]
       ^ TE3[(s1 >>  0) & 0xff]
       ^ K[p + 6];

    t3 = TE0[(s3 >> 24) & 0xff]
       ^ TE1[(s0 >> 16) & 0xff]
       ^ TE2[(s1 >>  8) & 0xff]
       ^ TE3[(s2 >>  0) & 0xff]
       ^ K[p + 7];

    p += 8;
    r -= 1;

    if (r == 0)
      break;

    s0 = TE0[(t0 >> 24) & 0xff]
       ^ TE1[(t1 >> 16) & 0xff]
       ^ TE2[(t2 >>  8) & 0xff]
       ^ TE3[(t3 >>  0) & 0xff]
       ^ K[p + 0];

    s1 = TE0[(t1 >> 24) & 0xff]
       ^ TE1[(t2 >> 16) & 0xff]
       ^ TE2[(t3 >>  8) & 0xff]
       ^ TE3[(t0 >>  0) & 0xff]
       ^ K[p + 1];

    s2 = TE0[(t2 >> 24) & 0xff]
       ^ TE1[(t3 >> 16) & 0xff]
       ^ TE2[(t0 >>  8) & 0xff]
       ^ TE3[(t1 >>  0) & 0xff]
       ^ K[p + 2];

    s3 = TE0[(t3 >> 24) & 0xff]
       ^ TE1[(t0 >> 16) & 0xff]
       ^ TE2[(t1 >>  8) & 0xff]
       ^ TE3[(t2 >>  0) & 0xff]
       ^ K[p + 3];
  }

  /* Apply last round and map cipher
     state to byte array block. */
  s0 = (TE2[(t0 >> 24) & 0xff] & 0xff000000)
     ^ (TE3[(t1 >> 16) & 0xff] & 0x00ff0000)
     ^ (TE0[(t2 >>  8) & 0xff] & 0x0000ff00)
     ^ (TE1[(t3 >>  0) & 0xff] & 0x000000ff)
     ^ K[p + 0];

  s1 = (TE2[(t1 >> 24) & 0xff] & 0xff000000)
     ^ (TE3[(t2 >> 16) & 0xff] & 0x00ff0000)
     ^ (TE0[(t3 >>  8) & 0xff] & 0x0000ff00)
     ^ (TE1[(t0 >>  0) & 0xff] & 0x000000ff)
     ^ K[p + 1];

  s2 = (TE2[(t2 >> 24) & 0xff] & 0xff000000)
     ^ (TE3[(t3 >> 16) & 0xff] & 0x00ff0000)
     ^ (TE0[(t0 >>  8) & 0xff] & 0x0000ff00)
     ^ (TE1[(t1 >>  0) & 0xff] & 0x000000ff)
     ^ K[p + 2];

  s3 = (TE2[(t3 >> 24) & 0xff] & 0xff000000)
     ^ (TE3[(t0 >> 16) & 0xff] & 0x00ff0000)
     ^ (TE0[(t1 >>  8) & 0xff] & 0x0000ff00)
     ^ (TE1[(t2 >>  0) & 0xff] & 0x000000ff)
     ^ K[p + 3];

  write32be(dst +  0, s0);
  write32be(dst +  4, s1);
  write32be(dst +  8, s2);
  write32be(dst + 12, s3);
}

void
aes_decrypt(const aes_t *ctx, unsigned char *dst, const unsigned char *src) {
  uint32_t s0, s1, s2, s3;
  uint32_t t0, t1, t2, t3;
  size_t r, p;

  /* Map byte array block to cipher
     state and add initial round key. */
  s0 = read32be(src +  0) ^ K[0];
  s1 = read32be(src +  4) ^ K[1];
  s2 = read32be(src +  8) ^ K[2];
  s3 = read32be(src + 12) ^ K[3];

  /* Nr - 1 full rounds */
  r = ctx->rounds >> 1;
  p = 0;

  for (;;) {
    t0 = TD0[(s0 >> 24) & 0xff]
       ^ TD1[(s3 >> 16) & 0xff]
       ^ TD2[(s2 >>  8) & 0xff]
       ^ TD3[(s1 >>  0) & 0xff]
       ^ K[p + 4];

    t1 = TD0[(s1 >> 24) & 0xff]
       ^ TD1[(s0 >> 16) & 0xff]
       ^ TD2[(s3 >>  8) & 0xff]
       ^ TD3[(s2 >>  0) & 0xff]
       ^ K[p + 5];

    t2 = TD0[(s2 >> 24) & 0xff]
       ^ TD1[(s1 >> 16) & 0xff]
       ^ TD2[(s0 >>  8) & 0xff]
       ^ TD3[(s3 >>  0) & 0xff]
       ^ K[p + 6];

    t3 = TD0[(s3 >> 24) & 0xff]
       ^ TD1[(s2 >> 16) & 0xff]
       ^ TD2[(s1 >>  8) & 0xff]
       ^ TD3[(s0 >>  0) & 0xff]
       ^ K[p + 7];

    p += 8;
    r -= 1;

    if (r == 0)
      break;

    s0 = TD0[(t0 >> 24) & 0xff]
       ^ TD1[(t3 >> 16) & 0xff]
       ^ TD2[(t2 >>  8) & 0xff]
       ^ TD3[(t1 >>  0) & 0xff]
       ^ K[p + 0];

    s1 = TD0[(t1 >> 24) & 0xff]
       ^ TD1[(t0 >> 16) & 0xff]
       ^ TD2[(t3 >>  8) & 0xff]
       ^ TD3[(t2 >>  0) & 0xff]
       ^ K[p + 1];

    s2 = TD0[(t2 >> 24) & 0xff]
       ^ TD1[(t1 >> 16) & 0xff]
       ^ TD2[(t0 >>  8) & 0xff]
       ^ TD3[(t3 >>  0) & 0xff]
       ^ K[p + 2];

    s3 = TD0[(t3 >> 24) & 0xff]
       ^ TD1[(t2 >> 16) & 0xff]
       ^ TD2[(t1 >>  8) & 0xff]
       ^ TD3[(t0 >>  0) & 0xff]
       ^ K[p + 3];
  }

  /* Apply last round and map cipher
     state to byte array block. */
  s0 = (TD4[(t0 >> 24) & 0xff] << 24)
     ^ (TD4[(t3 >> 16) & 0xff] << 16)
     ^ (TD4[(t2 >>  8) & 0xff] <<  8)
     ^ (TD4[(t1 >>  0) & 0xff] <<  0)
     ^ K[p + 0];

  s1 = (TD4[(t1 >> 24) & 0xff] << 24)
     ^ (TD4[(t0 >> 16) & 0xff] << 16)
     ^ (TD4[(t3 >>  8) & 0xff] <<  8)
     ^ (TD4[(t2 >>  0) & 0xff] <<  0)
     ^ K[p + 1];

  s2 = (TD4[(t2 >> 24) & 0xff] << 24)
     ^ (TD4[(t1 >> 16) & 0xff] << 16)
     ^ (TD4[(t0 >>  8) & 0xff] <<  8)
     ^ (TD4[(t3 >>  0) & 0xff] <<  0)
     ^ K[p + 2];

  s3 = (TD4[(t3 >> 24) & 0xff] << 24)
     ^ (TD4[(t2 >> 16) & 0xff] << 16)
     ^ (TD4[(t1 >>  8) & 0xff] <<  8)
     ^ (TD4[(t0 >>  0) & 0xff] <<  0)
     ^ K[p + 3];

  write32be(dst +  0, s0);
  write32be(dst +  4, s1);
  write32be(dst +  8, s2);
  write32be(dst + 12, s3);
}

#undef TE0
#undef TE1
#undef TE2
#undef TE3
#undef TD0
#undef TD1
#undef TD2
#undef TD3
#undef TD4
#undef RCON
#undef K

/**
 * Blowfish
 *
 * Resources:
 *   https://en.wikipedia.org/wiki/Blowfish_(cipher)
 *   https://www.schneier.com/blowfish.html
 *   https://github.com/joyent/node-bcrypt-pbkdf/blob/master/index.js
 */

static const uint32_t blowfish_S0[256] = {
  0xd1310ba6, 0x98dfb5ac, 0x2ffd72db, 0xd01adfb7,
  0xb8e1afed, 0x6a267e96, 0xba7c9045, 0xf12c7f99,
  0x24a19947, 0xb3916cf7, 0x0801f2e2, 0x858efc16,
  0x636920d8, 0x71574e69, 0xa458fea3, 0xf4933d7e,
  0x0d95748f, 0x728eb658, 0x718bcd58, 0x82154aee,
  0x7b54a41d, 0xc25a59b5, 0x9c30d539, 0x2af26013,
  0xc5d1b023, 0x286085f0, 0xca417918, 0xb8db38ef,
  0x8e79dcb0, 0x603a180e, 0x6c9e0e8b, 0xb01e8a3e,
  0xd71577c1, 0xbd314b27, 0x78af2fda, 0x55605c60,
  0xe65525f3, 0xaa55ab94, 0x57489862, 0x63e81440,
  0x55ca396a, 0x2aab10b6, 0xb4cc5c34, 0x1141e8ce,
  0xa15486af, 0x7c72e993, 0xb3ee1411, 0x636fbc2a,
  0x2ba9c55d, 0x741831f6, 0xce5c3e16, 0x9b87931e,
  0xafd6ba33, 0x6c24cf5c, 0x7a325381, 0x28958677,
  0x3b8f4898, 0x6b4bb9af, 0xc4bfe81b, 0x66282193,
  0x61d809cc, 0xfb21a991, 0x487cac60, 0x5dec8032,
  0xef845d5d, 0xe98575b1, 0xdc262302, 0xeb651b88,
  0x23893e81, 0xd396acc5, 0x0f6d6ff3, 0x83f44239,
  0x2e0b4482, 0xa4842004, 0x69c8f04a, 0x9e1f9b5e,
  0x21c66842, 0xf6e96c9a, 0x670c9c61, 0xabd388f0,
  0x6a51a0d2, 0xd8542f68, 0x960fa728, 0xab5133a3,
  0x6eef0b6c, 0x137a3be4, 0xba3bf050, 0x7efb2a98,
  0xa1f1651d, 0x39af0176, 0x66ca593e, 0x82430e88,
  0x8cee8619, 0x456f9fb4, 0x7d84a5c3, 0x3b8b5ebe,
  0xe06f75d8, 0x85c12073, 0x401a449f, 0x56c16aa6,
  0x4ed3aa62, 0x363f7706, 0x1bfedf72, 0x429b023d,
  0x37d0d724, 0xd00a1248, 0xdb0fead3, 0x49f1c09b,
  0x075372c9, 0x80991b7b, 0x25d479d8, 0xf6e8def7,
  0xe3fe501a, 0xb6794c3b, 0x976ce0bd, 0x04c006ba,
  0xc1a94fb6, 0x409f60c4, 0x5e5c9ec2, 0x196a2463,
  0x68fb6faf, 0x3e6c53b5, 0x1339b2eb, 0x3b52ec6f,
  0x6dfc511f, 0x9b30952c, 0xcc814544, 0xaf5ebd09,
  0xbee3d004, 0xde334afd, 0x660f2807, 0x192e4bb3,
  0xc0cba857, 0x45c8740f, 0xd20b5f39, 0xb9d3fbdb,
  0x5579c0bd, 0x1a60320a, 0xd6a100c6, 0x402c7279,
  0x679f25fe, 0xfb1fa3cc, 0x8ea5e9f8, 0xdb3222f8,
  0x3c7516df, 0xfd616b15, 0x2f501ec8, 0xad0552ab,
  0x323db5fa, 0xfd238760, 0x53317b48, 0x3e00df82,
  0x9e5c57bb, 0xca6f8ca0, 0x1a87562e, 0xdf1769db,
  0xd542a8f6, 0x287effc3, 0xac6732c6, 0x8c4f5573,
  0x695b27b0, 0xbbca58c8, 0xe1ffa35d, 0xb8f011a0,
  0x10fa3d98, 0xfd2183b8, 0x4afcb56c, 0x2dd1d35b,
  0x9a53e479, 0xb6f84565, 0xd28e49bc, 0x4bfb9790,
  0xe1ddf2da, 0xa4cb7e33, 0x62fb1341, 0xcee4c6e8,
  0xef20cada, 0x36774c01, 0xd07e9efe, 0x2bf11fb4,
  0x95dbda4d, 0xae909198, 0xeaad8e71, 0x6b93d5a0,
  0xd08ed1d0, 0xafc725e0, 0x8e3c5b2f, 0x8e7594b7,
  0x8ff6e2fb, 0xf2122b64, 0x8888b812, 0x900df01c,
  0x4fad5ea0, 0x688fc31c, 0xd1cff191, 0xb3a8c1ad,
  0x2f2f2218, 0xbe0e1777, 0xea752dfe, 0x8b021fa1,
  0xe5a0cc0f, 0xb56f74e8, 0x18acf3d6, 0xce89e299,
  0xb4a84fe0, 0xfd13e0b7, 0x7cc43b81, 0xd2ada8d9,
  0x165fa266, 0x80957705, 0x93cc7314, 0x211a1477,
  0xe6ad2065, 0x77b5fa86, 0xc75442f5, 0xfb9d35cf,
  0xebcdaf0c, 0x7b3e89a0, 0xd6411bd3, 0xae1e7e49,
  0x00250e2d, 0x2071b35e, 0x226800bb, 0x57b8e0af,
  0x2464369b, 0xf009b91e, 0x5563911d, 0x59dfa6aa,
  0x78c14389, 0xd95a537f, 0x207d5ba2, 0x02e5b9c5,
  0x83260376, 0x6295cfa9, 0x11c81968, 0x4e734a41,
  0xb3472dca, 0x7b14a94a, 0x1b510052, 0x9a532915,
  0xd60f573f, 0xbc9bc6e4, 0x2b60a476, 0x81e67400,
  0x08ba6fb5, 0x571be91f, 0xf296ec6b, 0x2a0dd915,
  0xb6636521, 0xe7b9f9b6, 0xff34052e, 0xc5855664,
  0x53b02d5d, 0xa99f8fa1, 0x08ba4799, 0x6e85076a
};

static const uint32_t blowfish_S1[256] = {
  0x4b7a70e9, 0xb5b32944, 0xdb75092e, 0xc4192623,
  0xad6ea6b0, 0x49a7df7d, 0x9cee60b8, 0x8fedb266,
  0xecaa8c71, 0x699a17ff, 0x5664526c, 0xc2b19ee1,
  0x193602a5, 0x75094c29, 0xa0591340, 0xe4183a3e,
  0x3f54989a, 0x5b429d65, 0x6b8fe4d6, 0x99f73fd6,
  0xa1d29c07, 0xefe830f5, 0x4d2d38e6, 0xf0255dc1,
  0x4cdd2086, 0x8470eb26, 0x6382e9c6, 0x021ecc5e,
  0x09686b3f, 0x3ebaefc9, 0x3c971814, 0x6b6a70a1,
  0x687f3584, 0x52a0e286, 0xb79c5305, 0xaa500737,
  0x3e07841c, 0x7fdeae5c, 0x8e7d44ec, 0x5716f2b8,
  0xb03ada37, 0xf0500c0d, 0xf01c1f04, 0x0200b3ff,
  0xae0cf51a, 0x3cb574b2, 0x25837a58, 0xdc0921bd,
  0xd19113f9, 0x7ca92ff6, 0x94324773, 0x22f54701,
  0x3ae5e581, 0x37c2dadc, 0xc8b57634, 0x9af3dda7,
  0xa9446146, 0x0fd0030e, 0xecc8c73e, 0xa4751e41,
  0xe238cd99, 0x3bea0e2f, 0x3280bba1, 0x183eb331,
  0x4e548b38, 0x4f6db908, 0x6f420d03, 0xf60a04bf,
  0x2cb81290, 0x24977c79, 0x5679b072, 0xbcaf89af,
  0xde9a771f, 0xd9930810, 0xb38bae12, 0xdccf3f2e,
  0x5512721f, 0x2e6b7124, 0x501adde6, 0x9f84cd87,
  0x7a584718, 0x7408da17, 0xbc9f9abc, 0xe94b7d8c,
  0xec7aec3a, 0xdb851dfa, 0x63094366, 0xc464c3d2,
  0xef1c1847, 0x3215d908, 0xdd433b37, 0x24c2ba16,
  0x12a14d43, 0x2a65c451, 0x50940002, 0x133ae4dd,
  0x71dff89e, 0x10314e55, 0x81ac77d6, 0x5f11199b,
  0x043556f1, 0xd7a3c76b, 0x3c11183b, 0x5924a509,
  0xf28fe6ed, 0x97f1fbfa, 0x9ebabf2c, 0x1e153c6e,
  0x86e34570, 0xeae96fb1, 0x860e5e0a, 0x5a3e2ab3,
  0x771fe71c, 0x4e3d06fa, 0x2965dcb9, 0x99e71d0f,
  0x803e89d6, 0x5266c825, 0x2e4cc978, 0x9c10b36a,
  0xc6150eba, 0x94e2ea78, 0xa5fc3c53, 0x1e0a2df4,
  0xf2f74ea7, 0x361d2b3d, 0x1939260f, 0x19c27960,
  0x5223a708, 0xf71312b6, 0xebadfe6e, 0xeac31f66,
  0xe3bc4595, 0xa67bc883, 0xb17f37d1, 0x018cff28,
  0xc332ddef, 0xbe6c5aa5, 0x65582185, 0x68ab9802,
  0xeecea50f, 0xdb2f953b, 0x2aef7dad, 0x5b6e2f84,
  0x1521b628, 0x29076170, 0xecdd4775, 0x619f1510,
  0x13cca830, 0xeb61bd96, 0x0334fe1e, 0xaa0363cf,
  0xb5735c90, 0x4c70a239, 0xd59e9e0b, 0xcbaade14,
  0xeecc86bc, 0x60622ca7, 0x9cab5cab, 0xb2f3846e,
  0x648b1eaf, 0x19bdf0ca, 0xa02369b9, 0x655abb50,
  0x40685a32, 0x3c2ab4b3, 0x319ee9d5, 0xc021b8f7,
  0x9b540b19, 0x875fa099, 0x95f7997e, 0x623d7da8,
  0xf837889a, 0x97e32d77, 0x11ed935f, 0x16681281,
  0x0e358829, 0xc7e61fd6, 0x96dedfa1, 0x7858ba99,
  0x57f584a5, 0x1b227263, 0x9b83c3ff, 0x1ac24696,
  0xcdb30aeb, 0x532e3054, 0x8fd948e4, 0x6dbc3128,
  0x58ebf2ef, 0x34c6ffea, 0xfe28ed61, 0xee7c3c73,
  0x5d4a14d9, 0xe864b7e3, 0x42105d14, 0x203e13e0,
  0x45eee2b6, 0xa3aaabea, 0xdb6c4f15, 0xfacb4fd0,
  0xc742f442, 0xef6abbb5, 0x654f3b1d, 0x41cd2105,
  0xd81e799e, 0x86854dc7, 0xe44b476a, 0x3d816250,
  0xcf62a1f2, 0x5b8d2646, 0xfc8883a0, 0xc1c7b6a3,
  0x7f1524c3, 0x69cb7492, 0x47848a0b, 0x5692b285,
  0x095bbf00, 0xad19489d, 0x1462b174, 0x23820e00,
  0x58428d2a, 0x0c55f5ea, 0x1dadf43e, 0x233f7061,
  0x3372f092, 0x8d937e41, 0xd65fecf1, 0x6c223bdb,
  0x7cde3759, 0xcbee7460, 0x4085f2a7, 0xce77326e,
  0xa6078084, 0x19f8509e, 0xe8efd855, 0x61d99735,
  0xa969a7aa, 0xc50c06c2, 0x5a04abfc, 0x800bcadc,
  0x9e447a2e, 0xc3453484, 0xfdd56705, 0x0e1e9ec9,
  0xdb73dbd3, 0x105588cd, 0x675fda79, 0xe3674340,
  0xc5c43465, 0x713e38d8, 0x3d28f89e, 0xf16dff20,
  0x153e21e7, 0x8fb03d4a, 0xe6e39f2b, 0xdb83adf7
};

static const uint32_t blowfish_S2[256] = {
  0xe93d5a68, 0x948140f7, 0xf64c261c, 0x94692934,
  0x411520f7, 0x7602d4f7, 0xbcf46b2e, 0xd4a20068,
  0xd4082471, 0x3320f46a, 0x43b7d4b7, 0x500061af,
  0x1e39f62e, 0x97244546, 0x14214f74, 0xbf8b8840,
  0x4d95fc1d, 0x96b591af, 0x70f4ddd3, 0x66a02f45,
  0xbfbc09ec, 0x03bd9785, 0x7fac6dd0, 0x31cb8504,
  0x96eb27b3, 0x55fd3941, 0xda2547e6, 0xabca0a9a,
  0x28507825, 0x530429f4, 0x0a2c86da, 0xe9b66dfb,
  0x68dc1462, 0xd7486900, 0x680ec0a4, 0x27a18dee,
  0x4f3ffea2, 0xe887ad8c, 0xb58ce006, 0x7af4d6b6,
  0xaace1e7c, 0xd3375fec, 0xce78a399, 0x406b2a42,
  0x20fe9e35, 0xd9f385b9, 0xee39d7ab, 0x3b124e8b,
  0x1dc9faf7, 0x4b6d1856, 0x26a36631, 0xeae397b2,
  0x3a6efa74, 0xdd5b4332, 0x6841e7f7, 0xca7820fb,
  0xfb0af54e, 0xd8feb397, 0x454056ac, 0xba489527,
  0x55533a3a, 0x20838d87, 0xfe6ba9b7, 0xd096954b,
  0x55a867bc, 0xa1159a58, 0xcca92963, 0x99e1db33,
  0xa62a4a56, 0x3f3125f9, 0x5ef47e1c, 0x9029317c,
  0xfdf8e802, 0x04272f70, 0x80bb155c, 0x05282ce3,
  0x95c11548, 0xe4c66d22, 0x48c1133f, 0xc70f86dc,
  0x07f9c9ee, 0x41041f0f, 0x404779a4, 0x5d886e17,
  0x325f51eb, 0xd59bc0d1, 0xf2bcc18f, 0x41113564,
  0x257b7834, 0x602a9c60, 0xdff8e8a3, 0x1f636c1b,
  0x0e12b4c2, 0x02e1329e, 0xaf664fd1, 0xcad18115,
  0x6b2395e0, 0x333e92e1, 0x3b240b62, 0xeebeb922,
  0x85b2a20e, 0xe6ba0d99, 0xde720c8c, 0x2da2f728,
  0xd0127845, 0x95b794fd, 0x647d0862, 0xe7ccf5f0,
  0x5449a36f, 0x877d48fa, 0xc39dfd27, 0xf33e8d1e,
  0x0a476341, 0x992eff74, 0x3a6f6eab, 0xf4f8fd37,
  0xa812dc60, 0xa1ebddf8, 0x991be14c, 0xdb6e6b0d,
  0xc67b5510, 0x6d672c37, 0x2765d43b, 0xdcd0e804,
  0xf1290dc7, 0xcc00ffa3, 0xb5390f92, 0x690fed0b,
  0x667b9ffb, 0xcedb7d9c, 0xa091cf0b, 0xd9155ea3,
  0xbb132f88, 0x515bad24, 0x7b9479bf, 0x763bd6eb,
  0x37392eb3, 0xcc115979, 0x8026e297, 0xf42e312d,
  0x6842ada7, 0xc66a2b3b, 0x12754ccc, 0x782ef11c,
  0x6a124237, 0xb79251e7, 0x06a1bbe6, 0x4bfb6350,
  0x1a6b1018, 0x11caedfa, 0x3d25bdd8, 0xe2e1c3c9,
  0x44421659, 0x0a121386, 0xd90cec6e, 0xd5abea2a,
  0x64af674e, 0xda86a85f, 0xbebfe988, 0x64e4c3fe,
  0x9dbc8057, 0xf0f7c086, 0x60787bf8, 0x6003604d,
  0xd1fd8346, 0xf6381fb0, 0x7745ae04, 0xd736fccc,
  0x83426b33, 0xf01eab71, 0xb0804187, 0x3c005e5f,
  0x77a057be, 0xbde8ae24, 0x55464299, 0xbf582e61,
  0x4e58f48f, 0xf2ddfda2, 0xf474ef38, 0x8789bdc2,
  0x5366f9c3, 0xc8b38e74, 0xb475f255, 0x46fcd9b9,
  0x7aeb2661, 0x8b1ddf84, 0x846a0e79, 0x915f95e2,
  0x466e598e, 0x20b45770, 0x8cd55591, 0xc902de4c,
  0xb90bace1, 0xbb8205d0, 0x11a86248, 0x7574a99e,
  0xb77f19b6, 0xe0a9dc09, 0x662d09a1, 0xc4324633,
  0xe85a1f02, 0x09f0be8c, 0x4a99a025, 0x1d6efe10,
  0x1ab93d1d, 0x0ba5a4df, 0xa186f20f, 0x2868f169,
  0xdcb7da83, 0x573906fe, 0xa1e2ce9b, 0x4fcd7f52,
  0x50115e01, 0xa70683fa, 0xa002b5c4, 0x0de6d027,
  0x9af88c27, 0x773f8641, 0xc3604c06, 0x61a806b5,
  0xf0177a28, 0xc0f586e0, 0x006058aa, 0x30dc7d62,
  0x11e69ed7, 0x2338ea63, 0x53c2dd94, 0xc2c21634,
  0xbbcbee56, 0x90bcb6de, 0xebfc7da1, 0xce591d76,
  0x6f05e409, 0x4b7c0188, 0x39720a3d, 0x7c927c24,
  0x86e3725f, 0x724d9db9, 0x1ac15bb4, 0xd39eb8fc,
  0xed545578, 0x08fca5b5, 0xd83d7cd3, 0x4dad0fc4,
  0x1e50ef5e, 0xb161e6f8, 0xa28514d9, 0x6c51133c,
  0x6fd5c7e7, 0x56e14ec4, 0x362abfce, 0xddc6c837,
  0xd79a3234, 0x92638212, 0x670efa8e, 0x406000e0
};

static const uint32_t blowfish_S3[256] = {
  0x3a39ce37, 0xd3faf5cf, 0xabc27737, 0x5ac52d1b,
  0x5cb0679e, 0x4fa33742, 0xd3822740, 0x99bc9bbe,
  0xd5118e9d, 0xbf0f7315, 0xd62d1c7e, 0xc700c47b,
  0xb78c1b6b, 0x21a19045, 0xb26eb1be, 0x6a366eb4,
  0x5748ab2f, 0xbc946e79, 0xc6a376d2, 0x6549c2c8,
  0x530ff8ee, 0x468dde7d, 0xd5730a1d, 0x4cd04dc6,
  0x2939bbdb, 0xa9ba4650, 0xac9526e8, 0xbe5ee304,
  0xa1fad5f0, 0x6a2d519a, 0x63ef8ce2, 0x9a86ee22,
  0xc089c2b8, 0x43242ef6, 0xa51e03aa, 0x9cf2d0a4,
  0x83c061ba, 0x9be96a4d, 0x8fe51550, 0xba645bd6,
  0x2826a2f9, 0xa73a3ae1, 0x4ba99586, 0xef5562e9,
  0xc72fefd3, 0xf752f7da, 0x3f046f69, 0x77fa0a59,
  0x80e4a915, 0x87b08601, 0x9b09e6ad, 0x3b3ee593,
  0xe990fd5a, 0x9e34d797, 0x2cf0b7d9, 0x022b8b51,
  0x96d5ac3a, 0x017da67d, 0xd1cf3ed6, 0x7c7d2d28,
  0x1f9f25cf, 0xadf2b89b, 0x5ad6b472, 0x5a88f54c,
  0xe029ac71, 0xe019a5e6, 0x47b0acfd, 0xed93fa9b,
  0xe8d3c48d, 0x283b57cc, 0xf8d56629, 0x79132e28,
  0x785f0191, 0xed756055, 0xf7960e44, 0xe3d35e8c,
  0x15056dd4, 0x88f46dba, 0x03a16125, 0x0564f0bd,
  0xc3eb9e15, 0x3c9057a2, 0x97271aec, 0xa93a072a,
  0x1b3f6d9b, 0x1e6321f5, 0xf59c66fb, 0x26dcf319,
  0x7533d928, 0xb155fdf5, 0x03563482, 0x8aba3cbb,
  0x28517711, 0xc20ad9f8, 0xabcc5167, 0xccad925f,
  0x4de81751, 0x3830dc8e, 0x379d5862, 0x9320f991,
  0xea7a90c2, 0xfb3e7bce, 0x5121ce64, 0x774fbe32,
  0xa8b6e37e, 0xc3293d46, 0x48de5369, 0x6413e680,
  0xa2ae0810, 0xdd6db224, 0x69852dfd, 0x09072166,
  0xb39a460a, 0x6445c0dd, 0x586cdecf, 0x1c20c8ae,
  0x5bbef7dd, 0x1b588d40, 0xccd2017f, 0x6bb4e3bb,
  0xdda26a7e, 0x3a59ff45, 0x3e350a44, 0xbcb4cdd5,
  0x72eacea8, 0xfa6484bb, 0x8d6612ae, 0xbf3c6f47,
  0xd29be463, 0x542f5d9e, 0xaec2771b, 0xf64e6370,
  0x740e0d8d, 0xe75b1357, 0xf8721671, 0xaf537d5d,
  0x4040cb08, 0x4eb4e2cc, 0x34d2466a, 0x0115af84,
  0xe1b00428, 0x95983a1d, 0x06b89fb4, 0xce6ea048,
  0x6f3f3b82, 0x3520ab82, 0x011a1d4b, 0x277227f8,
  0x611560b1, 0xe7933fdc, 0xbb3a792b, 0x344525bd,
  0xa08839e1, 0x51ce794b, 0x2f32c9b7, 0xa01fbac9,
  0xe01cc87e, 0xbcc7d1f6, 0xcf0111c3, 0xa1e8aac7,
  0x1a908749, 0xd44fbd9a, 0xd0dadecb, 0xd50ada38,
  0x0339c32a, 0xc6913667, 0x8df9317c, 0xe0b12b4f,
  0xf79e59b7, 0x43f5bb3a, 0xf2d519ff, 0x27d9459c,
  0xbf97222c, 0x15e6fc2a, 0x0f91fc71, 0x9b941525,
  0xfae59361, 0xceb69ceb, 0xc2a86459, 0x12baa8d1,
  0xb6c1075e, 0xe3056a0c, 0x10d25065, 0xcb03a442,
  0xe0ec6e0e, 0x1698db3b, 0x4c98a0be, 0x3278e964,
  0x9f1f9532, 0xe0d392df, 0xd3a0342b, 0x8971f21e,
  0x1b0a7441, 0x4ba3348c, 0xc5be7120, 0xc37632d8,
  0xdf359f8d, 0x9b992f2e, 0xe60b6f47, 0x0fe3f11d,
  0xe54cda54, 0x1edad891, 0xce6279cf, 0xcd3e7e6f,
  0x1618b166, 0xfd2c1d05, 0x848fd2c5, 0xf6fb2299,
  0xf523f357, 0xa6327623, 0x93a83531, 0x56cccd02,
  0xacf08162, 0x5a75ebb5, 0x6e163697, 0x88d273cc,
  0xde966292, 0x81b949d0, 0x4c50901b, 0x71c65614,
  0xe6c6c7bd, 0x327a140a, 0x45e1d006, 0xc3f27b9a,
  0xc9aa53fd, 0x62a80f00, 0xbb25bfe2, 0x35bdd2f6,
  0x71126905, 0xb2040222, 0xb6cbcf7c, 0xcd769c2b,
  0x53113ec0, 0x1640e3d3, 0x38abbd60, 0x2547adf0,
  0xba38209c, 0xf746ce76, 0x77afa1c5, 0x20756060,
  0x85cbfe4e, 0x8ae88dd8, 0x7aaaf9b0, 0x4cf9aa7e,
  0x1948c25c, 0x02fb8a8c, 0x01c36ae4, 0xd6ebe1f9,
  0x90d4f869, 0xa65cdea0, 0x3f09252d, 0xc208e69f,
  0xb74e6132, 0xce77e25b, 0x578fdfe3, 0x3ac372e6
};

static const uint32_t blowfish_P[18] = {
  0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344,
  0xa4093822, 0x299f31d0, 0x082efa98, 0xec4e6c89,
  0x452821e6, 0x38d01377, 0xbe5466cf, 0x34e90c6c,
  0xc0ac29b7, 0xc97c50dd, 0x3f84d5b5, 0xb5470917,
  0x9216d5d9, 0x8979fb1b
};

void
blowfish_init(blowfish_t *ctx,
              const unsigned char *key, size_t key_len,
              const unsigned char *salt, size_t salt_len) {
  if (key_len > 72)
    key_len = 72;

  if (salt_len > 1096)
    salt_len = 1096;

  memcpy(ctx->S[0], blowfish_S0, sizeof(blowfish_S0));
  memcpy(ctx->S[1], blowfish_S1, sizeof(blowfish_S1));
  memcpy(ctx->S[2], blowfish_S2, sizeof(blowfish_S2));
  memcpy(ctx->S[3], blowfish_S3, sizeof(blowfish_S3));
  memcpy(ctx->P, blowfish_P, sizeof(blowfish_P));

  if (salt_len > 0)
    blowfish_expandstate(ctx, key, key_len, salt, salt_len);
  else
    blowfish_expand0state(ctx, key, key_len);
}

/* Borrowed from nettle. */
#define F(c, x) \
  ((((c->S[0][(x >> 24) & 0xff] + c->S[1][(x >> 16) & 0xff]) \
    ^ c->S[2][(x >> 8) & 0xff]) + c->S[3][x & 0xff]) & 0xffffffff)

#define R(c, l, r, i) do { \
  l ^= c->P[i];            \
  r ^= F(c, l);            \
} while(0)

static void
blowfish_encipher(blowfish_t *ctx, uint32_t *x1, uint32_t *x2) {
  /* Borrowed from nettle. */
  uint32_t xl = *x1;
  uint32_t xr = *x2;

  R(ctx, xl, xr, 0);
  R(ctx, xr, xl, 1);
  R(ctx, xl, xr, 2);
  R(ctx, xr, xl, 3);
  R(ctx, xl, xr, 4);
  R(ctx, xr, xl, 5);
  R(ctx, xl, xr, 6);
  R(ctx, xr, xl, 7);
  R(ctx, xl, xr, 8);
  R(ctx, xr, xl, 9);
  R(ctx, xl, xr, 10);
  R(ctx, xr, xl, 11);
  R(ctx, xl, xr, 12);
  R(ctx, xr, xl, 13);
  R(ctx, xl, xr, 14);
  R(ctx, xr, xl, 15);

  xl ^= ctx->P[16];
  xr ^= ctx->P[17];

  *x1 = xr;
  *x2 = xl;
}

static void
blowfish_decipher(blowfish_t *ctx, uint32_t *x1, uint32_t *x2) {
  /* Borrowed from nettle. */
  uint32_t xl = *x1;
  uint32_t xr = *x2;

  R(ctx, xl, xr, 17);
  R(ctx, xr, xl, 16);
  R(ctx, xl, xr, 15);
  R(ctx, xr, xl, 14);
  R(ctx, xl, xr, 13);
  R(ctx, xr, xl, 12);
  R(ctx, xl, xr, 11);
  R(ctx, xr, xl, 10);
  R(ctx, xl, xr, 9);
  R(ctx, xr, xl, 8);
  R(ctx, xl, xr, 7);
  R(ctx, xr, xl, 6);
  R(ctx, xl, xr, 5);
  R(ctx, xr, xl, 4);
  R(ctx, xl, xr, 3);
  R(ctx, xr, xl, 2);

  xl ^= ctx->P[1];
  xr ^= ctx->P[0];

  *x1 = xr;
  *x2 = xl;
}

#undef F
#undef R

uint32_t
blowfish_stream2word(const unsigned char *data, size_t len, size_t *off) {
  uint32_t word;

  if (len == 0) {
    *off = 0;
    return 0;
  }

  word = ((uint32_t)data[(*off + 0) % len] << 24)
       | ((uint32_t)data[(*off + 1) % len] << 16)
       | ((uint32_t)data[(*off + 2) % len] << 8)
       | ((uint32_t)data[(*off + 3) % len] << 0);

  *off = (*off + 4) % len;

  return word;
}

void
blowfish_expand0state(blowfish_t *ctx,
                      const unsigned char *key,
                      size_t key_len) {
  uint32_t xl = 0;
  uint32_t xr = 0;
  size_t off = 0;
  size_t i, k;

  for (i = 0; i < 18; i++)
    ctx->P[i] ^= blowfish_stream2word(key, key_len, &off);

  for (i = 0; i < 18; i += 2) {
    blowfish_encipher(ctx, &xl, &xr);

    ctx->P[i + 0] = xl;
    ctx->P[i + 1] = xr;
  }

  for (i = 0; i < 4; i++) {
    for (k = 0; k < 256; k += 2) {
      blowfish_encipher(ctx, &xl, &xr);

      ctx->S[i][k + 0] = xl;
      ctx->S[i][k + 1] = xr;
    }
  }
}

void
blowfish_expandstate(blowfish_t *ctx,
                     const unsigned char *key, size_t key_len,
                     const unsigned char *data, size_t data_len) {
  uint32_t xl = 0;
  uint32_t xr = 0;
  size_t off = 0;
  size_t i, k;

  for (i = 0; i < 18; i++)
    ctx->P[i] ^= blowfish_stream2word(key, key_len, &off);

  off = 0;

  for (i = 0; i < 18; i += 2) {
    xl ^= blowfish_stream2word(data, data_len, &off);
    xr ^= blowfish_stream2word(data, data_len, &off);

    blowfish_encipher(ctx, &xl, &xr);

    ctx->P[i + 0] = xl;
    ctx->P[i + 1] = xr;
  }

  for (i = 0; i < 4; i++) {
    for (k = 0; k < 256; k += 2) {
      xl ^= blowfish_stream2word(data, data_len, &off);
      xr ^= blowfish_stream2word(data, data_len, &off);

      blowfish_encipher(ctx, &xl, &xr);

      ctx->S[i][k + 0] = xl;
      ctx->S[i][k + 1] = xr;
    }
  }
}

void
blowfish_enc(blowfish_t *ctx, uint32_t *data, size_t len) {
  size_t blocks = len / 2;

  while (blocks--) {
    blowfish_encipher(ctx, data + 0, data + 1);
    data += 2;
  }
}

void
blowfish_dec(blowfish_t *ctx, uint32_t *data, size_t len) {
  size_t blocks = len / 2;

  while (blocks--) {
    blowfish_decipher(ctx, data + 0, data + 1);
    data += 2;
  }
}

void
blowfish_encrypt(blowfish_t *ctx,
                 unsigned char *dst,
                 const unsigned char *src) {
  uint32_t xl = read32be(src + 0);
  uint32_t xr = read32be(src + 4);

  blowfish_encipher(ctx, &xl, &xr);

  write32be(dst + 0, xl);
  write32be(dst + 4, xr);
}

void
blowfish_decrypt(blowfish_t *ctx,
                 unsigned char *dst,
                 const unsigned char *src) {
  uint32_t xl = read32be(src + 0);
  uint32_t xr = read32be(src + 4);

  blowfish_decipher(ctx, &xl, &xr);

  write32be(dst + 0, xl);
  write32be(dst + 4, xr);
}

/*
 * Camellia
 *
 * Resources:
 *   https://en.wikipedia.org/wiki/Camellia_(cipher)
 *   https://tools.ietf.org/html/rfc3713
 *   https://github.com/aead/camellia/blob/master/camellia.go
 */

static const uint32_t camellia_SIGMA[12] = {
  0xa09e667f, 0x3bcc908b, 0xb67ae858, 0x4caa73b2,
  0xc6ef372f, 0xe94f82be, 0x54ff53a5, 0xf1d36f1c,
  0x10e527fa, 0xde682d1d, 0xb05688c2, 0xb3e6c1fd
};

static const uint32_t camellia_S1[256] = {
  0x70707000, 0x82828200, 0x2c2c2c00, 0xececec00,
  0xb3b3b300, 0x27272700, 0xc0c0c000, 0xe5e5e500,
  0xe4e4e400, 0x85858500, 0x57575700, 0x35353500,
  0xeaeaea00, 0x0c0c0c00, 0xaeaeae00, 0x41414100,
  0x23232300, 0xefefef00, 0x6b6b6b00, 0x93939300,
  0x45454500, 0x19191900, 0xa5a5a500, 0x21212100,
  0xededed00, 0x0e0e0e00, 0x4f4f4f00, 0x4e4e4e00,
  0x1d1d1d00, 0x65656500, 0x92929200, 0xbdbdbd00,
  0x86868600, 0xb8b8b800, 0xafafaf00, 0x8f8f8f00,
  0x7c7c7c00, 0xebebeb00, 0x1f1f1f00, 0xcecece00,
  0x3e3e3e00, 0x30303000, 0xdcdcdc00, 0x5f5f5f00,
  0x5e5e5e00, 0xc5c5c500, 0x0b0b0b00, 0x1a1a1a00,
  0xa6a6a600, 0xe1e1e100, 0x39393900, 0xcacaca00,
  0xd5d5d500, 0x47474700, 0x5d5d5d00, 0x3d3d3d00,
  0xd9d9d900, 0x01010100, 0x5a5a5a00, 0xd6d6d600,
  0x51515100, 0x56565600, 0x6c6c6c00, 0x4d4d4d00,
  0x8b8b8b00, 0x0d0d0d00, 0x9a9a9a00, 0x66666600,
  0xfbfbfb00, 0xcccccc00, 0xb0b0b000, 0x2d2d2d00,
  0x74747400, 0x12121200, 0x2b2b2b00, 0x20202000,
  0xf0f0f000, 0xb1b1b100, 0x84848400, 0x99999900,
  0xdfdfdf00, 0x4c4c4c00, 0xcbcbcb00, 0xc2c2c200,
  0x34343400, 0x7e7e7e00, 0x76767600, 0x05050500,
  0x6d6d6d00, 0xb7b7b700, 0xa9a9a900, 0x31313100,
  0xd1d1d100, 0x17171700, 0x04040400, 0xd7d7d700,
  0x14141400, 0x58585800, 0x3a3a3a00, 0x61616100,
  0xdedede00, 0x1b1b1b00, 0x11111100, 0x1c1c1c00,
  0x32323200, 0x0f0f0f00, 0x9c9c9c00, 0x16161600,
  0x53535300, 0x18181800, 0xf2f2f200, 0x22222200,
  0xfefefe00, 0x44444400, 0xcfcfcf00, 0xb2b2b200,
  0xc3c3c300, 0xb5b5b500, 0x7a7a7a00, 0x91919100,
  0x24242400, 0x08080800, 0xe8e8e800, 0xa8a8a800,
  0x60606000, 0xfcfcfc00, 0x69696900, 0x50505000,
  0xaaaaaa00, 0xd0d0d000, 0xa0a0a000, 0x7d7d7d00,
  0xa1a1a100, 0x89898900, 0x62626200, 0x97979700,
  0x54545400, 0x5b5b5b00, 0x1e1e1e00, 0x95959500,
  0xe0e0e000, 0xffffff00, 0x64646400, 0xd2d2d200,
  0x10101000, 0xc4c4c400, 0x00000000, 0x48484800,
  0xa3a3a300, 0xf7f7f700, 0x75757500, 0xdbdbdb00,
  0x8a8a8a00, 0x03030300, 0xe6e6e600, 0xdadada00,
  0x09090900, 0x3f3f3f00, 0xdddddd00, 0x94949400,
  0x87878700, 0x5c5c5c00, 0x83838300, 0x02020200,
  0xcdcdcd00, 0x4a4a4a00, 0x90909000, 0x33333300,
  0x73737300, 0x67676700, 0xf6f6f600, 0xf3f3f300,
  0x9d9d9d00, 0x7f7f7f00, 0xbfbfbf00, 0xe2e2e200,
  0x52525200, 0x9b9b9b00, 0xd8d8d800, 0x26262600,
  0xc8c8c800, 0x37373700, 0xc6c6c600, 0x3b3b3b00,
  0x81818100, 0x96969600, 0x6f6f6f00, 0x4b4b4b00,
  0x13131300, 0xbebebe00, 0x63636300, 0x2e2e2e00,
  0xe9e9e900, 0x79797900, 0xa7a7a700, 0x8c8c8c00,
  0x9f9f9f00, 0x6e6e6e00, 0xbcbcbc00, 0x8e8e8e00,
  0x29292900, 0xf5f5f500, 0xf9f9f900, 0xb6b6b600,
  0x2f2f2f00, 0xfdfdfd00, 0xb4b4b400, 0x59595900,
  0x78787800, 0x98989800, 0x06060600, 0x6a6a6a00,
  0xe7e7e700, 0x46464600, 0x71717100, 0xbababa00,
  0xd4d4d400, 0x25252500, 0xababab00, 0x42424200,
  0x88888800, 0xa2a2a200, 0x8d8d8d00, 0xfafafa00,
  0x72727200, 0x07070700, 0xb9b9b900, 0x55555500,
  0xf8f8f800, 0xeeeeee00, 0xacacac00, 0x0a0a0a00,
  0x36363600, 0x49494900, 0x2a2a2a00, 0x68686800,
  0x3c3c3c00, 0x38383800, 0xf1f1f100, 0xa4a4a400,
  0x40404000, 0x28282800, 0xd3d3d300, 0x7b7b7b00,
  0xbbbbbb00, 0xc9c9c900, 0x43434300, 0xc1c1c100,
  0x15151500, 0xe3e3e300, 0xadadad00, 0xf4f4f400,
  0x77777700, 0xc7c7c700, 0x80808000, 0x9e9e9e00
};

static const uint32_t camellia_S2[256] = {
  0x00e0e0e0, 0x00050505, 0x00585858, 0x00d9d9d9,
  0x00676767, 0x004e4e4e, 0x00818181, 0x00cbcbcb,
  0x00c9c9c9, 0x000b0b0b, 0x00aeaeae, 0x006a6a6a,
  0x00d5d5d5, 0x00181818, 0x005d5d5d, 0x00828282,
  0x00464646, 0x00dfdfdf, 0x00d6d6d6, 0x00272727,
  0x008a8a8a, 0x00323232, 0x004b4b4b, 0x00424242,
  0x00dbdbdb, 0x001c1c1c, 0x009e9e9e, 0x009c9c9c,
  0x003a3a3a, 0x00cacaca, 0x00252525, 0x007b7b7b,
  0x000d0d0d, 0x00717171, 0x005f5f5f, 0x001f1f1f,
  0x00f8f8f8, 0x00d7d7d7, 0x003e3e3e, 0x009d9d9d,
  0x007c7c7c, 0x00606060, 0x00b9b9b9, 0x00bebebe,
  0x00bcbcbc, 0x008b8b8b, 0x00161616, 0x00343434,
  0x004d4d4d, 0x00c3c3c3, 0x00727272, 0x00959595,
  0x00ababab, 0x008e8e8e, 0x00bababa, 0x007a7a7a,
  0x00b3b3b3, 0x00020202, 0x00b4b4b4, 0x00adadad,
  0x00a2a2a2, 0x00acacac, 0x00d8d8d8, 0x009a9a9a,
  0x00171717, 0x001a1a1a, 0x00353535, 0x00cccccc,
  0x00f7f7f7, 0x00999999, 0x00616161, 0x005a5a5a,
  0x00e8e8e8, 0x00242424, 0x00565656, 0x00404040,
  0x00e1e1e1, 0x00636363, 0x00090909, 0x00333333,
  0x00bfbfbf, 0x00989898, 0x00979797, 0x00858585,
  0x00686868, 0x00fcfcfc, 0x00ececec, 0x000a0a0a,
  0x00dadada, 0x006f6f6f, 0x00535353, 0x00626262,
  0x00a3a3a3, 0x002e2e2e, 0x00080808, 0x00afafaf,
  0x00282828, 0x00b0b0b0, 0x00747474, 0x00c2c2c2,
  0x00bdbdbd, 0x00363636, 0x00222222, 0x00383838,
  0x00646464, 0x001e1e1e, 0x00393939, 0x002c2c2c,
  0x00a6a6a6, 0x00303030, 0x00e5e5e5, 0x00444444,
  0x00fdfdfd, 0x00888888, 0x009f9f9f, 0x00656565,
  0x00878787, 0x006b6b6b, 0x00f4f4f4, 0x00232323,
  0x00484848, 0x00101010, 0x00d1d1d1, 0x00515151,
  0x00c0c0c0, 0x00f9f9f9, 0x00d2d2d2, 0x00a0a0a0,
  0x00555555, 0x00a1a1a1, 0x00414141, 0x00fafafa,
  0x00434343, 0x00131313, 0x00c4c4c4, 0x002f2f2f,
  0x00a8a8a8, 0x00b6b6b6, 0x003c3c3c, 0x002b2b2b,
  0x00c1c1c1, 0x00ffffff, 0x00c8c8c8, 0x00a5a5a5,
  0x00202020, 0x00898989, 0x00000000, 0x00909090,
  0x00474747, 0x00efefef, 0x00eaeaea, 0x00b7b7b7,
  0x00151515, 0x00060606, 0x00cdcdcd, 0x00b5b5b5,
  0x00121212, 0x007e7e7e, 0x00bbbbbb, 0x00292929,
  0x000f0f0f, 0x00b8b8b8, 0x00070707, 0x00040404,
  0x009b9b9b, 0x00949494, 0x00212121, 0x00666666,
  0x00e6e6e6, 0x00cecece, 0x00ededed, 0x00e7e7e7,
  0x003b3b3b, 0x00fefefe, 0x007f7f7f, 0x00c5c5c5,
  0x00a4a4a4, 0x00373737, 0x00b1b1b1, 0x004c4c4c,
  0x00919191, 0x006e6e6e, 0x008d8d8d, 0x00767676,
  0x00030303, 0x002d2d2d, 0x00dedede, 0x00969696,
  0x00262626, 0x007d7d7d, 0x00c6c6c6, 0x005c5c5c,
  0x00d3d3d3, 0x00f2f2f2, 0x004f4f4f, 0x00191919,
  0x003f3f3f, 0x00dcdcdc, 0x00797979, 0x001d1d1d,
  0x00525252, 0x00ebebeb, 0x00f3f3f3, 0x006d6d6d,
  0x005e5e5e, 0x00fbfbfb, 0x00696969, 0x00b2b2b2,
  0x00f0f0f0, 0x00313131, 0x000c0c0c, 0x00d4d4d4,
  0x00cfcfcf, 0x008c8c8c, 0x00e2e2e2, 0x00757575,
  0x00a9a9a9, 0x004a4a4a, 0x00575757, 0x00848484,
  0x00111111, 0x00454545, 0x001b1b1b, 0x00f5f5f5,
  0x00e4e4e4, 0x000e0e0e, 0x00737373, 0x00aaaaaa,
  0x00f1f1f1, 0x00dddddd, 0x00595959, 0x00141414,
  0x006c6c6c, 0x00929292, 0x00545454, 0x00d0d0d0,
  0x00787878, 0x00707070, 0x00e3e3e3, 0x00494949,
  0x00808080, 0x00505050, 0x00a7a7a7, 0x00f6f6f6,
  0x00777777, 0x00939393, 0x00868686, 0x00838383,
  0x002a2a2a, 0x00c7c7c7, 0x005b5b5b, 0x00e9e9e9,
  0x00eeeeee, 0x008f8f8f, 0x00010101, 0x003d3d3d
};

static const uint32_t camellia_S3[256] = {
  0x38003838, 0x41004141, 0x16001616, 0x76007676,
  0xd900d9d9, 0x93009393, 0x60006060, 0xf200f2f2,
  0x72007272, 0xc200c2c2, 0xab00abab, 0x9a009a9a,
  0x75007575, 0x06000606, 0x57005757, 0xa000a0a0,
  0x91009191, 0xf700f7f7, 0xb500b5b5, 0xc900c9c9,
  0xa200a2a2, 0x8c008c8c, 0xd200d2d2, 0x90009090,
  0xf600f6f6, 0x07000707, 0xa700a7a7, 0x27002727,
  0x8e008e8e, 0xb200b2b2, 0x49004949, 0xde00dede,
  0x43004343, 0x5c005c5c, 0xd700d7d7, 0xc700c7c7,
  0x3e003e3e, 0xf500f5f5, 0x8f008f8f, 0x67006767,
  0x1f001f1f, 0x18001818, 0x6e006e6e, 0xaf00afaf,
  0x2f002f2f, 0xe200e2e2, 0x85008585, 0x0d000d0d,
  0x53005353, 0xf000f0f0, 0x9c009c9c, 0x65006565,
  0xea00eaea, 0xa300a3a3, 0xae00aeae, 0x9e009e9e,
  0xec00ecec, 0x80008080, 0x2d002d2d, 0x6b006b6b,
  0xa800a8a8, 0x2b002b2b, 0x36003636, 0xa600a6a6,
  0xc500c5c5, 0x86008686, 0x4d004d4d, 0x33003333,
  0xfd00fdfd, 0x66006666, 0x58005858, 0x96009696,
  0x3a003a3a, 0x09000909, 0x95009595, 0x10001010,
  0x78007878, 0xd800d8d8, 0x42004242, 0xcc00cccc,
  0xef00efef, 0x26002626, 0xe500e5e5, 0x61006161,
  0x1a001a1a, 0x3f003f3f, 0x3b003b3b, 0x82008282,
  0xb600b6b6, 0xdb00dbdb, 0xd400d4d4, 0x98009898,
  0xe800e8e8, 0x8b008b8b, 0x02000202, 0xeb00ebeb,
  0x0a000a0a, 0x2c002c2c, 0x1d001d1d, 0xb000b0b0,
  0x6f006f6f, 0x8d008d8d, 0x88008888, 0x0e000e0e,
  0x19001919, 0x87008787, 0x4e004e4e, 0x0b000b0b,
  0xa900a9a9, 0x0c000c0c, 0x79007979, 0x11001111,
  0x7f007f7f, 0x22002222, 0xe700e7e7, 0x59005959,
  0xe100e1e1, 0xda00dada, 0x3d003d3d, 0xc800c8c8,
  0x12001212, 0x04000404, 0x74007474, 0x54005454,
  0x30003030, 0x7e007e7e, 0xb400b4b4, 0x28002828,
  0x55005555, 0x68006868, 0x50005050, 0xbe00bebe,
  0xd000d0d0, 0xc400c4c4, 0x31003131, 0xcb00cbcb,
  0x2a002a2a, 0xad00adad, 0x0f000f0f, 0xca00caca,
  0x70007070, 0xff00ffff, 0x32003232, 0x69006969,
  0x08000808, 0x62006262, 0x00000000, 0x24002424,
  0xd100d1d1, 0xfb00fbfb, 0xba00baba, 0xed00eded,
  0x45004545, 0x81008181, 0x73007373, 0x6d006d6d,
  0x84008484, 0x9f009f9f, 0xee00eeee, 0x4a004a4a,
  0xc300c3c3, 0x2e002e2e, 0xc100c1c1, 0x01000101,
  0xe600e6e6, 0x25002525, 0x48004848, 0x99009999,
  0xb900b9b9, 0xb300b3b3, 0x7b007b7b, 0xf900f9f9,
  0xce00cece, 0xbf00bfbf, 0xdf00dfdf, 0x71007171,
  0x29002929, 0xcd00cdcd, 0x6c006c6c, 0x13001313,
  0x64006464, 0x9b009b9b, 0x63006363, 0x9d009d9d,
  0xc000c0c0, 0x4b004b4b, 0xb700b7b7, 0xa500a5a5,
  0x89008989, 0x5f005f5f, 0xb100b1b1, 0x17001717,
  0xf400f4f4, 0xbc00bcbc, 0xd300d3d3, 0x46004646,
  0xcf00cfcf, 0x37003737, 0x5e005e5e, 0x47004747,
  0x94009494, 0xfa00fafa, 0xfc00fcfc, 0x5b005b5b,
  0x97009797, 0xfe00fefe, 0x5a005a5a, 0xac00acac,
  0x3c003c3c, 0x4c004c4c, 0x03000303, 0x35003535,
  0xf300f3f3, 0x23002323, 0xb800b8b8, 0x5d005d5d,
  0x6a006a6a, 0x92009292, 0xd500d5d5, 0x21002121,
  0x44004444, 0x51005151, 0xc600c6c6, 0x7d007d7d,
  0x39003939, 0x83008383, 0xdc00dcdc, 0xaa00aaaa,
  0x7c007c7c, 0x77007777, 0x56005656, 0x05000505,
  0x1b001b1b, 0xa400a4a4, 0x15001515, 0x34003434,
  0x1e001e1e, 0x1c001c1c, 0xf800f8f8, 0x52005252,
  0x20002020, 0x14001414, 0xe900e9e9, 0xbd00bdbd,
  0xdd00dddd, 0xe400e4e4, 0xa100a1a1, 0xe000e0e0,
  0x8a008a8a, 0xf100f1f1, 0xd600d6d6, 0x7a007a7a,
  0xbb00bbbb, 0xe300e3e3, 0x40004040, 0x4f004f4f
};

static const uint32_t camellia_S4[256] = {
  0x70700070, 0x2c2c002c, 0xb3b300b3, 0xc0c000c0,
  0xe4e400e4, 0x57570057, 0xeaea00ea, 0xaeae00ae,
  0x23230023, 0x6b6b006b, 0x45450045, 0xa5a500a5,
  0xeded00ed, 0x4f4f004f, 0x1d1d001d, 0x92920092,
  0x86860086, 0xafaf00af, 0x7c7c007c, 0x1f1f001f,
  0x3e3e003e, 0xdcdc00dc, 0x5e5e005e, 0x0b0b000b,
  0xa6a600a6, 0x39390039, 0xd5d500d5, 0x5d5d005d,
  0xd9d900d9, 0x5a5a005a, 0x51510051, 0x6c6c006c,
  0x8b8b008b, 0x9a9a009a, 0xfbfb00fb, 0xb0b000b0,
  0x74740074, 0x2b2b002b, 0xf0f000f0, 0x84840084,
  0xdfdf00df, 0xcbcb00cb, 0x34340034, 0x76760076,
  0x6d6d006d, 0xa9a900a9, 0xd1d100d1, 0x04040004,
  0x14140014, 0x3a3a003a, 0xdede00de, 0x11110011,
  0x32320032, 0x9c9c009c, 0x53530053, 0xf2f200f2,
  0xfefe00fe, 0xcfcf00cf, 0xc3c300c3, 0x7a7a007a,
  0x24240024, 0xe8e800e8, 0x60600060, 0x69690069,
  0xaaaa00aa, 0xa0a000a0, 0xa1a100a1, 0x62620062,
  0x54540054, 0x1e1e001e, 0xe0e000e0, 0x64640064,
  0x10100010, 0x00000000, 0xa3a300a3, 0x75750075,
  0x8a8a008a, 0xe6e600e6, 0x09090009, 0xdddd00dd,
  0x87870087, 0x83830083, 0xcdcd00cd, 0x90900090,
  0x73730073, 0xf6f600f6, 0x9d9d009d, 0xbfbf00bf,
  0x52520052, 0xd8d800d8, 0xc8c800c8, 0xc6c600c6,
  0x81810081, 0x6f6f006f, 0x13130013, 0x63630063,
  0xe9e900e9, 0xa7a700a7, 0x9f9f009f, 0xbcbc00bc,
  0x29290029, 0xf9f900f9, 0x2f2f002f, 0xb4b400b4,
  0x78780078, 0x06060006, 0xe7e700e7, 0x71710071,
  0xd4d400d4, 0xabab00ab, 0x88880088, 0x8d8d008d,
  0x72720072, 0xb9b900b9, 0xf8f800f8, 0xacac00ac,
  0x36360036, 0x2a2a002a, 0x3c3c003c, 0xf1f100f1,
  0x40400040, 0xd3d300d3, 0xbbbb00bb, 0x43430043,
  0x15150015, 0xadad00ad, 0x77770077, 0x80800080,
  0x82820082, 0xecec00ec, 0x27270027, 0xe5e500e5,
  0x85850085, 0x35350035, 0x0c0c000c, 0x41410041,
  0xefef00ef, 0x93930093, 0x19190019, 0x21210021,
  0x0e0e000e, 0x4e4e004e, 0x65650065, 0xbdbd00bd,
  0xb8b800b8, 0x8f8f008f, 0xebeb00eb, 0xcece00ce,
  0x30300030, 0x5f5f005f, 0xc5c500c5, 0x1a1a001a,
  0xe1e100e1, 0xcaca00ca, 0x47470047, 0x3d3d003d,
  0x01010001, 0xd6d600d6, 0x56560056, 0x4d4d004d,
  0x0d0d000d, 0x66660066, 0xcccc00cc, 0x2d2d002d,
  0x12120012, 0x20200020, 0xb1b100b1, 0x99990099,
  0x4c4c004c, 0xc2c200c2, 0x7e7e007e, 0x05050005,
  0xb7b700b7, 0x31310031, 0x17170017, 0xd7d700d7,
  0x58580058, 0x61610061, 0x1b1b001b, 0x1c1c001c,
  0x0f0f000f, 0x16160016, 0x18180018, 0x22220022,
  0x44440044, 0xb2b200b2, 0xb5b500b5, 0x91910091,
  0x08080008, 0xa8a800a8, 0xfcfc00fc, 0x50500050,
  0xd0d000d0, 0x7d7d007d, 0x89890089, 0x97970097,
  0x5b5b005b, 0x95950095, 0xffff00ff, 0xd2d200d2,
  0xc4c400c4, 0x48480048, 0xf7f700f7, 0xdbdb00db,
  0x03030003, 0xdada00da, 0x3f3f003f, 0x94940094,
  0x5c5c005c, 0x02020002, 0x4a4a004a, 0x33330033,
  0x67670067, 0xf3f300f3, 0x7f7f007f, 0xe2e200e2,
  0x9b9b009b, 0x26260026, 0x37370037, 0x3b3b003b,
  0x96960096, 0x4b4b004b, 0xbebe00be, 0x2e2e002e,
  0x79790079, 0x8c8c008c, 0x6e6e006e, 0x8e8e008e,
  0xf5f500f5, 0xb6b600b6, 0xfdfd00fd, 0x59590059,
  0x98980098, 0x6a6a006a, 0x46460046, 0xbaba00ba,
  0x25250025, 0x42420042, 0xa2a200a2, 0xfafa00fa,
  0x07070007, 0x55550055, 0xeeee00ee, 0x0a0a000a,
  0x49490049, 0x68680068, 0x38380038, 0xa4a400a4,
  0x28280028, 0x7b7b007b, 0xc9c900c9, 0xc1c100c1,
  0xe3e300e3, 0xf4f400f4, 0xc7c700c7, 0x9e9e009e
};

#define SIGMA camellia_SIGMA
#define S1 camellia_S1
#define S2 camellia_S2
#define S3 camellia_S3
#define S4 camellia_S4

#define FEIS(r, i0, i1, i2, i3, k0, k1) do { \
  uint32_t t0, t1, z;                        \
                                             \
  t0 = k0 ^ r[i0];                           \
  t1 = k1 ^ r[i1];                           \
                                             \
  z = S4[(t0 >>  0) & 0xff]                  \
    ^ S3[(t0 >>  8) & 0xff]                  \
    ^ S2[(t0 >> 16) & 0xff]                  \
    ^ S1[(t0 >> 24) & 0xff];                 \
                                             \
  r[i3] ^= (z >> 8) | (z << (32 - 8));       \
                                             \
  t0 = z ^ S1[(t1 >>  0) & 0xff]             \
         ^ S4[(t1 >>  8) & 0xff]             \
         ^ S3[(t1 >> 16) & 0xff]             \
         ^ S2[(t1 >> 24) & 0xff];            \
                                             \
  r[i2] ^= t0;                               \
  r[i3] ^= t0;                               \
} while (0)

#define ROTL(r, i0, i1, i2, i3, n) do {       \
  uint32_t t = r[i0] >> (32 - n);             \
                                              \
  r[i0] = (r[i0] << n) | (r[i1] >> (32 - n)); \
  r[i1] = (r[i1] << n) | (r[i2] >> (32 - n)); \
  r[i2] = (r[i2] << n) | (r[i3] >> (32 - n)); \
  r[i3] = (r[i3] << n) | t;                   \
} while (0)

#define SET2(x, x0, x1, y, y0, y1) do { \
  x[x0] = y[y0];                        \
  x[x1] = y[y1];                        \
} while (0)

#define SET4(x, x0, x1, x2, x3, y, y0, y1, y2, y3) do { \
  x[x0] = y[y0];                                        \
  x[x1] = y[y1];                                        \
  x[x2] = y[y2];                                        \
  x[x3] = y[y3];                                        \
} while (0)

#define XOR4(x, x0, x1, x2, x3, y, y0, y1, y2, y3) do { \
  x[x0] ^= y[y0];                                       \
  x[x1] ^= y[y1];                                       \
  x[x2] ^= y[y2];                                       \
  x[x3] ^= y[y3];                                       \
} while (0)

static void
camellia128_init(camellia_t *ctx, const unsigned char *key) {
  uint32_t *k = ctx->key;
  uint32_t s[4];

  memset(ctx, 0, sizeof(*ctx));

  s[0] = read32be(key +  0);
  s[1] = read32be(key +  4);
  s[2] = read32be(key +  8);
  s[3] = read32be(key + 12);

  SET4(k, 0, 1, 2, 3, s, 0, 1, 2, 3);

  FEIS(s, 0, 1, 2, 3, SIGMA[0], SIGMA[1]);
  FEIS(s, 2, 3, 0, 1, SIGMA[2], SIGMA[3]);

  XOR4(s, 0, 1, 2, 3, k, 0, 1, 2, 3);

  FEIS(s, 0, 1, 2, 3, SIGMA[4], SIGMA[5]);
  FEIS(s, 2, 3, 0, 1, SIGMA[6], SIGMA[7]);
  SET4(k, 4, 5, 6, 7, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 15); /* KA << 15 */
  SET4(k, 12, 13, 14, 15, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 15); /* KA << 30 */
  SET4(k, 16, 17, 18, 19, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 15); /* KA << 45 */
  SET2(k, 24, 25, s, 0, 1);

  ROTL(s, 0, 1, 2, 3, 15); /* KA << 60 */
  SET4(k, 28, 29, 30, 31, s, 0, 1, 2, 3);

  ROTL(s, 1, 2, 3, 0, 2); /* KA << 94 */
  SET4(k, 40, 41, 42, 43, s, 1, 2, 3, 0);

  ROTL(s, 1, 2, 3, 0, 17); /* KA << 111 */
  SET4(k, 48, 49, 50, 51, s, 1, 2, 3, 0);
  SET4(s, 0, 1, 2, 3, k, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 15); /* KL << 15 */
  SET4(k, 8, 9, 10, 11, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 30); /* KL << 45 */
  SET4(k, 20, 21, 22, 23, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 15); /* KL << 60 */
  SET2(k, 26, 27, s, 2, 3);

  ROTL(s, 0, 1, 2, 3, 17); /* KL << 77 */
  SET4(k, 32, 33, 34, 35, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 17); /* KL << 94 */
  SET4(k, 36, 37, 38, 39, s, 0, 1, 2, 3);

  ROTL(s, 0, 1, 2, 3, 17); /* KL << 111 */
  SET4(k, 44, 45, 46, 47, s, 0, 1, 2, 3);
}

static void
camellia128_encrypt(camellia_t *ctx,
                    unsigned char *dst,
                    const unsigned char *src) {
  uint32_t *k = ctx->key;
  uint32_t r[4];
  uint32_t t;

  r[0] = read32be(src +  0);
  r[1] = read32be(src +  4);
  r[2] = read32be(src +  8);
  r[3] = read32be(src + 12);

  XOR4(r, 0, 1, 2, 3, k, 0, 1, 2, 3);

  FEIS(r, 0, 1, 2, 3, k[4], k[5]);
  FEIS(r, 2, 3, 0, 1, k[6], k[7]);
  FEIS(r, 0, 1, 2, 3, k[8], k[9]);
  FEIS(r, 2, 3, 0, 1, k[10], k[11]);
  FEIS(r, 0, 1, 2, 3, k[12], k[13]);
  FEIS(r, 2, 3, 0, 1, k[14], k[15]);

  t = r[0] & k[16];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[19];
  r[0] ^= r[1] | k[17];
  t = r[2] & k[18];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[20], k[21]);
  FEIS(r, 2, 3, 0, 1, k[22], k[23]);
  FEIS(r, 0, 1, 2, 3, k[24], k[25]);
  FEIS(r, 2, 3, 0, 1, k[26], k[27]);
  FEIS(r, 0, 1, 2, 3, k[28], k[29]);
  FEIS(r, 2, 3, 0, 1, k[30], k[31]);

  t = r[0] & k[32];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[35];
  r[0] ^= r[1] | k[33];
  t = r[2] & k[34];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[36], k[37]);
  FEIS(r, 2, 3, 0, 1, k[38], k[39]);
  FEIS(r, 0, 1, 2, 3, k[40], k[41]);
  FEIS(r, 2, 3, 0, 1, k[42], k[43]);
  FEIS(r, 0, 1, 2, 3, k[44], k[45]);
  FEIS(r, 2, 3, 0, 1, k[46], k[47]);

  XOR4(r, 2, 3, 0, 1, k, 48, 49, 50, 51);

  write32be(dst +  0, r[2]);
  write32be(dst +  4, r[3]);
  write32be(dst +  8, r[0]);
  write32be(dst + 12, r[1]);
}

static void
camellia128_decrypt(camellia_t *ctx,
                    unsigned char *dst,
                    const unsigned char *src) {
  uint32_t *k = ctx->key;
  uint32_t r[4];
  uint32_t t;

  r[0] = read32be(src +  0);
  r[1] = read32be(src +  4);
  r[2] = read32be(src +  8);
  r[3] = read32be(src + 12);

  XOR4(r, 3, 2, 1, 0, k, 51, 50, 49, 48);

  FEIS(r, 0, 1, 2, 3, k[46], k[47]);
  FEIS(r, 2, 3, 0, 1, k[44], k[45]);
  FEIS(r, 0, 1, 2, 3, k[42], k[43]);
  FEIS(r, 2, 3, 0, 1, k[40], k[41]);
  FEIS(r, 0, 1, 2, 3, k[38], k[39]);
  FEIS(r, 2, 3, 0, 1, k[36], k[37]);

  t = r[0] & k[34];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[33];
  r[0] ^= r[1] | k[35];
  t = r[2] & k[32];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[30], k[31]);
  FEIS(r, 2, 3, 0, 1, k[28], k[29]);
  FEIS(r, 0, 1, 2, 3, k[26], k[27]);
  FEIS(r, 2, 3, 0, 1, k[24], k[25]);
  FEIS(r, 0, 1, 2, 3, k[22], k[23]);
  FEIS(r, 2, 3, 0, 1, k[20], k[21]);

  t = r[0] & k[18];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[17];
  r[0] ^= r[1] | k[19];
  t = r[2] & k[16];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[14], k[15]);
  FEIS(r, 2, 3, 0, 1, k[12], k[13]);
  FEIS(r, 0, 1, 2, 3, k[10], k[11]);
  FEIS(r, 2, 3, 0, 1, k[8], k[9]);
  FEIS(r, 0, 1, 2, 3, k[6], k[7]);
  FEIS(r, 2, 3, 0, 1, k[4], k[5]);

  XOR4(r, 1, 0, 3, 2, k, 3, 2, 1, 0);

  write32be(dst +  0, r[2]);
  write32be(dst +  4, r[3]);
  write32be(dst +  8, r[0]);
  write32be(dst + 12, r[1]);
}

static void
camellia256_init(camellia_t *ctx, const unsigned char *key, size_t key_len) {
  uint32_t *k = ctx->key;
  uint32_t s[4];

  k[0] = read32be(key + 0);
  k[1] = read32be(key + 4);
  k[2] = read32be(key + 8);
  k[3] = read32be(key + 12);

  k[8] = read32be(key + 16);
  k[9] = read32be(key + 20);

  if (key_len == 24) {
    k[10] = ~k[8];
    k[11] = ~k[9];
  } else if (key_len == 32) {
    k[10] = read32be(key + 24);
    k[11] = read32be(key + 28);
  } else {
    ASSERT(0);
  }

  s[0] = k[8] ^ k[0];
  s[1] = k[9] ^ k[1];
  s[2] = k[10] ^ k[2];
  s[3] = k[11] ^ k[3];

  FEIS(s, 0, 1, 2, 3, SIGMA[0], SIGMA[1]);
  FEIS(s, 2, 3, 0, 1, SIGMA[2], SIGMA[3]);

  XOR4(s, 0, 1, 2, 3, k, 0, 1, 2, 3);
  FEIS(s, 0, 1, 2, 3, SIGMA[4], SIGMA[5]);
  FEIS(s, 2, 3, 0, 1, SIGMA[6], SIGMA[7]);

  SET4(k, 12, 13, 14, 15, s, 0, 1, 2, 3);

  XOR4(s, 0, 1, 2, 3, k, 8, 9, 10, 11);
  FEIS(s, 0, 1, 2, 3, SIGMA[8], SIGMA[9]);
  FEIS(s, 2, 3, 0, 1, SIGMA[10], SIGMA[11]);

  SET4(k, 4, 5, 6, 7, s, 0, 1, 2, 3);
  ROTL(s, 0, 1, 2, 3, 30); /* KB << 30 */
  SET4(k, 20, 21, 22, 23, s, 0, 1, 2, 3);
  ROTL(s, 0, 1, 2, 3, 30); /* KB << 60 */
  SET4(k, 40, 41, 42, 43, s, 0, 1, 2, 3);
  ROTL(s, 1, 2, 3, 0, 19); /* KB << 111 */
  SET4(k, 64, 65, 66, 67, s, 1, 2, 3, 0);

  SET4(s, 0, 1, 2, 3, k, 8, 9, 10, 11);
  ROTL(s, 0, 1, 2, 3, 15); /* KR << 15 */
  SET4(k, 8, 9, 10, 11, s, 0, 1, 2, 3);
  ROTL(s, 0, 1, 2, 3, 15); /* KR << 30 */
  SET4(k, 16, 17, 18, 19, s, 0, 1, 2, 3);
  ROTL(s, 0, 1, 2, 3, 30); /* KR << 60 */
  SET4(k, 36, 37, 38, 39, s, 0, 1, 2, 3);
  ROTL(s, 1, 2, 3, 0, 2); /* KR << 94 */
  SET4(k, 52, 53, 54, 55, s, 1, 2, 3, 0);

  SET4(s, 0, 1, 2, 3, k, 12, 13, 14, 15);
  ROTL(s, 0, 1, 2, 3, 15); /* KA << 15 */
  SET4(k, 12, 13, 14, 15, s, 0, 1, 2, 3);
  ROTL(s, 0, 1, 2, 3, 30); /* KA << 45 */
  SET4(k, 28, 29, 30, 31, s, 0, 1, 2, 3);
  /* KA << 77 */
  SET4(k, 48, 49, 50, 51, s, 1, 2, 3, 0);
  ROTL(s, 1, 2, 3, 0, 17); /* KA << 94 */
  SET4(k, 56, 57, 58, 59, s, 1, 2, 3, 0);

  SET4(s, 0, 1, 2, 3, k, 0, 1, 2, 3);
  ROTL(s, 1, 2, 3, 0, 13); /* KL << 45 */
  SET4(k, 24, 25, 26, 27, s, 1, 2, 3, 0);
  ROTL(s, 1, 2, 3, 0, 15); /* KL << 60 */
  SET4(k, 32, 33, 34, 35, s, 1, 2, 3, 0);
  ROTL(s, 1, 2, 3, 0, 17); /* KL << 77 */
  SET4(k, 44, 45, 46, 47, s, 1, 2, 3, 0);
  ROTL(s, 2, 3, 0, 1, 2); /* KL << 111 */
  SET4(k, 60, 61, 62, 63, s, 2, 3, 0, 1);
}

static void
camellia256_encrypt(camellia_t *ctx,
                    unsigned char *dst,
                    const unsigned char *src) {
  uint32_t *k = ctx->key;
  uint32_t r[4];
  uint32_t t;

  r[0] = read32be(src +  0);
  r[1] = read32be(src +  4);
  r[2] = read32be(src +  8);
  r[3] = read32be(src + 12);

  XOR4(r, 0, 1, 2, 3, k, 0, 1, 2, 3);

  FEIS(r, 0, 1, 2, 3, k[4], k[5]);
  FEIS(r, 2, 3, 0, 1, k[6], k[7]);
  FEIS(r, 0, 1, 2, 3, k[8], k[9]);
  FEIS(r, 2, 3, 0, 1, k[10], k[11]);
  FEIS(r, 0, 1, 2, 3, k[12], k[13]);
  FEIS(r, 2, 3, 0, 1, k[14], k[15]);

  t = r[0] & k[16];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[19];
  r[0] ^= r[1] | k[17];
  t = r[2] & k[18];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[20], k[21]);
  FEIS(r, 2, 3, 0, 1, k[22], k[23]);
  FEIS(r, 0, 1, 2, 3, k[24], k[25]);
  FEIS(r, 2, 3, 0, 1, k[26], k[27]);
  FEIS(r, 0, 1, 2, 3, k[28], k[29]);
  FEIS(r, 2, 3, 0, 1, k[30], k[31]);

  t = r[0] & k[32];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[35];
  r[0] ^= r[1] | k[33];
  t = r[2] & k[34];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[36], k[37]);
  FEIS(r, 2, 3, 0, 1, k[38], k[39]);
  FEIS(r, 0, 1, 2, 3, k[40], k[41]);
  FEIS(r, 2, 3, 0, 1, k[42], k[43]);
  FEIS(r, 0, 1, 2, 3, k[44], k[45]);
  FEIS(r, 2, 3, 0, 1, k[46], k[47]);

  t = r[0] & k[48];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[51];
  r[0] ^= r[1] | k[49];
  t = r[2] & k[50];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[52], k[53]);
  FEIS(r, 2, 3, 0, 1, k[54], k[55]);
  FEIS(r, 0, 1, 2, 3, k[56], k[57]);
  FEIS(r, 2, 3, 0, 1, k[58], k[59]);
  FEIS(r, 0, 1, 2, 3, k[60], k[61]);
  FEIS(r, 2, 3, 0, 1, k[62], k[63]);

  XOR4(r, 2, 3, 0, 1, k, 64, 65, 66, 67);

  write32be(dst +  0, r[2]);
  write32be(dst +  4, r[3]);
  write32be(dst +  8, r[0]);
  write32be(dst + 12, r[1]);
}

static void
camellia256_decrypt(camellia_t *ctx,
                    unsigned char *dst,
                    const unsigned char *src) {
  uint32_t *k = ctx->key;
  uint32_t r[4];
  uint32_t t;

  r[0] = read32be(src +  0);
  r[1] = read32be(src +  4);
  r[2] = read32be(src +  8);
  r[3] = read32be(src + 12);

  XOR4(r, 3, 2, 1, 0, k, 67, 66, 65, 64);

  FEIS(r, 0, 1, 2, 3, k[62], k[63]);
  FEIS(r, 2, 3, 0, 1, k[60], k[61]);
  FEIS(r, 0, 1, 2, 3, k[58], k[59]);
  FEIS(r, 2, 3, 0, 1, k[56], k[57]);
  FEIS(r, 0, 1, 2, 3, k[54], k[55]);
  FEIS(r, 2, 3, 0, 1, k[52], k[53]);

  t = r[0] & k[50];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[49];
  r[0] ^= r[1] | k[51];
  t = r[2] & k[48];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[46], k[47]);
  FEIS(r, 2, 3, 0, 1, k[44], k[45]);
  FEIS(r, 0, 1, 2, 3, k[42], k[43]);
  FEIS(r, 2, 3, 0, 1, k[40], k[41]);
  FEIS(r, 0, 1, 2, 3, k[38], k[39]);
  FEIS(r, 2, 3, 0, 1, k[36], k[37]);

  t = r[0] & k[34];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[33];
  r[0] ^= r[1] | k[35];
  t = r[2] & k[32];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[30], k[31]);
  FEIS(r, 2, 3, 0, 1, k[28], k[29]);
  FEIS(r, 0, 1, 2, 3, k[26], k[27]);
  FEIS(r, 2, 3, 0, 1, k[24], k[25]);
  FEIS(r, 0, 1, 2, 3, k[22], k[23]);
  FEIS(r, 2, 3, 0, 1, k[20], k[21]);

  t = r[0] & k[18];
  r[1] ^= (t << 1) | (t >> (32 - 1));
  r[2] ^= r[3] | k[17];
  r[0] ^= r[1] | k[19];
  t = r[2] & k[16];
  r[3] ^= (t << 1) | (t >> (32 - 1));

  FEIS(r, 0, 1, 2, 3, k[14], k[15]);
  FEIS(r, 2, 3, 0, 1, k[12], k[13]);
  FEIS(r, 0, 1, 2, 3, k[10], k[11]);
  FEIS(r, 2, 3, 0, 1, k[8], k[9]);
  FEIS(r, 0, 1, 2, 3, k[6], k[7]);
  FEIS(r, 2, 3, 0, 1, k[4], k[5]);

  XOR4(r, 1, 0, 3, 2, k, 3, 2, 1, 0);

  write32be(dst +  0, r[2]);
  write32be(dst +  4, r[3]);
  write32be(dst +  8, r[0]);
  write32be(dst + 12, r[1]);
}

void
camellia_init(camellia_t *ctx, unsigned int bits, const unsigned char *key) {
  ctx->bits = bits;

  switch (ctx->bits) {
    case 128:
      camellia128_init(ctx, key);
      break;
    case 192:
      camellia256_init(ctx, key, 24);
      break;
    case 256:
      camellia256_init(ctx, key, 32);
      break;
    default:
      ASSERT(0);
      break;
  }
}

void
camellia_encrypt(camellia_t *ctx,
                 unsigned char *dst,
                 const unsigned char *src) {
  if (ctx->bits == 128)
    camellia128_encrypt(ctx, dst, src);
  else
    camellia256_encrypt(ctx, dst, src);
}

void
camellia_decrypt(camellia_t *ctx,
                 unsigned char *dst,
                 const unsigned char *src) {
  if (ctx->bits == 128)
    camellia128_decrypt(ctx, dst, src);
  else
    camellia256_decrypt(ctx, dst, src);
}

#undef SIGMA
#undef S1
#undef S2
#undef S3
#undef S4
#undef FEIS
#undef ROTL
#undef SET2
#undef SET4
#undef XOR4

/*
 * Cipher
 */

static int
cipher_ctx_init(cipher_t *ctx, const unsigned char *key, size_t key_len) {
  switch (ctx->type) {
    case CIPHER_AES128: {
      if (key_len != 16)
        return 0;

      ctx->block_size = 16;

      if (ctx->encrypt)
        aes_init_encrypt(&ctx->ctx.aes, 128, key);
      else
        aes_init_decrypt(&ctx->ctx.aes, 128, key);

      break;
    }

    case CIPHER_AES192: {
      if (key_len != 24)
        return 0;

      ctx->block_size = 16;

      if (ctx->encrypt)
        aes_init_encrypt(&ctx->ctx.aes, 192, key);
      else
        aes_init_decrypt(&ctx->ctx.aes, 192, key);

      break;
    }

    case CIPHER_AES256: {
      if (key_len != 32)
        return 0;

      ctx->block_size = 16;

      if (ctx->encrypt)
        aes_init_encrypt(&ctx->ctx.aes, 256, key);
      else
        aes_init_decrypt(&ctx->ctx.aes, 256, key);

      break;
    }

    case CIPHER_BLOWFISH: {
      if (key_len < 1 || key_len > 72)
        return 0;

      ctx->block_size = 8;

      blowfish_init(&ctx->ctx.blowfish, key, key_len, NULL, 0);

      break;
    }

    case CIPHER_CAMELLIA128: {
      if (key_len != 16)
        return 0;

      ctx->block_size = 16;

      camellia_init(&ctx->ctx.camellia, 128, key);

      break;
    }

    case CIPHER_CAMELLIA192: {
      if (key_len != 24)
        return 0;

      ctx->block_size = 16;

      camellia_init(&ctx->ctx.camellia, 192, key);

      break;
    }

    case CIPHER_CAMELLIA256: {
      if (key_len != 32)
        return 0;

      ctx->block_size = 16;

      camellia_init(&ctx->ctx.camellia, 256, key);

      break;
    }

    default: {
      return 0;
    }
  }

  return 1;
}

static void
cipher_ctx_encrypt(cipher_t *ctx,
                   unsigned char *dst,
                   const unsigned char *src) {
  switch (ctx->type) {
    case CIPHER_AES128:
    case CIPHER_AES192:
    case CIPHER_AES256:
      aes_encrypt(&ctx->ctx.aes, dst, src);
      break;
    case CIPHER_BLOWFISH:
      blowfish_encrypt(&ctx->ctx.blowfish, dst, src);
      break;
    case CIPHER_CAMELLIA128:
    case CIPHER_CAMELLIA192:
    case CIPHER_CAMELLIA256:
      camellia_encrypt(&ctx->ctx.camellia, dst, src);
      break;
    default:
      ASSERT(0);
      break;
  }
}

static void
cipher_ctx_decrypt(cipher_t *ctx,
                   unsigned char *dst,
                   const unsigned char *src) {
  switch (ctx->type) {
    case CIPHER_AES128:
    case CIPHER_AES192:
    case CIPHER_AES256:
      aes_decrypt(&ctx->ctx.aes, dst, src);
      break;
    case CIPHER_BLOWFISH:
      blowfish_decrypt(&ctx->ctx.blowfish, dst, src);
      break;
    case CIPHER_CAMELLIA128:
    case CIPHER_CAMELLIA192:
    case CIPHER_CAMELLIA256:
      camellia_decrypt(&ctx->ctx.camellia, dst, src);
      break;
    default:
      ASSERT(0);
      break;
  }
}

static int
cipher_mode_init(cipher_t *ctx, const unsigned char *iv, size_t iv_len) {
  switch (ctx->mode) {
    case CIPHER_MODE_ECB: {
      if (iv_len != 0)
        return 0;

      break;
    }

    case CIPHER_MODE_CBC: {
      if (iv_len != ctx->block_size)
        return 0;

      memcpy(ctx->prev, iv, iv_len);

      break;
    }

    case CIPHER_MODE_CTR: {
      if (iv_len != ctx->block_size)
        return 0;

      memcpy(ctx->ctr, iv, iv_len);

      break;
    }

    case CIPHER_MODE_CFB: {
      if (iv_len != ctx->block_size)
        return 0;

      memcpy(ctx->prev, iv, iv_len);

      break;
    }

    case CIPHER_MODE_OFB: {
      if (iv_len != ctx->block_size)
        return 0;

      memcpy(ctx->state, iv, iv_len);

      break;
    }

    case CIPHER_MODE_GCM: {
      return 0;
    }

    default: {
      return 0;
    }
  }

  return 1;
}

static void
cipher_mode_update(cipher_t *ctx,
                   unsigned char *dst,
                   const unsigned char *src) {
  size_t i;

  switch (ctx->mode) {
    case CIPHER_MODE_ECB: {
      if (ctx->encrypt)
        cipher_ctx_encrypt(ctx, dst, src);
      else
        cipher_ctx_encrypt(ctx, dst, src);

      break;
    }

    case CIPHER_MODE_CBC: {
      if (ctx->encrypt) {
        for (i = 0; i < ctx->block_size; i++)
          dst[i] = src[i] ^ ctx->prev[i];

        cipher_ctx_encrypt(ctx, dst, dst);

        memcpy(ctx->prev, dst, ctx->block_size);
      } else {
        cipher_ctx_decrypt(ctx, dst, src);

        for (i = 0; i < ctx->block_size; i++)
          dst[i] ^= ctx->prev[i];

        memcpy(ctx->prev, src, ctx->block_size);
      }

      break;
    }

    case CIPHER_MODE_CTR: {
      cipher_ctx_encrypt(ctx, ctx->state, ctx->ctr);

      for (i = ctx->block_size; i-- > 0;) {
        if (++ctx->ctr[i] != 0x00)
          break;
      }

      for (i = 0; i < ctx->block_size; i++)
        dst[i] = src[i] ^ ctx->state[i];

      break;
    }

    case CIPHER_MODE_CFB: {
      cipher_ctx_encrypt(ctx, ctx->state, ctx->prev);

      for (i = 0; i < ctx->block_size; i++)
        dst[i] = src[i] ^ ctx->state[i];

      if (ctx->encrypt)
        memcpy(ctx->prev, dst, ctx->block_size);
      else
        memcpy(ctx->prev, src, ctx->block_size);

      break;
    }

    case CIPHER_MODE_OFB: {
      cipher_ctx_encrypt(ctx, ctx->state, ctx->state);

      for (i = 0; i < ctx->block_size; i++)
        dst[i] = src[i] ^ ctx->state[i];

      break;
    }

    case CIPHER_MODE_GCM: {
      ASSERT(0);
      break;
    }

    default: {
      ASSERT(0);
      break;
    }
  }
}

static int
cipher_mode_final(cipher_t *ctx, unsigned char *out, size_t *out_len) {
  size_t i;

  switch (ctx->mode) {
    case CIPHER_MODE_ECB:
    case CIPHER_MODE_CBC: {
      size_t left, end;

      if (ctx->encrypt) {
        left = ctx->block_size - ctx->block_pos;

        memcpy(out, ctx->block, ctx->block_pos);

        for (i = ctx->block_pos; i < ctx->block_size; i++)
          out[i] = left;

        cipher_mode_update(ctx, out, out);

        *out_len = ctx->block_size;
      } else {
        if (ctx->block_pos != 0)
          return 0;

        if (ctx->last_size == 0)
          return 0;

        ASSERT(ctx->last_size == ctx->block_size);

        memcpy(out, ctx->last, ctx->block_size);

        left = out[ctx->block_size - 1];

        if (left == 0 || left > ctx->block_size)
          return 0;

        end = ctx->block_size - left;

        for (i = end; i < ctx->block_size; i++) {
          if (out[i] != left)
            return 0;
        }

        *out_len = end;
      }

      break;
    }

    case CIPHER_MODE_CTR: {
      cipher_ctx_encrypt(ctx, ctx->state, ctx->ctr);

      for (i = 0; i < ctx->block_pos; i++)
        out[i] = ctx->block[i] ^ ctx->state[i];

      *out_len = ctx->block_pos;

      break;
    }

    case CIPHER_MODE_CFB: {
      cipher_ctx_encrypt(ctx, ctx->state, ctx->prev);

      for (i = 0; i < ctx->block_pos; i++)
        out[i] = ctx->block[i] ^ ctx->state[i];

      *out_len = ctx->block_pos;

      break;
    }

    case CIPHER_MODE_OFB: {
      cipher_ctx_encrypt(ctx, ctx->state, ctx->state);

      for (i = 0; i < ctx->block_pos; i++)
        out[i] = ctx->block[i] ^ ctx->state[i];

      *out_len = ctx->block_pos;

      break;
    }

    case CIPHER_MODE_GCM: {
      return 0;
    }

    default: {
      return 0;
    }
  }

  return 1;
}

int
cipher_init(cipher_t *ctx, int type, int mode, int encrypt,
            const unsigned char *key, size_t key_len,
            const unsigned char *iv, size_t iv_len) {
  memset(ctx, 0, sizeof(*ctx));

  ctx->type = type;
  ctx->mode = mode;
  ctx->encrypt = encrypt;

  if (!cipher_ctx_init(ctx, key, key_len))
    return 0;

  if (!cipher_mode_init(ctx, iv, iv_len))
    return 0;

  return 1;
}

void
cipher_update(cipher_t *ctx,
              unsigned char *output, size_t *output_len,
              const unsigned char *input, size_t input_len) {
  int has_padding = ctx->mode <= CIPHER_MODE_CBC;
  size_t bsize = ctx->block_size;
  size_t bpos = ctx->block_pos;
  size_t ilen = input_len;
  size_t olen = 0;
  size_t ipos = 0;
  size_t opos = 0;

  if (ilen == 0) {
    *output_len = 0;
    return;
  }

  ctx->block_pos = (ctx->block_pos + ilen) % bsize;

  if (has_padding)
    olen += ctx->last_size;

  if (bpos > 0) {
    size_t want = bsize - bpos;

    if (want > ilen)
      want = ilen;

    memcpy(ctx->block + bpos, input + ipos, want);

    bpos += want;
    ilen -= want;
    ipos += want;

    if (bpos < bsize) {
      *output_len = 0;
      return;
    }

    olen += bsize;
  }

  olen += ilen - (ilen % bsize);

  *output_len = olen;

  if (has_padding) {
    memcpy(output + opos, ctx->last, ctx->last_size);
    opos += ctx->last_size;
  }

  if (bpos > 0) {
    cipher_mode_update(ctx, output + opos, ctx->block);
    opos += bsize;
  }

  while (ilen >= bsize) {
    cipher_mode_update(ctx, output + opos, input + ipos);
    opos += bsize;
    ipos += bsize;
    ilen -= bsize;
  }

  if (ilen > 0)
    memcpy(ctx->block, input + ipos, ilen);

  if (has_padding && olen > 0) {
    memcpy(ctx->last, output + olen - bsize, bsize);

    ctx->last_size = bsize;

    *output_len = olen - bsize;
  }
}

int
cipher_final(cipher_t *ctx, unsigned char *output, size_t *output_len) {
  return cipher_mode_final(ctx, output, output_len);
}
