typedef struct schnorr_vector_s {
  const char *priv;
  const char *pub;
  const char *aux;
  const char *msg;
  const char *sig;
  int result;
  const char *comment;
} schnorr_vector_t;

static const schnorr_vector_t schnorr_vectors[] = {
  {
    "0000000000000000000000000000000000000000000000000000000000000003",
    "f9308a019258c31049344f85f89d5229b531c845836f99b08601f113bce036f9",
    "0000000000000000000000000000000000000000000000000000000000000000",
    "0000000000000000000000000000000000000000000000000000000000000000",
    "067e337ad551b2276ec705e43f0920926a9ce08ac68159f9d258c9bba412781c9f059fcdf4824f13b3d7c1305316f956704bb3fea2c26142e18acd90a90c947e",
    1,
    ""
  },
  {
    "b7e151628aed2a6abf7158809cf4f3c762e7160f38b4da56a784d9045190cfef",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "0000000000000000000000000000000000000000000000000000000000000001",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "0e12b8c520948a776753a96f21abd7fdc2d7d0c0ddc90851be17b04e75ef86a47ef0da46c4dc4d0d1bcb8668c2ce16c54c7c23a6716ede303af86774917cf928",
    1,
    ""
  },
  {
    "c90fdaa22168c234c4c6628b80dc1cd129024e088a67cc74020bbea63b14e5c9",
    "dd308afec5777e13121fa72b9cc1b7cc0139715309b086c960e18fd969774eb8",
    "c87aa53824b4d7ae2eb035a2b5bbbccc080e76cdc6d1692c4b0b62d798e6d906",
    "7e2d58d8b3bcdf1abadec7829054f90dda9805aab56c77333024b9d0a508b75c",
    "fc012f9fb8fe00a358f51ef93dce0dc0c895f6e9a87c6c4905bc820b0c3677616b8737d14e703af8e16e22e5b8f26227d41e5128f82d86f747244cc289c74d1d",
    1,
    ""
  },
  {
    "0b432b2677937381aef05bb02a66ecd012773062cf3fa2549e44f58ed2401710",
    "25d1dff95105f5253c4022f628a996ad3a0d95fbf21d468a1b33f8c160d8f517",
    "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
    "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
    "fc132d4e426dff535aec0fa7083ac5118bc1d5fffd848abd8290c23f271ca0dd11aedcea3f55da9bd677fe29c9dda0cf878bce43fde0e313d69d1af7a5ae8369",
    1,
    "test fails if msg is reduced modulo p or n"
  },
  {
    "",
    "d69c3509bb99e412e68b0fe8544e72837dfa30746d8be2aa65975f29d22dc7b9",
    "",
    "4df3c3f68fcc83b27e9d42c90431a72499f17875c81a599b566c9889b9696703",
    "00000000000000000000003b78ce563f89a0ed9414f5aa28ad0d96d6795f9c630ec50e5363e227acac6f542ce1c0b186657e0e0d1a6ffe283a33438de4738419",
    1,
    ""
  },
  {
    "",
    "eefdea4cdb677750a420fee807eacf21eb9898ae79b9768766e4faa04a2d4a34",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "7036d6bfe1837ae919631039a2cf652a295dfac9a8bbb0806014b2f48dd7c807941607b563abba414287f374a332ba3636de009ee1ef551a17796b72b68b8a24",
    0,
    "public key not on the curve"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "f9308a019258c31049344f85f89d5229b531c845836f99b08601f113bce036f995a579da959fa739fce39e8bd16fecb5cdcf97060b2c73cde60e87abca1aa5d9",
    0,
    "has_square_y(R) is 0"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "f8704654f4687b7365ed32e796de92761390a3bcc495179bfe073817b7ed32824e76b987f7c1f9a751ef5c343f7645d3cffc7d570b9a7192ebf1898e1344e3bf",
    0,
    "negated message"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "7036d6bfe1837ae919631039a2cf652a295dfac9a8bbb0806014b2f48dd7c8076be9f84a9c5445bebd780c8b5ccd45c883d0dc47cd594b21a858f31a19aab71d",
    0,
    "negated s value"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "00000000000000000000000000000000000000000000000000000000000000009915ee59f07f9dbbaedc31bfcc9b34ad49de669cd24773bced77dda36d073ec8",
    0,
    "sG - eP is infinite. Test fails in single verification if has_square_y(inf) is defined as 1 and x(inf) as 0"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "0000000000000000000000000000000000000000000000000000000000000001c7ec918b2b9cf34071bb54bed7eb4bb6bab148e9a7e36e6b228f95dfa08b43ec",
    0,
    "sG - eP is infinite. Test fails in single verification if has_square_y(inf) is defined as 1 and x(inf) as 1"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "4a298dacae57395a15d0795ddbfd1dcb564da82b0f269bc70a74f8220429ba1d941607b563abba414287f374a332ba3636de009ee1ef551a17796b72b68b8a24",
    0,
    "sig{0:32} is not an X coordinate on the curve"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f941607b563abba414287f374a332ba3636de009ee1ef551a17796b72b68b8a24",
    0,
    "sig{0:32} is equal to field size"
  },
  {
    "",
    "dff1d77f2a671c5f36183726db2341be58feae1da2deced843240f7b502ba659",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "7036d6bfe1837ae919631039a2cf652a295dfac9a8bbb0806014b2f48dd7c807fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",
    0,
    "sig{32:64} is equal to curve order"
  },
  {
    "",
    "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc30",
    "",
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89",
    "7036d6bfe1837ae919631039a2cf652a295dfac9a8bbb0806014b2f48dd7c807941607b563abba414287f374a332ba3636de009ee1ef551a17796b72b68b8a24",
    0,
    "public key is not a valid X coordinate because it exceeds the field size"
  }
};
