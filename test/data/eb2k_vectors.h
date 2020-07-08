typedef struct eb2k_vector_s {
  const char *pass;
  const char *salt;
  size_t key_len;
  size_t iv_len;
  const char *key;
  const char *iv;
} eb2k_vector_t;

static const eb2k_vector_t eb2k_vectors[] = {
  {
    "1234567890",
    "286e37f0b031c0d2c72c674b116e5342",
    16,
    16,
    "1fc6b89794d36b3ce33087de2c107fb6",
    "8df8d21d44bfe94ab64cde7d6496857f"
  },
  {
    "1234567890",
    "113bdee3c90ee160997b2e373660368e",
    16,
    16,
    "85560be2bd49c607ca036024999c3a30",
    "daddbffe81ee02330c98d492c8f430bd"
  },
  {
    "foo",
    "113bdee3c90ee160997b2e373660368e",
    16,
    16,
    "b25272c2fa948e94d01242e297762876",
    "562dedfd27b08b383254ffe1efc7b60d"
  },
  {
    "1234567890",
    "bdb867fd73611c0586852f406b934243",
    16,
    16,
    "85dad080e8a280f6f9ef18808afea519",
    "a3ada9edf9fb2c4433a6efbbc210505a"
  },
  {
    "1234567890",
    "0eca5c2e3a32fa90e57d685621eb9786",
    16,
    16,
    "0a9b842e6bbf13a26a5e99e5ec625857",
    "c48bd7cfd5803964a8132326692ef6d8"
  },
  {
    "foo",
    "0eca5c2e3a32fa90e57d685621eb9786",
    16,
    16,
    "df157b6c6d652f934a5fee158c5cc0e7",
    "9e434c07188328c657638b42b9b2040d"
  },
  {
    "1234567890",
    "f27df2f6f6856e9bb4e2fe4c5269162a",
    16,
    16,
    "654a6043e4e021375d0954f3b1a66fb4",
    "361d6661e43465d2ae203380f0a30c74"
  },
  {
    "1234567890",
    "acf9985e9f755604b48dc03f48c8f584",
    16,
    16,
    "7b8b628d62ca33fa5133bed68e5d18b9",
    "218e80364ecac6bf3fee95781baa98e6"
  },
  {
    "foo",
    "acf9985e9f755604b48dc03f48c8f584",
    16,
    16,
    "040b57f89272efc3e225672799a9d2b7",
    "7b611924696718c6fbfb600965ee587e"
  }
};
