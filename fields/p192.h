#ifdef BCRYPTO_EC_64BIT
typedef uint64_t p192_fe_word_t;
#define P192_FIELD_WORDS 4
#include "p192_64.h"
#else
typedef uint32_t p192_fe_word_t;
#define P192_FIELD_WORDS 7
#include "p192_32.h"
#endif
