/*!
 * hrtime.c - high-resolution time for libtorsion benchmarks
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#define TORSION_ENTROPY_H

#if defined(__unix) || defined(__unix__)
#  define TORSION_UNIX
#elif defined(__APPLE__) && defined(__MACH__)
/* unix not defined on GCC/Clang. */
#  define TORSION_UNIX
#elif defined(__linux__) || defined(_AIX)
/* unix not defined on IBM XL C. */
#  define TORSION_UNIX
#elif defined(__HAIKU__)
/* unix not defined on GCC (Haiku Patch).
   Note that it _is_ defined on Clang. */
#  define TORSION_UNIX
#elif defined(__QNX__)
/* unix not defined on GCC. */
#  define TORSION_UNIX
#elif defined(__MVS__)
/* unix not defined on Clang or XL C. */
#  define TORSION_UNIX
#endif

#include "../src/entropy/hrt.c"
