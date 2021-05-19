/*!
 * posix.h - posix compat for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef TORSION_POSIX_H
#define TORSION_POSIX_H

/* We need to do this because:
 *
 * 1. A lot of feature testing code in C standard libraries
 *    disables almost every POSIX feature if __STRICT_ANSI__
 *    (or __STDC__ == 1) is defined.
 *
 * 2. Some OSes won't even work at all without _POSIX_SOURCE
 *    defined.
 *
 * 3. Some OSes just have really bad defaults.
 *
 * For reference:
 *
 *   - glibc: https://github.com/bminor/glibc/blob/6c57d32/include/features.h
 *   - bionic: https://github.com/aosp-mirror/platform_bionic/blob/a1112fd/libc/include/sys/cdefs.h#L162
 *   - musl: https://github.com/ifduyue/musl/blob/1febd21/include/features.h
 *   - cygwin: https://github.com/cygwin/cygwin/blob/8050ef2/newlib/libc/include/sys/features.h
 *   - darwin: https://github.com/apple/darwin-xnu/blob/d4061fb/bsd/sys/cdefs.h#L792
 *   - freebsd: https://github.com/freebsd/freebsd-src/blob/cfad8bd/sys/sys/cdefs.h#L638
 *   - openbsd: https://github.com/openbsd/src/blob/43b1a0f/sys/sys/cdefs.h#L261
 *   - netbsd: https://github.com/NetBSD/src/blob/ee650a6/sys/sys/featuretest.h
 *   - dragonfly: https://github.com/DragonFlyBSD/DragonFlyBSD/blob/c97dc9d/sys/sys/cdefs.h#L597
 *   - sun: https://docs.oracle.com/cd/E19253-01/816-5175/standards-5/index.html
 *          https://github.com/illumos/illumos-gate/blob/9ecd05b/usr/src/uts/common/sys/feature_tests.h
 *          https://github.com/illumos/illumos-gate/blob/9ecd05b/usr/src/uts/common/sys/unistd.h
 *   - hpux: https://nixdoc.net/man-pages/HP-UX/man5/stdsyms.5.html
 *   - aix: ??
 *   - z/os: https://www.ibm.com/docs/en/zos/2.3.0?topic=files-feature-test-macros
 *   - qnx: http://www.qnx.com/developers/docs/6.5.0_sp1/topic/com.qnx.doc.neutrino_prog/devel.html
            https://github.com/vocho/openqnx/blob/cc95df3/trunk/lib/c/public/sys/platform.h
 *   - haiku: https://github.com/haiku/haiku/blob/144f45a/headers/compatibility/bsd/features.h
 *            https://github.com/haiku/haiku/blob/144f45a/headers/posix/unistd.h
 *   - minix: https://github.com/Stichting-MINIX-Research-Foundation/minix/blob/4db99f4/sys/sys/featuretest.h
 *   - redox: https://github.com/redox-os/relibc/blob/9790289/include/bits/unistd.h
 *   - vxworks: ??
 *   - cloudabi: https://github.com/NuxiNL/cloudlibc/blob/7e5c649/src/include/unistd.h
 *   - wasi: https://github.com/WebAssembly/wasi-libc/blob/575e157/libc-top-half/musl/include/features.h
 *   - emscripten: https://github.com/emscripten-core/emscripten/blob/dee59ba/system/lib/libc/musl/include/features.h
 *   - posix: https://pubs.opengroup.org/onlinepubs/7908799/xsh/compilation.html
 *            https://pubs.opengroup.org/onlinepubs/007904975/functions/xsh_chap02_02.html
 */

#if defined(__linux__) || defined(__CYGWIN__)
#  undef _GNU_SOURCE
#  define _GNU_SOURCE
#elif defined(__gnu_hurd__) || (defined(__GNU__) && defined(__MACH__))
#  undef _GNU_SOURCE
#  define _GNU_SOURCE
#elif defined(__FreeBSD_kernel__) && !defined(__FreeBSD__)
#  undef _GNU_SOURCE
#  define _GNU_SOURCE
#elif defined(__sun) && defined(__SVR4)
#  undef _XOPEN_SOURCE
#  undef __EXTENSIONS__
#  undef _REENTRANT
#  define _XOPEN_SOURCE 500
#  define __EXTENSIONS__
#  define _REENTRANT
#elif defined(__hpux)
#  undef _HPUX_SOURCE
#  undef _REENTRANT
#  define _HPUX_SOURCE
#  define _REENTRANT
#elif defined(_AIX)
#  undef _ALL_SOURCE
#  undef _THREAD_SAFE
#  define _ALL_SOURCE
#  define _THREAD_SAFE
#elif defined(__MVS__)
#  undef _ALL_SOURCE
#  define _ALL_SOURCE
#elif defined(__QNX__)
#  undef _QNX_SOURCE
#  define _QNX_SOURCE
#elif defined(__HAIKU__)
#  undef _BSD_SOURCE
#  define _BSD_SOURCE
#elif defined(__wasi__) || defined(__EMSCRIPTEN__)
#  undef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

/*
 * POSIX Detection
 */

#if defined(__unix) || defined(__unix__)
#  define TORSION_POSIX
#elif defined(__APPLE__) && defined(__MACH__)
/* unix not defined on GCC/Clang. */
#  define TORSION_POSIX
#elif defined(__linux__) || defined(_AIX)
/* unix not defined on IBM XL C. */
#  define TORSION_POSIX
#elif defined(__MVS__)
/* unix not defined on Clang or XL C. */
#  define TORSION_POSIX
#elif defined(__QNX__)
/* unix not defined on GCC. */
#  define TORSION_POSIX
#elif defined(__HAIKU__)
/* unix not defined on GCC (Haiku Patch).
   Note that it _is_ defined on Clang. */
#  define TORSION_POSIX
#elif defined(__vxworks)
/* unix not defined on GCC (?). */
#  define TORSION_POSIX
#elif defined(__CloudABI__)
/* unix not defined on Clang. */
#  define TORSION_POSIX
#elif defined(__wasi__)
/* unix not defined on Clang. */
#  define TORSION_POSIX
#elif defined(__EMSCRIPTEN__)
/* unix not defined on (raw) Clang, but
   it is defined by emcc since 2016. */
#  define TORSION_POSIX
#endif

#endif /* TORSION_POSIX_H */
