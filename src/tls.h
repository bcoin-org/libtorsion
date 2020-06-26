/*!
 * tls.h - thread-local storage for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef _TORSION_TLS_H
#define _TORSION_TLS_H

/* TLS Compiler Support
 *
 * GCC:
 *
 * - TLS first implemented in GCC 3.3.1[1].
 * - Full TLS emulation by GCC 4.3[2].
 * - Note that qtbase abides by 4.8 for some reason[3].
 *
 * Clang:
 *
 * - TLS first implemented in Clang 2.6 (?)[4].
 * - Support for Mac OSX added in Clang 3.0[5].
 * - Support for x86{,-64} Windows added in Clang 3.5[6].
 * - Support for x86 Cygwin added in Clang 3.8[7].
 * - Support for Apple iOS added in Clang 3.8[8].
 * - Support for ARM Windows added in Clang 3.9[9].
 * - Support for OpenBSD added in Clang 5.0[10].
 * - Support for iOS Simulator added in Clang 7.0[11].
 * - Support still lacking for ARM/x86-64 Cygwin as of Clang 10.0[12].
 * - Support still lacking for Haiku as of Clang 10.0[13].
 *
 * Intel C:
 *
 * - Unknown when TLS was implemented, but this post[14]
 *   suggests it existed at least as far back as 2006 (9.0/10.0).
 * - Curiously, qtbase requires version 15.0[15] (perhaps
 *   because 15.0 is when C++ and Apple gained support)[16].
 *
 * MSVC:
 *
 * - Unknown when TLS was implemented, but this repo[17]
 *   from 2009 suggests it existed in VS .NET 2002 (7.0).
 *
 * Sun Pro/Solaris Studio:
 *
 * - TLS first implemented in 5.9[18].
 * - qtbase does not use Sun TLS at all[19].
 *
 * IBM XL C:
 *
 * - TLS first implemented for Linux in XL C 8.0[20].
 * - Support added for AIX in XL C 10.1[20].
 * - Note that -qtls must be passed on the command line.
 *
 * C++ Builder/Borland:
 *
 * - This seems to be rumored by wikipedia[21], so everyone does it.
 *
 * Digital Mars C/C++:
 *
 * - Ditto[21].
 *
 * Green Hill C/C++:
 *
 * - TLS first implemented in 20.15.5, used by qtbase[22].
 *
 * ARM C:
 *
 * - Mentioned on the GNU mailing list[23][24].
 *
 * [1] https://gcc.gnu.org/onlinedocs/gcc-3.3.1/gcc/Thread-Local.html
 * [2] grep the 4.3 tree for "thread-local storage not supported"
 * [3] https://github.com/qt/qtbase/blob/44da43e/src/corelib/global/qcompilerdetection.h#L898
 * [4] https://github.com/llvm/llvm-project/blob/release/2.6.x/clang/lib/Basic/Targets.cpp
 * [5] https://github.com/llvm/llvm-project/blob/release/3.0.x/clang/lib/Basic/Targets.cpp#L202
 * [6] https://github.com/llvm/llvm-project/blob/release/3.5.x/clang/lib/Basic/Targets.cpp#L3110
 * [7] https://github.com/llvm/llvm-project/blob/release/3.8.x/clang/lib/Basic/Targets.cpp#L3855
 * [8] https://github.com/llvm/llvm-project/blob/release/3.8.x/clang/lib/Basic/Targets.cpp#L197
 * [9] https://github.com/llvm/llvm-project/blob/release/3.9.x/clang/lib/Basic/Targets.cpp#L5508
 * [10] https://github.com/llvm/llvm-project/blob/release/5.x/clang/lib/Basic/Targets.cpp#L555
 * [11] https://github.com/llvm/llvm-project/blob/release/7.x/clang/lib/Basic/Targets/OSTargets.h#L103
 * [12] https://github.com/llvm/llvm-project/blob/release/10.x/clang/lib/Basic/Targets/X86.h#L819
 * [13] https://github.com/llvm/llvm-project/blob/release/10.x/clang/lib/Basic/Targets/OSTargets.h#L310
 * [14] https://community.intel.com/t5/Intel-C-Compiler/thread-local-storage-linking-problems/td-p/932631
 * [15] https://github.com/qt/qtbase/blob/44da43e/src/corelib/global/qcompilerdetection.h#L646
 * [16] https://community.intel.com/t5/Intel-C-Compiler/Mach-O-thread-local-storage/td-p/948267
 * [17] https://github.com/snaewe/loki-lib/commit/7d8e59abc8f48785d564ddabab5ba3f01cd24444
 * [18] lost my source on this one
 * [19] https://github.com/qt/qtbase/blob/44da43e/src/corelib/global/qcompilerdetection.h#L434
 * [20] https://www.ibm.com/support/pages/node/318521#6
 * [21] https://en.wikipedia.org/wiki/Thread_local_storage#C_and_C++
 * [22] https://github.com/qt/qtbase/blob/44da43e/src/corelib/global/qcompilerdetection.h#L343
 * [23] https://lists.gnu.org/archive/html/bug-gnulib/2019-06/msg00063.html
 * [24] https://developer.arm.com/docs/dui0472/latest/compiler-specific-features/__declspecthread
 */

/* Now to stop and talk about Apple for a second...
 *
 * The TLS situation on Apple is mind-bogglingly stupid.
 * This results from the fact that the Mach-O executable
 * format is inferior to both ELF and PE, and also sheer
 * incompetence on Apple's part in altering Clang features
 * and intentially disabling TLS when it was supported.
 *
 * Note that as we go through the explanation below, we
 * have to distinguish Apple Clang from Real Clang.
 *
 * The Facts:
 *
 *   - Apple enabled TLS in Xcode 8.0[1] (Real Clang 3.9.0[2]).
 *   - Intel enabled TLS in ICC 15.0.0[3].
 *   - Real Clang 3.0 only supported OSX TLS[4].
 *   - Real Clang 3.8 added support for iOS TLS[5].
 *   - Real Clang 7.0 added support for iOS Simulator TLS[6].
 *   - Apple Clang 5.0 (Xcode 5.0.0) included Real Clang 3.3[2] (no TLS[1]).
 *   - Apple Clang 8.0.0 (Xcode 8.0) included Real Clang 3.9.0[1][2].
 *   - Apple Clang 10.0.1 (Xcode 10.2) included Real Clang 7.0.0[2].
 *   - The minimum OSX version required must be >=10.7 for TLS[7].
 *   - The minimum iOS 64-bit version required must be >=8 for TLS[7].
 *   - The minimum iOS 32-bit version required must be >=9 for TLS[7].
 *   - The minimum iOS simulator version required must be >=10 for TLS[7].
 *
 * [1] https://stackoverflow.com/a/29929949
 * [2] https://en.wikipedia.org/wiki/Xcode#Xcode_7.0_-_11.x_(since_Free_On-Device_Development)
 * [3] https://community.intel.com/t5/Intel-C-Compiler/Mach-O-thread-local-storage/td-p/948267
 * [4] https://github.com/llvm/llvm-project/blob/release/3.0.x/clang/lib/Basic/Targets.cpp#L202
 * [5] https://github.com/llvm/llvm-project/blob/release/3.8.x/clang/lib/Basic/Targets.cpp#L197
 * [6] https://github.com/llvm/llvm-project/blob/release/7.x/clang/lib/Basic/Targets/OSTargets.h#L95
*/

/* And now on to Android quirks.
 *
 * Android apparently has issues compiling with NDK < r12[1].
 * This appears to be fixed in NDK r12 according to SO[2].
 * Note that NDK r12 upgrades Clang to 8.0[3].
 *
 * [1] https://stackoverflow.com/questions/27191214
 * [2] https://stackoverflow.com/a/27195324
 * [3] https://developer.android.com/ndk/downloads/revision_history
 */

#if defined(__APPLE__) && defined(__MACH__)
#  if defined(__ENVIRONMENT_IPHONE_OS_VERSION_MIN_REQUIRED__)
/* Note: Clang serializes the short form version prior to 3.9.1. */
#    if __ENVIRONMENT_IPHONE_OS_VERSION_MIN_REQUIRED__ >= 100000 /* 10.0 */
#      define TORSION__HAVE_IPHONE_OS
#      define TORSION__HAVE_APPLE_OS
#    endif
#  elif defined(__ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__)
#    if __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ >= 1070 /* 10.7 */
#      define TORSION__HAVE_APPLE_OS
#    endif
#  endif
#endif

#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  define TORSION_HAVE_TLS
#elif defined(__BORLANDC__)
#  define TORSION_HAVE_TLS
#elif defined(__DMC__)
#  define TORSION_HAVE_TLS
#elif defined(__ghs)
#  if __GHS_VERSION_NUMBER >= 201505 /* 20.15.5 */
#    define TORSION_HAVE_TLS
#  endif
#elif defined(__CC_ARM)
#  define TORSION_HAVE_TLS
#elif defined(__SUNPRO_C)
#  if __SUNPRO_C >= 0x590 /* 5.9 */
#    define TORSION_HAVE_TLS
#  endif
#elif defined(__INTEL_COMPILER)
#  if defined(__APPLE__) && defined(__MACH__)
#    if defined(TORSION__HAVE_APPLE_OS) && __INTEL_COMPILER >= 1500 /* 15.0.0 */
#      define TORSION_HAVE_TLS
#    endif
#  else
#    if __INTEL_COMPILER >= 1000 /* 10.0.0 */
#      define TORSION_HAVE_TLS
#    endif
#  endif
#elif defined(__clang__) && defined(__clang_major__) && defined(__clang_minor__)
#  if defined(_WIN32) || defined(__MINGW32__)
#    ifdef _M_ARM
#      if ((__clang_major__ << 16) + __clang_minor__ >= 0x30009) /* 3.9 */
#        define TORSION_HAVE_TLS
#      endif
#    else
#      if ((__clang_major__ << 16) + __clang_minor__ >= 0x30005) /* 3.5 */
#        define TORSION_HAVE_TLS
#      endif
#    endif
#  elif defined(__CYGWIN__)
#    if defined(_X86_)
#      if ((__clang_major__ << 16) + __clang_minor__ >= 0x30008) /* 3.8 */
#        define TORSION_HAVE_TLS
#      endif
#    elif defined(__x86_64__)
/* Still no support as of Clang 10.0. */
#    elif defined(_ARM_)
/* Still no support as of Clang 10.0. */
#    endif
#  elif defined(__APPLE__) && defined(__MACH__)
#    ifdef TORSION__HAVE_APPLE_OS
#      if defined(__apple_build_version__) /* >= 8000038 */
#        if ((__clang_major__ << 16) + __clang_minor__ >= 0x80000) /* 8.0 */
#          define TORSION_HAVE_TLS
#        endif
#      elif defined(TORSION__HAVE_IPHONE_OS)
#        if ((__clang_major__ << 16) + __clang_minor__ >= 0x30008) /* 3.8 */
#          define TORSION_HAVE_TLS
#        endif
#      else
#        if ((__clang_major__ << 16) + __clang_minor__ >= 0x30000) /* 3.0 */
#          define TORSION_HAVE_TLS
#        endif
#      endif
#    endif
#  elif defined(__OpenBSD__)
#    if ((__clang_major__ << 16) + __clang_minor__ >= 0x50000) /* 5.0 */
#      define TORSION_HAVE_TLS
#    endif
#  elif defined(__HAIKU__)
/* Still unsupported as of Clang 10.0. */
#  elif defined(__ANDROID__)
#    if ((__clang_major__ << 16) + __clang_minor__ >= 0x30008) /* 3.8 */
#      define TORSION_HAVE_TLS
#    endif
#  else
#    if ((__clang_major__ << 16) + __clang_minor__ >= 0x20006) /* 2.6 */
#      define TORSION_HAVE_TLS
#    endif
#  endif
#elif defined(__ibmxl__)
#  error "XL C did not define __clang__."
#elif defined(__xlC__)
#  if defined(__linux__)
#    if __xlC__ >= 0x0800 /* 8.0 */
#      define TORSION_HAVE_TLS
#    endif
#  elif defined(_AIX)
#    if __xlC__ >= 0x1001 /* 10.1 */
#      define TORSION_HAVE_TLS
#    endif
#  endif
#elif defined(_MSC_VER)
#  if _MSC_VER >= 1300 /* VS .NET 2002 (7.0) */
#    define TORSION_HAVE_TLS
#  endif
#elif defined(__GNUC__) && defined(__GNUC_MINOR__)
#  ifdef __linux__
#    if ((__GNUC__ << 16) + __GNUC_MINOR__ >= 0x30003) /* 3.3 */
#      define TORSION_HAVE_TLS
#    endif
#  else
#    if ((__GNUC__ << 16) + __GNUC_MINOR__ >= 0x40003) /* 4.3 */
#      define TORSION_HAVE_TLS
#    endif
#  endif
#endif

#ifdef TORSION_HAVE_TLS
/* See: https://software.intel.com/en-us/forums/intel-c-compiler/topic/721059 */
#  if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#    if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1800 /* 18.0.0 */
#      include <limits.h>
#      ifndef __STDC_NO_THREADS__
#        define TORSION_STDC_THREADS
#      endif
#    endif
#  endif
#  if defined(__EMSCRIPTEN__) || defined(__wasm__)
#    define TORSION_TLS
#  elif defined(TORSION_STDC_THREADS)
#    define TORSION_TLS _Thread_local
#  elif defined(_WIN32) && !defined(__MINGW32__)
#    define TORSION_TLS __declspec(thread)
#  else
#    define TORSION_TLS __thread
#  endif
#else
#  define TORSION_TLS
#endif

#endif /* _TORSION_TLS_H */
