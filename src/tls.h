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
 * - Supports TLS via __thread[1].
 * - TLS first implemented in GCC 3.3[1].
 * - TLS emulation added in GCC 4.3[2].
 *
 * Clang:
 *
 * - Supports TLS via __thread, __declspec(thread)[3].
 * - TLS first implemented in Clang 2.0[4].
 * - TLS emulation added in Clang 3.8[5].
 * - TLS "modernized" in Clang 2.6[6].
 * - Support for Mac OSX added in Clang 2.8[7].
 * - Support for Mac OSX "modernized" in Clang 3.0[8].
 * - Support for System Z added in Clang 3.3[9].
 * - Support for __has_extension(c_thread_local) added in 3.4[10].
 * - Support for x86 Windows added in Clang 3.5[11].
 * - Support for x86-64 Windows added in Clang 3.5[12].
 * - Support for x86 Cygwin added in Clang 3.8[13].
 * - Support for Apple iOS added in Clang 3.8[14].
 * - Support for ARM Windows added in Clang 3.9[15].
 * - Support for OpenBSD added in Clang 5.0[16].
 * - Support for iOS Simulator added in Clang 7.0[17].
 * - Support for RISC-V added in Clang 9.0[18].
 * - No support for x86-64 Cygwin as of Clang 10.0[19].
 * - No support for ARM Cygwin as of Clang 10.0[20].
 * - No support for Haiku as of Clang 10.0[21].
 *
 * Intel C:
 *
 * - Supports TLS via __thread[22], __declspec(thread)[23].
 * - Intel mentions __thread in documentation from 2004[22].
 *   This would suggest support from version 8.x onward.
 * - Furthermore, this post[24] suggests __thread existed
 *   at least as far back as 2006 (version 9.0 or 10.0).
 * - Apple and Mach-O support implemented in 15.0[25].
 *
 * MSVC:
 *
 * - Supports TLS via __declspec(thread)[26].
 * - Unknown when TLS was implemented, but this repo[27]
 *   from 2009 suggests it existed in VS .NET 2002 (7.0).
 *
 * Sun Pro C / Sun Studio / Solaris Studio:
 *
 * - Supports TLS via __thread[28].
 * - First mentioned in documentation for Sun Studio 12[29].
 *   This would suggest support from 5.9 onward.
 *
 * IBM XL C:
 *
 * - Supports TLS via __thread[30].
 * - Support added for Linux in XL C 8.0[31].
 * - Support added for AIX in XL C 10.1[31].
 * - Note that -qtls must be passed on the command line.
 *
 * C++ Builder:
 *
 * - Supports TLS via __thread[32], __declspec(thread)[33].
 * - Mentioned in C++ Builder 2009 documentation[32].
 *
 * Digital Mars C/C++:
 *
 * - Supports TLS via __declspec(thread) (32 bit only)[34].
 * - TLS supported since at least 2005 (8.42n)[35].
 *
 * ARM CC:
 *
 * - Supports TLS via __thread, __declspec(thread)[36][37].
 * - Mentioned on the gnulib mailing list[38].
 *
 * HP ANSI C:
 *
 * - Supports TLS via __thread[39].
 * - The release notes suggest this has been the
 *   case since at least version A.05.55.02.
 *
 * Watcom C:
 *
 * - TLS supported via __declspec(thread)[40].
 * - TLS supported since at least version 11.0c[40].
 *   Notable as this predates Open Watcom.
 *
 * Wind River Compiler (Diab C):
 *
 * - TLS supported via __thread[41][42].
 * - TLS supported since at least 2007 (5.6)[42].
 *
 * NWCC:
 *
 * - TLS supported via __thread[43].
 *
 * CompCert:
 *
 * - TLS not yet supported[44].
 *
 * C11:
 *
 * - C11 specifies support for _Thread_local[45].
 * - Support can be tested by checking both:
 *
 *     __STDC_VERSION__ >= 201112L
 *     !defined(__STDC_NO_THREADS__)
 *
 *   However, some compilers do not define STDC_NO_THREADS
 *   or do not define it directly (in particular, Intel C
 *   versions less than 18.0.0[46]).
 *
 * [1] https://gcc.gnu.org/onlinedocs/gcc-3.3.1/gcc/Thread-Local.html
 * [2] https://github.com/gcc-mirror/gcc/commit/8893239dc4ed32bd3bb4e00d6e43b859554ab82a
 * [3] https://clang.llvm.org/docs/AttributeReference.html#thread
 * [4] https://releases.llvm.org/2.0/docs/ReleaseNotes.html
 * [5] https://releases.llvm.org/3.8.0/docs/ReleaseNotes.html
 * [6] https://github.com/llvm/llvm-project/blob/llvmorg-2.6.0/clang/lib/Basic/Targets.cpp
 * [7] https://github.com/llvm/llvm-project/blob/llvmorg-2.8.0/clang/lib/Basic/Targets.cpp#L153
 * [8] https://github.com/llvm/llvm-project/blob/llvmorg-3.0.0/clang/lib/Basic/Targets.cpp#L202
 * [9] https://github.com/llvm/llvm-project/blob/llvmorg-3.3.0/clang/lib/Basic/Targets.cpp#L4352
 * [10] https://github.com/llvm/llvm-project/blob/llvmorg-3.4.0/clang/lib/Lex/PPMacroExpansion.cpp#L998
 * [11] https://github.com/llvm/llvm-project/blob/llvmorg-3.5.0/clang/lib/Basic/Targets.cpp#L3110
 * [12] https://github.com/llvm/llvm-project/blob/llvmorg-3.8.0/clang/lib/Basic/Targets.cpp#L4133
 * [13] https://github.com/llvm/llvm-project/blob/llvmorg-3.8.0/clang/lib/Basic/Targets.cpp#L3855
 * [14] https://github.com/llvm/llvm-project/blob/llvmorg-3.8.0/clang/lib/Basic/Targets.cpp#L232
 * [15] https://github.com/llvm/llvm-project/blob/llvmorg-3.9.0/clang/lib/Basic/Targets.cpp#L5495
 * [16] https://github.com/llvm/llvm-project/blob/llvmorg-5.0.0/clang/lib/Basic/Targets.cpp#L555
 * [17] https://github.com/llvm/llvm-project/blob/llvmorg-7.0.0/clang/lib/Basic/Targets/OSTargets.h#L103
 * [18] https://github.com/llvm/llvm-project/blob/llvmorg-9.0.0/clang/lib/Basic/Targets/RISCV.h#L24
 * [19] https://github.com/llvm/llvm-project/blob/llvmorg-10.0.0/clang/lib/Basic/Targets/X86.h#L819
 * [20] https://github.com/llvm/llvm-project/blob/llvmorg-10.0.0/clang/lib/Basic/Targets/ARM.cpp#L1208
 * [21] https://github.com/llvm/llvm-project/blob/llvmorg-10.0.0/clang/lib/Basic/Targets/OSTargets.h#L310
 * [22] https://software.intel.com/sites/default/files/ae/4f/6320
 * [23] https://community.intel.com/t5/Intel-C-Compiler/Thread-local-storage-support-on-Windows/td-p/949321
 * [24] https://community.intel.com/t5/Intel-C-Compiler/thread-local-storage-linking-problems/td-p/932631
 * [25] https://community.intel.com/t5/Intel-C-Compiler/Mach-O-thread-local-storage/td-p/948267
 * [26] https://docs.microsoft.com/en-us/cpp/c-language/thread-local-storage?view=vs-2019
 * [27] https://github.com/snaewe/loki-lib/commit/7d8e59abc8f48785d564ddabab5ba3f01cd24444
 * [28] https://docs.oracle.com/cd/E18659_01/html/821-1383/bkaeg.html
 * [29] https://docs.oracle.com/cd/E19205-01/819-5267/bkaeg/index.html
 * [30] https://www.ibm.com/support/knowledgecenter/en/SSXVZZ_13.1.3/com.ibm.xlcpp1313.lelinux.doc/language_ref/thread.html
 * [31] https://www.ibm.com/support/pages/node/318521#6
 * [32] http://docs.embarcadero.com/products/rad_studio/delphiAndcpp2009/HelpUpdate2/EN/html/devwin32/threadsusingthreadlocalvariables_xml.html
 * [33] http://docwiki.embarcadero.com/RADStudio/Sydney/en/Declspec(thread)
 * [34] https://digitalmars.com/ctg/ctgLanguageImplementation.html#declspec
 * [35] https://www.digitalmars.com/d/archives/c++/setjmp_longjmp_code_crashing_5923.html
 * [36] https://developer.arm.com/docs/dui0472/latest/compiler-specific-features/__declspecthread
 * [37] https://developer.arm.com/docs/dui0491/g/compiler-specific-features/__declspec-attributes
 * [38] https://lists.gnu.org/archive/html/bug-gnulib/2019-06/msg00063.html
 * [39] http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.172.9698&rep=rep1&type=pdf
 * [40] http://www.os2site.com/sw/dev/watcom/11.0c/c_readme.txt
 * [41] https://community.synopsys.com/s/article/Multiple-parse-warnings-for-defined-variables
 * [42] http://read.pudn.com/downloads259/doc/1193608/wr_compiler_error_messages_reference_5.6.pdf
 * [43] http://nwcc.sourceforge.net/features.html
 * [44] https://github.com/AbsInt/CompCert/issues/268
 * [45] https://en.cppreference.com/w/c/keyword/_Thread_local
 * [46] https://software.intel.com/en-us/forums/intel-c-compiler/topic/721059
 */

/* Apple Quirks
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
 * A useful SO answer points out that Apple Clang's
 * changes affect __has_extension(c_thread_local)[7].
 *
 * [1] https://stackoverflow.com/a/29929949
 * [2] https://en.wikipedia.org/wiki/Xcode#Xcode_7.0_-_11.x_(since_Free_On-Device_Development)
 * [3] https://community.intel.com/t5/Intel-C-Compiler/Mach-O-thread-local-storage/td-p/948267
 * [4] https://github.com/llvm/llvm-project/blob/llvmorg-3.0.0/clang/lib/Basic/Targets.cpp#L202
 * [5] https://github.com/llvm/llvm-project/blob/llvmorg-3.8.0/clang/lib/Basic/Targets.cpp#L232
 * [6] https://github.com/llvm/llvm-project/blob/llvmorg-7.0.0/clang/lib/Basic/Targets/OSTargets.h#L103
 * [7] https://stackoverflow.com/a/23850891
 */

/* Android Quirks
 *
 * Android apparently has issues compiling with NDK < r12[1].
 * According to SO, this appears to be fixed with NDK r12[2],
 * which upgrades Clang to 3.8[3][4].
 *
 * Note that the Android NDK does not provide a #define for
 * us to check[5]. We could check the Clang version, but NDK
 * r11 _also_ used Clang 3.8[6]. Instead, we must check for
 * NDK r15, which was upgraded to Clang 5.0[7].
 *
 * [1] https://stackoverflow.com/questions/27191214
 * [2] https://stackoverflow.com/a/27195324
 * [3] https://developer.android.com/ndk/downloads/revision_history
 * [4] https://github.com/android/ndk/blob/master/Changelogs/Changelog-r12.md#clang
 * [5] https://github.com/android/ndk/issues/407
 * [6] https://github.com/android/ndk/blob/master/Changelogs/Changelog-r11.md#clang
 * [7] https://github.com/android/ndk/blob/master/Changelogs/Changelog-r15.md#clang
 */

#ifndef TORSION_HAVE_CONFIG

#ifndef __has_extension
#  define __has_extension(x) 0
#endif

/* Some hackery to get apple versions. */
#if defined(__APPLE__) && defined(__MACH__)
#  if defined(__ENVIRONMENT_IPHONE_OS_VERSION_MIN_REQUIRED__)
/* Note: Clang serializes the short form version prior to 3.9.1. */
#    if __ENVIRONMENT_IPHONE_OS_VERSION_MIN_REQUIRED__ >= 100000 /* 10.0 */
#      define TORSION__APPLE_OS
#    endif
#  elif defined(__ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__)
#    if __ENVIRONMENT_MAC_OS_X_VERSION_MIN_REQUIRED__ >= 1070 /* 10.7 */
#      define TORSION__APPLE_OS
#    endif
#  endif
#endif

/* Detect TLS support. */
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  define TORSION_TLS_NONE
#elif defined(__DMC__)
#  if defined(_M_IX86) && __DMC__ >= 0x843 /* 8.43 */
#    define TORSION_TLS_MSVC
#  endif
#elif defined(__HP_aCC)
#  if _HP_aCC >= 55502 /* A.05.55.02 */
#    define TORSION_TLS_GNUC
#  endif
#elif defined(__HP_cc)
#  define TORSION_TLS_GNUC
#elif defined(__WATCOMC__)
#  if __WATCOMC__ >= 1200 /* 12.0 */
#    define TORSION_TLS_MSVC
#  endif
#elif defined(__DCC__)
#  if __VERSION_NUMBER__ >= 5600 /* 5.6 */
#    define TORSION_TLS_GNUC
#  endif
#elif defined(__NWCC__)
#  define TORSION_TLS_GNUC
#elif defined(__SUNPRO_C)
#  if __SUNPRO_C >= 0x590 /* 5.9 */
#    define TORSION_TLS_GNUC
#  endif
#elif defined(__INTEL_COMPILER)
#  if defined(__APPLE__) && defined(__MACH__)
#    if defined(TORSION__APPLE_OS) && __INTEL_COMPILER >= 1500 /* 15.0.0 */
#      define TORSION_TLS_GNUC
#    endif
#  elif __INTEL_COMPILER >= 800 /* 8.0.0 */
#    define TORSION_TLS_BOTH
#  endif
#elif defined(__clang__)
#  if defined(__apple_build_version__)
#    if defined(TORSION__APPLE_OS) && __apple_build_version__ >= 8000038 /* 800.0.38 */
#      define TORSION_TLS_GNUC
#    endif
#  elif __has_extension(c_thread_local)
#    if defined(__APPLE__) && defined(__MACH__)
#      ifdef TORSION__APPLE_OS
#        define TORSION_TLS_GNUC
#      endif
#    elif defined(__ANDROID__)
#      if ((__clang_major__ << 16) + __clang_minor__ >= 0x50000) /* 5.0 */
#        define TORSION_TLS_GNUC
#      endif
#    else
#      define TORSION_TLS_BOTH
#    endif
#  endif
#elif defined(__xlC__)
/* Newer XL C versions are based on clang and should be caught above. */
#  if defined(__linux__)
#    if __xlC__ >= 0x0800 /* 8.0.0 */
#      define TORSION_TLS_GNUC
#    endif
#  elif defined(_AIX)
#    if __xlC__ >= 0x0A01 /* 10.1.0 */
#      define TORSION_TLS_GNUC
#    endif
#  endif
#elif defined(__CC_ARM)
/* Newer ARM CC versions are based on clang and should be caught above. */
#  if __ARMCC_VERSION >= 510000 /* 5.1 */
#    define TORSION_TLS_BOTH
#  endif
#elif defined(__BORLANDC__)
/* Newer C++ Builder versions are based on clang and should be caught above. */
#  if __BORLANDC__ >= 0x613 /* C++ Builder 2009 */
#    define TORSION_TLS_BOTH
#  endif
#elif defined(_MSC_VER)
#  if _MSC_VER >= 1300 /* VS .NET 2002 (7.0) */
#    define TORSION_TLS_MSVC
#  endif
#elif defined(__GNUC__) && defined(__GNUC_MINOR__)
#  if defined(__ELF__) && (defined(__alpha__) || defined(__i386__)  \
                        || defined(__x86_64__) || defined(__ia64__) \
                        || defined(__s390__) || defined(__s390x__))
#    if ((__GNUC__ << 16) + __GNUC_MINOR__ >= 0x30003) /* 3.3 */
#      define TORSION_TLS_GNUC
#    endif
#  else
#    if ((__GNUC__ << 16) + __GNUC_MINOR__ >= 0x40003) /* 4.3 */
#      define TORSION_TLS_GNUC
#    endif
#  endif
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  include <limits.h> /* <stdc-predef.h> */
#  ifndef __STDC_NO_THREADS__
#    define TORSION_TLS_STDC
#  endif
#endif

#ifdef TORSION_TLS_BOTH
#  ifdef _WIN32
#    define TORSION_TLS_MSVC
#  else
#    define TORSION_TLS_GNUC
#  endif
#endif

/* Pick thread-local keyword. */
#if defined(TORSION_TLS_NONE)
#  define TORSION_HAVE_TLS
#  define TORSION_TLS
#elif defined(TORSION_TLS_MSVC)
#  define TORSION_HAVE_TLS
#  define TORSION_TLS __declspec(thread)
#elif defined(TORSION_TLS_GNUC)
#  define TORSION_HAVE_TLS
#  define TORSION_TLS __thread
#elif defined(TORSION_TLS_STDC)
#  define TORSION_HAVE_TLS
#  define TORSION_TLS _Thread_local
#else
#  define TORSION_TLS
#endif

#else /* TORSION_HAVE_CONFIG */

/* Pick thread-local keyword. */
#ifdef TORSION_HAVE_TLS
#  if defined(__EMSCRIPTEN__) || defined(__wasm__)
#    define TORSION_TLS
#  elif defined(_WIN32) && !defined(__MINGW32__)
#    define TORSION_TLS __declspec(thread)
#  else
#    define TORSION_TLS __thread
#  endif
#else
#  define TORSION_TLS
#endif

#endif /* TORSION_HAVE_CONFIG */

/* Configure fallback. */
#ifndef TORSION_HAVE_TLS
#  if defined(__EMSCRIPTEN__) || defined(__wasm__)
/* No threads support. */
#elif (defined(__unix) || defined(__unix__))    \
   || (defined(__APPLE__) && defined(__MACH__))
#    ifdef _REENTRANT
#      define TORSION_HAVE_PTHREAD
#    endif
#  endif
#endif

/* Allow some overrides (for testing). */
#ifdef TORSION_NO_TLS
#  undef TORSION_HAVE_TLS
#endif

#ifdef TORSION_NO_PTHREAD
#  undef TORSION_HAVE_PTHREAD
#endif

#endif /* _TORSION_TLS_H */
