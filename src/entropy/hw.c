/*!
 * hw.c - hardware entropy for libtorsion
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 *
 * Resources:
 *   https://en.wikipedia.org/wiki/Time_Stamp_Counter
 *   https://en.wikipedia.org/wiki/CPUID
 *   https://en.wikipedia.org/wiki/RDRAND
 *
 * Windows (x86 / x64):
 *   https://docs.microsoft.com/en-us/cpp/intrinsics/rdtsc
 *   https://docs.microsoft.com/en-us/cpp/intrinsics/cpuid-cpuidex
 *   https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_rdrand32_step
 *   https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_rdrand64_step
 *   https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_rdseed32_step
 *   https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_rdseed64_step
 *
 * Windows (arm64):
 *   https://docs.microsoft.com/en-us/cpp/intrinsics/arm64-intrinsics
 *
 * x86{,-64}:
 *   https://www.felixcloutier.com/x86/rdtsc
 *   https://www.felixcloutier.com/x86/cpuid
 *   https://www.felixcloutier.com/x86/rdrand
 *   https://www.felixcloutier.com/x86/rdseed
 *
 * ARMv8.5-A (rndr, rndrrs, pmccntr_el0):
 *   https://developer.arm.com/documentation/dui0068/b/ARM-Instruction-Reference/Miscellaneous-ARM-instructions/MRS
 *   https://developer.arm.com/documentation/ddi0595/2021-03/AArch64-Registers/RNDR--Random-Number
 *   https://developer.arm.com/documentation/ddi0595/2021-03/AArch64-Registers/RNDRRS--Reseeded-Random-Number
 *   https://developer.arm.com/documentation/ddi0595/2021-03/AArch64-Registers/PMCCNTR-EL0--Performance-Monitors-Cycle-Count-Register
 *   https://developer.arm.com/documentation/ddi0595/2021-03/AArch64-Registers/ID-AA64ISAR0-EL1--AArch64-Instruction-Set-Attribute-Register-0
 *   https://developer.arm.com/documentation/ddi0595/2021-03/AArch64-Registers/ID-DFR0-EL1--AArch32-Debug-Feature-Register-0
 *
 * POWER9/POWER10 (darn, mftb):
 *   https://www.docdroid.net/tWT7hjD/powerisa-v30-pdf
 *   https://openpowerfoundation.org/?resource_lib=power-isa-version-3-0
 */

/**
 * Hardware Entropy
 *
 * One simple source of hardware entropy is the current cycle
 * count. This is accomplished via RDTSC on x86 CPUs. We only
 * use RDTSC if there is an instrinsic for it (win32) or if
 * the compiler supports inline ASM (gcc/clang).
 *
 * For ARM and PPC, we use PMCCNTR_EL0 and MFTB respectively.
 *
 * For other hardware, we fallback to whatever system clocks
 * are available. See `hrt.c` for a list of functions.
 *
 * The CPUID instruction can serve as good source of "static"
 * entropy for seeding (see env.c).
 *
 * x86{,-64} also offers hardware entropy in the form of RDRAND
 * and RDSEED. There are concerns that these instructions may
 * be backdoored in some way. This is not an issue as we only
 * use hardware entropy to supplement our full entropy pool.
 *
 * On POWER9 and POWER10, the `darn` (Deliver A Random Number)
 * instruction is available. We have `torsion_rdrand` as well
 * as `torsion_rdseed` return the output of `darn` if this is
 * the case.
 *
 * ARMv8.5-A provides new system registers (RNDR and RNDRRS)
 * to be used with the MRS instruction. Similar to `darn`, we
 * have `torsion_{rdrand,rdseed}` output the proper values.
 *
 * For other hardware, torsion_rdrand and torsion_rdseed are
 * no-ops returning zero. torsion_has_rd{rand,seed} MUST be
 * checked before calling torsion_rd{rand,seed}.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "entropy.h"

#undef HAVE_RDTSC
#undef HAVE_CPUIDEX
#undef HAVE_RDRAND
#undef HAVE_RDRAND32
#undef HAVE_RDRAND64
#undef HAVE_RDSEED
#undef HAVE_RDSEED32
#undef HAVE_RDSEED64
#undef HAVE_ASM_RDTSC
#undef HAVE_READSTATUSREG
#undef HAVE_ASM_INTEL
#undef HAVE_ASM_X86
#undef HAVE_ASM_X64
#undef HAVE_ASM_ARM64
#undef HAVE_ASM_PPC64
#undef HAVE_GETAUXVAL
#undef HAVE_ELF_AUX_INFO
#undef HAVE_POWER_SET
#undef HAVE_AUXVAL
#undef HAVE_PERFMON /* Define if ARM FEAT_PMUv3 is supported. */

/* Detect intrinsic and ASM support. */
#if defined(_MSC_VER)
#  if defined(_M_IX86) || defined(_M_AMD64) || defined(_M_X64)
#    if _MSC_VER >= 1400 /* VS 2005 */
#      include <intrin.h> /* __cpuidex, __rdtsc */
#      pragma intrinsic(__rdtsc)
#      define HAVE_RDTSC
#    endif
#    if _MSC_VER >= 1600 /* VS 2010 */
#      pragma intrinsic(__cpuidex)
#      define HAVE_CPUIDEX
#    endif
#    if _MSC_VER >= 1700 /* VS 2012 */
#      include <immintrin.h> /* _rd{rand,seed}{32,64}_step */
#      define HAVE_RDRAND
#      if defined(_M_AMD64) || defined(_M_X64)
#        define HAVE_RDRAND64
#      else
#        define HAVE_RDRAND32
#      endif
#    endif
#    if _MSC_VER >= 1800 /* VS 2013 */
#      define HAVE_RDSEED
#      if defined(_M_AMD64) || defined(_M_X64)
#        define HAVE_RDSEED64
#      else
#        define HAVE_RDSEED32
#      endif
#    endif
#    if !defined(HAVE_RDTSC) && defined(_M_IX86)
#      define HAVE_ASM_RDTSC /* fallback to _asm */
#    endif
#  elif defined(_M_ARM64)
#    include <intrin.h> /* _ReadStatusReg */
#    define HAVE_READSTATUSREG
#  endif
#elif (defined(__GNUC__) && __GNUC__ >= 4) || defined(__IBM_GCC_ASM)
#  if defined(__amd64__) || defined(__x86_64__)
#    define HAVE_ASM_INTEL
#    define HAVE_ASM_X64
#  elif defined(__i386__)
#    define HAVE_ASM_INTEL
#    define HAVE_ASM_X86
#  elif defined(__aarch64__)
#    define HAVE_ASM_ARM64
#  elif defined(__powerpc64__) || defined(_ARCH_PPC64) || defined(__PPC64__)
#    define HAVE_ASM_PPC64
#  endif
#endif

/* Some insanity to detect features at runtime. */
#if defined(HAVE_ASM_ARM64) || defined(HAVE_ASM_PPC64)
#  if defined(__GLIBC__) && defined(__GLIBC_PREREQ)
#    define TORSION_GLIBC_PREREQ __GLIBC_PREREQ
#  else
#    define TORSION_GLIBC_PREREQ(maj, min) 0
#  endif
#  if TORSION_GLIBC_PREREQ(2, 16)
#    include <sys/auxv.h> /* getauxval */
#    define HAVE_GETAUXVAL
#    define HAVE_AUXVAL
#  elif defined(__FreeBSD__)
#    include <sys/param.h>
#    if defined(__FreeBSD_version) && __FreeBSD_version >= 1200000 /* 12.0 */
#      include <sys/auxv.h> /* elf_aux_info */
#      define HAVE_ELF_AUX_INFO
#      define HAVE_AUXVAL
#    endif
#  elif defined(HAVE_ASM_PPC64) && defined(_AIX53) && !defined(__PASE__)
#    include <sys/systemcfg.h> /* __power_set */
#    ifndef __power_set
#      define __power_set(x) (_system_configuration.implementation & (x))
#    endif
#    define HAVE_POWER_SET
#  endif
#  ifdef HAVE_AUXVAL
#    ifndef AT_HWCAP
#      define AT_HWCAP 16
#    endif
#    ifndef AT_HWCAP2
#      define AT_HWCAP2 26
#    endif
#  endif
#endif

/*
 * Auxiliary Value
 */

#if defined(HAVE_GETAUXVAL)
#  define torsion_auxval getauxval
#elif defined(HAVE_ELF_AUX_INFO)
__attribute__((unused)) static unsigned long
torsion_auxval(unsigned long type) {
  unsigned long val;

  if (elf_aux_info(type, &val, sizeof(val)) != 0)
    return 0;

  return val;
}
#endif

/*
 * Timestamp Counter
 */

uint64_t
torsion_rdtsc(void) {
#if defined(HAVE_RDTSC)
  return __rdtsc();
#elif defined(HAVE_ASM_RDTSC)
  _asm rdtsc
#elif defined(HAVE_READSTATUSREG) && defined(HAVE_PERFMON)
  return _ReadStatusReg(0x5ce8); /* ARM64_PMCCNTR_EL0 */
#elif defined(HAVE_ASM_X86)
  uint64_t ts;

  __asm__ __volatile__ (
    "rdtsc\n"
    : "=A" (ts)
  );

  return ts;
#elif defined(HAVE_ASM_X64)
  uint64_t lo, hi;

  __asm__ __volatile__ (
    "rdtsc\n"
    : "=a" (lo),
      "=d" (hi)
  );

  return (hi << 32) | lo;
#elif defined(HAVE_ASM_ARM64) && defined(HAVE_PERFMON)
  uint64_t ts;

  /* Note that `mrs %0, pmccntr_el0` can be
   * spelled out as:
   *
   *   .inst (0xd5200000 | 0x1b9d00 | %0)
   *              |            |       |
   *             mrs        sysreg    reg
   *
   * Requires FEAT_PMUv3. We _could_ check:
   *
   *   ((ID_DFR0_EL1 >> 24) & 15) >= 3
   *
   * But that is an EL1 register and the
   * kernel doesn't emulate the debug
   * feature registers (yet).
   */
  __asm__ __volatile__ (
    "mrs %0, s3_3_c9_c13_0\n" /* PMCCNTR_EL0 */
    : "=r" (ts)
  );

  return ts;
#elif defined(HAVE_ASM_PPC64)
  uint64_t ts;

  /* mftb  %0 = mftb  %0, 268
   *          = mfspr %0, 268
   * mftbu %0 = mftb  %0, 269
   *          = mfspr %0, 269
   *
   * mfspr available since 2.01
   * (excludes 601 and POWER3).
   */
  __asm__ __volatile__ (
    "mfspr %0, 268\n" /* mftb %0 */
    : "=r" (ts)
  );

  return ts;
#else
  /* Fall back to high-resolution time. */
  return torsion_hrtime();
#endif
}

/*
 * CPUID
 */

int
torsion_has_cpuid(void) {
#if defined(HAVE_CPUIDEX)
  return 1;
#elif defined(HAVE_ASM_X86)
  uint32_t ax, bx;

  __asm__ __volatile__ (
    "pushfl\n"
    "pushfl\n"
    "popl %k0\n"
    "movl %k0, %k1\n"
    "xorl $0x00200000, %k0\n"
    "pushl %k0\n"
    "popfl\n"
    "pushfl\n"
    "popl %k0\n"
    "popfl\n"
    : "=&r" (ax),
      "=&r" (bx)
    :
    : "cc"
  );

  return ((ax ^ bx) >> 21) & 1;
#elif defined(HAVE_ASM_X64)
  return 1;
#else
  return 0;
#endif
}

void
torsion_cpuid(uint32_t *a,
              uint32_t *b,
              uint32_t *c,
              uint32_t *d,
              uint32_t leaf,
              uint32_t subleaf) {
#if defined(HAVE_CPUIDEX)
  unsigned int regs[4];

  __cpuidex((int *)regs, leaf, subleaf);

  *a = regs[0];
  *b = regs[1];
  *c = regs[2];
  *d = regs[3];
#elif defined(HAVE_ASM_X86)
  if (torsion_has_cpuid()) {
    __asm__ __volatile__ (
      "xchgl %%ebx, %k1\n"
      "cpuid\n"
      "xchgl %%ebx, %k1\n"
      : "=a" (*a), "=&r" (*b), "=c" (*c), "=d" (*d)
      : "0" (leaf), "2" (subleaf)
    );
  } else {
    *a = 0;
    *b = 0;
    *c = 0;
    *d = 0;
  }
#elif defined(HAVE_ASM_X64)
  __asm__ __volatile__ (
    "cpuid\n"
    : "=a" (*a), "=b" (*b), "=c" (*c), "=d" (*d)
    : "0" (leaf), "2" (subleaf)
  );
#else
  (void)leaf;
  (void)subleaf;

  *a = 0;
  *b = 0;
  *c = 0;
  *d = 0;
#endif
}

/*
 * RDRAND/RDSEED
 */

int
torsion_has_rdrand(void) {
#if defined(HAVE_ASM_INTEL) && defined(__RDRND__)
  /* Explicitly built with RDRAND support (-mrdrnd). */
  return 1;
#elif defined(HAVE_RDRAND) || defined(HAVE_ASM_INTEL)
  uint32_t eax, ebx, ecx, edx;

  torsion_cpuid(&eax, &ebx, &ecx, &edx, 0, 0);

  if (eax < 1)
    return 0;

  torsion_cpuid(&eax, &ebx, &ecx, &edx, 1, 0);

  return (ecx >> 30) & 1;
#elif defined(HAVE_ASM_ARM64) && defined(__ARM_FEATURE_RNG)
  /* Explicitly built with ARM RNG support (-march=armv8.5-a+rng). */
  return 1;
#elif defined(HAVE_ASM_ARM64) && defined(HAVE_AUXVAL)
  /* Bit 16 = RNG support (HWCAP2_RNG) */
  if ((torsion_auxval(AT_HWCAP2) >> 16) & 1)
    return 1;

  /* Bit 11 = MRS emulation (HWCAP_CPUID) */
  /* https://www.kernel.org/doc/html/latest/arm64/cpu-feature-registers.html */
  if ((torsion_auxval(AT_HWCAP) >> 11) & 1) {
    uint64_t isar0;

    /* Note that `mrs %0, id_aa64isar0_el1` can be
     * spelled out as:
     *
     *   .inst (0xd5200000 | 0x180600 | %0)
     *              |            |       |
     *             mrs        sysreg    reg
     */
    __asm__ __volatile__ (
      "mrs %0, s3_0_c0_c6_0\n" /* ID_AA64ISAR0_EL1 */
      : "=r" (isar0)
    );

    /* Bits 63-60 = RNDR (0b0001) */
    return (isar0 >> 60) >= 1;
  }

  return 0;
#elif defined(HAVE_ASM_PPC64) && (defined(_ARCH_PWR9) || defined(_ARCH_PWR10))
  /* Explicitly built for POWER9 (-mcpu=power9 or -mpower9-vector). */
  return 1;
#elif defined(HAVE_ASM_PPC64) && defined(HAVE_AUXVAL)
  /* Bit 21 = DARN support (PPC_FEATURE2_DARN) */
  return (torsion_auxval(AT_HWCAP2) >> 21) & 1;
#elif defined(HAVE_ASM_PPC64) && defined(HAVE_POWER_SET)
  /* Check for POWER9 or greater. */
  return __power_set(0xffffffffU << 17) != 0;
#else
  return 0;
#endif
}

int
torsion_has_rdseed(void) {
#if defined(HAVE_ASM_INTEL) && defined(__RDSEED__)
  /* Explicitly built with RDSEED support (-mrdseed). */
  return 1;
#elif defined(HAVE_RDSEED) || defined(HAVE_ASM_INTEL)
  uint32_t eax, ebx, ecx, edx;

  torsion_cpuid(&eax, &ebx, &ecx, &edx, 0, 0);

  if (eax < 7)
    return 0;

  torsion_cpuid(&eax, &ebx, &ecx, &edx, 7, 0);

  return (ebx >> 18) & 1;
#elif defined(HAVE_ASM_ARM64)
  return torsion_has_rdrand();
#elif defined(HAVE_ASM_PPC64)
  return torsion_has_rdrand();
#else
  return 0;
#endif
}

uint64_t
torsion_rdrand(void) {
#if defined(HAVE_RDRAND32)
  unsigned int lo, hi;
  int i;

  for (i = 0; i < 10; i++) {
    if (_rdrand32_step(&lo))
      break;
  }

  for (i = 0; i < 10; i++) {
    if (_rdrand32_step(&hi))
      break;
  }

  return ((uint64_t)hi << 32) | lo;
#elif defined(HAVE_RDRAND64)
  unsigned __int64 r;
  int i;

  for (i = 0; i < 10; i++) {
    if (_rdrand64_step(&r))
      break;
  }

  return r;
#elif defined(HAVE_ASM_X86)
  /* Borrowed from Bitcoin Core. */
  uint32_t lo, hi;
  uint8_t ok;
  int i;

  for (i = 0; i < 10; i++) {
    __asm__ __volatile__ (
      ".byte 0x0f, 0xc7, 0xf0\n" /* rdrand %eax */
      "setc %b1\n"
      : "=a" (lo), "=q" (ok)
      :
      : "cc"
    );

    if (ok)
      break;
  }

  for (i = 0; i < 10; i++) {
    __asm__ __volatile__ (
      ".byte 0x0f, 0xc7, 0xf0\n" /* rdrand %eax */
      "setc %b1\n"
      : "=a" (hi), "=q" (ok)
      :
      : "cc"
    );

    if (ok)
      break;
  }

  return ((uint64_t)hi << 32) | lo;
#elif defined(HAVE_ASM_X64)
  /* Borrowed from Bitcoin Core. */
  uint64_t r;
  uint8_t ok;
  int i;

  for (i = 0; i < 10; i++) {
    __asm__ __volatile__ (
      ".byte 0x48, 0x0f, 0xc7, 0xf0\n" /* rdrand %rax */
      "setc %b1\n"
      : "=a" (r), "=q" (ok)
      :
      : "cc"
    );

    if (ok)
      break;
  }

  return r;
#elif defined(HAVE_ASM_ARM64)
  uint64_t r;
  uint32_t ok;
  int i;

  for (i = 0; i < 10; i++) {
    /* Note that `mrs %0, rndr` can be spelled out as:
     *
     *   .inst (0xd5200000 | 0x1b2400 | %0)
     *              |            |       |
     *             mrs        sysreg    reg
     *
     * Though, this presents a difficulty in that %0
     * will be expanded to `x<n>` instead of `<n>`.
     * That is to say, %0 becomes x3 instead of 3.
     *
     * We can solve this with some crazy macros like
     * the linux kernel does, but it's probably not
     * worth the effort.
     */
    __asm__ __volatile__ (
      "mrs %0, s3_3_c2_c4_0\n" /* RNDR */
      "cset %w1, ne\n"
      : "=r" (r), "=r" (ok)
      :
      : "cc"
    );

    if (ok)
      break;
  }

  return r;
#elif defined(HAVE_ASM_PPC64)
  uint64_t r;
  int i;

  for (i = 0; i < 10; i++) {
    /* Darn modes:
     *
     *   0 = 32 bit (conditioned)
     *   1 = 64 bit (conditioned)
     *   2 = 64 bit (raw)
     *   3 = reserved
     *
     * Spelling below was taken from the linux kernel
     * (after stripping out a load of preprocessor).
     */
    __asm__ __volatile__ (
      ".long (0x7c0005e6 | (%0 << 21) | (1 << 16))\n" /* darn %0, 1 */
      : "=r" (r)
    );

    if (r != UINT64_MAX)
      break;

    r = 0;
  }

  return r;
#else
  return 0;
#endif
}

uint64_t
torsion_rdseed(void) {
#if defined(HAVE_RDSEED32)
  unsigned int lo, hi;

  for (;;) {
    if (_rdseed32_step(&lo))
      break;

    _asm { rep nop }
  }

  for (;;) {
    if (_rdseed32_step(&hi))
      break;

    _asm { rep nop }
  }

  return ((uint64_t)hi << 32) | lo;
#elif defined(HAVE_RDSEED64)
  unsigned __int64 r;

  for (;;) {
    if (_rdseed64_step(&r))
      break;

    _mm_pause();
  }

  return r;
#elif defined(HAVE_ASM_X86)
  /* Borrowed from Bitcoin Core. */
  uint32_t lo, hi;
  uint8_t ok;

  for (;;) {
    __asm__ __volatile__ (
      ".byte 0x0f, 0xc7, 0xf8\n" /* rdseed %eax */
      "setc %b1\n"
      : "=a" (lo), "=q" (ok)
      :
      : "cc"
    );

    if (ok)
      break;

    __asm__ __volatile__ (
      ".byte 0xf3, 0x90\n" /* pause */
    );
  }

  for (;;) {
    __asm__ __volatile__ (
      ".byte 0x0f, 0xc7, 0xf8\n" /* rdseed %eax */
      "setc %b1\n"
      : "=a" (hi), "=q" (ok)
      :
      : "cc"
    );

    if (ok)
      break;

    __asm__ __volatile__ (
      ".byte 0xf3, 0x90\n" /* pause */
    );
  }

  return ((uint64_t)hi << 32) | lo;
#elif defined(HAVE_ASM_X64)
  /* Borrowed from Bitcoin Core. */
  uint64_t r;
  uint8_t ok;

  for (;;) {
    __asm__ __volatile__ (
      ".byte 0x48, 0x0f, 0xc7, 0xf8\n" /* rdseed %rax */
      "setc %b1\n"
      : "=a" (r), "=q" (ok)
      :
      : "cc"
    );

    if (ok)
      break;

    __asm__ __volatile__ ("pause\n");
  }

  return r;
#elif defined(HAVE_ASM_ARM64)
  uint64_t r;
  uint32_t ok;

  for (;;) {
    /* Note that `mrs %0, rndrrs` can be spelled out as:
     *
     *   .inst (0xd5200000 | 0x1b2420 | %0)
     *              |            |       |
     *             mrs        sysreg    reg
     */
    __asm__ __volatile__ (
      "mrs %0, s3_3_c2_c4_1\n" /* RNDRRS */
      "cset %w1, ne\n"
      : "=r" (r), "=r" (ok)
      :
      : "cc"
    );

    if (ok)
      break;

    __asm__ __volatile__ ("yield\n");
  }

  return r;
#elif defined(HAVE_ASM_PPC64)
  uint64_t r;

  for (;;) {
    __asm__ __volatile__ (
      ".long (0x7c0005e6 | (%0 << 21) | (2 << 16))\n" /* darn %0, 2 */
      : "=r" (r)
    );

    if (r != UINT64_MAX)
      break;

    /* https://stackoverflow.com/questions/5425506 */
    __asm__ __volatile__ ("or 27, 27, 27\n" ::: "cc");
  }

  return r;
#else
  return 0;
#endif
}

/*
 * Hardware Entropy
 */

int
torsion_hwrand(void *dst, size_t size) {
  unsigned char *data = (unsigned char *)dst;
  int has_rdrand = torsion_has_rdrand();
  int has_rdseed = torsion_has_rdseed();
  uint64_t x;
  int i;

  if (!has_rdrand && !has_rdseed) {
    if (size > 0)
      memset(dst, 0, size);

    return 0;
  }

  while (size > 0) {
    if (has_rdseed) {
      x = torsion_rdseed();
    } else {
      x = 0;

      /* Idea from Bitcoin Core: force rdrand to reseed. */
      for (i = 0; i < 1024; i++)
        x ^= torsion_rdrand();
    }

    if (size < sizeof(x)) {
      memcpy(data, &x, size);
      break;
    }

    memcpy(data, &x, sizeof(x));

    data += sizeof(x);
    size -= sizeof(x);
  }

  return 1;
}
