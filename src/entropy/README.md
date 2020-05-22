# entropy/

This directory contains various entropy sources for seeding an RNG.

- `env.c` - Manual entropy gathering (`/proc`, `sysctl`,
  `HKEY_PERFORMANCE_DATA`).
- `hw.c` - Hardware entropy sources (`cpuid`, `rdtsc`, `rdrand`, `rdseed`).
- `sys.c` - OS/System entropy sources (`getrandom`, `getentropy`,
  `sysctl(kern.arandom)`, `/dev/{u,}random`, etc).

All three combined can ensure that a seeded RNG is not compromised in any way,
as an adversary would have to backdoor all of our entropy sources (very
unlikely). Seeding an RNG in this way results in _very_ strong randomness.

We support windows, unix¹, and other oddball OSes like fuchsia and vxworks.

Emscripten (WASM, asm.js), WASI and CloudABI are also supported.

See each file's comments for a more in-depth description of features.

## Footnotes

1. Including (but not limited to): Linux, OSX, iOS, OpenBSD, FreeBSD, NetBSD,
DragonFly BSD, Solaris/Illumos.
