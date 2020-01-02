#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#ifdef _WIN32
/* For SecureZeroMemory (actually defined in winbase.h). */
#include <windows.h>
#endif

static void
cleanse(void *ptr, size_t len) {
#if defined(_WIN32)
  /* https://github.com/jedisct1/libsodium/blob/3b26a5c/src/libsodium/sodium/utils.c#L112 */
  SecureZeroMemory(ptr, len);
#elif defined(__GNUC__)
  /* https://github.com/torvalds/linux/blob/37d4e84/include/linux/string.h#L233 */
  /* https://github.com/torvalds/linux/blob/37d4e84/include/linux/compiler-gcc.h#L21 */
  memset(ptr, 0, len);
  __asm__ __volatile__("": :"r"(ptr) :"memory");
#else
  /* http://www.daemonology.net/blog/2014-09-04-how-to-zero-a-buffer.html */
  static void *(*const volatile memset_ptr)(void *, int, size_t) = memset;
  (memset_ptr)(ptr, 0, len);
#endif
}

static uint32_t
is_zero(const unsigned char *a, size_t size) {
  size_t i = 0;
  uint32_t z = 0;

  for (i = 0; i < size; i++)
    z |= (uint32_t)a[i];

  return (z - 1) >> 31;
}

static uint32_t
less_than(const unsigned char *a,
          const unsigned char *b,
          int size,
          int endian) {
  int le = (endian == -1);
  int32_t eq = -1;
  int32_t gt = 0;
  int i = le ? size - 1 : 0;

  for (; le ? i >= 0 : i < size; le ? i-- : i++) {
    int32_t x = a[i];
    int32_t y = b[i];

    gt = (~eq & gt) | (eq & ((x - y) >> 31));
    eq = eq & (((x ^ y) - 1) >> 31);
  }

  return (uint32_t)(~eq & 1 & gt);
}
