cmake_minimum_required(VERSION 3.4)
project(libtorsion LANGUAGES C)

include(GNUInstallDirs)
include(CheckCCompilerFlag)
include(CheckIncludeFile)

set(CMAKE_C_VISIBILITY_PRESET default)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 90)

option(TORSION_ENABLE_64BIT "Use 64-bit backends if available" ON)
option(TORSION_ENABLE_ASM "Use inline x86-64 assembly if available" ON)
option(TORSION_ENABLE_INT128 "Use __int128 if available" ON)
option(TORSION_ENABLE_ASSERT "Enable assertions" ON)
option(TORSION_ENABLE_CRT "Enable chinese remainder theorem for RSA" OFF)
option(TORSION_ENABLE_CTIME "Enable valgrind constant-time test" OFF)
option(TORSION_ENABLE_LIBSECP256K1 "Enable libsecp256k1 field element backend" OFF)
option(TORSION_ENABLE_RNG "Enable RNG" ON)
option(TORSION_ENABLE_TESTS "Enable tests" ON)
option(TORSION_ENABLE_VERIFY "Enable scalar bounds checks" OFF)

if(TORSION_ENABLE_CTIME)
  set(TORSION_ENABLE_ASSERT OFF)
  set(TORSION_ENABLE_VERIFY OFF)
endif()

if(MSVC)
check_c_compiler_flag(/W3 TORSION_HAS_W3)
check_c_compiler_flag(/wd4146 TORSION_HAS_WMINUS_UNSIGNED)
check_c_compiler_flag(/wd4244 TORSION_HAS_WCONV_DATA_LOSS1)
check_c_compiler_flag(/wd4267 TORSION_HAS_WCONV_DATA_LOSS2)
check_c_compiler_flag(/wd4477 TORSION_HAS_WPRINTF_FORMAT)
check_c_compiler_flag(/O2 TORSION_HAS_O2)
else()
check_c_compiler_flag(-g TORSION_HAS_G)
check_c_compiler_flag(-pedantic TORSION_HAS_PEDANTIC)
check_c_compiler_flag(-Wall TORSION_HAS_WALL)
check_c_compiler_flag(-Wextra TORSION_HAS_WEXTRA)
check_c_compiler_flag(-Wshadow TORSION_HAS_WSHADOW)
check_c_compiler_flag(-Wcast-align TORSION_HAS_WCAST_ALIGN)
check_c_compiler_flag(-Wno-implicit-fallthrough TORSION_HAS_WFALLTHROUGH)
check_c_compiler_flag(-Wno-declaration-after-statement TORSION_HAS_WDECL)
check_c_compiler_flag(-Wno-overlength-strings TORSION_HAS_WOVERLENGTH)
check_c_compiler_flag(-Wno-long-long TORSION_HAS_WLONGLONG)
check_c_compiler_flag(-O3 TORSION_HAS_O3)
endif()

if("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8")
  set(TORSION_HAS_64BIT TRUE)
endif()

check_c_source_compiles("
#ifndef __EMSCRIPTEN__
#  error \"compiler is not emscripten\"
#endif
int main(void) { return 0; }
" TORSION_EMSCRIPTEN)

check_c_source_compiles("
#include <stdint.h>
#if !defined(__amd64__) && !defined(__x86_64__)
#  error \"not an x64 platform\"
#endif
#if defined(__CYGWIN__) || defined(__MINGW32__)
#  error \"no asm on unix lookalikes\"
#endif
int main(void) {
  uint32_t stream[8], state[8];
  __asm__ __volatile__(
    \"movups (%%rsi), %%xmm0;\"
    \"movups 16(%%rsi), %%xmm1;\"
    \"pxor %%xmm1, %%xmm0;\"
    \"movups %%xmm0,(%%rdi);\"
    \"movups %%xmm1,16(%%rdi);\"
    :
    : \"D\" (stream), \"S\" (state)
    : \"xmm0\", \"xmm1\", \"cc\", \"memory\"
  );
  return 0;
}
" TORSION_HAS_ASM)

check_c_source_compiles("
typedef signed __int128 int128_t;
typedef unsigned __int128 uint128_t;
int main(void) {
  uint128_t r;
  r = ((uint128_t)1) << 64;
  r *= 113;
  r >>= 65;
  return r & 1;
}
" TORSION_HAS_INT128)

check_include_file(valgrind/memcheck.h TORSION_HAS_VALGRIND)

list(APPEND torsion_defines TORSION_HAVE_CONFIG)

if(TORSION_ENABLE_64BIT AND TORSION_HAS_64BIT AND NOT TORSION_EMSCRIPTEN)
  list(APPEND torsion_defines TORSION_HAVE_64BIT)
  if(TORSION_ENABLE_ASM AND TORSION_HAS_ASM)
    list(APPEND torsion_defines TORSION_HAVE_ASM_X64)
  endif()
  if(TORSION_ENABLE_INT128 AND TORSION_HAS_INT128)
    list(APPEND torsion_defines TORSION_HAVE_INT128)
  endif()
endif()

if(NOT TORSION_ENABLE_ASSERT)
  list(APPEND torsion_defines TORSION_NO_ASSERT)
endif()

if(TORSION_ENABLE_CRT)
  list(APPEND torsion_defines TORSION_USE_CRT)
endif()

if(TORSION_ENABLE_CTIME)
  if(NOT TORSION_HAS_VALGRIND)
    message(FATAL_ERROR "valgrind must be installed for constant time test")
  endif()
  list(APPEND torsion_defines TORSION_VALGRIND)
endif()

if(TORSION_ENABLE_LIBSECP256K1)
  list(APPEND torsion_defines TORSION_USE_LIBSECP256K1)
endif()

if(TORSION_ENABLE_RNG)
  list(APPEND torsion_defines TORSION_HAVE_RNG)
endif()

if(TORSION_ENABLE_TESTS)
  list(APPEND torsion_defines TORSION_TEST)
  if(TORSION_HAS_G)
    list(APPEND torsion_cflags -g)
  endif()
endif()

if(TORSION_ENABLE_VERIFY)
  if(NOT TORSION_ENABLE_TESTS)
    message(FATAL_ERROR "scalar bounds verification for testing only")
  endif()
  list(APPEND torsion_defines TORSION_VERIFY)
endif()

if(MSVC)
  if(TORSION_HAS_W3)
    list(APPEND torsion_cflags /W3)
  endif()
  if(TORSION_HAS_WMINUS_UNSIGNED)
    list(APPEND torsion_cflags /wd4146)
  endif()
  if(TORSION_HAS_WCONV_DATA_LOSS1)
    list(APPEND torsion_cflags /wd4244)
  endif()
  if(TORSION_HAS_WCONV_DATA_LOSS2)
    list(APPEND torsion_cflags /wd4267)
  endif()
  if(TORSION_HAS_WPRINTF_FORMAT)
    list(APPEND torsion_cflags /wd4477)
  endif()
  if(TORSION_HAS_O2)
    list(APPEND torsion_cflags /O2)
  endif()
else()
  if(TORSION_HAS_PEDANTIC)
    list(APPEND torsion_cflags -pedantic)
  endif()
  if(TORSION_HAS_WALL)
    list(APPEND torsion_cflags -Wall)
  endif()
  if(TORSION_HAS_WEXTRA)
    list(APPEND torsion_cflags -Wextra)
  endif()
  if(TORSION_HAS_WSHADOW)
    list(APPEND torsion_cflags -Wshadow)
  endif()
  if(TORSION_HAS_WCAST_ALIGN)
    list(APPEND torsion_cflags -Wcast-align)
  endif()
  if(TORSION_HAS_WFALLTHROUGH)
    list(APPEND torsion_cflags -Wno-implicit-fallthrough)
  endif()
  if(TORSION_HAS_WDECL)
    list(APPEND torsion_cflags -Wno-declaration-after-statement)
  endif()
  if(TORSION_HAS_WOVERLENGTH)
    list(APPEND torsion_cflags -Wno-overlength-strings)
  endif()
  if(TORSION_HAS_WLONGLONG)
    list(APPEND torsion_cflags -Wno-long-long)
  endif()
  if(TORSION_HAS_O3)
    list(APPEND torsion_cflags -O3)
  endif()
endif()

set(torsion_sources
    src/aead.c
    src/asn1.c
    src/chacha20.c
    src/cipher.c
    src/ecc.c
    src/encoding.c
    src/drbg.c
    src/dsa.c
    src/hash.c
    src/internal.c
    src/kdf.c
    src/mpi.c
    src/poly1305.c
    src/rsa.c
    src/salsa20.c
    src/secretbox.c
    src/siphash.c
    src/util.c)

if(TORSION_ENABLE_RNG)
  list(APPEND torsion_sources src/entropy/env.c)
  list(APPEND torsion_sources src/entropy/hw.c)
  list(APPEND torsion_sources src/entropy/sys.c)
  list(APPEND torsion_sources src/rand.c)
endif()

add_library(torsion SHARED ${torsion_sources})
target_compile_definitions(torsion PRIVATE ${torsion_defines})
target_compile_options(torsion PRIVATE ${torsion_cflags})
target_include_directories(torsion PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_library(torsion_a STATIC ${torsion_sources})
target_compile_definitions(torsion_a PRIVATE ${torsion_defines})
target_compile_options(torsion_a PRIVATE ${torsion_cflags})
target_include_directories(torsion_a PUBLIC ${PROJECT_SOURCE_DIR}/include)

if(TORSION_ENABLE_TESTS)
  add_executable(tests test/test.c test/hrtime.c)
  target_compile_definitions(tests PRIVATE ${torsion_defines})
  target_compile_options(tests PRIVATE ${torsion_cflags})
  target_link_libraries(tests torsion)
endif()

if(TORSION_ENABLE_CTIME)
  add_executable(ctime-test test/ctime-test.c)
  target_compile_definitions(ctime-test PRIVATE ${torsion_defines})
  target_compile_options(ctime-test PRIVATE ${torsion_cflags})
  target_link_libraries(ctime-test torsion)
endif()

if(UNIX)
  set(LIBS "")
  set(PACKAGE_VERSION "0.0.0")
  set_target_properties(torsion PROPERTIES VERSION 0.0.0 SOVERSION 0)
  set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
  set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
  set(prefix ${CMAKE_INSTALL_PREFIX})
  configure_file(libtorsion.pc.in libtorsion.pc @ONLY)

  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DOCDIR})
  install(FILES ${PROJECT_BINARY_DIR}/libtorsion.pc
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
  install(TARGETS torsion LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(TARGETS torsion_a ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(WIN32)
  install(DIRECTORY include/ DESTINATION include)
  install(FILES LICENSE DESTINATION .)
  install(TARGETS torsion torsion_a
          RUNTIME DESTINATION lib/$<CONFIG>
          ARCHIVE DESTINATION lib/$<CONFIG>)
endif()

message(STATUS "Build Options:

   CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}
   CMAKE_SYSTEM_NAME:    ${CMAKE_SYSTEM_NAME}
   CMAKE_C_COMPILER:     ${CMAKE_C_COMPILER}
   CMAKE_C_FLAGS:        ${CMAKE_C_FLAGS_${_build_type}} ${CMAKE_C_FLAGS}
   torsion_cflags:       ${torsion_cflags}
   torsion_defines:      ${torsion_defines}
")
