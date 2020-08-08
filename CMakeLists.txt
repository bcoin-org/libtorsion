# CMakeLists.txt - cmake build for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

cmake_minimum_required(VERSION 3.4)
project(libtorsion LANGUAGES C)

set(TORSION_PKG_VERSION 0.0.0)
set(TORSION_LIB_VERSION 0:0:0)

set(CPACK_PACKAGE_VERSION_MAJOR 0)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_VENDOR bcoin-org)
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE)
set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/README.md)
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${TORSION_PKG_VERSION}")

include(cmake/AppendCCompilerFlag.cmake)
include(cmake/TargetLinkOptions.cmake)
include(CheckIncludeFile)
include(CheckLibraryExists)
include(CTest)
include(GNUInstallDirs)

set(TORSION_WASM 0)

if(EMSCRIPTEN OR WASI_SDK_PREFIX)
  if(CMAKE_C_FLAGS STREQUAL "-v")
    set(CMAKE_C_FLAGS "")
  endif()
  set(TORSION_WASM 1)
endif()

if(NOT TORSION_WASM)
  include(CPack)
endif()

set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 90)

if(NOT TORSION_WASM)
  set(CMAKE_C_VISIBILITY_PRESET hidden)
  set(CMAKE_OSX_DEPLOYMENT_TARGET 10.7)
endif()

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

option(TORSION_ENABLE_ASM "Use inline x86-64 assembly if available" ON)
option(TORSION_ENABLE_COVERAGE "Enable coverage" OFF)
option(TORSION_ENABLE_CRT "Enable chinese remainder theorem for RSA" OFF)
option(TORSION_ENABLE_DEBUG "Enable debug build" OFF)
option(TORSION_ENABLE_INT128 "Use __int128 if available" ON)
option(TORSION_ENABLE_LIBSECP256K1 "Use libsecp256k1 field element backend" OFF)
option(TORSION_ENABLE_PTHREAD "Use pthread as a fallback for TLS" ON)
option(TORSION_ENABLE_RNG "Enable RNG" ON)
option(TORSION_ENABLE_TLS "Enable TLS" ON)
option(TORSION_ENABLE_VERIFY "Enable scalar bounds checks" OFF)

if(TORSION_WASM)
  set(TORSION_INITIAL_MEMORY "16777216" CACHE STRING "WASM initial memory")
  set(TORSION_MAX_MEMORY "2147483648" CACHE STRING "WASM maximum memory")
  set(TORSION_STACK_SIZE "5242880" CACHE STRING "WASM stack size")
endif()

if(EMSCRIPTEN)
  set(TORSION_ENVIRONMENT "node" CACHE STRING "Emscripten environment")
endif()

set(torsion_cflags)
set(torsion_ldflags)

if(MSVC)
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /DTORSION_DEBUG")

  append_c_compiler_flag(torsion_cflags /wd4146
                                        /wd4244
                                        /wd4267
                                        /wd4334)
else()
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DTORSION_DEBUG")

  append_c_compiler_flag(torsion_cflags -pedantic
                                        -Wall
                                        -Wextra
                                        -Wcast-align
                                        -Wno-implicit-fallthrough
                                        -Wno-long-long
                                        -Wno-overlength-strings
                                        -Wshadow)
endif()

check_c_source_compiles([=[
int main(void) {
  unsigned char ptr[32];
  int x = 0;
  __asm__ ("" : "+r" (x) ::);
  __asm__ __volatile__("" :: "r" (ptr) : "memory");
  return x;
}
]=] TORSION_HAS_ASM)

check_c_source_compiles([=[
#if !defined(__amd64__) && !defined(__x86_64__)
#  error "not an x64 platform"
#endif
#include <stdint.h>
int main(void) {
  uint32_t stream[8], state[8];
  __asm__ __volatile__(
    "movups (%%rsi), %%xmm0\\n"
    "movups 16(%%rsi), %%xmm1\\n"
    "pxor %%xmm1, %%xmm0\\n"
    "movups %%xmm0,(%%rdi)\\n"
    "movups %%xmm1,16(%%rdi)\\n"
    :
    : "D" (stream), "S" (state)
    : "xmm0", "xmm1", "cc", "memory"
  );
  return 0;
}
]=] TORSION_HAS_ASM_X64)

check_c_source_compiles([=[
#ifndef __wasm__
typedef char check_64bit_t[sizeof(void *) == 8 ? 1 : -1];
#endif
typedef signed __int128 int128_t;
typedef unsigned __int128 uint128_t;
int main(void) {
  uint128_t r = 1;
  r <<= 64;
  r *= 113;
  r >>= 65;
  return r & 1;
}
]=] TORSION_HAS_INT128)

check_c_source_compiles([=[
#if defined(_WIN32) && !defined(__MINGW32__)
#  define TLS __declspec(thread)
#else
#  define TLS __thread
#endif
static TLS int value;
int main(void) {
  value = 1;
  return 0;
}
]=] TORSION_HAS_TLS)

set(TORSION_HAS_THREADS OFF)

if(NOT TORSION_WASM)
  if(WIN32)
    set(TORSION_HAS_THREADS ON)
  else()
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads)
    if(CMAKE_USE_PTHREADS_INIT)
      set(TORSION_HAS_THREADS ON)
    endif()
  endif()
endif()

check_c_source_compiles([=[
#if defined(__EMSCRIPTEN__) || defined(__wasm__)
#  error "cannot use fork with wasm"
#endif
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
int main(void) {
  unsigned char data[32];
  int pfds[2];
  int status;
  pid_t pid;
  pipe(pfds);
  pid = fork();
  close(pfds[1]);
  read(pfds[0], data, 32);
  close(pfds[0]);
  waitpid(pid, &status, 0);
  WIFEXITED(status);
  WEXITSTATUS(status);
  return 0;
}
]=] TORSION_HAS_FORK)

set(TORSION_HAS_VALGRIND)
set(TORSION_HAS_ZLIB)

if(NOT TORSION_WASM AND NOT WIN32)
  check_include_file(valgrind/memcheck.h TORSION_HAS_VALGRIND)
  if(TORSION_ENABLE_RNG)
    check_include_file(zlib.h TORSION_HAS_ZLIB_H)
    if (TORSION_HAS_ZLIB_H)
      check_library_exists(z compress2 "" TORSION_HAS_ZLIB)
    endif()
  endif()
endif()

set(torsion_defines)

list(APPEND torsion_defines TORSION_HAVE_CONFIG)

if(WIN32 AND NOT TORSION_WASM)
  list(APPEND torsion_defines _WIN32_WINNT=0x0600)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND NOT TORSION_WASM)
  list(APPEND torsion_defines _POSIX_C_SOURCE=200112)
endif()

if(TORSION_ENABLE_ASM AND TORSION_HAS_ASM)
  list(APPEND torsion_defines TORSION_HAVE_ASM)
endif()

if(TORSION_ENABLE_ASM AND TORSION_HAS_ASM_X64)
  list(APPEND torsion_defines TORSION_HAVE_ASM_X64)
endif()

if(TORSION_ENABLE_COVERAGE)
  if(TORSION_WASM)
    message(FATAL_ERROR "coverage not supported with wasm")
  endif()
  list(APPEND torsion_cflags -O0 --coverage)
  list(APPEND torsion_ldflags --coverage)
  list(APPEND torsion_defines TORSION_COVERAGE)
endif()

if(TORSION_ENABLE_CRT)
  list(APPEND torsion_defines TORSION_USE_CRT)
endif()

if(TORSION_ENABLE_DEBUG)
  if(MSVC)
    append_c_compiler_flag(torsion_cflags /Zi)
  else()
    append_c_compiler_flag(torsion_cflags -g)
  endif()
  list(APPEND torsion_defines TORSION_DEBUG)
endif()

if(TORSION_HAS_FORK)
  list(APPEND torsion_defines TORSION_HAVE_FORK)
endif()

if(TORSION_ENABLE_INT128 AND TORSION_HAS_INT128)
  list(APPEND torsion_defines TORSION_HAVE_INT128)
endif()

if(TORSION_ENABLE_LIBSECP256K1)
  list(APPEND torsion_defines TORSION_USE_LIBSECP256K1)
endif()

if(TORSION_ENABLE_RNG)
  list(APPEND torsion_defines TORSION_HAVE_RNG)
endif()

if(TORSION_HAS_THREADS)
  list(APPEND torsion_defines TORSION_HAVE_THREADS)
endif()

if(TORSION_ENABLE_TLS AND TORSION_HAS_TLS)
  list(APPEND torsion_defines TORSION_HAVE_TLS)
endif()

if(TORSION_ENABLE_VERIFY)
  list(APPEND torsion_defines TORSION_VERIFY)
endif()

if(TORSION_HAS_ZLIB)
  list(APPEND torsion_defines TORSION_HAVE_ZLIB)
endif()

set(torsion_sources src/aead.c
                    src/asn1.c
                    src/cipher.c
                    src/ecc.c
                    src/encoding.c
                    src/drbg.c
                    src/dsa.c
                    src/hash.c
                    src/ies.c
                    src/internal.c
                    src/kdf.c
                    src/mac.c
                    src/mpi.c
                    src/rsa.c
                    src/stream.c
                    src/util.c)

if(TORSION_ENABLE_RNG)
  list(APPEND torsion_sources src/entropy/env.c
                              src/entropy/hw.c
                              src/entropy/sys.c
                              src/rand.c)
endif()

set(torsion_includes ${PROJECT_SOURCE_DIR}/include)

if(TORSION_WASM)
  if(EMSCRIPTEN)
    # CMAKE_CROSSCOMPILING_EMULATOR is mistakenly quoted by emcmake.
    string(REPLACE "\"" "" CMAKE_CROSSCOMPILING_EMULATOR
           "${CMAKE_CROSSCOMPILING_EMULATOR}")

    set(torsion_ldflags "SHELL:-s SINGLE_FILE=1"
                        "SHELL:-s ASSERTIONS=0"
                        "SHELL:-s NODEJS_CATCH_EXIT=0"
                        "SHELL:-s NODEJS_CATCH_REJECTION=0"
                        "SHELL:-s ALLOW_MEMORY_GROWTH=1"
                        "SHELL:-s INITIAL_MEMORY=${TORSION_INITIAL_MEMORY}"
                        "SHELL:-s MAXIMUM_MEMORY=${TORSION_MAX_MEMORY}"
                        "SHELL:-s TOTAL_STACK=${TORSION_STACK_SIZE}"
                        "SHELL:-s ENVIRONMENT=${TORSION_ENVIRONMENT}")

    set(torsion_sources_lib "")
    set(torsion_ldflags_lib "SHELL:-s EXPORTED_FUNCTIONS=@etc/exports.json")
  else()
    set(CMAKE_EXECUTABLE_SUFFIX ".wasm")

    set(torsion_ldflags -Wl,--allow-undefined
                        -Wl,--initial-memory=${TORSION_INITIAL_MEMORY}
                        -Wl,--max-memory=${TORSION_MAX_MEMORY}
                        -Wl,-z -Wl,stack-size=${TORSION_STACK_SIZE}
                        -Wl,--stack-first)

    # Note: We have to manually choose crt1-reactor.o
    # until reactor support is widely available via:
    #
    #   -mexec-model=reactor
    #
    # See: https://github.com/llvm/llvm-project/commit/d496437
    #      https://reviews.llvm.org/D62922
    set(torsion_syslib ${CMAKE_SYSROOT}/lib/${CMAKE_C_COMPILER_TARGET})
    set(torsion_sources_lib ${torsion_syslib}/crt1-reactor.o)
    set(torsion_ldflags_lib -nostartfiles
                            -Wl,--entry=_initialize
                            -Wl,--export-dynamic
                            -Wl,--export=malloc
                            -Wl,--export=free)
  endif()

  add_library(torsion_o OBJECT ${torsion_sources})
  target_compile_definitions(torsion_o PUBLIC ${torsion_defines}
                                       PRIVATE TORSION_BUILD)
  target_compile_options(torsion_o PUBLIC ${torsion_cflags})
  target_include_directories(torsion_o PUBLIC ${torsion_includes})
  target_link_options(torsion_o INTERFACE ${torsion_ldflags})

  add_executable(torsion ${torsion_sources_lib})
  target_link_options(torsion PRIVATE ${torsion_ldflags_lib})
  target_link_libraries(torsion PRIVATE torsion_o)

  add_executable(torsion_bench test/bench.c test/hrtime.c)
  target_link_libraries(torsion_bench PRIVATE torsion_o)

  add_executable(torsion_test test/test.c)
  target_link_libraries(torsion_test PRIVATE torsion_o)
  add_test(NAME test_wasm COMMAND torsion_test)
else()
  set(torsion_libs)
  set(pkgconf_libs)
  set(test_sources test/test.c)
  set(test_libs)

  if(TORSION_ENABLE_PTHREAD AND TORSION_HAS_THREADS AND NOT WIN32)
    if(NOT TORSION_ENABLE_TLS OR NOT TORSION_HAS_TLS)
      list(APPEND torsion_defines TORSION_HAVE_PTHREAD)
      list(APPEND torsion_libs Threads::Threads)
      list(APPEND pkgconf_libs ${CMAKE_THREAD_LIBS_INIT})
    endif()
  endif()

  if(TORSION_HAS_THREADS)
    list(APPEND test_sources test/thread.c)
    if(NOT WIN32)
      list(APPEND test_libs Threads::Threads)
    endif()
  endif()

  if(TORSION_HAS_ZLIB)
    list(APPEND test_libs z)
  endif()

  add_library(torsion_o OBJECT ${torsion_sources})
  target_compile_definitions(torsion_o PRIVATE ${torsion_defines} TORSION_BUILD)
  target_compile_options(torsion_o PRIVATE ${torsion_cflags})
  target_include_directories(torsion_o PRIVATE ${torsion_includes})
  set_property(TARGET torsion_o PROPERTY POSITION_INDEPENDENT_CODE ON)

  add_library(torsion SHARED $<TARGET_OBJECTS:torsion_o>)
  target_compile_definitions(torsion INTERFACE ${torsion_defines})
  target_compile_options(torsion INTERFACE ${torsion_cflags})
  target_include_directories(torsion INTERFACE ${torsion_includes})
  target_link_options(torsion PUBLIC ${torsion_ldflags})
  target_link_libraries(torsion PRIVATE ${torsion_libs})

  add_library(torsion_static STATIC $<TARGET_OBJECTS:torsion_o>)
  target_compile_definitions(torsion_static INTERFACE ${torsion_defines})
  target_compile_options(torsion_static INTERFACE ${torsion_cflags})
  target_include_directories(torsion_static INTERFACE ${torsion_includes})
  target_link_options(torsion_static INTERFACE ${torsion_ldflags})
  target_link_libraries(torsion_static INTERFACE ${torsion_libs})

  if(UNIX)
    string(REPLACE ":" "." TORSION_VERSION "${TORSION_LIB_VERSION}")
    string(REGEX MATCH "^[0-9]+" TORSION_SOVERSION "${TORSION_LIB_VERSION}")
    set_target_properties(torsion PROPERTIES VERSION ${TORSION_VERSION}
                                             SOVERSION ${TORSION_SOVERSION})
    set_target_properties(torsion_static PROPERTIES OUTPUT_NAME torsion)
  endif()

  add_executable(torsion_bench test/bench.c test/hrtime.c)
  target_link_libraries(torsion_bench PRIVATE torsion)

  add_executable(torsion_test ${test_sources})
  target_link_libraries(torsion_test PRIVATE torsion ${test_libs})
  add_test(NAME test_shared COMMAND torsion_test)

  add_executable(torsion_test_static ${test_sources})
  target_link_libraries(torsion_test_static PRIVATE torsion_static ${test_libs})
  add_test(NAME test_static COMMAND torsion_test_static)

  if(TORSION_HAS_VALGRIND)
    add_executable(torsion_ctgrind test/ctgrind.c)
    target_link_libraries(torsion_ctgrind PRIVATE torsion)
  endif()

  if(UNIX)
    set(prefix ${CMAKE_INSTALL_PREFIX})
    set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
    set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
    set(PACKAGE_VERSION "${TORSION_PKG_VERSION}")
    string(REPLACE ";" " " LIBS "${pkgconf_libs}")
    configure_file(libtorsion.pc.in libtorsion.pc @ONLY)
  endif()

  if(UNIX)
    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES LICENSE
            DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME})
    install(FILES README.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
    install(FILES ${PROJECT_BINARY_DIR}/libtorsion.pc
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
    install(TARGETS torsion LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(TARGETS torsion_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
  endif()

  if(WIN32)
    install(DIRECTORY include/ DESTINATION include)
    install(FILES LICENSE DESTINATION .)
    install(FILES README.md DESTINATION .)
    install(TARGETS torsion torsion_static
            RUNTIME DESTINATION lib/$<CONFIG>
            ARCHIVE DESTINATION lib/$<CONFIG>)
  endif()
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type)

message(STATUS "Build Options:

   TORSION_PKG_VERSION:  ${TORSION_PKG_VERSION}
   TORSION_LIB_VERSION:  ${TORSION_LIB_VERSION}
   CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}
   CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}
   CMAKE_SYSTEM_NAME:    ${CMAKE_SYSTEM_NAME}
   CMAKE_C_COMPILER:     ${CMAKE_C_COMPILER}
   CMAKE_C_FLAGS:        ${CMAKE_C_FLAGS_${_build_type}} ${CMAKE_C_FLAGS}
   torsion_cflags:       ${torsion_cflags}
   torsion_defines:      ${torsion_defines}
   torsion_libs:         ${torsion_libs}
")
