# CMakeLists.txt - cmake build for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

#
# Initialization
#

cmake_minimum_required(VERSION 3.4)
project(libtorsion VERSION 0.0.0 LANGUAGES C)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(WASI OR EMSCRIPTEN)
  set(TORSION_WASM 1)
else()
  set(TORSION_WASM 0)
endif()

#
# Toolchain/Platform Fixes
#

if(EMSCRIPTEN)
  # CMAKE_CROSSCOMPILING_EMULATOR is mistakenly quoted by emcmake.
  string(REPLACE "\"" "" CMAKE_CROSSCOMPILING_EMULATOR
         "${CMAKE_CROSSCOMPILING_EMULATOR}")
endif()

if(WASI)
  # No idea why this isn't set.
  set(CMAKE_EXECUTABLE_SUFFIX ".wasm")
endif()

#
# Includes
#

include(cmake/AppendCCompilerFlag.cmake)
include(cmake/CheckCThreadLocalStorage.cmake)
include(cmake/TargetLinkOptions.cmake)
include(CheckCSourceCompiles)
include(CheckIncludeFile)
include(CheckLibraryExists)
include(CTest)
include(GNUInstallDirs)

#
# Options
#

option(TORSION_ENABLE_ASM "Use inline assembly if available" ON)
option(TORSION_ENABLE_COVERAGE "Enable coverage" OFF)
option(TORSION_ENABLE_DEBUG "Enable debug build" OFF)
option(TORSION_ENABLE_INT128 "Use __int128 if available" ON)
option(TORSION_ENABLE_RNG "Enable RNG" ON)
option(TORSION_ENABLE_TLS "Enable thread-local storage" ON)
option(TORSION_ENABLE_VERIFY "Enable scalar bounds checks" OFF)

if(TORSION_WASM)
  set(TORSION_INITIAL_MEMORY "16777216" CACHE STRING "WASM initial memory")
  set(TORSION_MAX_MEMORY "2147483648" CACHE STRING "WASM maximum memory")
  set(TORSION_STACK_SIZE "5242880" CACHE STRING "WASM stack size")
endif()

if(EMSCRIPTEN)
  set(TORSION_ENVIRONMENT "node" CACHE STRING "Emscripten environment")
endif()

#
# Global Flags
#

if(NOT PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME)
  unset(CMAKE_C_STANDARD_REQUIRED)
  unset(CMAKE_C_EXTENSIONS)
  unset(CMAKE_C_STANDARD)
  unset(CMAKE_C_VISIBILITY_PRESET)
  unset(CMAKE_OSX_DEPLOYMENT_TARGET)
endif()

if(NOT TORSION_WASM)
  set(CMAKE_C_VISIBILITY_PRESET hidden)
endif()

if(APPLE)
  set(CMAKE_OSX_DEPLOYMENT_TARGET 10.7)
endif()

#
# Feature Testing
#

if(TORSION_ENABLE_ASM)
  check_c_source_compiles([=[
#   if defined(__native_client__) || defined(__IBMC__) || defined(__IBMCPP__)
#     error "inline asm not supported"
#   endif
    int main(void) {
      unsigned long z = 953;
      unsigned long x = 109;
      unsigned long y = 577;
      unsigned long c;
      __asm__ __volatile__ (
        ""
        : "+r" (z), "=&r" (c)
#   if defined(__TINYC__)
        : "rm" (x), "rm" (y)
#   else
        : "%rm" (x), "rm" (y)
#   endif
        : "cc", "memory"
      );
      return z & 0x7f;
    }
  ]=] TORSION_HAS_ASM)
else()
  set(TORSION_HAS_ASM)
endif()

if(BUILD_TESTING AND NOT TORSION_WASM)
  check_c_source_compiles([=[
#   include <sys/types.h>
#   include <sys/wait.h>
#   include <unistd.h>
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
      return status;
    }
  ]=] TORSION_HAS_FORK)
else()
  set(TORSION_HAS_FORK)
endif()

if(TORSION_ENABLE_INT128)
  check_c_source_compiles([=[
    typedef signed __int128 xint128_t;
    typedef unsigned __int128 xuint128_t;
    typedef char check_voidptr_t[sizeof(void *) >= 8 ? 1 : -1];
    typedef char check_int128_t[sizeof(xint128_t) == 16 ? 1 : -1];
    typedef char check_uint128_t[sizeof(xuint128_t) == 16 ? 1 : -1];
    int main(int argc, char **argv) {
      xint128_t c = argv[0][0];
      xuint128_t r = argc + c;
      while (argc--) r *= r;
      return r >> 121;
    }
  ]=] TORSION_HAS_INT128)
else()
  set(TORSION_HAS_INT128)
endif()

if(TORSION_ENABLE_TLS)
  check_c_thread_local_storage(TORSION_TLS)
else()
  set(TORSION_TLS)
endif()

set(TORSION_HAS_MEMCHECK 0)
set(TORSION_HAS_ZLIB 0)

if(BUILD_TESTING AND NOT TORSION_WASM AND NOT WIN32)
  check_include_file(valgrind/memcheck.h TORSION_HAS_MEMCHECK_H)

  if(TORSION_HAS_MEMCHECK_H)
    set(TORSION_HAS_MEMCHECK 1)
  endif()

  if (TORSION_ENABLE_RNG)
    check_include_file(zlib.h TORSION_HAS_ZLIB_H)
    check_library_exists(z compress2 "" TORSION_HAS_ZLIB_A)

    if(TORSION_HAS_ZLIB_H AND TORSION_HAS_ZLIB_A)
      set(TORSION_HAS_ZLIB 1)
    endif()
  endif()
endif()

#
# Flags
#

set(torsion_cflags)
set(torsion_ldflags)

if(MSVC)
  append_c_compiler_flag(torsion_cflags /wd4146
                                        /wd4244
                                        /wd4267
                                        /wd4334)
else()
  append_c_compiler_flag(torsion_cflags -pedantic
                                        -Wall
                                        -Wextra
                                        -Wcast-align
                                        -Wno-implicit-fallthrough
                                        -Wno-overlength-strings
                                        -Wshadow)
endif()

if(TORSION_ENABLE_COVERAGE)
  list(APPEND torsion_cflags -O0 --coverage)
  list(APPEND torsion_ldflags --coverage)
endif()

if(TORSION_ENABLE_DEBUG)
  if(MSVC)
    append_c_compiler_flag(torsion_cflags /Zi)
  else()
    append_c_compiler_flag(torsion_cflags -g)
  endif()
endif()

if(MINGW)
  # No support for #pragma comment(lib).
  list(APPEND torsion_cflags -Wno-unknown-pragmas)
  # Ensure we are redistributable on windows.
  list(APPEND torsion_ldflags -static-libgcc)
endif()

#
# Defines
#

set(torsion_defines)

list(APPEND torsion_defines TORSION_HAVE_CONFIG)

if(TORSION_HAS_ASM)
  list(APPEND torsion_defines TORSION_HAVE_ASM)
endif()

if(TORSION_ENABLE_COVERAGE)
  list(APPEND torsion_defines TORSION_COVERAGE)
endif()

if(TORSION_ENABLE_DEBUG)
  list(APPEND torsion_defines TORSION_DEBUG)
else()
  list(APPEND torsion_defines $<$<CONFIG:Debug>:TORSION_DEBUG>)
endif()

if(TORSION_HAS_FORK)
  list(APPEND torsion_defines TORSION_HAVE_FORK)
endif()

if(TORSION_HAS_INT128)
  list(APPEND torsion_defines TORSION_HAVE_INT128)
endif()

if(TORSION_ENABLE_RNG)
  list(APPEND torsion_defines TORSION_HAVE_RNG)
endif()

if(TORSION_TLS)
  list(APPEND torsion_defines TORSION_TLS=${TORSION_TLS})
endif()

if(TORSION_ENABLE_VERIFY)
  list(APPEND torsion_defines TORSION_VERIFY)
endif()

if(TORSION_HAS_ZLIB)
  list(APPEND torsion_defines TORSION_HAVE_ZLIB)
endif()

#
# Feature Test Macros
#

if(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
  list(APPEND torsion_defines _TS_ERRNO)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "AIX")
  list(APPEND torsion_defines _THREAD_SAFE_ERRNO)
endif()

#
# Libraries
#

set(torsion_libs)
set(pkgconf_libs)
set(test_libs)

if(TORSION_TLS)
  list(APPEND torsion_libs Threads::TLS)
  get_property(pkgconf_libs TARGET Threads::TLS
                            PROPERTY INTERFACE_LINK_LIBRARIES)
endif()

if(TORSION_HAS_ZLIB)
  list(APPEND test_libs z)
endif()

if(MINGW)
  # No support for #pragma comment(lib).
  list(APPEND torsion_libs -liphlpapi -lpsapi -lbcrypt)
endif()

#
# Includes
#

set(torsion_includes ${PROJECT_SOURCE_DIR}/include)

#
# Sources
#

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

set(rng_sources src/entropy/env.c
                src/entropy/hrt.c
                src/entropy/hw.c
                src/entropy/sys.c
                src/rand.c)

set(test_sources test/test.c test/utils.c)
set(bench_sources test/bench.c test/hrtime.c test/utils.c)
set(ctgrind_sources test/ctgrind.c test/utils.c)

if(TORSION_ENABLE_RNG)
  list(APPEND torsion_sources ${rng_sources})
endif()

#
# Targets
#

if(NOT PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME)
  # Bulding a subproject. Keep it simple.
  if(TORSION_WASM)
    add_library(torsion OBJECT ${torsion_sources})
    target_compile_definitions(torsion PRIVATE ${torsion_defines})
    target_compile_options(torsion PRIVATE ${torsion_cflags})
    target_include_directories(torsion PUBLIC ${torsion_includes})
  else()
    if(COMMAND add_node_library)
      add_node_library(torsion STATIC ${torsion_sources})
    else()
      add_library(torsion STATIC ${torsion_sources})
      set_property(TARGET torsion PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif()

    target_compile_definitions(torsion PRIVATE ${torsion_defines})
    target_compile_options(torsion PRIVATE ${torsion_cflags})
    target_include_directories(torsion PUBLIC ${torsion_includes})
    target_link_options(torsion INTERFACE ${torsion_ldflags})
    target_link_libraries(torsion INTERFACE ${torsion_libs})
  endif()
elseif(TORSION_WASM)
  if(EMSCRIPTEN)
    set(torsion_ldflags "SHELL:-s SINGLE_FILE=1"
                        "SHELL:-s ASSERTIONS=0"
                        "SHELL:-s NODEJS_CATCH_EXIT=0"
                        "SHELL:-s NODEJS_CATCH_REJECTION=0"
                        "SHELL:-s ALLOW_MEMORY_GROWTH=1"
                        "SHELL:-s INITIAL_MEMORY=${TORSION_INITIAL_MEMORY}"
                        "SHELL:-s MAXIMUM_MEMORY=${TORSION_MAX_MEMORY}"
                        "SHELL:-s TOTAL_STACK=${TORSION_STACK_SIZE}"
                        "SHELL:-s ENVIRONMENT=${TORSION_ENVIRONMENT}")

    set(torsion_ldflags_lib "SHELL:-s EXPORTED_FUNCTIONS=@etc/exports.json")
  else()
    set(torsion_ldflags -Wl,--allow-undefined
                        -Wl,--initial-memory=${TORSION_INITIAL_MEMORY}
                        -Wl,--max-memory=${TORSION_MAX_MEMORY}
                        -Wl,-z,stack-size=${TORSION_STACK_SIZE}
                        -Wl,--stack-first)

    set(torsion_ldflags_lib -mexec-model=reactor
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

  add_executable(torsion "")
  target_link_options(torsion PRIVATE ${torsion_ldflags_lib})
  target_link_libraries(torsion PRIVATE torsion_o)

  if(BUILD_TESTING)
    add_executable(torsion_bench ${bench_sources})
    target_link_libraries(torsion_bench PRIVATE torsion_o)

    add_executable(torsion_test ${test_sources})
    target_link_libraries(torsion_test PRIVATE torsion_o)
    add_test(NAME test_wasm COMMAND torsion_test)

    if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
      add_custom_target(bench COMMAND torsion_bench)
      add_custom_target(check COMMAND ${CMAKE_MAKE_PROGRAM} test)
    endif()
  endif()
else()
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

  if(MINGW)
    # Useful for creating more-compatible import libraries on windows:
    #
    #   lib.exe /OUT:[name].lib /DEF:[name].def /NAME:[name].dll /MACHINE:X64
    #
    # Or:
    #
    #   llvm-dlltool -l [name].lib -d [name].def -D [name].dll -m i386:x86-64
    #
    # More-compatible static libraries can already be created with:
    #
    #   llvm-lib /out:[name].lib [name].a
    target_link_options(torsion PRIVATE -Wl,--output-def,libtorsion.def)
  endif()

  if(UNIX)
    set_target_properties(torsion PROPERTIES VERSION 0.0.0 SOVERSION 0)
    set_target_properties(torsion_static PROPERTIES OUTPUT_NAME torsion)

    configure_file(libtorsion.pc.in libtorsion.pc @ONLY)

    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES LICENSE
            DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME})
    install(FILES README.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
    install(FILES ${PROJECT_BINARY_DIR}/libtorsion.pc
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
    install(TARGETS torsion LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(TARGETS torsion_static ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()

  if(WIN32)
    install(DIRECTORY include/ DESTINATION include)
    install(FILES LICENSE README.md DESTINATION .)
    install(TARGETS torsion torsion_static
            RUNTIME DESTINATION lib/$<CONFIG>
            ARCHIVE DESTINATION lib/$<CONFIG>)
  endif()

  if(BUILD_TESTING)
    add_executable(torsion_bench ${bench_sources})
    target_link_libraries(torsion_bench PRIVATE torsion)

    add_executable(torsion_test ${test_sources})
    target_link_libraries(torsion_test PRIVATE torsion ${test_libs})
    add_test(NAME test_shared COMMAND torsion_test)

    add_executable(torsion_test_static ${test_sources})
    target_link_libraries(torsion_test_static PRIVATE torsion_static
                                                      ${test_libs})
    add_test(NAME test_static COMMAND torsion_test_static)

    if(TORSION_HAS_MEMCHECK)
      add_executable(torsion_ctgrind ${ctgrind_sources})
      target_link_libraries(torsion_ctgrind PRIVATE torsion)
    endif()

    # Not sure why cmake fails to do this on tinycc.
    if(CMAKE_C_COMPILER_ID STREQUAL "TinyCC")
      set(rpath -Wl,-rpath,${PROJECT_SOURCE_DIR})

      target_link_options(torsion_test PRIVATE ${rpath})
      target_link_options(torsion_bench PRIVATE ${rpath})

      if(TORSION_HAS_MEMCHECK)
        target_link_options(torsion_ctgrind PRIVATE ${rpath})
      endif()
    endif()

    find_program(TORSION_VALGRIND valgrind)

    if(TORSION_VALGRIND)
      if(TORSION_HAS_MEMCHECK AND NOT TORSION_ENABLE_DEBUG)
        add_test(NAME test_ctgrind COMMAND ${TORSION_VALGRIND}
                                           --leak-check=full
                                           --error-exitcode=1
                                           $<TARGET_FILE:torsion_ctgrind>)
      endif()

      add_test(NAME test_valgrind COMMAND ${TORSION_VALGRIND}
                                          --leak-check=full
                                          --error-exitcode=1
                                          $<TARGET_FILE:torsion_test>)
    endif()

    if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
      add_custom_target(bench COMMAND torsion_bench)
      add_custom_target(check COMMAND ${CMAKE_MAKE_PROGRAM} test)

      if(TORSION_VALGRIND)
        if(TORSION_HAS_MEMCHECK AND NOT TORSION_ENABLE_DEBUG)
          add_custom_target(ctgrind COMMAND ${TORSION_VALGRIND}
                                            --leak-check=full
                                            --error-exitcode=1
                                            $<TARGET_FILE:torsion_ctgrind>
                                    DEPENDS torsion_ctgrind)
        endif()

        add_custom_target(valgrind COMMAND ${TORSION_VALGRIND}
                                           --leak-check=full
                                           --error-exitcode=1
                                           $<TARGET_FILE:torsion_test>
                                   DEPENDS torsion_test)
      endif()
    endif()
  endif()
endif()

#
# Output
#

string(TOUPPER "${CMAKE_BUILD_TYPE}" torsion_build_type)

message(STATUS "Build Options:

   CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}
   CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}
   CMAKE_SYSTEM_NAME:    ${CMAKE_SYSTEM_NAME}
   CMAKE_C_COMPILER:     ${CMAKE_C_COMPILER}
   CMAKE_C_FLAGS:        ${CMAKE_C_FLAGS_${torsion_build_type}} ${CMAKE_C_FLAGS}
   torsion_cflags:       ${torsion_cflags}
   torsion_ldflags:      ${torsion_ldflags}
   torsion_defines:      ${torsion_defines}
   torsion_libs:         ${torsion_libs}
")
