cmake_minimum_required(VERSION 3.4)
project(libtorsion LANGUAGES C)

include(GNUInstallDirs)
include(CheckCCompilerFlag)
include(CheckIncludeFile)

set(CMAKE_C_VISIBILITY_PRESET default)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 90)

option(TORSION_ENABLE_ASSERT "Enable assertions" ON)
option(TORSION_ENABLE_CRT "Enable chinese remainder theorem for RSA" OFF)
option(TORSION_ENABLE_CTIME "Enable valgrind constant-time test" OFF)
option(TORSION_ENABLE_LIBSECP256K1 "Enable libsecp256k1 field element backend" OFF)
option(TORSION_ENABLE_RNG "Enable RNG" OFF)
option(TORSION_ENABLE_TESTS "Enable tests" OFF)
option(TORSION_ENABLE_VERIFY "Enable scalar bounds checks" OFF)

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
check_include_file(valgrind/memcheck.h TORSION_HAS_VALGRIND)

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
  if(TORSION_ENABLE_ASSERT)
    message(FATAL_ERROR "cannot run constant-time test with assertions")
  endif()
  if(TORSION_ENABLE_VERIFY)
    message(FATAL_ERROR "cannot run constant-time test with bounds verification")
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
# todo: disable msvc warnings
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
  list(APPEND src/entropy/env.c)
  list(APPEND src/entropy/hw.c)
  list(APPEND src/entropy/sys.c)
  list(APPEND src/rand.c)
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
  set(torsion_test_sources test/test.c test/hrtime.c)

  add_executable(torsion_run_tests ${torsion_test_sources})
  target_compile_definitions(torsion_run_tests PRIVATE ${torsion_defines})
  target_compile_options(torsion_run_tests PRIVATE ${torsion_cflags})
  target_link_libraries(torsion_run_tests torsion)
  add_test(NAME torsion_test
           COMMAND torsion_run_tests
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(torsion_run_tests_a ${torsion_test_sources})
  target_compile_definitions(torsion_run_tests_a PRIVATE ${torsion_defines})
  target_compile_options(torsion_run_tests_a PRIVATE ${torsion_cflags})
  target_link_libraries(torsion_run_tests_a torsion_a)
  add_test(NAME torsion_test_a
           COMMAND torsion_run_tests_a
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(TORSION_ENABLE_CTIME)
  set(torsion_ctime_sources test/ctime-test.c)

  add_executable(torsion_run_ctime_test ${torsion_ctime_sources})
  target_compile_definitions(torsion_run_ctime_test PRIVATE ${torsion_defines})
  target_compile_options(torsion_run_ctime_test PRIVATE ${torsion_cflags})
  target_link_libraries(torsion_run_ctime_test torsion)
  add_test(NAME torsion_ctime_test
           COMMAND torsion_run_ctime_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(torsion_run_ctime_test_a ${torsion_ctime_sources})
  target_compile_definitions(torsion_run_ctime_test_a PRIVATE ${torsion_defines})
  target_compile_options(torsion_run_ctime_test_a PRIVATE ${torsion_cflags})
  target_link_libraries(torsion_run_ctime_test_a torsion_a)
  add_test(NAME torsion_ctime_test_a
           COMMAND torsion_run_ctime_test_a
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(UNIX)
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
