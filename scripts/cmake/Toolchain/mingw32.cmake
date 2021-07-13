# mingw32.cmake - mingw32 toolchain for cmake
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_CROSSCOMPILING_EMULATOR wine64)

set(MINGW32_ARCH "x86_64" CACHE STRING "mingw32 arch")

# https://github.com/ruslo/polly/blob/d5fb153/compiler/gcc-cross-compile.cmake
set(CMAKE_AR               "${MINGW32_ARCH}-w64-mingw32-ar"       CACHE PATH "Archiver")
set(CMAKE_ASM_COMPILER     "${MINGW32_ARCH}-w64-mingw32-as"       CACHE PATH "Assembler")
set(CMAKE_C_COMPILER       "${MINGW32_ARCH}-w64-mingw32-cc"       CACHE PATH "C Compiler")
set(CMAKE_C_PREPROCESSOR   "${MINGW32_ARCH}-w64-mingw32-cpp"      CACHE PATH "C Preprocessor")
set(CMAKE_CXX_COMPILER     "${MINGW32_ARCH}-w64-mingw32-c++"      CACHE PATH "C++ Compiler")
set(CMAKE_DLLTOOL          "${MINGW32_ARCH}-w64-mingw32-dlltool"  CACHE PATH "dlltool")
set(CMAKE_Fortran_COMPILER "${MINGW32_ARCH}-w64-mingw32-gfortran" CACHE PATH "Fortran Compiler")
set(CMAKE_LINKER           "${MINGW32_ARCH}-w64-mingw32-ld"       CACHE PATH "Linker")
set(CMAKE_NM               "${MINGW32_ARCH}-w64-mingw32-nm"       CACHE PATH "nm")
set(CMAKE_OBJCOPY          "${MINGW32_ARCH}-w64-mingw32-objcopy"  CACHE PATH "objcopy")
set(CMAKE_OBJDUMP          "${MINGW32_ARCH}-w64-mingw32-objdump"  CACHE PATH "objdump")
set(CMAKE_RANLIB           "${MINGW32_ARCH}-w64-mingw32-ranlib"   CACHE PATH "ranlib")
set(CMAKE_RC_COMPILER      "${MINGW32_ARCH}-w64-mingw32-windres"  CACHE PATH "Resource Compiler")
set(CMAKE_STRIP            "${MINGW32_ARCH}-w64-mingw32-strip"    CACHE PATH "strip")

set(CMAKE_FIND_ROOT_PATH /usr/${MINGW32_ARCH}-w64-mingw32)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
