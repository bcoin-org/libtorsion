# mingw32.cmake - mingw32 toolchain for cmake
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_CROSSCOMPILING_EMULATOR wine64)

set(MINGW_ARCH "x86_64" CACHE STRING "mingw32 arch")

# https://github.com/ruslo/polly/blob/d5fb153/compiler/gcc-cross-compile.cmake
set(CMAKE_C_COMPILER       "${MINGW_ARCH}-w64-mingw32-cc"       CACHE PATH "C compiler")
set(CMAKE_CXX_COMPILER     "${MINGW_ARCH}-w64-mingw32-c++"      CACHE PATH "C++ compiler")
set(CMAKE_ASM_COMPILER     "${MINGW_ARCH}-w64-mingw32-as"       CACHE PATH "Assembler")
set(CMAKE_C_PREPROCESSOR   "${MINGW_ARCH}-w64-mingw32-cpp"      CACHE PATH "Preprocessor")
set(CMAKE_STRIP            "${MINGW_ARCH}-w64-mingw32-strip"    CACHE PATH "strip")
set(CMAKE_AR               "${MINGW_ARCH}-w64-mingw32-ar"       CACHE PATH "Archiver")
set(CMAKE_LINKER           "${MINGW_ARCH}-w64-mingw32-ld"       CACHE PATH "Linker")
set(CMAKE_NM               "${MINGW_ARCH}-w64-mingw32-nm"       CACHE PATH "nm")
set(CMAKE_OBJCOPY          "${MINGW_ARCH}-w64-mingw32-objcopy"  CACHE PATH "objcopy")
set(CMAKE_OBJDUMP          "${MINGW_ARCH}-w64-mingw32-objdump"  CACHE PATH "objdump")
set(CMAKE_RANLIB           "${MINGW_ARCH}-w64-mingw32-ranlib"   CACHE PATH "ranlib")
set(CMAKE_DLLTOOL          "${MINGW_ARCH}-w64-mingw32-dlltool"  CACHE PATH "dlltool")
set(CMAKE_RC_COMPILER      "${MINGW_ARCH}-w64-mingw32-windres"  CACHE PATH "WindowsRC")

set(CMAKE_FIND_ROOT_PATH /usr/${MINGW_ARCH}-w64-mingw32)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
