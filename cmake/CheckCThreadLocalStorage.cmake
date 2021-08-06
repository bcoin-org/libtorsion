# CheckCThreadLocalStorage.cmake - tls check for c
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND check_c_thread_local_storage)
  return()
endif()

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)
include(CheckCSourceRuns)

function(_check_c_emutls name code flags result)
  string(TOUPPER "${result}_${name}" var)

  set(dir ${CMAKE_BINARY_DIR}/CMakeFiles/CheckEmuTLS)
  set(src ${dir}/${name}.c)
  set(bin ${dir}/${name}${CMAKE_EXECUTABLE_SUFFIX})
  set(found 0)

  file(MAKE_DIRECTORY ${dir})
  file(WRITE ${src} "${code}\n")

  try_compile(${var} ${CMAKE_BINARY_DIR} ${src} COPY_FILE ${bin}
              CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${flags})

  if(${var} AND EXISTS "${bin}")
    # There is evidence that some non-GNU platforms also do TLS
    # emulation. It's possible this includes 32-bit AIX, but I
    # cannot confirm this.
    #
    # TODO: Find other platforms with emulated TLS and figure
    #       out how to detect it.
    file(STRINGS ${bin} emutls LIMIT_COUNT 1 REGEX "__emutls_get_address")

    if(emutls)
      set(found 1)
    endif()
  endif()

  file(REMOVE_RECURSE "${dir}")

  set(${result} ${found} PARENT_SCOPE)
endfunction()

function(check_c_thread_local_storage result_keyword result_flags result_emutls)
  if(DEFINED "${result_keyword}" AND DEFINED "${result_flags}"
                                 AND DEFINED "${result_emutls}")
    return()
  endif()

  if(NOT CMAKE_REQUIRED_QUIET AND NOT CMAKE_VERSION VERSION_LESS 3.17)
    message(CHECK_START "Checking for thread-local storage")
    set(verbose 1)
  else()
    set(verbose 0)
  endif()

  set(CMAKE_REQUIRED_FLAGS "")
  set(CMAKE_REQUIRED_QUIET 1)

  # XL requires a special flag. Don't ask me why.
  # Note that CMake handles -qthreaded for us.
  if(CMAKE_C_COMPILER_ID MATCHES "^XL")
    check_c_compiler_flag(-qtls CMAKE_C_HAVE_FLAG_QTLS)
    if(CMAKE_C_HAVE_FLAG_QTLS)
      set(CMAKE_REQUIRED_FLAGS "-qtls")
    endif()
  endif()

  # Various TLS keywords.
  #
  # The last keyword is not widely known, but there is evidence
  # that Compaq C for Tru64 UNIX supported it at one point.
  set(keywords __thread "__declspec(thread)" "__declspec(__thread)")

  # Prepend or append _Thread_local according to the C standard.
  if (DEFINED CMAKE_C_STANDARD AND CMAKE_C_STANDARD LESS 11)
    list(APPEND keywords _Thread_local)
  else()
    list(INSERT keywords 0 _Thread_local)
  endif()

  # We try to run the executable when not cross compiling. There
  # are far too many instances of TLS code successfully building
  # but not running.
  set(tls "")
  set(flags "")
  set(emutls 0)

  foreach(keyword ${keywords})
    string(REGEX REPLACE "[^0-9A-Za-z]" "_" name "${keyword}")
    string(TOUPPER "${result_keyword}_${name}" var)

    # The thread-local variable must have external linkage otherwise
    # the optimizer may remove the TLS code. GCC and Clang refuse to
    # optimize the below code (even with -O3 enabled).
    set(code "${keyword} int x; int main(void) { x = 1; return !x; }")

    if(CMAKE_CROSSCOMPILING)
      check_c_source_compiles("${code}" ${var})
    else()
      check_c_source_runs("${code}" ${var})
    endif()

    if(${var})
      set(tls "${keyword}")
      set(flags "${CMAKE_REQUIRED_FLAGS}")
      _check_c_emutls(${name} "${code}" "${flags}" emutls)
      break()
    endif()
  endforeach()

  set(${result_keyword} "${tls}" CACHE INTERNAL "TLS keyword")
  set(${result_flags} "${flags}" CACHE INTERNAL "TLS flags")
  set(${result_emutls} ${emutls} CACHE INTERNAL "TLS emulation")

  if(verbose)
    if(tls AND flags)
      message(CHECK_PASS "${tls} (flags=${flags}, emulated=${emutls})")
    elseif(tls)
      message(CHECK_PASS "${tls} (emulated=${emutls})")
    else()
      message(CHECK_FAIL "Failed")
    endif()
  endif()
endfunction()
