# CheckCThreadLocalStorage.cmake - tls check for c
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND check_c_thread_local_storage)
  return()
endif()

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)
include(CheckCSourceRuns)

function(check_c_thread_local_storage name flag)
  set(${name} "" PARENT_SCOPE)
  set(${flag} "" PARENT_SCOPE)

  set(CMAKE_REQUIRED_FLAGS "")

  # XL requires a special flag. Don't ask me why.
  # Note that CMake handles -qthreaded for us.
  if(CMAKE_C_COMPILER_ID MATCHES "^XL")
    check_c_compiler_flag(-qtls CMAKE_HAVE_XL_TLS)
    if(CMAKE_HAVE_XL_TLS)
      set(CMAKE_REQUIRED_FLAGS "-qtls")
    endif()
  endif()

  # Various TLS keywords. We prepend or append
  # _Thread_local depending on the C standard.
  # The last keyword, __declspec(__thread), is
  # not widely known, but there is evidence
  # that Compaq C for Tru64 UNIX supported it
  # at one point.
  set(keywords __thread "__declspec(thread)" "__declspec(__thread)")

  if (DEFINED CMAKE_C_STANDARD AND CMAKE_C_STANDARD LESS 11)
    list(APPEND keywords _Thread_local)
  else()
    list(INSERT keywords 0 _Thread_local)
  endif()

  foreach(keyword ${keywords})
    string(TOUPPER "CMAKE_HAVE_C_TLS_${keyword}" varname)
    string(REGEX REPLACE "[^A-Z0-9]" "_" varname "${varname}")

    set(code "static ${keyword} int x; int main(void) { x = 1; return !x; }")

    if(CMAKE_CROSSCOMPILING)
      check_c_source_compiles("${code}" ${varname})
    else()
      check_c_source_runs("${code}" ${varname})
    endif()

    if(${varname})
      set(${name} ${keyword} PARENT_SCOPE)
      if(CMAKE_REQUIRED_FLAGS)
        set(${flag} ${CMAKE_REQUIRED_FLAGS} PARENT_SCOPE)
      endif()
      break()
    endif()
  endforeach()
endfunction()
