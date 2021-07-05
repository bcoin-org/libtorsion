# CheckCThreadLocalStorage.cmake - tls check for c
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND check_c_thread_local_storage)
  return()
endif()

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)

function(check_c_thread_local_storage name)
  set(${name} "" PARENT_SCOPE)
  set(CMAKE_TLS_LIBS_INIT "" PARENT_SCOPE)

  set(CMAKE_REQUIRED_FLAGS "")
  set(CMAKE_REQUIRED_LIBRARIES "")

  # XL requires a special flag. Don't ask me why.
  # Note that CMake handles -qthreaded for us.
  if(CMAKE_C_COMPILER_ID MATCHES "^XL")
    check_c_compiler_flag(-qtls CMAKE_HAVE_XL_TLS)
    if(CMAKE_HAVE_XL_TLS)
      set(CMAKE_REQUIRED_FLAGS "-qtls")
    endif()
  endif()

  # Compilers which do TLS emulation sometimes
  # need a threads library. This includes 32-bit
  # AIX compilers (GCC and XL).
  if(CMAKE_SYSTEM_NAME STREQUAL "AIX")
    set(CMAKE_REQUIRED_LIBRARIES "-lpthread")
  endif()

  # Various TLS keywords. We prepend or append
  # _Thread_local depending on the C standard.
  # The last keyword, __declspec(__thread), is
  # not widely known, but there is evidence
  # that Compaq C for Tru64 UNIX supported it
  # at one point.
  set(keywords "__declspec(thread)" __thread "__declspec(__thread)")

  if (DEFINED CMAKE_C_STANDARD AND CMAKE_C_STANDARD LESS 11)
    list(APPEND keywords _Thread_local)
  else()
    list(INSERT keywords 0 _Thread_local)
  endif()

  foreach(keyword ${keywords})
    string(TOUPPER "CMAKE_HAVE_C_TLS_${keyword}" varname)
    string(REGEX REPLACE "[^A-Z0-9]" "_" varname "${varname}")

    check_c_source_compiles("
      static ${keyword} int value;
      int main(void) {
        value = 1;
        return 0;
      }
    " ${varname})

    if(${varname})
      set(${name} ${keyword} PARENT_SCOPE)

      if(CMAKE_REQUIRED_LIBRARIES)
        set(CMAKE_TLS_LIBS_INIT "${CMAKE_REQUIRED_LIBRARIES}" PARENT_SCOPE)
      endif()

      if(NOT TARGET Threads::TLS)
        add_library(Threads::TLS INTERFACE IMPORTED)

        if(CMAKE_REQUIRED_FLAGS)
          set_property(TARGET Threads::TLS
                       PROPERTY INTERFACE_COMPILE_OPTIONS
                                "${CMAKE_REQUIRED_FLAGS}")
        endif()

        if(CMAKE_REQUIRED_LIBRARIES)
          set_property(TARGET Threads::TLS
                       PROPERTY INTERFACE_LINK_LIBRARIES
                                "${CMAKE_REQUIRED_LIBRARIES}")
        endif()
      endif()

      break()
    endif()
  endforeach()
endfunction()
