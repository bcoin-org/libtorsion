# CheckCThreadLocalStorage.cmake - tls check for c
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND check_c_thread_local_storage)
  return()
endif()

include(CheckCCompilerFlag)
include(CheckCSourceCompiles)

function(check_c_thread_local_storage ident_name flags_name libs_name)
  set(${ident_name} "" PARENT_SCOPE)
  set(${flags_name} "" PARENT_SCOPE)
  set(${libs_name} "" PARENT_SCOPE)

  set(CMAKE_REQUIRED_FLAGS "")
  set(CMAKE_REQUIRED_LIBRARIES "")

  if(CMAKE_C_COMPILER_ID MATCHES "^XL")
    # CMake handles -qthreaded for us.
    check_c_compiler_flag(-qtls CMAKE_HAVE_XL_TLS)
    if(CMAKE_HAVE_XL_TLS)
      set(CMAKE_REQUIRED_FLAGS "-qtls")
    endif()
  endif()

  if(CMAKE_SYSTEM_NAME STREQUAL "AIX")
    set(CMAKE_REQUIRED_LIBRARIES "-lpthread")
  endif()

  set(keywords "__declspec(thread)" __thread)

  if (DEFINED CMAKE_C_STANDARD AND CMAKE_C_STANDARD LESS "11")
    list(APPEND keywords _Thread_local)
  else()
    list(PREPEND keywords _Thread_local)
  endif()

  foreach(keyword IN LISTS keywords)
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
      set(${ident_name} ${keyword} PARENT_SCOPE)
      set(${flags_name} "${CMAKE_REQUIRED_FLAGS}" PARENT_SCOPE)
      set(${libs_name} "${CMAKE_REQUIRED_LIBRARIES}" PARENT_SCOPE)
      break()
    endif()
  endforeach()
endfunction()
