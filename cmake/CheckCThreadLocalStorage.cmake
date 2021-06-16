# CheckCThreadLocalStorage.cmake - tls check for c
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND check_c_thread_local_storage)
  return()
endif()

include(CheckCSourceCompiles)

function(check_c_thread_local_storage name)
  set(keywords "__declspec(thread)" __thread)

  if (DEFINED CMAKE_C_STANDARD AND CMAKE_C_STANDARD LESS "11")
    list(APPEND keywords _Thread_local)
  else()
    list(PREPEND keywords _Thread_local)
  endif()

  set(${name} "" PARENT_SCOPE)

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
      set(${name} ${keyword} PARENT_SCOPE)
      break()
    endif()
  endforeach()
endfunction()
