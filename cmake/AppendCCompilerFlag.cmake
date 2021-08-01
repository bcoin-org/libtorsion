# AppendCCompilerFlag.cmake - checked c flags appending
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND append_c_compiler_flag)
  return()
endif()

include(CheckCCompilerFlag)

function(append_c_compiler_flag list)
  check_c_compiler_flag(-Werror=unknown-warning-option
                        CMAKE_HAVE_UNKNOWN_WARNING_OPTION)

  foreach(flag ${ARGN})
    string(TOUPPER "CMAKE_HAVE_C_FLAG${flag}" name)
    string(REGEX REPLACE "[^A-Z0-9]" "_" name "${name}")

    if(CMAKE_HAVE_UNKNOWN_WARNING_OPTION)
      check_c_compiler_flag("-Werror=unknown-warning-option ${flag}" ${name})
    else()
      check_c_compiler_flag(${flag} ${name})
    endif()

    if(${name})
      list(APPEND ${list} ${flag})
    endif()
  endforeach()

  set(${list} ${${list}} PARENT_SCOPE)
endfunction()
