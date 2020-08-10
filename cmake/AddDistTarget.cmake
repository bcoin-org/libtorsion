# AddDistTarget.cmake - dist target for cmake
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(COMMAND add_dist_target)
  return()
endif()

function(add_dist_target)
  if(NOT PROJECT_NAME)
    message(FATAL_ERROR "add_dist_target(): project must be initialized")
  endif()

  if(NOT UNIX OR NOT CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    return()
  endif()

  if(NOT PROJECT_NAME STREQUAL CMAKE_PROJECT_NAME)
    return()
  endif()

  if(PROJECT_VERSION)
    set(distname ${PROJECT_NAME}-${PROJECT_VERSION})
  else()
    set(distname ${PROJECT_NAME})
  endif()

  set(tarname ${distname}.tar)
  set(distfiles)

  foreach(file IN LISTS ARGN)
    if(file STREQUAL "" OR file MATCHES "^/")
      message(FATAL_ERROR "add_dist_target(): invalid filename")
    endif()
    list(APPEND distfiles ${PROJECT_SOURCE_DIR}/${file})
  endforeach()

  add_custom_target(dist
    COMMAND rm -rf ${distname} ${tarname} ${tarname}.gz
    COMMAND mkdir ${distname}
    COMMAND cp -r ${distfiles} ${distname}/
    COMMAND tar chof ${tarname} ${distname}
    COMMAND gzip --best ${tarname}
    COMMAND rm -rf ${distname}
  )

  add_custom_target(distclean
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND rm -rf ${distname} ${tarname} ${tarname}.gz
  )
endfunction()
