# LibtoolEmulator.cmake - libtool versioning for cmake
# Copyright (c) 2021, Christopher Jeffrey (MIT License).
# https://github.com/chjj

if(DEFINED __LIBTOOL_EMU)
  return()
endif()

set(__LIBTOOL_EMU 1)

if(CMAKE_C_COMPILER_LOADED)
  include(CheckSymbolExists)
elseif(CMAKE_CXX_COMPILER_LOADED)
  include(CheckCXXSymbolExists)
endif()

#
# Private Functions
#

function(_libtool_ld_is_gnu result)
  if(CMAKE_C_COMPILER_LOADED)
    set(cc ${CMAKE_C_COMPILER})
  elseif(CMAKE_CXX_COMPILER_LOADED)
    set(cc ${CMAKE_CXX_COMPILER})
  else()
    set(${result} 0 PARENT_SCOPE)
    return()
  endif()

  execute_process(COMMAND ${cc} -Wl,-v /dev/null
                  OUTPUT_VARIABLE stdout
                  ERROR_VARIABLE stderr)

  if("${stdout}" MATCHES "^GNU ld")
    set(${result} 1 PARENT_SCOPE)
  else()
    set(${result} 0 PARENT_SCOPE)
  endif()
endfunction()

function(_libtool_has_elf result)
  if(CMAKE_C_COMPILER_LOADED)
    check_symbol_exists(__ELF__ "" __LIBTOOL_HAS_ELF)
  elseif(CMAKE_CXX_COMPILER_LOADED)
    check_cxx_symbol_exists(__ELF__ "" __LIBTOOL_HAS_ELF)
  else()
    set(${result} 1 PARENT_SCOPE)
    return()
  endif()

  if(__LIBTOOL_HAS_ELF)
    set(${result} 1 PARENT_SCOPE)
  else()
    set(${result} 0 PARENT_SCOPE)
  endif()
endfunction()

function(_libtool_macho_versions target compat_version current_version)
  if(CMAKE_VERSION VERSION_LESS "3.17")
    set(ldflags -Wl,-compatibility_version,${compat_version}
                -Wl,-current_version,${current_version})

    if(COMMAND target_link_options)
      target_link_options(${target} PRIVATE ${ldflags})
    else()
      target_link_libraries(${target} PRIVATE ${ldflags})
    endif()
  else()
    set_target_properties(${target} PROPERTIES
                          MACHO_COMPATIBILITY_VERSION ${compat_version}
                          MACHO_CURRENT_VERSION ${current_version})
  endif()
endfunction()

function(_libtool_version target scheme info isnum)
  # https://github.com/autotools-mirror/libtool/blob/544fc0e/build-aux/ltmain.in#L6898
  # https://github.com/autotools-mirror/libtool/blob/544fc0e/build-aux/ltmain.in#L6973
  if (NOT info MATCHES "^[0-9]+(:[0-9]+(:[0-9]+)?)?$")
    message(FATAL_ERROR "'${info}' is not valid version information.")
  endif()

  string(REPLACE ":" ";" parts ${info})

  list(GET parts 0 current)
  list(GET parts 1 revision)
  list(GET parts 2 age)

  if(NOT revision)
    set(revision 0)
  endif()

  if(NOT age)
    set(age 0)
  endif()

  set(irix_inc 1)

  if(isnum)
    set(major ${current})
    set(minor ${revision})
    set(patch ${age})

    if(scheme MATCHES "darwin|freebsd-elf|linux|osf|windows")
      math(EXPR current "${major} + ${minor}")
      set(age ${minor})
      set(revision ${patch})
    elseif(scheme MATCHES "freebsd-aout|qnx|sunos")
      set(current ${major})
      set(revision ${minor})
      set(age 0)
    elseif(scheme MATCHES "irix|nonstopux")
      math(EXPR current "${major} + ${minor}")
      set(age ${minor})
      set(revision ${minor})
      set(irix_inc 0)
    endif()
  endif()

  if(age GREATER current)
    message(WARNING "AGE '${age}' is greater than the "
                    "current interface number '${current}'.")
    message(FATAL_ERROR "'${info}' is not valid version information.")
  endif()

  set(major)
  set(version)
  set(compat_version)
  set(current_version)

  if(scheme STREQUAL "darwin")
    math(EXPR major "${current} - ${age}")
    set(version "${major}.${age}.${revision}")
    math(EXPR compat_version "${current} + 1")
    set(current_version "${compat_version}.${revision}")
  elseif(scheme STREQUAL "freebsd-aout")
    set(major ${current})
    set(version "${current}.${revision}")
  elseif(scheme STREQUAL "freebsd-elf")
    math(EXPR major "${current} - ${age}")
    set(version "${major}.${age}.${revision}")
  elseif(scheme MATCHES "irix|nonstopux")
    math(EXPR major "${current} - ${age} + ${irix_inc}")
    set(version "${major}.${revision}")
  elseif(scheme STREQUAL "linux")
    math(EXPR major "${current} - ${age}")
    set(version "${major}.${age}.${revision}")
  elseif(scheme STREQUAL "osf")
    math(EXPR major "${current} - ${age}")
    set(version "${current}.${age}.${revision}")
  elseif(scheme STREQUAL "qnx")
    set(major ${current})
    set(version ${current})
  elseif(scheme STREQUAL "sco")
    set(major ${current})
    set(version ${current})
  elseif(scheme MATCHES "sunos")
    set(major ${current})
    set(version "${current}.${revision}")
  elseif(scheme STREQUAL "windows")
    math(EXPR major "${current} - ${age}")
    set(version ${major})
  else()
    message(FATAL_ERROR "Invalid versioning scheme.")
  endif()

  if(scheme STREQUAL "darwin")
    set_target_properties(${target} PROPERTIES VERSION ${major}
                                               SOVERSION ${major})
    _libtool_macho_versions(${target} ${compat_version} ${current_version})
  elseif(scheme STREQUAL "windows" AND MINGW)
    # Libtool does _not_ set any version information. Also, CMake
    # notably does not add the release suffix to mingw libraries.
    get_property(name TARGET ${target} PROPERTY OUTPUT_NAME)
    set_property(TARGET ${target} PROPERTY OUTPUT_NAME "${name}-${version}")
  elseif(scheme MATCHES "-aout$")
    set_target_properties(${target} PROPERTIES VERSION ${version}
                                               SOVERSION ${version})
  else()
    set_target_properties(${target} PROPERTIES VERSION ${version}
                                               SOVERSION ${major})
  endif()
endfunction()

# The CMake docs claim that CMAKE_SYSTEM_NAME is simply
# the result of `uname -s`. This is not true. It makes
# the following substitutions:
#
# OS             uname -s           substitution
# ----------------------------------------------
# Windows        [none]             Windows
# Windows CE     [none]             WindowsCE
# Windows Phone  [none]             WindowsPhone
# Windows Store  [none]             WindowsStore
# Android        Linux              Android
# BSD/OS         BSD/OS             BSDOS
# kFreeBSD       GNU/kFreeBSD       kFreeBSD
# Cygwin         CYGWIN_NT-*        CYGWIN
# Msys           MSYS_NT-*          MSYS
# z/OS           OS/390             OS390 (??)
# OS/2           OS/2               OS2 (??)
# Tru64 UNIX     OSF1               Tru64 (??)
# unknown        [falsey value]     UnknownOS
function(_libtool_set_version target info isnum)
  # https://github.com/autotools-mirror/libtool/blob/544fc0e/m4/libtool.m4#L2407
  # https://github.com/autotools-mirror/autoconf/blob/378351d/build-aux/config.guess#L179
  if(CMAKE_SYSTEM_NAME STREQUAL "AIX")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "[Aa]miga[Oo][Ss]"
      OR CMAKE_SYSTEM_PROCESSOR MATCHES "^Amiga")
    # Don't know what type this is. Use linux for now.
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Android" OR ANDROID)
    # Android doesn't support versioning.
  elseif(CMAKE_SYSTEM_NAME STREQUAL "BeOS" OR BEOS)
    # BeOS doesn't support versioning.
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Bitrig") # OpenBSD
    _libtool_version(${target} sunos-aout ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "BSDOS")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^CYGWIN" OR CYGWIN)
    _libtool_version(${target} windows ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" OR APPLE OR IOS)
    _libtool_version(${target} darwin ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "dgux")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "DragonFly")
    _libtool_version(${target} freebsd-elf ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "DYNIX/ptx")
    _libtool_version(${target} linux ${info} ${isnum}) # sysv4
  elseif(CMAKE_SYSTEM_NAME STREQUAL "ekkoBSD") # OpenBSD
    _libtool_version(${target} sunos-aout ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    if(CMAKE_SYSTEM_VERSION MATCHES "^[23]")
      _libtool_version(${target} freebsd-aout ${info} ${isnum})
    else()
      _libtool_version(${target} freebsd-elf ${info} ${isnum})
    endif()
  elseif(CMAKE_SYSTEM_NAME MATCHES "^GNU")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Haiku" OR HAIKU)
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "HP-UX")
    _libtool_version(${target} sunos-elf ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Interix")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^IRIX$|^IRIX64$")
    _libtool_ld_is_gnu(gnu_ld)
    if(gnu_ld)
      _libtool_version(${target} linux ${info} ${isnum})
    else()
      _libtool_version(${target} irix ${info} ${isnum})
    endif()
  elseif(CMAKE_SYSTEM_NAME STREQUAL "kFreeBSD")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "LibertyBSD") # OpenBSD
    _libtool_version(${target} sunos-aout ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "MidnightBSD") # FreeBSD
    _libtool_version(${target} freebsd-elf ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^MINGW" OR MINGW)
    _libtool_version(${target} windows ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "MirBSD") # OpenBSD
    _libtool_version(${target} sunos-aout ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "MP-RAS") # ??
    _libtool_version(${target} linux ${info} ${isnum}) # sysv4
  elseif(CMAKE_SYSTEM_NAME MATCHES "^MSYS" OR MSYS)
    _libtool_version(${target} windows ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "NetBSD")
    _libtool_has_elf(has_elf)
    if(has_elf)
      _libtool_version(${target} sunos-elf ${info} ${isnum})
    else()
      _libtool_version(${target} sunos-aout ${info} ${isnum})
    endif()
  elseif(CMAKE_SYSTEM_NAME STREQUAL "NEWS-OS")
    _libtool_version(${target} linux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^NONSTOP_KERNEL$|^NonStop-UX$")
    _libtool_version(${target} nonstopux ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "OpenBSD")
    _libtool_version(${target} sunos-aout ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^OS/?2$" OR OS2)
    _libtool_version(${target} windows ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^OSF1$|^Tru64$")
    _libtool_version(${target} osf ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^PW") # PW32
    _libtool_version(${target} windows ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "QNX" OR QNXNTO)
    _libtool_version(${target} qnx ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "rdos")
    # No dynamic linker.
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Redox")
    # Unknown.
  elseif(CMAKE_SYSTEM_NAME STREQUAL "ReliantUNIX")
    _libtool_version(${target} linux ${info} ${isnum}) # sysv4
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Rhapsody")
    _libtool_version(${target} darwin ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^SCO_SV$|^UnixWare$")
    _libtool_version(${target} sco ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "SINIX")
    _libtool_version(${target} linux ${info} ${isnum}) # sysv4
  elseif(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
    if(CMAKE_SYSTEM_VERSION MATCHES "^[0-4]")
      _libtool_version(${target} sunos-aout ${info} ${isnum})
    else()
      _libtool_version(${target} linux ${info} ${isnum})
    endif()
  elseif(CMAKE_SYSTEM_NAME STREQUAL "syllable")
    # Unknown.
  elseif(CMAKE_SYSTEM_NAME STREQUAL "TPF")
    # Cross-compile only. Assume linux.
    _libtool_version(${target} linux ${info})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^UNIX_S(ystem_)?V$"
     AND CMAKE_SYSTEM_VERSION MATCHES "^4")
    _libtool_version(${target} linux ${info} ${isnum}) # sysv4
  elseif(CMAKE_SYSTEM_NAME STREQUAL "UTS4") ## ??
    _libtool_version(${target} linux ${info})
  elseif(CMAKE_SYSTEM_NAME MATCHES "^Windows"
      OR WIN32 OR WINCE OR WINDOWS_PHONE OR WINDOWS_STORE)
    _libtool_version(${target} windows ${info} ${isnum})
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Zircon" OR FUCHSIA)
    # Probably linux given that fuchsia is gnu-like.
    _libtool_version(${target} linux ${info} ${isnum})
  else()
    # No dynamic linker.
  endif()
endfunction()

#
# Public Functions
#

function(set_target_version_info target info)
  _libtool_set_version(${target} ${info} 0)
endfunction()

function(set_target_version_number target number)
  _libtool_set_version(${target} ${number} 1)
endfunction()

function(set_target_release target version)
  if (NOT version MATCHES "^[0-9]+(\.[0-9]+(\.[0-9]+)?)?$")
    message(FATAL_ERROR "'${version}' is not a valid release number.")
  endif()

  get_property(name TARGET ${target} PROPERTY OUTPUT_NAME)

  string(REGEX REPLACE "-[0-9]+(\.[0-9]+(\.[0-9]+)?)?$" "" name "${name}")

  set_property(TARGET ${target} PROPERTY OUTPUT_NAME "${name}-${version}")
endfunction()
