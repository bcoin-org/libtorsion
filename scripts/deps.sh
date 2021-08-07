#!/bin/sh

# deps.sh - dependency generator for libtorsion
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/bcoin-org/libtorsion

set -e

global_host=unix
global_oflag='o'
global_cflag='-c'

deps_get() {
  if gcc --version > /dev/null 2>& 1; then
    get_cc=gcc
  elif clang --version > /dev/null 2>& 1; then
    get_cc=clang
  else
    echo 'No compiler available.' >& 2
    return 1
  fi

  $get_cc -I./include -MF - -MM "$1" \
    | sed -e 's/^ *//'               \
          -e 's/^.*: .*\.c //'       \
          -e 's/\\$//g'              \
          -e 's/ *$//g'              \
    | tr ' ' '\n'                    \
    | sort -u                        \
    | tr '\n' ' '
}

deps_format() {
  fmt_src=$1
  fmt_obj=`basename $fmt_src | sed -e "s/c$/$global_oflag/"`
  fmt_rule="$fmt_obj: $fmt_src "
  fmt_space=`echo "$fmt_rule" | sed -e 's/./ /g'`
  fmt_deps=`deps_get $fmt_src`
  fmt_len=`echo "$fmt_deps" | wc -w`
  fmt_last=`expr $fmt_len - 1 || true`
  fmt_max=0
  fmt_cur=0

  for fmt_dep in $fmt_deps; do
    fmt_cnt=`echo "$fmt_space$fmt_dep" | wc -c`

    if test $fmt_cnt -gt $fmt_max; then
      fmt_max=$fmt_cnt
    fi
  done

  for fmt_dep in $fmt_deps; do
    fmt_dep="$fmt_rule$fmt_dep"
    fmt_cnt=`echo "$fmt_dep" | wc -c`

    if test $fmt_cur -eq $fmt_last; then
      echo "$fmt_dep"
      break
    fi

    while test $fmt_cnt -lt $fmt_max; do
      fmt_dep="$fmt_dep "
      fmt_cnt=`expr $fmt_cnt + 1`
    done

    echo "$fmt_dep \\"

    fmt_rule="$fmt_space"
    fmt_cur=`expr $fmt_cur + 1`
  done

  echo '	$(CC) '$global_cflag' $(LIB_CFLAGS) '$fmt_src
  echo ''
}

deps_print() {
  echo '#'
  echo '# Library Objects'
  echo '#'
  echo ''

  deps_format src/aead.c
  deps_format src/asn1.c
  deps_format src/cipher.c
  deps_format src/ecc.c
  deps_format src/encoding.c
  if test x"$global_host" != x'dmc'; then
    deps_format src/entropy/hw.c
    deps_format src/entropy/sys.c
  fi
  deps_format src/drbg.c
  deps_format src/dsa.c
  deps_format src/hash.c
  deps_format src/ies.c
  deps_format src/internal.c
  deps_format src/kdf.c
  deps_format src/mac.c
  deps_format src/mpi.c
  if test x"$global_host" != x'dmc'; then
    deps_format src/rand.c
  fi
  deps_format src/rsa.c
  deps_format src/stream.c
  deps_format src/util.c

  echo '#'
  echo '# Test Objects'
  echo '#'
  echo ''

  deps_format test/bench.c
  if test x"$global_host" = x'unix'; then
    deps_format test/ctgrind.c
  fi
  deps_format test/hrtime.c
  deps_format test/test.c
  deps_format test/utils.c
}

deps_replace() {
  rep_top=0
  rep_bot=0
  rep_cur=0

  while IFS= read -r rep_line; do
    if test x"$rep_line" = x'# Library Objects'; then
      rep_top=`expr $rep_cur - 1`
    fi
    if echo x"$rep_line" | grep '^x	$(CC) '$global_cflag' ' > /dev/null; then
      rep_bot=`expr $rep_cur + 1`
    fi
    rep_cur=`expr $rep_cur + 1`
  done < "$1"

  if test $rep_top -eq 0 -o $rep_bot -eq 0; then
    echo 'Could not find objects.' >& 2
    return 1
  fi

  rep_cur=0

  while IFS= read -r rep_line; do
    if test $rep_cur -eq $rep_top; then
      deps_print
    elif test $rep_cur -ge $rep_top -a $rep_cur -le $rep_bot; then
      :
    else
      echo x"$rep_line" | sed -e 's/^x//'
    fi
    rep_cur=`expr $rep_cur + 1`
  done < "$1"
}

main() {
  main_action=print
  main_file=/dev/null

  while test $# != 0; do
    case "$1" in
      --unix|-u)
        global_host=unix
        global_oflag='o'
        global_cflag='-c'
      ;;
      --msvc|-m)
        global_host=msvc
        global_oflag='obj'
        global_cflag='/c'
      ;;
      --dmc|-d)
        global_host=dmc
        global_oflag='obj'
        global_cflag='-c'
      ;;
      --print|-p)
        main_action=print
      ;;
      --replace|-i)
        main_action=replace
        main_file="$2"
        shift
      ;;
      --version|-v)
        main_action=version
      ;;
      --help|-h)
        main_action=help
      ;;
      *)
        echo "Unknown option '$main_arg'." >& 2
        return 1
      ;;
    esac
    shift
  done

  if ! test -d .git; then
    echo 'Invalid working directory.' >& 2
    return 1
  fi

  if ! test -f include/torsion/ecc.h; then
    echo 'Invalid working directory.' >& 2
    return 1
  fi

  case "$main_action" in
    print)
      deps_print
    ;;
    replace)
      main_buf=`deps_replace "$main_file"`
      echo "$main_buf" > "$main_file"
    ;;
    version)
      echo 'deps.sh 0.0.0'
    ;;
    help)
      echo '  Usage: deps.sh [options]'
      echo ''
      echo '  Options:'
      echo ''
      echo '    -u, --unix            use unix arguments'
      echo '    -m, --msvc            use msvc arguments'
      echo '    -d, --dmc             use dmc arguments'
      echo '    -p, --print           print makefile rules'
      echo '    -i, --replace <file>  modify file in place'
      echo '    -v, --version         output version number'
      echo '    -h, --help            output usage information'
      echo ''
    ;;
    *)
      echo "Unknown action '$main_action'." >& 2
      return 1
    ;;
  esac
}

main "$@"
