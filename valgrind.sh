#!/bin/sh

set -ex
export LD_LIBRARY_PATH=./.libs
exec valgrind \
  --track-origins=yes \
  --leak-check=full \
  --show-leak-kinds=all \
  ./tests "$@"
