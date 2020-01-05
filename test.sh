#!/bin/sh

./build.sh
exec ./test "$@"
