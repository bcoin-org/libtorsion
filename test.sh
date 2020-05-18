#!/bin/sh

set -ex
export LD_LIBRARY_PATH=./.libs
exec ./tests "$@"
