#!/bin/sh

set -ex
exec valgrind --leak-check=full --show-leak-kinds=all ./test
