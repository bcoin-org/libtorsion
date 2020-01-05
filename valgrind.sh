#!/bin/sh

set -ex
exec valgrind --track-origins=yes --leak-check=full --show-leak-kinds=all ./test
