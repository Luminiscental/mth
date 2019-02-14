#!/bin/sh

mkdir -p build

pushd build

rm -f CMakeCache.txt
cmake ..
make
./mth_test

popd
