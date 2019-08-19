#!/bin/sh

mkdir -p build

pushd build

rm -f CMakeCache.txt
cmake .. -DTESTS=ON
make
./mth_test

popd
