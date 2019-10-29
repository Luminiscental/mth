#!/bin/sh

mkdir -p build

clang-format -i -style=file ./src/* ./include/**/*

pushd build

rm -f CMakeCache.txt
cmake .. -DTESTS=ON
make
./mth_test

popd
