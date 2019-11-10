#!/bin/sh

# pre-build actions:

clang-format -i -style=file ./src/* ./include/**/* # auto-format
export CXX=/usr/bin/clang++ # g++ has appalling errors

mkdir -p build
pushd build

# build:

rm -f CMakeCache.txt
cmake .. -DTESTS=ON
make

# test:

./mth_test

popd
