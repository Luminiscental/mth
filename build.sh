#!/bin/sh

pushd build

flag_cmake=false
flag_run=false

while getopts ":cr" opt; do
  case $opt in
    c)
        flag_cmake=true
        ;;
    r)
        flag_run=true
        ;;
  esac
done

if [ "$flag_cmake" = true ] ; then
    cmake ..
fi

make

if [ "$flag_run" = true ] ; then
    ./testing
fi

popd
