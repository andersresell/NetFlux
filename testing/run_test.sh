#!/bin/bash

clear

FLAGS="-O2 -Wall -DNDEBUG" 
FLAGS="-g -Wall -std=c++20"

g++ -o test test.cpp $FLAGS
build_status=$?
if [ $build_status == 0 ]; then
    ./test
fi
