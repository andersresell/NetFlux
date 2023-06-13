#!/bin/bash

clear

FLAGS="-O2 -Wall -DNDEBUG" 
FLAGS="-g -Wall -std=c++20"
FLAGS="-g -Wall"

g++ -o test test.cpp -lyaml-cpp $FLAGS
build_status=$?
if [ $build_status == 0 ]; then
    ./test
fi
