#!/bin/bash

clear

FLAGS="-g -Wall -lyaml-cpp -lboost_serialization -D_GLIBCXX_DEBUG"


g++ -o test test1.cpp $FLAGS
build_status=$?
if [ $build_status == 0 ]; then
    #mpirun -n 1 ./test
    ./test
    rm test
fi
