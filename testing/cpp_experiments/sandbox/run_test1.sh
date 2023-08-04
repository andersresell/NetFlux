#!/bin/bash

clear

FLAGS="-g -Wall -std=c++20  -I/usr/include/x86_64-linux-gnu/openmpi -L/usr/lib/x86_64-linux-gnu/openmpi/lib/ -lmpi -lyaml-cpp -lboost_serialization"


g++ -o test test1.cpp $FLAGS
build_status=$?
if [ $build_status == 0 ]; then
    #mpirun -n 1 ./test
    ./test
    rm test
fi
