#!/bin/bash

clear

FLAGS="-g -Wall -std=c++20  -lyaml-cpp -lboost_serialization"


mpicxx -o mpidbg mpidbg.cpp $FLAGS
build_status=$?
if [ $build_status == 0 ]; then
    mpirun -n 3 ./test
    rm test
fi
