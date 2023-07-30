#!/bin/bash

clear

FLAGS="-g -Wall -std=c++20  -I/usr/include/x86_64-linux-gnu/openmpi -L/usr/lib/x86_64-linux-gnu/openmpi/lib/  -lyaml-cpp -lboost_serialization"


mpicxx -o test mpi.cpp $FLAGS
build_status=$?
if [ $build_status == 0 ]; then
    mpirun -n 3 ./test
    rm test
fi
