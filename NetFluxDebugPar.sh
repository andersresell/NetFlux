#! /bin/bash

sim_dir=$1

clear 
cd ./src
echo building debug
make debug
exit_status=$?
if [ $exit_status -eq 0 ]; then
    cd ../build_debug
    echo "Build successful. Running..."
    mpirun -np 2 ./NetFlux "../"$sim_dir
fi