#! /bin/bash

sim_dir=$1

clear 
cd ./src
make NetFlux
exit_status=$?
if [ $exit_status -eq 0 ]; then
    cd ../build
    echo "Build successful. Running..."
    ./NetFlux "../"$sim_dir
fi