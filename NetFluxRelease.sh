#! /bin/bash

sim_dir=$1

clear 
cd ./src
echo building release
make release
exit_status=$?
if [ $exit_status -eq 0 ]; then
    cd ../build_release
    echo "Build successful. Running..."
    ./NetFlux "../"$sim_dir
fi