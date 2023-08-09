#! /bin/bash

sim_dir=$1

clear 
cd ./src
echo building release
make release
cd ..
exit_status=$?
if [ $exit_status -eq 0 ]; then
    echo "Build successful. Running..."
    build_release/NetFlux $sim_dir
fi