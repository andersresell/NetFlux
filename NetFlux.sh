#! /bin/bash

config_file=$1

clear 
cd ./src
make NetFlux
exit_status=$?
if [ $exit_status -eq 0 ]; then
    cd ../build

    echo "Build successful. Running..."
    ./NetFlux "../"$config_file
fi