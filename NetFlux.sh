#! /bin/bash

clear 
cd ./src
make NetFlux
exit_status=$?
if [ $exit_status -eq 0 ]; then
    echo "Build successful. Running..."
    ../build/NetFlux
fi