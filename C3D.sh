#! /bin/bash

clear 
cd ./src
make C3D
exit_status=$?
if [ $exit_status -eq 0 ]; then
    echo "Build successful. Running..."
    ../build/C3D
fi