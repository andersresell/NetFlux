#!/bin/bash

FLAGS="-O2 -Wall -DNDEBUG" 
FLAGS="-g -Wall"

g++ -o test test.cpp $FLAGS

./test
