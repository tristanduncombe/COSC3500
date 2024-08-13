#!/bin/sh
gcc -o electron_simulation electron_simulation.cpp -lstdc++ -o -lhdf5 -lhdf5_cpp -lm -std=c++11 
./electron_simulation