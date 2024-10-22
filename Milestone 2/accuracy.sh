#!/bin/sh
g++ -o electron_simulation electron_simulation.cpp -lhdf5_cpp -lhdf5 -std=c++11

# Submit the batch script
sbatch accuracy-run.sh