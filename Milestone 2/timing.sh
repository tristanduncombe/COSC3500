#!/bin/sh

# Compile the C++ code with OpenMP and AVX support
g++ -o electron_simulation electron_simulation.cpp -std=c++11 -fopenmp -mavx -lhdf5_cpp -lhdf5

# Submit the batch script
sbatch timing_run.sh