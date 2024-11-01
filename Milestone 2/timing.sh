#!/bin/bash

# Compile the C++ code with MPI, OpenMP, and AVX support using g++
g++ -std=c++11 -O2 -mavx -fopenmp \
    -I/usr/include/openmpi-x86_64 \
    -L/usr/lib64/openmpi/lib \
    -I/usr/include/hdf5/serial \
    -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
    electron_simulation.cpp -o electron_simulation \
    -lmpi_cxx -lmpi \
    -lhdf5_cpp -lhdf5

# Submit the batch script to Slurm
sbatch timing_run.sh