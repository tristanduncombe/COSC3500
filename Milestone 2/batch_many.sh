#!/bin/bash

# Compile the C++ code with OpenMP, AVX, and MPI support
g++ -std=c++11 -O2 -mavx -fopenmp \
    -I/usr/include/openmpi-x86_64 \
    -L/usr/lib64/openmpi/lib \
    -I/usr/include/hdf5/serial \
    -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
    electron_simulation.cpp -o electron_simulation \
    -lmpi_cxx -lmpi \
    -lhdf5_cpp -lhdf5

# Function to submit a job with sbatch
submit_job() {
    local numElectrons=$1
    local numFrames=$2
    local jobName="electron_${numElectrons}_${numFrames}"

    sbatch --job-name="$jobName" \
           --partition=cosc3500 \
           --account=cosc3500 \
           --nodes=3 \
           --ntasks=3 \
           --cpus-per-task=4 \
           --time=0-00:08:00 \
           --output=output_${numElectrons}_${numFrames}_%j.txt \
           --error=error_${numElectrons}_${numFrames}_%j.txt \
           --wrap="module load compiler-rt/latest mkl/latest mpi/openmpi-x86_64 && \
                   export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK && \
                   cd \"\$SLURM_SUBMIT_DIR\" && \
                   ./electron_simulation -timing -electrons $numElectrons -frames $numFrames -no-hdf5"
}

# Test numElectrons from 1000 to 10000 with 1000 increments
for numElectrons in {1000..10000..1000}; do
    submit_job $numElectrons 1000
done

# Test numFrames from 1000 to 10000 with 1000 increments
for numFrames in {1000..10000..1000}; do
    submit_job 1000 $numFrames
done