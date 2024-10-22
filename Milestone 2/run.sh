#!/bin/sh
#SBATCH --job-name="electron"
#SBATCH --partition=a100
#SBATCH --gres=gpu:1
module load cuda/11.1

# Compile the project with CUDA enabled
nvcc -o electron_simulation electron_simulation.cpp -lhdf5_cpp -lhdf5 -std=c++11 -L$(/opt/local/stow/cuda-11.1/lib64) -lcudart -lcublas

# Batch the program
sbatch ./batch.sh