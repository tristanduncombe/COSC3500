#!/bin/sh
#SBATCH --job-name="electron"
#SBATCH --partition=a100
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

export PATH=/usr/lib64/openmpi/bin:$PATH

# Run the compiled program
./electron_simulation