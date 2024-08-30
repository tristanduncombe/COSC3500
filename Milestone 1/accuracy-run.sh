#!/bin/sh
#SBATCH --job-name="electron"
#SBATCH --partition=a100

# Run the program with the -timing argument
./electron_simulation -accuracy