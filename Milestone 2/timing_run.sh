#!/bin/sh
#SBATCH --job-name="electron"
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1M
#SBATCH --time=0-00:08:00

# Run the program with the -timing argument
./electron_simulation -timing