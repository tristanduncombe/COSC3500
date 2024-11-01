#!/bin/sh
#SBATCH --job-name="electron"
#SBATCH --partition=cosc3500
#SBATCH --account=cosc3500
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1M
#SBATCH --time=0-00:08:00

# Load necessary modules
module load compiler-rt/latest
module add mkl/latest
module add mpi/openmpi-x86_64
module load cuda/11.1

export PATH=/opt/local/stow/cuda-11.1/bin:$PATH
export PATH=/usr/lib64/openmpi/bin:$PATH

# Run the simulation using srun (recommended with Slurm)
mpiexec ./electron_simulation -timing