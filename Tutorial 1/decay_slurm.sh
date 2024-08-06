#!/bin/bash
#SBATCH --job-name=decay
#SBATCH --partition=cosc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G # memory (MB)
#SBATCH --time=0-00:01 # time (D-HH:MM)

time ./decay
