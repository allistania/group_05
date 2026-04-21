#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=cpp_debug
#SBATCH --time=0:10:00
#SBATCH --output=./stdout/slurm-%j.out

./Shem
