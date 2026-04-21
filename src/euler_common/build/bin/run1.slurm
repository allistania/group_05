#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=python_job
#SBATCH --time=0:10:00
#SBATCH --output=./stdout/slurm-%j.out



python3 graph.py