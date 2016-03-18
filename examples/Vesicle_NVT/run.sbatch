#!/bin/bash
#SBATCH --job-name="myMPI"
#SBATCH --output="myMPI.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH -t 02:00:00

# Comet sbatch script 

pwd

time ./run.sh
