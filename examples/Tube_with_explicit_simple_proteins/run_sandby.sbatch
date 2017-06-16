#!/bin/bash
#SBATCH --job-name=em2_ibar
#SBATCH --output=bd_example_sbatch.out
#SBATCH --error=bd_example_sbatch.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCG --ntasks-per-node=1
#SBATCH --constraint=ib
#SBATCH --partition=sandyb

let NCORE=SLURM_NNODES*SLURM_CPUS_ON_NODE

echo "Total number of core:", $NCORE

mpirun -np 1 ~/lammps_Nov_2016/src/lmp_rcc < em2.in
