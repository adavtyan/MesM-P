#!/bin/bash
#SBATCH --job-name=em2_ibar
#SBATCH --output=bd_example_sbatch.out
#SBATCH --error=bd_example_sbatch.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --constraint=ib
#SBATCH --partition=gavoth-ivb
#SBATCH --qos=gavoth

##gavoth-sdb
##gavoth-ivb
##gavoth-brw

let NCORE=SLURM_NNODES*SLURM_CPUS_ON_NODE

echo "Total number of core:", $NCORE

mpirun -np $NCORE ~/lammps_Nov_2016/src/lmp_rcc < em2.in
