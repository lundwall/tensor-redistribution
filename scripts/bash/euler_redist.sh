#!/bin/bash

#SBATCH --ntasks=12
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=12
#SBATCH --mem-per-cpu=25G
#SBATCH --constraint=EPYC_7H12
#SBATCH --output=out_redist.txt

export OMP_NUM_THREADS=8
srun --cpu-bind=verbose,cores /cluster/home/mlundwall/dphpc/build/redistribute/redistribute