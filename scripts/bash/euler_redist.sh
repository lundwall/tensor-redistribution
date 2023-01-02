#!/bin/bash
#SBATCH --ntasks=12
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=12
#SBATCH --constraint=EPYC_7763&ibfabric7
#SBATCH --output=out_redist.txt
#SBATCH --time=24:00:00

export OMP_NUM_THREADS=8
srun --cpu-bind=cores /cluster/home/mlundwall/dphpc/build/redistribute/redistribute