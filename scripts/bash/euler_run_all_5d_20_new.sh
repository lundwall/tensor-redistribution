#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-socket=1
#SBATCH --nodes=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=out_mt_5d.txt
#SBATCH --constraint=EPYC_7H12
export OMP_NUM_THREADS=8
srun --cpu-bind=cores ./build/measurement/5d_single_20