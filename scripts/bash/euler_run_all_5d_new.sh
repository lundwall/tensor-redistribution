#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-socket=1
#SBATCH --nodes=2
#SBATCH --mem-per-cpu=16G
#SBATCH --output=out_mt_5d.txt
#SBATCH --constraint=EPYC_7763&ibfabric7
#SBATCH --time=5:00:00
export OMP_NUM_THREADS=4
srun --cpu-bind=cores ./build/measurement/5d_single
