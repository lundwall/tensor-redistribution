#!/bin/bash

#SBATCH --ntasks=12
#SBATCH --ntasks-per-socket=12
#SBATCH --nodes=1
#SBATCH --array=1-8

for mode in manual datatype; do
    if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
        #SBATCH --cpus-per-task=1
        echo 1
    elif [ $SLURM_ARRAY_TASK_ID -eq 2 ]; then
        #SBATCH --cpus-per-task=2
        echo 2
    elif [ $SLURM_ARRAY_TASK_ID -eq 3 ]; then
        #SBATCH --cpus-per-task=3
        echo 3
    elif [ $SLURM_ARRAY_TASK_ID -eq 4 ]; then
        #SBATCH --cpus-per-task=4
        echo 4
    elif [ $SLURM_ARRAY_TASK_ID -eq 5 ]; then
        #SBATCH --cpus-per-task=5
        echo 5
    elif [ $SLURM_ARRAY_TASK_ID -eq 6 ]; then
        #SBATCH --cpus-per-task=6
        echo 6
    elif [ $SLURM_ARRAY_TASK_ID -eq 7 ]; then
        #SBATCH --cpus-per-task=7
        echo 7
    elif [ $SLURM_ARRAY_TASK_ID -eq 8 ]; then
        #SBATCH --cpus-per-task=8
        echo 8
    fi
    export OMP_NUM_THREADS=$SLURM_ARRAY_TASK_ID
    srun --cpu-bind=cores /cluster/home/mlundwall/dphpc/build/redistribute/redistribute "$mode"
done