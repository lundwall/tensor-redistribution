#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-socket=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000
#SBATCH --output=2d_s2.out

export OMP_NUM_THREADS=4
for twos in {0..4}
do
    for threes in {0..3}
    do
        for fives in {0..2}
        do
            srun --cpu-bind=cores ./build/manual/2dtransmit_manual $((24300000 / ((2 ** twos) * (3 ** threes) * (5 ** fives))))
        done
    done
done