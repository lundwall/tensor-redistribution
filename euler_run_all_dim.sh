#!/usr/bin/env sh

for num_thread in {1..8}
do
    for num_dim in {2..5}
    do
    export OMP_NUM_THREADS=${num_thread}
    sbatch --output=out_mt_${num_dim}d_${num_thread}.txt --mem-per-cpu=2G --ntasks=2*${num_thread} --ntasks-per-node=1 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=${num_thread} ./build/manual/${num_dim}dtransmit_manual"
    done
done