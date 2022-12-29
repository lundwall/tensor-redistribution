#!/usr/bin/env sh

sbatch --output=out_mt_5d.txt --mem-per-cpu=8G --time=5:00:00 --ntasks=16 --ntasks-per-node=8 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=8 ../../build/measurement/5d_single_30"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

