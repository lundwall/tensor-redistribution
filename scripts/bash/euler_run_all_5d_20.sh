#!/usr/bin/env sh

export OMP_NUM_THREADS=1
sbatch --output=out_mt_5d_1t.txt --mem-per-cpu=2G --ntasks=2 --ntasks-per-node=1 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=1 ../../build/measurement/5d_single_20"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

#: << 'COMMENT'
export OMP_NUM_THREADS=2
sbatch --output=out_mt_5d_2t.txt --mem-per-cpu=2G --ntasks=4 --ntasks-per-node=2 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=2 ../../build/measurement/5d_single_20"
sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

export OMP_NUM_THREADS=3
sbatch --output=out_mt_5d_3t.txt --mem-per-cpu=2G --ntasks=6 --ntasks-per-node=3 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=3 ../../build/measurement/5d_single_20"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

#: << 'COMMENT'
export OMP_NUM_THREADS=4
sbatch --output=out_mt_5d_4t.txt --mem-per-cpu=2G --ntasks=8 --ntasks-per-node=4 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=4 ../../build/measurement/5d_single_20"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

export OMP_NUM_THREADS=5
sbatch --output=out_mt_5d_5t.txt --mem-per-cpu=2G --ntasks=10 --ntasks-per-node=5 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=5 ../../build/measurement/5d_single_20"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

export OMP_NUM_THREADS=6
sbatch --output=out_mt_5d_6t.txt --mem-per-cpu=2G --ntasks=12 --ntasks-per-node=6 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=6 ../../build/measurement/5d_single_20"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

export OMP_NUM_THREADS=7
sbatch --output=out_mt_5d_7t.txt --mem-per-cpu=2G --ntasks=14 --ntasks-per-node=7 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=7 ../../build/measurement/5d_single_20"

sleep 1
while squeue | grep -m 1 "wrap"; do sleep 1 ; done

export OMP_NUM_THREADS=8
sbatch --output=out_mt_5d_8t.txt --mem-per-cpu=2G --ntasks=16 --ntasks-per-node=8 --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=8 ../../build/measurement/5d_single_20"
#COMMENT
