#!/usr/bin/env sh

#rm out_*.txt
#rm lsb.*
timestamp=$(date +%Y%m%d-%H%M%S)
mkdir ${timestamp}
cd ${timestamp}

for num_dim_datatype in {2..5}
do
    sbatch --output=out_datatype_${num_dim_datatype}d.txt --mem-per-cpu=4G --ntasks=2 --ntasks-per-node=1 --wrap="mpirun ./../../build/datatype/${num_dim_datatype}dtransmit"
    sleep 1
    while squeue | grep -m 1 "wrap"; do sleep 1; done
done

for num_thread in {1..8}
do
    for num_dim in {2..5}
    do
    export OMP_NUM_THREADS=${num_thread}
    sbatch --output=out_mt_${num_dim}d_${num_thread}.txt --mem-per-cpu=4G --ntasks=`expr 2 \* ${num_thread}` --ntasks-per-node=${num_thread} --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n 2 --map-by node:PE=${num_thread} ./../../build/manual/${num_dim}dtransmit_manual"
    sleep 1
    while squeue | grep -m 1 "wrap"; do sleep 1; done
    done
done
