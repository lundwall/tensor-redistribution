# DPHPC

## Environment set up
The package and their versions we are using:

gcc/8.2.0 openmpi/4.0.2 cmake/3.20.3


## How to build

First make sure the working directory is at dphpc/ directory (root directory of this project).
Then

```
mkdir build
cd build
cmake ..
make
```
We use the files ```5d_single_10```, ```5d_single_20```, ```5d_single_30``` to run the 5d single block tramsmission experiments. Each would
run with sub block size = 10**5 ints, 20 ** 5 ints, 30 ** 5 ints, respectively.

```scripts/bash/euler_redist.sh``` performs all the redistribution experiments on euler, and ```scripts/python/plot_redist.sh``` does the plots in the report.

## How to run
We run the experiments on Euler. For each experiment we use command like this:

```
export OMP_NUM_THREADS=1
sbatch --output=[output_log_file] --mem-per-cpu=[mem_size] --ntasks=[thread_number] --ntasks-per-node=[thread_num_per_node] --wrap="unset LSB_AFFINITY_HOSTFILE; mpirun -n [mpi_process_number] --map-by node:PE=[thread_per_mpi_process] [execution_file]"
```

We have a script to handle settings of different threads. dphpc/scripts/bash/euler_run_all_5d_[10/20/30].sh would run the 5d single transmission experiment with specified sub block size for thread = 1 to thread = 8, and the lsb files would be generated.
