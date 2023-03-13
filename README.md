# Accelerating Redistribution of Tensors

## Abstract

The performance of high-dimensional multilinear algebra kernels on massively parallel system is often dominated by data movement. Recent studies on improving the kernel execution performance introduced a data movement-optimal distribution schedule. After the computation of intermediaries, the distribution often changes, making redistribution of tensors necessary. This work focuses on improving the performance of tensor redistribution through optimising data communication. 

In this project, we propose and compare various methods of accelerating the redistribution of multi-dimensional tensors with MPI. The methods experimented with different MPI communication modes, acceleration with multithreading, and pipelining. In addition, a C++ library for point-to-point tensor transmission and redistribution was implemented. When compared to MPI’s inbuilt custom datatypes, our proposed method achieved an improvement of up to 1.98x.

## Environment set up
We are using the following package versions:

gcc/8.2.0 openmpi/4.0.2 cmake/3.20.3

On Euler, simply run:

```
./configure.sh
```


## How to build

First, make sure the working directory is at the root of this project.
Then:

```
mkdir build
cd build
cmake ..
make
```

## How to run
We use ```5d_single.cpp``` to run the 5D single block tramsmission experiments. 
```scripts/bash/euler_run_all_5d_new.sh``` helps run the experiments on euler, and ```scripts/python/plot_single_block_transmission.py``` does the plots.

```scripts/bash/euler_redist.sh``` performs all the redistribution experiments on euler, and ```scripts/python/plot_redist.sh``` does the plots in the report.
