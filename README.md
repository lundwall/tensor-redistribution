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

## How to run
We use the files ```5d_single.cpp``` to run the 5d single block tramsmission experiments. 
```scripts/bash/euler_run_all_5d_new.sh``` helps run the experiments on euler, and ```scripts/python/plot_single_block_transmission.py``` does the plots.

```scripts/bash/euler_redist.sh``` performs all the redistribution experiments on euler, and ```scripts/python/plot_redist.sh``` does the plots in the report.
