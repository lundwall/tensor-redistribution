#!/bin/bash

# euler loading
module load gcc openmpi cmake

# build lsb if needed
if [ ! -z "$1" -a "$1" = "--get-lsb" ]; then
    git clone https://github.com/spcl/liblsb.git
    cd liblsb
    ./configure --prefix=$HOME/.local
    make -j2
    make install
fi

# # path variables
export LD_LIBRARY_PATH=$HOME/.local/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$HOME/.local/lib:$LIBRARY_PATH
export CPATH=$HOME/.local/include:$CPATH

# make
rm -rf build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RELEASE
make -j2
