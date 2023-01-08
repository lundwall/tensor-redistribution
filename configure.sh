#!/bin/bash

# euler loading
env2lmod
module load gcc openmpi cmake

# build lsb if needed
if [ ! -z "$1" -a "$1" = "--get-lsb" ]; then
    git clone https://github.com/spcl/liblsb.git
    cd liblsb
    sed -i '480s/.*/    return 0;/' sync/hca_sync.cpp
    sed -i '29s/^/\/\//' sync/hca_sync.h
    sed -i '61s/^/\/\//' sync/hca_sync.h
    ./configure --prefix=$HOME/.local --enable-sync
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
