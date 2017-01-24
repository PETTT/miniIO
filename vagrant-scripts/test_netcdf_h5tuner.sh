#!/bin/bash

# Utility script to run netcdf-c tests using libautotuner.
#
# Currently a number of paths are hard coded.  This will
# change as I have a chance to parameterize.
#
# For now, assume that libautotuner.so and config.xml
# exist in $HOME

set -e

if [ -d netcdf-c ]; then
    echo "Using existing netcdf-c directory."
else
    git clone https://github.com/Unidata/netcdf-c
fi

cd netcdf-c
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=$(which mpicc)
make -j 4

LD_PRELOAD=$HOME/libautotuner.so H5TUNER_VERBOSE=4 H5TUNER_CONFIG_FILE=$HOME/config.xml make test ARGS="-V"
