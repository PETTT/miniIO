#!/bin/bash

echo "Running HDF5 tests."
mpirun -np 8 ./amr --tasks 2 2 2 --size 9 9 9 --levels 3 --tsteps 1 --hdf5
echo ""
echo "Runing netCDF tests."
mpirun -np 8 ./amr --tasks 2 2 2 --size 9 9 9 --levels 3 --tsteps 1 --nc
echo ""
echo "Finished."
