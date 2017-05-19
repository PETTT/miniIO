#!/bin/bash

echo "Running HDF5 tests."
mpirun -n 8 ./unstruct --points 24 --hdf5
echo ""
echo "Runing netCDF tests."
mpirun -n 8 ./unstruct --points 24 --nc
echo ""
echo "Finished."
