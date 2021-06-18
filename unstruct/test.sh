#!/bin/bash

# Script for testing unstruct

PTS=512


NT=2

output="--hdf5"
#mpirun -np 8 ./unstruct --pointspertask $PTS --tsteps $NT $output

#n2=`grep -oP 'nelems2: \K\w+' out`
#n3=`grep -oP 'nelems3: \K\w+' out`

output="--hdf5 --hdf5_chunk 1 1 1 --hdf5_compress gzip,9"
mpirun -np 8 ./unstruct --pointspertask $PTS --tsteps $NT $output


