#!/bin/bash

# Script for testing unstruct

PTS=512


NT=2

output="--hdf5"
mpirun -np 8 ./unstruct --pointspertask $PTS --tsteps $NT $output

output="--hdf5 --hdf5_chunk 18 --hdf5_compress gzip,9"
mpirun -np 8 ./unstruct --pointspertask $PTS --tsteps $NT $output


