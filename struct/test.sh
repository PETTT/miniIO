#!/bin/bash

# Script for testing struct

INP=2
JNP=4

NI=64
NJ=128
NK=256

NT=2

CHUNK_Y=$(($NJ/$JNP))
CHUNK_Z=32

output="--hdf5"
mpirun -np 8 ./struct --tasks $INP $JNP --size $NI $NJ $NK --tsteps $NT $output

output="--hdf5 --hdf5_chunk $CHUNK_Y $CHUNK_Z --hdf5_compress gzip,9"
mpirun -np 8 ./struct --tasks $INP $JNP --size $NI $NJ $NK --tsteps $NT $output


