#!/bin/bash

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.2 0.2 0.2 \
    --centertask --freq 8 8 8 --tsteps 20 --sin2gauss

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.2 0.2 0.2 \
    --centertask --freq 8 8 8 --tsteps 80 --tstart 20 --gaussmove

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.2 0.2 0.2 \
    --centertask --freq 8 8 8 --tsteps 20 --tstart 100 --gaussresize \
    --backward --sigmaend 0.4 0.4 0.4

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.4 0.4 0.4 \
    --centertask --freq 8 8 8 --tsteps 80 --tstart 120 --gaussmove --backward

