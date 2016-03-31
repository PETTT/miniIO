#!/bin/bash

output="--pvtp"
#output="--pvti"
#output="--pvti --pvtp"

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.2 0.2 0.2 \
    --centertask --freq 8 8 8 --tsteps 20 --sin2gauss $output

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.2 0.2 0.2 \
    --centertask --freq 8 8 8 --tsteps 80 --tstart 20 --gaussmove $output

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.2 0.2 0.2 \
    --centertask --freq 8 8 8 --tsteps 20 --tstart 100 --gaussresize \
    --backward --sigmaend 0.4 0.4 0.4 $output

mpirun -np 4 ./cartiso --tasks 2 2 1 --size 256 256 256 --sigma 0.4 0.4 0.4 \
    --centertask --freq 8 8 8 --tsteps 80 --tstart 120 --gaussmove --backward $output

