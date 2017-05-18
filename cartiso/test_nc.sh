#!/bin/bash

set -e

output="--nci --ncp"
#output="--nci"
#output="--ncp"

mpirun -np 8 ./cartiso --tasks 2 2 2 --size 256 256 256 --sigma 0.29 0.29 0.29 --centertask --freq 2 2 2 --tsteps 20 --sin2gauss $output
