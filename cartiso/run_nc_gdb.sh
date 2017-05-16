#!/bin/bash
# Utility script to easily run a netCDF file.

set -e

mpirun --np 8 xterm -bg white -fg black -e gdb -x gdb.txt ./cartiso
#mpirun --np 4 xterm -e "gdb ./struct --tasks 2 2 --size 2 2 2 --tsteps 1 --nc"
#mpirun --np 4 ./struct --tasks 2 2 --size 16 12 2 --tsteps 1 --nc
