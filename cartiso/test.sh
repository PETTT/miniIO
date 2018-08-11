#!/bin/bash

# Script to demonstrate a range of cartiso features & cycle through all
# four time step modes

# First script parameter determines output type
# Both full and isosurface output types (don't do that for benchmarks)
if [ "$1" == vtk ]; then     # VTK XML Output options
    output="--pvti --pvtp"
elif [ "$1" == adios ]; then    # ADIOS options
    output="--adiosfull MPI --adiosiso MPI"
elif [ "$1" == nc ]; then     # NetCDF
    output="--nci --ncp"
elif [ "$1" == hdf5 ];then   # HDF5 with chunking parameters too
    output="--hdf5i --hdf5p --hdf5i_chunk 256 256 256 --hdf5p_chunk 100"
else
    echo "usage: test.sh <vtk or adios or hdf5 or nc>"
fi

tsk="2 2 2"    # 8 ranks with 2x2x2 decomposition
fq="2 2 2"     # 2x2x2 sinusoid frequencies
s=256          # Grid will be 256^3
sg=0.29        # Initial Gaussian sigma
sge=0.55       # Ending Gaussian sigma

# Phase 1: 20 time steps of sin2gauss mode
mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 20 --sin2gauss $output

# Phase 2: 80 time steps of gaussmove mode
mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 80 --tstart 20 --gaussmove $output

# Phase 3: 20 time steps of gaussresize mode
mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 20 --tstart 100 --gaussresize \
    --backward --sigmaend $sge $sge $sge $output

# Phase 4: 80 time steps of backward gaussmove mode
mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sge $sge $sge \
    --centertask --freq $fq --tsteps 80 --tstart 120 --gaussmove --backward $output

