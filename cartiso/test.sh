#!/bin/bash

if [ "$1" == adios ];then
output="--pvtp"
output="--pvti"
output="--pvti --pvtp"
output="--adiosfull MPI"
else
 if [ "$1" == hdf5 ];then
    output="--hdf5i --hdf5p --hdf5i_chunk 256 256 256 --hdf5p_chunk 100"
 else
    echo "usage: test.sh <adios or hdf5>"
 fi
fi

tsk="2 2 2"
fq="2 2 2"
s=256
sg=0.29
sge=0.55

mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 20 --sin2gauss $output

#mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
#    --centertask --freq $fq --tsteps 80 --tstart 20 --gaussmove $output
#
#mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
#    --centertask --freq $fq --tsteps 20 --tstart 100 --gaussresize \
#    --backward --sigmaend $sge $sge $sge $output
#
#mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sge $sge $sge \
#    --centertask --freq $fq --tsteps 80 --tstart 120 --gaussmove --backward $output

