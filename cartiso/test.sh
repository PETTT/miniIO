#!/bin/bash

output="--pvtp"
#output="--pvti"
#output="--pvti --pvtp"

tsk="2 2 2"
fq="2 2 2"
s=256
sg=0.29
sge=0.55

mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 20 --sin2gauss $output

mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 80 --tstart 20 --gaussmove $output

mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sg $sg $sg \
    --centertask --freq $fq --tsteps 20 --tstart 100 --gaussresize \
    --backward --sigmaend $sge $sge $sge $output

mpirun -np 8 ./cartiso --tasks $tsk --size $s $s $s --sigma $sge $sge $sge \
    --centertask --freq $fq --tsteps 80 --tstart 120 --gaussmove --backward $output

