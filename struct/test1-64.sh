#!/bin/bash

## Struct test script for 1 - 64 ranks

# Set up work directory
w=/mnt/ssd/teststruct1-64.$BASHPID
mkdir -p $w
cp ./struct ./processlog.py $w
cd $w

# Parameter spaces to cover
ps=("1,1" "2,1" "2,2" "4,2" "4,4" "8,4" "11,4" "8,8")    # "i,j" ranks
outputs=( "--adios POSIX" "--adios MPI_LUSTRE" )    # output types
comps=( "" "--adios_transform zlib" "--adios_transform szip" "--adios_transform zfp:accuracy=0.0001" "--adios_transform sz:absolute=0.0001" )
balances=( "" "--balance" )      # to balance or not to balance
cleanup=1     # Remove data directories after each run: 0 or 1
progress=1    # Print progress: 0 or 1
mpi=mpirun      # MPI launcher
ni=128x         # i points (with x: per i rank)
nj=128x         # j points (with x: per j rank)
nk=1024            # k levels stays constant
noiseij=10x       # frequency of spatial noise (with x: per i,j rank)
noisek=10         # frequency of spatial noise k
ts=11     # time steps

for p in "${ps[@]}"; do
    IFS="," read -a pp <<< "${p}"    # Convert i,j ranks to array
    pi=${pp[0]}     # i ranks
    pj=${pp[1]}     # j ranks
    nc=$(($pi * $pj))    # total ranks
    tsk="$pi $pj"    # to pass to command line
    i=0

    for output in "${outputs[@]}"; do
      for comp in "${comps[@]}"; do
        for balance in "${balances[@]}"; do
            name=$(printf "test.%02d.%02d" $nc $i)
            log=../${name}.log
            dir=${name}.d
            mkdir $dir
            cd $dir
            echo "tsk=$tsk ni=$ni nj=$nj nk=$nk ts=$ts output=$output $comp $balance hints=$MPICH_MPIIO_HINTS" > $log
            echo "$mpi -n $nc ../struct --tasks $tsk --size $ni $nj $nk --tsteps $ts --noisespacefreq $noiseij $noiseij $noisek $output $comp $balance" >> $log
	    if [ $progress -ne 0 ]; then echo "$tsk $output $comp $balance ..."; fi

	    # Real work happens here
            $mpi -n $nc ../struct --tasks $tsk --size $ni $nj $nk --tsteps $ts --noisespacefreq $noiseij $noiseij $noisek $output $comp $balance 2>&1 >> $log

            ls -l >> $log
            if [ $(ls -1d *.dir 2>/dev/null | wc -l) -gt 0 ]; then  # list ADIOS .dirs if exist
                ls -l *.dir >> $log
            fi

	    if [ $progress -ne 0 ]; then ../processlog.py -a $log; fi

            cd ..
            if [ $cleanup -ne 0 ]; then
                rm -r $dir
            fi
            ((i++))
        done
      done
    done
done

