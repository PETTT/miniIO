#!/bin/sh
# source me (sh/ksh/bash)
# For known systems, to set environment variables before running make

### All DSRC systems
if [ -n "$BC_HOST" ]; then

    export CC=mpicc
    export OPT="-O3 -xhost -ip"     # Assumes intel compiler

    ## All Cray's using the CrayPE compiler environment
    if [ -n "$CRAYPE_DIR" ]; then

        echo "Using DSRC Cray with CrayPE configuration"
        prgenv=`echo $LOADEDMODULES | grep -Po 'PrgEnv-.*?/'`
        echo "   Current environment: $prgenv"
        export CC=cc     # CrayPE compiler wrapper
        if [ "$prgenv" == "PrgEnv-intel/" ]; then
            export OPT="-O3 -xhost -ip"
        elif [ "$prgenv" == "PrgEnv-gnu/" ]; then
            export OPT="-O3"
        else
            export OPT="-O2"
        fi
        echo "   Set OPT=$OPT"

    ## Topaz/SGI
    elif [ "$BC_HOST" == "topaz" ]; then
        echo "Using DSRC Topaz configuration"
        export CC=icc    # Assume intel compiler for now
        export LIBS="-lmpi"

    ## Other DSRC
    else
        echo "Using default DSRC configuration"
    fi

### All POD systems (all seem to define the following variable)
elif [ -n "$POD_PBSSERVERS" ]; then

    echo "Using default POD configuration"
    export CC=mpicc
    export OPT="-O3 -xhost -ip"     # Assumes intel compiler

### Nothing else to try
else 
    echo "No special configuration detected"
    if [ -n "`which mpicc`" ]; then
        export CC=mpicc
        echo "   Found mpicc - setting CC=$CC"
    fi
fi

