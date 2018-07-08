#!/bin/sh
# source me (sh/ksh/bash)
# For known systems, to set environment variables before running make

### All DSRC systems
if [ -n "$BC_HOST" ]; then

    export CC=mpicc

    ## All Cray's using the CrayPE compiler environment
    if [ -n "$CRAYPE_DIR" ]; then

        echo "Using DSRC Cray with CrayPE configuration"
        echo "   Current environment: `echo $LOADEDMODULES | grep -Po 'PrgEnv-.*?/'`"
        export CC=cc     # CrayPE compiler wrapper

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

### Nothing else to try
else 
    echo "No special configuration detected - doing nothing."
fi

