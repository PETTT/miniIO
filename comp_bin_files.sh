#!/bin/bash
#
# Quick, temporary utility script to compare a netcdf and an hdf5 file.
#

if [ "x$1" == "x" ]; then
    echo "Usage: $0 [identifer]"
    exit 1
fi

NCFILE=$1.nc
HFILE=$1.h5

echo "Dumping $1"
ncdump $NCFILE > tmpnc.txt

echo "Dumping $2"
ncdump $HFILE > tmphdf5.txt

echo "Comparing"
diff tmpnc.txt tmphdf5.txt

echo ""
echo "Finished"
