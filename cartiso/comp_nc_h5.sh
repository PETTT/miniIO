#!/bin/bash
#
# Quick, temporary utility script to compare a netcdf and an hdf5 file.
#

if [ "x$1" == "x" ]; then
    echo "Usage: $0 [file 1] [file 2]"
    exit 1
fi

if [ "x$2" == "x" ]; then
    echo "Usage: $0 [file 1] [file 2]"
    exit 1
fi

echo "Dumping $1"
ncdump $1 > tmpa.txt

echo "Dumping $2"
ncdump $2 > tmpb.txt

echo "Comparing"
diff tmpa.txt tmpb.txt

echo ""
echo "Finished"
