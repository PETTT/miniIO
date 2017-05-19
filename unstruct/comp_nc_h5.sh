#!/bin/bash
#
# Quick, temporary utility script to compare a netcdf and an hdf5 file.
#

if [ "x$1" == "x" ]; then
    echo "Usage: $0 [file id]"
    exit 1
fi

NCFILE="unstruct.nc.d/$1.d/r.nc"
H5FILE="unstruct.hdf5.d/$1.d/r.h5"

echo "Dumping $NCFILE"
ncdump $NCFILE > tmpnc.txt

echo "Dumping $H5FILE"
ncdump $H5FILE > tmph5.txt

echo "Comparing"
diff tmpnc.txt tmph5.txt
echo ""
echo "Finished"
