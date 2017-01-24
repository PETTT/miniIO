#!/bin/bash
#
# Utility script for downloading and installing h5tuner

set -e

echo "==================="
echo "Cloning from HDFGroup/H5Tuner"
echo "==================="

if [ $# -lt 1 ]; then
	echo "You must specify a target directory to copy the config and autotuner libraries to."
	exit 0
fi

rm -rf H5Tuner

git clone https://bitbucket.hdfgroup.org/scm/engility/h5tuner.git
cd h5tuner

autoreconf -if

CC=`which mpicc` ./configure && make -j 4
cp src/libautotuner.so $1
cp examples/config.xml $1
