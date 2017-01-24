#!/bin/bash
#
# Roughly follows:
#
# https://github.com/HDFGroup/H5Tuner/wiki/Example

set -e

echo "Installing a few things via package manager, just to be sure."
sudo apt-get install -y libmxml-dev

echo "Installing zlib"
wget http://zlib.net/zlib-1.2.8.tar.gz
tar xvzf zlib-1.2.8.tar.gz
pushd zlib-1.2.8
./configure --prefix=/usr
make -j 4 && sudo make install
popd


echo "Installing hdf5"
wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.17/src/hdf5-1.8.17.tar.gz
tar xvzf hdf5-1.8.17.tar.gz
pushd hdf5-1.8.17
CC=mpicc FC=mpif90 CXX=mpic++ ./configure --prefix=/usr --enable-parallel --with-zlib=/usr --enable-hl --enable-shared --disable-static --enable-using-memchecker --enable-debug=all --enable-codestack
make -j 4 && sudo make install
popd
