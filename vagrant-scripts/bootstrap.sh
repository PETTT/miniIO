#!/bin/bash

set -u

###
# System stuff.
###

apt-get update
apt-get upgrade -y

apt-get -y install emacs ubuntu-dev-tools mpich libmpich-dev libcurl4-openssl-dev autoconf libtool libmxml-dev libtool cmake

###
# Update system language.
###
/usr/sbin/locale-gen
echo 'LANG="en_US"' > /etc/default/locale
echo 'LANGUAGE="en_US:"' >> /etc/default/locale
/usr/bin/localectl set-locale LANG="en_US.UTF-8"

###
# Run dependency, h5tuner installation scripts.
###

cp /vagrant/vagrant-scripts/install_hdf5tune*.sh /home/vagrant/Desktop
cd /home/vagrant/Desktop

su -c '/home/vagrant/Desktop/install_hdf5tune_deps.sh' -s /bin/sh vagrant
su -c '/home/vagrant/Desktop/install_hdf5tune.sh' -s /bin/sh vagrant

git clone http://github.com/Unidata/netcdf-c /home/vagrant/Desktop/netcdf-c

cd /home/vagrant/Desktop/netcdf-c
git checkout h5tuner
cp AUTOTUNER_README.md /home/vagrant/Desktop


chown -R vagrant:vagrant /home/vagrant
