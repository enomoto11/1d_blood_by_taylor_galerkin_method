#!/bin/sh
mkdir build
cd build
cmake -D enable_GLOG=OFF  -D EIGEN_DIR=/Users/t-enomoto/research/1d_blood/lib/eigen -D EIGEN_INCLUDE_DIR=/Users/t-enomoto/research/1d_blood/lib/eigen ..
make && make install
