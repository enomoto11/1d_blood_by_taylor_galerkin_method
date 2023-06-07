#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
      -D enable_GLOG=ON \
      ..

      # -D CMAKE_INSTALL_PREFIX=/home/totani/bin/BDIM-cpp \
      # -D TP_DIR=/home/totani/lib/TextParser-1.8.5 \
      # -D enable_GLOG=ON \
      # -D GLOG_DIR=/home/totani/lib/glog \
make && make install