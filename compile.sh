#!/bin/bash

NCORES=12
unamestr=`uname`
ENABLE_TIME_MEASUREMENTS=${1:-OFF}

rm -rf deploy
mkdir -p build
cd build
cmake .. -DENABLE_TIME_MEASUREMENTS=${ENABLE_TIME_MEASUREMENTS}
make -j $NCORES
cd ..

mkdir deploy
cp ./build/heistream deploy/
cp ./build/heistream_edge deploy/
#rm -r build
