#!/bin/bash

NCORES=12
unamestr=`uname`

rm -rf deploy
rm -rf build
mkdir build
cd build 
cmake ../
make -j $NCORES
cd ..

mkdir deploy
cp ./build/heistream deploy/
rm -r build
