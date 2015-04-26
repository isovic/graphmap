#! /bin/sh

# cd libdivsufsort-2.0.1
# make clean
# rm -rf build
# mkdir build
# cd build
# cmake ..
# cd ..
# ./configure
# mkdir build
# cd build
# cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="/usr/local" ..
# make

cd libdivsufsort-2.0.1
make clean
rm -rf build
rm -rf lib_install
mkdir lib_install
mkdir build
./configure --prefix=$PWD/../lib_install
cd build
cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF .. # -DCMAKE_INSTALL_PREFIX="$../lib_install" ..
make
