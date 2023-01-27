#!/bin/sh

git submodule update --init
cd ./lib/kahip/KaHIP/
rm -rf installed build

mkdir installed
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../installed
make "$@"
make install
