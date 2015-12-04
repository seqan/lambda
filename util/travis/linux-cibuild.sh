#!/bin/sh

## build
cd build
cmake ..
make -j 2

## test
ctest .
