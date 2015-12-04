#!/bin/sh

## print commands and exit on error
set -e
set -x

## build
mkdir -p build
cd build
cmake ..
make -j 2

## test
ctest .
