#!/bin/bash
# This script builds the executable using cmake
# Usage: ./build.sh

mkdir renders
mkdir build

cd build
cmake ..
cmake --build .
