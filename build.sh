#!/bin/bash
# This script builds the executable using cmake
# Usage: ./build.sh

mkdir renders

cmake .
cmake --build .
