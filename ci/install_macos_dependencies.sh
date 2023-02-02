#!/bin/sh -e

brew install llvm cmake libomp python3 coreutils
export CMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS -Xpreprocessor -fopenmp"
