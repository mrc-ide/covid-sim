#!/bin/sh -e

brew install llvm cmake libomp python3 coreutils
export CMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS -Xpreprocessor -fopenmp"
export OpenMP_C_LIB_NAMES="libomp"
export OpenMP_C_FLAGS="-fopenmp"
