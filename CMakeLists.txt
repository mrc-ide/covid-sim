# CMakeLists.txt

# CMake setup
cmake_minimum_required (VERSION 3.8)

# Project initialisation
project("CovidSim")

# Work around some policy behaviours
if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif()

# Set a default build type if none was specified
set(default_build_type "RelWithDebInfo")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(USE_OPENMP "Compile with OpenMP parallelism enabled" ON)

# Country configuration for compile-time properties
set(COUNTRY "UK" CACHE STRING "The country being modeled.")
set_property(CACHE COUNTRY PROPERTY STRINGS UK US)

# Packages used
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Python3 needed for testing
find_package(Python3 REQUIRED)

enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(src)
add_subdirectory(tests)
