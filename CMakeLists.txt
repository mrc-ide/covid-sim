# CMakeLists.txt

# CMake setup
cmake_minimum_required (VERSION 3.8)

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

# Project initialisation
project("CovidSim")

option(USE_OPENMP "Compile with OpenMP parallelism enabled" ON)
option(ENABLE_REGRESSION_TEST "Enables regression testing targets" ON)

# Packages used
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)

# CovidSim target
add_subdirectory(src)
if(ENABLE_REGRESSION_TEST)
  enable_testing()
  add_subdirectory(tests)
endif()
