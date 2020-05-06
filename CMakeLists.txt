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

# Do we use the Win32 BMP routines or the generic ones?
option(USE_WIN32_BMP
  "Use Windows BMP specific output routines instead of generic ones.  If ON disables testing."
  OFF)

# Packages used
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Python3 needed for testing
if (CMAKE_VERSION VERSION_LESS 3.12)
  set(Python_ADDITIONAL_VERSIONS 3 3.9 3.8 3.7 3.6)
  find_package(PythonInterp REQUIRED)
  set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}")
else()
  find_package(Python3 REQUIRED)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)

add_subdirectory(src)

if(WIN32 AND USE_WIN32_BMP)
  message(WARNING "USE_WIN32_BMP is enabled.  Therefore testing is disabled.")
  message(WARNING "To enable testing invoke with `cmake -DUSE_WIN32_BMP=NO`")
else()
enable_testing()
add_subdirectory(tests)
endif()
