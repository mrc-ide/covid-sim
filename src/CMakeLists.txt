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

# Packages used
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)

# CovidSim target
add_executable(CovidSim CovidSim.cpp CovidSim.h binio.cpp binio.h Rand.cpp Rand.h Constants.h Country.h MachineDefines.h Error.cpp Error.h Dist.cpp Dist.h Kernels.cpp Kernels.h Bitmap.cpp Bitmap.h Model.h Param.h SetupModel.cpp SetupModel.h SharedFuncs.h ModelMacros.h InfStat.h CalcInfSusc.cpp CalcInfSusc.h Sweep.cpp Sweep.h Update.cpp Update.h)
target_include_directories(CovidSim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(USE_OPENMP)
  target_link_libraries(CovidSim PUBLIC OpenMP::OpenMP_CXX)
endif()
if(WIN32)
  target_link_libraries(CovidSim PUBLIC Gdiplus.lib Vfw32.lib)
  target_compile_definitions(CovidSim PUBLIC  "_CRT_SECURE_NO_WARNINGS")
elseif(UNIX)
  target_compile_definitions(CovidSim PUBLIC UNIX)
else()
  message(FATAL_ERROR "Unknown operating system type: ${CMAKE_SYSTEM_NAME}")
endif()
