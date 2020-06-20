# CMakeLists.txt for src directory

# Set up the IDE
set(MAIN_SRC_FILES CovidSim.cpp BinIO.cpp Rand.cpp Error.cpp Dist.cpp
  Kernels.cpp Bitmap.cpp SetupModel.cpp CalcInfSusc.cpp Sweep.cpp Update.cpp
  Param.cpp Direction.cpp InverseCdf.cpp Memory.cpp CLI.cpp)
set(MAIN_HDR_FILES CovidSim.h BinIO.h Rand.h Constants.h Country.h Error.h
  Dist.h Kernels.h Bitmap.h Model.h Param.h SetupModel.h ModelMacros.h
  InfStat.h CalcInfSusc.h Sweep.h Update.h MicroCellPosition.hpp Direction.hpp
  InverseCdf.h Memory.h CLI.h)
source_group(covidsim\\main FILES ${MAIN_SRC_FILES} ${MAIN_HDR_FILES})

# CovidSim target
add_executable(CovidSim ${MAIN_SRC_FILES} ${MAIN_HDR_FILES})
target_include_directories(CovidSim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Geometry)
add_subdirectory(Models)

target_link_libraries(CovidSim PUBLIC geometrylib)
if(USE_OPENMP)
  target_link_libraries(CovidSim PUBLIC OpenMP::OpenMP_CXX)
endif()
if(WIN32)
  target_link_libraries(CovidSim PUBLIC Gdiplus.lib Vfw32.lib)
  target_compile_definitions(CovidSim PUBLIC  "_CRT_SECURE_NO_WARNINGS")
  add_compile_options("fp:strict")
elseif(UNIX)
  target_compile_definitions(CovidSim PUBLIC UNIX)
else()
  message(FATAL_ERROR "Unknown operating system type: ${CMAKE_SYSTEM_NAME}")
endif()
