# CMakeLists.txt for src directory

# CovidSim target
add_executable(CovidSim CovidSim.cpp CovidSim.h binio.cpp binio.h Rand.cpp Rand.h Constants.h Country.h MachineDefines.h Error.cpp Error.h Dist.cpp Dist.h Kernels.cpp Kernels.h Bitmap.cpp Bitmap.h Model.h Param.h SetupModel.cpp SetupModel.h SharedFuncs.h ModelMacros.h InfStat.h CalcInfSusc.cpp CalcInfSusc.h Sweep.cpp Sweep.h Update.cpp Update.h)
target_include_directories(CovidSim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if(USE_OPENMP)
  target_link_libraries(CovidSim PUBLIC OpenMP::OpenMP_CXX)
endif()
if(WIN32)
  target_link_libraries(CovidSim PUBLIC Gdiplus.lib Vfw32.lib)
  target_compile_definitions(CovidSim PUBLIC  "_CRT_SECURE_NO_WARNINGS")

  # Configure using BMP or PNG images.
  if(NOT USE_WIN32_BMP)
    add_compile_definitions("NO_WIN32_BM")
  endif()
elseif(UNIX)
  target_compile_definitions(CovidSim PUBLIC UNIX)
else()
  message(FATAL_ERROR "Unknown operating system type: ${CMAKE_SYSTEM_NAME}")
endif()
