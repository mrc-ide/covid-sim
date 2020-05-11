include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()
include(GoogleTest)

FetchContent_Declare(Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.3.7
)

FetchContent_GetProperties(Eigen3)
if(NOT eigen3_POPULATED)
  FetchContent_Populate(Eigen3)
  add_library(eigen INTERFACE)
  target_include_directories(eigen INTERFACE ${eigen3_SOURCE_DIR})
  add_library(Eigen3::Eigen ALIAS eigen)
  set(eigen3_SOURCE_DIR ${eigen3_SOURCE_DIR} CACHE PATH "Eigen3 path")
endif()
