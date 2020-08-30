cmake_policy(SET CMP0079 NEW)

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

include_directories("${gtest_SOURCE_DIR}/include")

FIND_PACKAGE(Threads REQUIRED)

function(add_unit_tests)
  # Parse options
  set(_options "")
  set(_one_value TARGET)
  set(_multi_value SOURCES LIBRARIES)

  cmake_parse_arguments(_aut "${_options}" "${_one_value}" "${_multi_value}" ${ARGN})

  # Check parameters
  if(NOT _aut_TARGET)
    message(FATAL "No target specified to add_unit_tests")
  endif()
  if(NOT _aut_SOURCES)
    message(FATAL "Must specify sources in add_unit_tests for test ${_aut_TARGET}")
  endif()

  # Create target and wire it up to ctest
  add_executable(${_aut_TARGET} ${_aut_SOURCES})
  target_link_libraries(${_aut_TARGET} ${_aut_LIBRARIES} gtest_main ${CMAKE_THREAD_LIBS_INIT})
  target_include_directories(${_aut_TARGET} PUBLIC "${CMAKE_SOURCE_DIR}/src")
  gtest_discover_tests(${_aut_TARGET})
endfunction()

ADD_SUBDIRECTORY(geometry)
