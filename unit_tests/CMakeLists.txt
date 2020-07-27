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

# from list of files we'll create tests test_name.cpp -> test_name
function(add_unit_tests TEST_SRC_FILES TEST_LIBRARIES)
  foreach(_test_file ${TEST_SRC_FILES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    target_link_libraries(${_test_name} ${TEST_LIBRARIES} gtest_main ${CMAKE_THREAD_LIBS_INIT})
    gtest_discover_tests(${_test_name})
  endforeach()
endfunction()

ADD_SUBDIRECTORY(geometry)
