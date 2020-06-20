include_directories(../../src/Geometry)

set(TEST_GEOMETRY_SRC_FILES BoundingBox.cpp)

source_group(test\\geometry FILES ${TEST_GEOMETRY_SRC_FILES})

FIND_PACKAGE(Threads REQUIRED)

# from list of files we'll create tests test_name.cpp -> test_name
foreach(_test_file ${TEST_GEOMETRY_SRC_FILES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    target_link_libraries(${_test_name} geometrylib gtest_main ${CMAKE_THREAD_LIBS_INIT})
    gtest_discover_tests(${_test_name})
endforeach()
