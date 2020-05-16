set(TEST_NODDY_SRC_FILES Holder.cpp Noddy.cpp)

source_group(test\\noddy FILES ${TEST_SRC_FILES})

FIND_PACKAGE(Threads REQUIRED)

add_executable(test1 Noddy.cpp)
target_link_libraries(test1 gtest_main ${CMAKE_THREAD_LIBS_INIT})
add_test(test1 test1)
add_executable(test2 Holder.cpp)
target_link_libraries(test2 gtest_main ${CMAKE_THREAD_LIBS_INIT})
add_test(test2 test2)

add_custom_target(test2_run
  COMMAND test2
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_dependencies(test1 test2_run)
