include_directories("${gtest_SOURCE_DIR}/include" "${eigen3_SOURCE_DIR}")

enable_testing()

ADD_SUBDIRECTORY(noddy)
ADD_SUBDIRECTORY( geometry )
