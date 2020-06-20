set(TEST_EXAMPLE_SRC_FILES Example.cpp)

source_group(test\\example FILES ${TEST_EXAMPLE_SRC_FILES})

add_unit_tests(${TEST_EXAMPLE_SRC_FILES} "")
