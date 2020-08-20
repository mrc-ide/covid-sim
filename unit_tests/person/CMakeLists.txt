set(TEST_PERSON_SRC_FILES Person.cpp)

source_group(test\\person FILES ${TEST_PERSON_SRC_FILES})

add_unit_tests(${TEST_PERSON_SRC_FILES} person)
