find_package(PythonInterp)

if (${PYTHONINTERP_FOUND} AND ${PYTHON_VERSION_MAJOR} LESS 3)
    message (FATAL_ERROR "Python interpreter must be version 3 or higher, please supply a valid interpreter via the PYTHON_EXECUTABLE variable")
endif ()

if (${PYTHONINTERP_FOUND})
    add_test(
        NAME
            regressiontest_UK_100th
        COMMAND
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/regressiontest_UK_100th.py
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_LIST_DIR}
    )
    set_tests_properties (
        regressiontest_UK_100th
        PROPERTIES
            PASS_REGULAR_EXPRESSION "SUCCESS"
            FAIL_REGULAR_EXPRESSION "FAILURE;TEST FAILED"
    )
    add_test(
        NAME
            regressiontest_US_based
        COMMAND
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/regressiontest_US_based.py
        WORKING_DIRECTORY
            ${CMAKE_CURRENT_LIST_DIR}
    )
    set_tests_properties (
        regressiontest_US_based
        PROPERTIES
            PASS_REGULAR_EXPRESSION "SUCCESS"
            FAIL_REGULAR_EXPRESSION "FAILURE;TEST FAILED"
    )
endif ()