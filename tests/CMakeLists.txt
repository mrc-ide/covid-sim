# Add integration tests for CovidSim

# The "test-accept" target runs all tests and accepts their output.
add_custom_target("test-accept")

function(add_integration_test TESTNAME INPUT OUTPUT_ROOT POPFILE R)
  foreach(j RANGE 1 2)
    set(test_name "${TESTNAME}-j${j}")
    set(accept_name "${test_name}-accept")
    add_test(NAME "${test_name}"
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
              "--input" "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}"
              "--output" "${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_ROOT}-j${j}"
              "--checksums" "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}/results-j${j}.cksum"
              "--covidsim" "$<TARGET_FILE:CovidSim>"
              "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/${POPFILE}"
              "--r" "${R}"
              "-j" "${j}")
    set_property(TEST "${test_name}" PROPERTY PROCESSORS ${j})
    add_custom_target("${accept_name}"
      DEPENDS CovidSim
      COMMAND "${Python3_EXECUTABLE}"
              "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
              "--input" "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}"
              "--output" "${CMAKE_CURRENT_BINARY_DIR}/${OUTPUT_ROOT}-j${j}-accept"
              "--checksums" "${CMAKE_CURRENT_SOURCE_DIR}/${INPUT}/results-j${j}.cksum"
              "--covidsim" "$<TARGET_FILE:CovidSim>"
              "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/${POPFILE}"
              "--r" "${R}"
              "-j" "${j}"
              "--accept")
    add_dependencies(test-accept "${accept_name}")
  endforeach()
endfunction()

add_integration_test("inttest-uk-based" "uk-input" "uk-output" "wpop_eur.txt.gz" "1.1")
add_integration_test("inttest-us-based" "us-input" "us-output" "wpop_usacan.txt.gz" "1.5")
