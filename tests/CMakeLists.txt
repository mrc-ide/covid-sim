# Add integration tests for CovidSim

add_test(
  NAME "inttest-uk-based-1thread"
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
          "--input" "${CMAKE_CURRENT_SOURCE_DIR}/uk-input"
          "--output" "${CMAKE_CURRENT_BINARY_DIR}/uk-j1-output"
          "--checksums" "${CMAKE_CURRENT_SOURCE_DIR}/uk-input/results-j1.cksum"
          "--covidsim" "$<TARGET_FILE:CovidSim>"
          "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/wpop_eur.txt.gz"
          "--r" "1.1"
          "-j" "1")

add_test(
  NAME "inttest-uk-based-2threads"
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
          "--input" "${CMAKE_CURRENT_SOURCE_DIR}/uk-input"
          "--output" "${CMAKE_CURRENT_BINARY_DIR}/uk-j2-output"
          "--checksums" "${CMAKE_CURRENT_SOURCE_DIR}/uk-input/results-j2.cksum"
          "--covidsim" "$<TARGET_FILE:CovidSim>"
          "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/wpop_eur.txt.gz"
          "--r" "1.1"
          "-j" "2")
set_property(TEST "inttest-uk-based-2threads" PROPERTY PROCESSORS 2)

add_test(
  NAME "inttest-us-based-1thread"
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
          "--input" "${CMAKE_CURRENT_SOURCE_DIR}/us-input"
          "--output" "${CMAKE_CURRENT_BINARY_DIR}/us-j1-output"
          "--checksums" "${CMAKE_CURRENT_SOURCE_DIR}/us-input/results-j1.cksum"
          "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/wpop_usacan.txt.gz"
          "--covidsim" "$<TARGET_FILE:CovidSim>"
          "-j" "1")

add_test(
  NAME "inttest-us-based-2threads"
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
          "--input" "${CMAKE_CURRENT_SOURCE_DIR}/us-input"
          "--output" "${CMAKE_CURRENT_BINARY_DIR}/us-j2-output"
          "--checksums" "${CMAKE_CURRENT_SOURCE_DIR}/us-input/results-j2.cksum"
          "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/wpop_usacan.txt.gz"
          "--covidsim" "$<TARGET_FILE:CovidSim>"
          "-j" "2")
set_property(TEST "inttest-us-based-2threads" PROPERTY PROCESSORS 2)
