# Add integration tests for CovidSim

add_test(
  NAME "inttest-uk-based"
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
          "--input" "${CMAKE_CURRENT_SOURCE_DIR}/uk-input"
          "--output" "${CMAKE_CURRENT_BINARY_DIR}/uk-output"
          "--covidsim" "$<TARGET_FILE:CovidSim>"
          "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/wpop_eur.txt.gz"
          "--r" "1.1")

add_test(
  NAME "inttest-us-based"
  COMMAND "${Python3_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/integration-test.py"
          "--input" "${CMAKE_CURRENT_SOURCE_DIR}/us-input"
          "--output" "${CMAKE_CURRENT_BINARY_DIR}/us-output"
          "--popfile" "${CMAKE_SOURCE_DIR}/data/populations/wpop_usacan.txt.gz"
          "--covidsim" "$<TARGET_FILE:CovidSim>")
