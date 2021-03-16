#include <gtest/gtest.h>
#include "Error.h"

TEST(ErrorDeathTest, ErrorCritical) {
  ASSERT_DEATH({
    ErrorCritical("Test Error", "test-error.cpp", 7);
  }, "test-error.cpp line \\d+ \\] Test Error");
}

TEST(ErrorDeathTest, err_critical_fmt) {
  ASSERT_DEATH({
    ERR_CRITICAL_FMT("Test %s\n", "TEST_STRING");
  }, "test-error.cpp line \\d+ \\] Test TEST_STRING");
}

TEST(ErrorDeathTest, err_critical) {
  ASSERT_DEATH({
    ERR_CRITICAL("Test ERR_CRITICAL");
  }, "Test ERR_CRITICAL");
}
