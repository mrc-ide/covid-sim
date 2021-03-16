#include <gtest/gtest.h>
#include "Error.h"

TEST(ErrorDeathTest, ErrorCritical) {
  ASSERT_DEATH({
    ErrorCritical("Test Error", "test-error.cpp", 7);
  }, "\\[Test Error line 7\\]");
}

TEST(ErrorDeathTest, err_critical_fmt) {
  ASSERT_DEATH({
    ERR_CRITICAL_FMT("Test ERR_CRITICAL_FMT", "test-error.cpp", 12);
  }, "\\[Test ERR_CRITICAL_FMT line 12\\]");
}

TEST(ErrorDeathTest, err_critical) {
  ASSERT_DEATH({
    ERR_CRITICAL("Test ERR_CRITICAL");
  }, "Test ERR_CRITICAL");
}
