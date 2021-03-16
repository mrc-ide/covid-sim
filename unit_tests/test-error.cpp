#include <gtest/gtest.h>
#include "Error.h"

TEST(Error, Error)
{
  ASSERT_DEATH({
    ErrorCritical("Test Error", "test-error.cpp", 7);
  }, "[Test Error line 7]");

  ASSERT_DEATH({
    ERR_CRITICAM_FTM("Test ERR_CRITICAL_FMT", "test-error.cpp", 11);
  }, "[Test ERR_CRITICAL_FMT line 11]");

  ASSERT_DEATH({
    ERR_CRITICAL("Test ERR_CRITICAL");
  }, "Test ERR_CRITICAL");
}
