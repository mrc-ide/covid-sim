#include <gtest/gtest.h>
#include "Error.h"

TEST(Error, Error)
{
  ASSERT_DEATH({
    ErrorCritical("Test Error", "test-error.cpp", 7)
  }, "[Test Error line 7]");
}

