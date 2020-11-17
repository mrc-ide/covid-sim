# Writing Google Test tests for CovidSim

The unit test framework chosen is called Google Test and is a fairly widely
used framework.

The basic structure for adding a test involves creating a test executable, and
adding it to the list of unit tests in the appropriate CMakeLists.txt

## Overview

### Initial Test

Start with a basic test file.  Somewhere in the `unit_tests` tree create a test
source file.  In this example we'll use `test-basic.cpp` as an example.

In the test file put the basic scaffolding needed to use GoogleTest:

```c++
#include <gtest/gtest.h>

TEST(CovidSimBasicTest, BasicTest)
{
    ASSERT_EQ(1, 1);
}
```

### Add to build system.

At the bottom of `./unit_tests/CMakeLists.txt` (or an appropriate CMakeLists.txt in a subdirectory) add a call to create the test:

```
add_unit_tests(TARGET test-basic SOURCES test-basic.cpp)
```

Now rerun cmake and run the tests.  From the command line, something like:

```sh
cd <build directory>
cmake . # Update config to know about new test executable
cmake --build . # Build CovidSim & Tests
ctest # Run all tests
```

See the (build documentation)[./build.md] for further details.

### Writing actual tests

Tests go in the test source file (`test-basic.cpp` in this example).

For full documentation on GoogleTest
(read their docs)[https://github.com/google/googletest/blob/v1.10.x/googletest/docs/primer.md].

Note that we use v0.10 of Google Test and the trunk development has moved on with changing some nomenculture, so check you are reading the correct docs version.

However, in summary tests look like:

```c++
TEST(<TESTSUITENAME>, <TEST>)
{
    // Do stuff - set code up, call functions, etc.
    ASSERT_EQ(X, Y); /* Checks X == Y, for everything except C-Style strings */
    ASSERT_STREQ(A, B); /* Checks that the C-Style strings A & B are the same. */
}
```

To test behaviours of particular CovidSim functions, include the appropriate header, and then add the source file to the SOURCES list for that test in CMakeLists.txt.

For example if testing functions in SetupModel.cpp, you would `#include "SetupModel.h"` in your sources and then have your CMakeLists.txt entry for the test look like:

```
add_unit_tests(TARGET test-basic SOURCES test-basic.cpp ${CMAKE_SOURCE_DIR}/src/SetupModel.cpp)
```

Note: The limitation here is that we cannot test functions in CovidSim.cpp because that has a `main()` function which will get the test framework confused.  The solution here is to move the functions under test out of CovidSim.cpp

Note: You have to specify a full path to the SetupModel.cpp code.

Note: That CMake uses `/` as a directory separator internally even on Windows.

### Debugging Unit Tests

The unit tests are stand alone programs that can be run by hand.  For example on UNIX.

```sh
cd build/unit-tests
./test-basic
```

The GoogleTest framework provides some extra command line options pass `--help` to see them.  Otherwise use your standard debugger to test the executable.
