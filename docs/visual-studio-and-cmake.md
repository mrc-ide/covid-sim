# Using Visual Studio and CMake

## Introduction

This document walks through using the [CMake](https://www.cmake.org) based
build-system within Visual Studio to build & develop CovidSim

## Requirements

This document assumes you are using Visual Stduio 2019 version 16.6 (release
20 May 2020).

You should also have installed
[Python 3.8](https://www.microsoft.com/store/productId/9MSSZTT1N39L).

It also assumes you have a clean checkout of mrc-ide/covid-sim.

## Basic Usage.

Open the project in Visual Studio by choosing the directory or the
[`CMakeLists.txt`](../CMakeLists.txt) file in the root directory.

There are default configurations provided:

 * "Debug" - builds debug build

 * "Release" - builds a release build (with Debug Information included)

These can be chosen from the "Configuration" drop-down list in the middle
of the Standard toolbar.

To build the project go to the Build menu, and choose "Build All".

If the "Build All" option is not available go to the "Project" menu and choose
"Rescan Solution".

## Running Tests

The integration tests can be run from the "Test" menu by choosing the
"Run CTests for CovidSim" menu option.

Note that these don't display much information whilst running and take about
30 minutes to run on my machine.

## Building from the Command Line

:note: The most important piece of information in this section is do not do
manual builds within the source checkout.  Doing so will confuse Visual Studio.

The easiest way to build from the command-line is to go to the
"Solution Explorer - Folder View" - by default on the left of the screen.
Right click on the top level "covid-sim" folder and choose "Open Developer
Command Prompt".

Change to the parent directory: `cd ..`

Make a build directory: `md covid-sim-build`

Change to that directory: `cd covid-sim-build`

Configure cmake: `cmake ..\covid-sim -G Ninja`

Build CovidSim: `cmake --build .`

Run tests (method 1): `cmake --build . --target tests`

Run tests (method 2): `ctest -j 2 -V`

The second method of running tests may not work depending on your system
config - but if it does it will run the tests in parallel and output more
information whilst the tests are run.

## Accepting changes to test outputs

To accept the changes to test output.  Open a Developer Command prompt in the
root covid-sim source directory and do the following:

```sh
cd ..
md test-output-build
cd test-output-build
cmake ..\covid-sim -G Ninja
cmake --build .
..\covid-sim\tests\regressiontest_UK_100th.py --accept --output uk-output --covidsim src\CovidSim.exe
..\covid-sim\tests\regressiontest_US_based.py --accept --output us-output --covidsim src\CovidSim.exe
```

These will update the expected checksums which should be committed to the git repo in the normal way.
