# Overview of the build system

The COVID-19 CovidSim model uses the [CMake](www.cmake.org) build tool to
generate build files for other build systems. Currently, building using clang
and gcc with Makefiles and MSVC with Visual Studio are supported.

## Building from the Command Line (all platforms)

### Build requirements

To build `CovidSim` from the command line you need to have the following tools
and libraries installed:

 * [git](https://git-scm.com/) Source code control.
 * [Cmake >=3.8](https://www.cmake.org): Controls the build
 * Make system (GNU Make, Ninja Build, MSBuild all have been tested)
 * Compiler and related tools (GCC >=7, Clang >=10, MSVC >=16.6 have all been
   tested)
 * OpenMP libraries: Strongly recommended, enables multi-threading.
 * [Python >=3.6](https://www.python.org): Used to drive the tests.
 * [Doxygen](https://www.doxygen.nl/index.html): Optional, used to generate
   documentation.
 * [GraphViz](https://graphviz.org/): Optional, used by doxygen to generate
   images in the documentation.

### Cloning the repo

From the command line:

```sh
git clone https://github.com/mrc-ide/covid-sim
```

### Configuring

Configuring is slightly different depending on your platform

#### Unix-like (Linux, macOS etc...)

```sh
cd covid-sim
mkdir build
cd build
cmake ..
```

#### Windows

On Windows it is recommended to use Ninja Build as the make system instead of
Visual Studio.  This is to simplify the instructions following, and because
CMake's support for Ninja Build is more feature rich.  Ninja Build is
installed with CMake in Visual Studio.

```sh
cd covid-sim
md build
cd build
cmake .. -G Ninja
```

#### Configure options

Additional configuration variables can be provided in the `cmake` invocation.

##### OpenMP

`USE_OPENMP` determines whether the model is compiled with parallelization
using OpenMP. This option defaults to on, but can be disabled by passing
`-DUSE_OPENMP=OFF`. The simulation is designed to be run on multi-core systems.
Performance improvements are approximately linear up to 24 to 32 cores,
depending on memory performance.

##### Build type

For Makefile builds, use `-DCMAKE_BUILD_TYPE=` to specify the output format:

- `Debug`
- `MinSizeRel`
- `Release`
- `RelWithDebInfo`

By default, Makefile builds will use `RelWithDebInfo`.  To build a Debug build
you may invoke `cmake` like the following:

```sh
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

##### Different build system

You may use a different build system to the default by using the `-G
<generator>` option to CMake.

In particular for Windows it may be useful to use the Ninja generator instead
of the Windows one.  For example:

```sh
cmake .. -G Ninja
```

### Building

```sh
cmake --build . --target all
```

If using the default Windows build system (MSBuild) you may also need to
specify the config in use with the `--config <cfg>` command line option.

To specify how many cores the build system can use when building the tools use
the `-j N` option, where `N` is the number of CPUs to use.  For example:

```sh
cmake --build . --target all -j 6
```

### Testing

To run the tests do:

```sh
ctest
```

The tests can take a while, and may produce no user-visible output for >10
minutes.  Therefore, it may be worthwhile to increase the verbosity of the
tests:

```sh
ctest -V
```

The tests may be run in parallel by specifying `-j N` to the command line.  For
example:

```sh
ctest -V -j 6
```

### Accepting changes in test results

If the output of `CovidSim` has changed such that the tests start to fail and
you are happy that the changes in output are acceptable then the make target
`test-accept` will update the expected test results.

To accept the test changes do something like the following in the build directory:

```sh
cmake --build . --target test-accept -j 6
git add tests/*-input/results-j*.cksum
git commit -m"Update expected results."
```

### Generating Doxygen docs

If `doxygen` is installed there will be a `doxygen` build target that builds
the documentation.  Build the documentation as follows:

```sh
cmake --build . --target doxygen
```

## Building with Visual Studio project files from Cmake

Visual Studio 2019 supports using CMake to manage the build directly by
selecting File -> Open -> Cmake... and opening `CMakeLists.txt`. Then
Visual Studio's normal build shortcuts will update the CMake configuration
as well as building the project.

CMake options are configured using the `CMakeSettings.json` file, which
Visual Studio will generate when `CMakeLists.txt` is opened.

### Testing

To enable Visual Studio to pick up the tests added by CMake:

 * Open the `Tests` menu in the Visual Studio menu bar

 * Choose `Run CTests`.

Because of the way Visual Studio runs the tests they can only be run serially.
It is recommended that you use the command line for running tests.

### Accepting changes in test results

If the output of `CovidSim` has changed such that the tests start to fail and
you are happy that the changes in output are acceptable then the make target
`test-accept` will update the expected test results.  To do this in Visual
Studio do:

 * From the menu bar "View" menu choose "Solution Explorer"

 * At the top of the "Solution Explorer" window is a drop down menu icon
   (window and folder with circling arrows).  From this menu choose "CMake
   Targets View".

 * In the tree that appears expand: "covid-sim", "CovidSim (executable)".

 * Right click on "test-accept (utility target)".

 * In the menu that appears click "Build".

## Visual Studio solution

A manually created VS-2019 solution and project is included for convenience,
but it should not be considered the source of truth for the project.

### Testing

The regression tests are not supported in the Visual Studio stand-alone
solution.

## Doxygen

The CMake custom target `doxygen` builds the Doxygen documentation
(`make doxygen`). Some IDEs do not expose this so alternatively setting the
`BUILD_DOC` CMake option to `ON` will cause the documentation to be built
every time the software is built regardless of whether anything has changed.

To generate any diagrams [GraphViz](https://graphviz.org/) must be installed.

### Replaying warnings in Visual Studio

The Doxygen warnings can be navigated to in Visual Studio:

1. From the `Tools` menu select `External Tools...` to launch the
   `External Tools` dialogue box.

2. Press the `Add` button.

3. For `Title` enter `CovidSim`.

4. For `Command` enter `cmd.exe`.

5. For `Arguments` enter `/c type \path\to\warnings.txt`.

6. Check `Use Output window`.

7. Click `Apply`.

8. From the `Tools` menu select `CovidSim` to replay the file.

9. Double click on a line in the `Output` window which begins with a filename
   to jump to that file.
