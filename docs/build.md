# Overview of the build system

The COVID-19 CovidSim model uses the [CMake](www.cmake.org) build tool to
generate build files for other build systems. Currently, building using clang
and gcc with Makefiles and MSVC with Visual Studio are supported.

## Building with Makefiles

From the command line inside a git clone, run the following:

```sh
mkdir build
cd build
cmake ..
make
```

### Testing

Once `make` has completed use:

```sh
make test
```

to run the regressions tests.

The tests can take a while, and may produce no output for >10 minutes.
Therefore using:

```sh
make test ARGS="-V"
```

or

```sh
ctest -V
```

May be more reassuring that something is happening.

## Building with Visual Studio project files from Cmake

From the command line inside a git clone, run the following:

```sh
mkdir build
cd build
cmake ..
```

This will create project files inside the `build` directory that can be opened
in Visual Studio. Modifications to the CMake configuration may require
regenerating the Visual Studio projects.

### Testing

To enable Visual Studio to pick up the tests added by CMake:

 * Open the `Tests` menu in the Visual Studio menu bar

 * Choose `Run CTests`.

## Building directly with CMake in Visual Studio 2019

Visual Studio 2019 supports using CMake to manage the build directly by
selecting File -> Open -> Cmake... and opening `src/CMakeLists.txt`. Then
Visual Studio's normal build shortcuts will update the CMake configuration
as well as building the project.

CMake options are configured using the `CMakeSettings.json` file, which
Visual Studio will generate when `CMakeLists.txt` is opened.

### Testing

To enable Visual Studio to pick up the tests added by CMake:

 * Open the `Tests` menu in the Visual Studio menu bar

 * Choose `Run CTests`.

### Build options

Additional configuration variables can be provided in the `cmake` invocation.

- `USE_OPENMP` determines whether the model is compiled with parallelization
using OpenMP. This option defaults to on, but can be disabled by passing
`-DUSE_OPENMP=OFF`. The simulation is designed to be run on multi-core systems.
Performance improvements are approximately linear up to 24 to 32 cores,
depending on memory performance.

For Makefile builds, use `-DCMAKE_BUILD_TYPE=` to specify the output format:

- `Debug`
- `MinSizeRel`
- `Release`
- `RelWithDebInfo`

By default, Makefile builds will use `RelWithDebInfo`.

## VisualStudio solution

A manually created VS-2019 solution and project is included for convenience,
but it should not be considered the source of truth for the project.

## Testing

The regression tests are not supported in the Visual Studio stand-alone
solution.
