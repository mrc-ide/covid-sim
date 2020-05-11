# Visual Studio Solutions

`covid-sim.sln` exists because certain Visual Studio tools do not work with
CMake builds. Third party git repositories are used by `covid-sim` and are,
unfortunately, only downloaded as part of a CMake build. Hence it is
required to open the parent directory in Visual Studio and perform a
CMake release build before building `covid-sim.sln`.
