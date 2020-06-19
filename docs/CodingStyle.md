# Coding Style

This document describes the coding style for new additions to the CovidSim project.

Coding Styles are always a controversial topic, and there is no "one size"
fits all.  So the following are only guidelines, and sometimes breaking them
is justified.

In particular these apply to new code and not to pre-existing code.

## Enforcement

As far as possible the coding style is automatically applied by
[.editorconfig](../.editorconfig) and [.clang-format](../.clang-format)
formatting files.

Many editors will use these either automatically or through plugins.

The CI loop runs `clang-format-9` over the code and CI will fail if it finds
any changes to make.

## C++ Coding Style

Basic style is C++14, Allman, mostly Microsoft Visual Studio defaults.  See
the [.clang-format](../.clang-format) file for what we use.

Whitespace: lines are 100 characters long, with tabs used for indents,
tabstops at 2 spaces.

## Commenting

All added structures, classes, functions should be documented using
[doxygen](https://www.doxygen.nl/index.html) formats.

In particular use:

 * `/** ... */` for multiline documentation comments.  Each continuation line
   to be preceded by a `*`.

 * `///< ...` for single line documentation comments after declarations.

 * `/// ...` for single line documentation comments.

 * In documentation use `\` as the prefix for commands not `@`.  (so `\brief`
   not `@brief`.

 * Comments that aren't intended for the doxygen documentation should follow
   similar rules but use `/* ... */` and `// ...`

 * All comments should be full sentences ending with a full-stop.

## Naming Conventions

 * Header files: \*.h
 * Source files: \*.cpp
 * Namespaces: CamelCase
 * Classes, Structures, Class Enums, Enum Members: CamelCase
 * Template Parameters: Short CAPS (So T, or IT not ITERATOR or Iterator)
 * Function parameters and variables: snake\_case
 * Class, structure member variables: snake\_case\_trailing\_underscore\_

## Language Features

 * C++14 + OpenMP standard features only.

 * Follow the
   [C++ Core Guidelines](http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
   as far as possible.

 * Header files should contain a block of functionality.  In particular this
   means that in general more than one-class per header file.  (Header files
   are similar granularity to namespaces).

 * Header files should include the minimum number of `#include`s necessary,
   and use forward declarations wherever possible.

 * Source files should include the appropriate headers, and not depend on
   headers including other headers.

 * All new classes/functions should be placed in a namespace.

 * `using namespace` should not be used to import symbols into the global
   namespace.

 * No exceptions in code that might be run in a multi-threaded environment.
   Excpetions and OpenMP are hard to get to play well with each other.

 * Don't use explicit `this` unless required.

 * Functions and methods should be marked:

   * noexcept when they don't throw exceptions

   * (Methods) const when they don't modify state.
