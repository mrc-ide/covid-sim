#ifndef COVIDSIM_COVIDSIM_H_INCLUDED_
#define COVIDSIM_COVIDSIM_H_INCLUDED_

// Compiling with Clang from VS requires libomp.lib to be included in
// Project, Properties, Input, Additional Dependencies, however
// that causes VS/Intel compilers to fail to build if libomp.lib is
// not found. The following should include the lib file only on Clang.

#ifdef _WIN32
  #ifdef __clang__
#pragma comment(lib, "libomp.lib")
  #endif
#endif

#include "MachineDefines.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <time.h>

#include "Country.h"
#include "Constants.h"

#endif // COVIDSIM_COVIDSIM_H_INCLUDED_
