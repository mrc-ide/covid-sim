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

#include <climits>
#include <cstdlib>
#include <cstddef>

#include <cmath>
#include <cstring>
#include <csignal>
#include <ctime>

#include "Country.h"
#include "Constants.h"
#include "Dist.h"
#include "Files.h"
#include "ReadParams.h"

int ChooseTriggerVariableAndValue(int AdUnit);
double ChooseThreshold(int AdUnit, double WhichThreshold);
void UpdateCurrentInterventionParams(double t);
void UpdateCFRs(double t_CalTime);
void RecordAdminAgeBreakdowns(int t_int);
void RecordQuarNotInfected(int n, unsigned short int TimeStepNow);

#endif // COVIDSIM_COVIDSIM_H_INCLUDED_
