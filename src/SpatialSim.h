#pragma once

#ifndef SPATIALSIM_SPATIALSIM_H_INCLUDED_
#define SPATIALSIM_SPATIALSIM_H_INCLUDED_

#include "MachineDefines.h"

#define KMP_LIBRARY throughput

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

  /*
	#define HOST_TREATED(x) (0)
	#define HOST_TO_BE_TREATED(x) (0)
	#define PLACE_TREATED(x,y) (0)
  */
#define FAST_US

//// place type codes
#define PlaceType_PrimarySchool			0
#define PlaceType_SecondarySchool		1
#define PlaceType_University			2
#define PlaceType_Office				3

#endif // SPATIALSIM_SPATIALSIM_H_INCLUDED_
