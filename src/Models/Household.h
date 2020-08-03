#ifndef COVIDSIM_MODELS_HOUSEHOLD_H_INCLUDED_
#define COVIDSIM_MODELS_HOUSEHOLD_H_INCLUDED_

#include "geometry/Vector2.h"

struct Household
{
	int FirstPerson;
	unsigned short int nh; // number people in household
	unsigned short int nhr;
	CovidSim::Geometry::Vector2f loc;
};

#endif
