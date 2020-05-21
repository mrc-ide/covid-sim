#pragma once

#include "Coordinates/Location.hpp"

struct Household
{
	int FirstPerson;
	unsigned short int nh; // number people in household
	Location loc;
	unsigned short int nhr;

	float distance_to(Household other);
	float distance_squared_to(Household other);
};