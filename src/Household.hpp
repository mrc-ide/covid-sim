#pragma once

#include "Location.hpp"

struct Household
{
	int FirstPerson;
	unsigned short int nh; // number people in household
	Location loc;
	unsigned short int nhr;
};