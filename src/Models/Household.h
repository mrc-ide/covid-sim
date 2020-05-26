#pragma once

#include "Place.h"

struct Household
{
	int FirstPerson;
	unsigned short int nh; // number people in household
	Vector2<float> loc;
	unsigned short int nhr;

	float distance_to(const Household &other);
	float distance_squared_to(const Household &other);

	float distance_to(const Place &other);
	float distance_squared_to(const Place &other);
};