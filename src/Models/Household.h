#pragma once

struct Household
{
	int FirstPerson;
	unsigned short int nh; // number people in household
	unsigned short int nhr;
	Geometry::Vector2<float> loc;
};