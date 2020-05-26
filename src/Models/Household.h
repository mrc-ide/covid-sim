#pragma once

#include "Place.h"

namespace Models {
	struct Household {
		int FirstPerson;
		unsigned short int nh; // number people in household
		Geometry::Vector2<float> loc;
		unsigned short int nhr;

		float distance_to(const Household &other);

		float distance_to_squared(const Household &other);

		float distance_to(const Place &other);

		float distance_to_squared(const Place &other);
	};
}