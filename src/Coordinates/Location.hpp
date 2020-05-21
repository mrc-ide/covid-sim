#pragma once

#include "Vector2.hpp"

struct Location : Vector2<float> {
	Location(float x, float y) : Vector2<float>(x, y){};

	float distance_to(Location other) const;
	float distance_squared_to(Location other) const;
};

