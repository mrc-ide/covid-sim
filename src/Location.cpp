#include <cmath>

#include "Location.hpp"
#include "Dist.h"

Location::Location() :
	x(0),
	y(0)
{}

Location::Location(float x, float y) :
	x(x),
	y(y)
{}

Location::Location(double x, double y) :
	x((float)x),
	y((float)y)
{}

float Location::distance_to(Location other) const {
	return std::sqrt(this->distance_squared_to(other));
}

float Location::distance_squared_to(Location other) const {
	return (float)dist2_raw(this->x, this->y, other.x, other.y);
}
