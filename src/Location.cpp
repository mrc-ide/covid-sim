#include <cmath>

#include "Location.hpp"
#include "Dist.h"

float Location::distance_to(Location other) const {
	return std::sqrt(this->distance_squared_to(other));
}

float Location::distance_squared_to(Location other) const {
	return dist2_raw(this->x, this->y, other.x, other.y);
}
