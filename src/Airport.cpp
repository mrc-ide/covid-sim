#include "Airport.h"

double Airport::distance_to(Airport *other) const {
	return this->loc.distance_to(other->loc);
}

double Airport::distance_squared_to(Airport* other) const {
	return this->loc.distance_squared_to(other->loc);
}
