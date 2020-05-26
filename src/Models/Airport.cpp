#include "Airport.h"

double Airport::distance_to(const Airport &other) const {
	return this->loc.distance_to(other.loc);
}

double Airport::distance_to_squared(const Airport &other) const {
	return this->loc.distance_to_squared(other.loc);
}

double Airport::distance_to(const Place &other) const {
	return this->loc.distance_to(other.loc);
}
double Airport::distance_to_squared(const Place &other) const {
	return this->loc.distance_to_squared(other.loc);
}

