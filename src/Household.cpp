#include "Household.hpp"

float Household::distance_to(Household other) {
	return this->loc.distance_to(other.loc);
}

float Household::distance_squared_to(Household other) {
	return this->loc.distance_squared_to(other.loc);
}
