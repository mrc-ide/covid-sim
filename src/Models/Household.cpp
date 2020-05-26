#include "Household.h"

float Household::distance_to(const Household &other) {
	return this->loc.distance_to(other.loc);
}

float Household::distance_squared_to(const Household &other) {
	return this->loc.distance_squared_to(other.loc);
}

float Household::distance_to(const Place &other) {
	return this->loc.distance_to(other.loc);
}

float Household::distance_squared_to(const Place &other) {
	return this->loc.distance_squared_to(other.loc);
}