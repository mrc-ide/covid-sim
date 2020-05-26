#include "Household.h"

using namespace Models;

float Household::distance_to(const Household &other) {
	return this->loc.distance_to(other.loc);
}

float Household::distance_to_squared(const Household &other) {
	return this->loc.distance_to_squared(other.loc);
}

float Household::distance_to(const Place &other) {
	return this->loc.distance_to(other.loc);
}

float Household::distance_to_squared(const Place &other) {
	return this->loc.distance_to_squared(other.loc);
}