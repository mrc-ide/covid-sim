#include "Person.h"

#include "../Model.h"

double Person::distance_to(const Person &other) const {
	return std::sqrt(this->distance_to_squared(other));
}

double Person::distance_to_squared(const Person &other) const {
	return Households[this->hh].distance_to_squared(Households[other.hh]);
}
