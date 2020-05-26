#include "Person.h"

#include "../Model.h"

double Person::distance_to(const Person &other) const {
	return std::sqrt(this->distance_squared_to(other));
}

double Person::distance_squared_to(const Person &other) const {
	return Households[this->hh].distance_squared_to(Households[other.hh]);
}
