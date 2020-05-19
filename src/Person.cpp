#include "Model.h"
#include <cmath>

bool Person::isAlive() const {
	return !this->isDead();
}

bool Person::isDead() const {
	return std::abs(this->inf) == InfStat_Dead;
}