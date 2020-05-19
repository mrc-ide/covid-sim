#include "Model.h"
#include <cstdlib>

bool Person::isAlive() const {
	return !this->isDead();
}

bool Person::isDead() const {
	return std::abs(this->inf) == InfStat_Dead;
}