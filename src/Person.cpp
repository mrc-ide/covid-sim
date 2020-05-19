#include "Model.h"
#include <cstdlib>

bool Person::is_alive() const {
	return !this->is_dead();
}

bool Person::is_dead() const {
	return std::abs(this->inf) == InfStat_Dead;
}