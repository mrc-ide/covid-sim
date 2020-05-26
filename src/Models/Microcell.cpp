#include <cstdlib>
#include "Microcell.h"
#include "Param.h"
#include "Model.h"

using namespace Models;
using namespace Geometry;

int Microcell::index() const {
	return (int)(this - Mcells);
}

Vector2<double> Microcell::position() const {
	int x = std::abs(this->index() / P.number_of_micro_cells().height);
	int y = std::abs(this->index() % P.number_of_micro_cells().height);
	return Vector2<double>(x, y) * (Vector2<double>)P.in_microcells_;
}

double Microcell::distance_to_squared(Microcell *other) const {
	return this->position().distance_to_squared(other->position());
}
