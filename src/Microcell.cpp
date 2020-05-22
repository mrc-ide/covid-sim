#include <cmath>
#include "Microcell.h"
#include "Param.h"
#include "Model.h"

int Microcell::index() const {
	return (int)(this - Mcells);
}

Vector2<double> Microcell::position() const {
	int x = std::abs(this->index() / P.get_number_of_micro_cells_high());
	int y = std::abs(this->index() % P.get_number_of_micro_cells_high());
	return Vector2<double>(x, y) * (Vector2<double>)P.in_microcells_;
}

double Microcell::distance_squared_to(Microcell *other) const {
	return this->position().distance_squared_to(other->position());
}
