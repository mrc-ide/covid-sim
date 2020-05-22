#include "Cell.hpp"
#include "Model.h"
#include "Param.h"
#include <cmath>

int Cell::index() const {
	return (int)(this - Cells);
}

Vector2<double> Cell::position() const {
	int x = std::abs(this->index() / P.number_of_cells.height);
	int y = std::abs(this->index() % P.number_of_cells.height);
	return Vector2<double>(x, y) * (Vector2<double>)P.in_cells_;
}
