#include "Cell.h"
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


double dist2_cc(Cell* a, Cell* b)
{
	if (P.DoUTM_coords) {
		return a->position().distance_squared_to(b->position());
	}
	else
	{
		Vector2<double> delta = a->position() - b->position();
		if (P.DoPeriodicBoundaries)
		{
			if (delta.x > P.in_degrees_.width * 0.5) delta.x = P.in_degrees_.width - delta.x;
			if (delta.y > P.in_degrees_.height * 0.5) delta.y = P.in_degrees_.height - delta.y;
		}
		return delta.length_squared();
	}
}