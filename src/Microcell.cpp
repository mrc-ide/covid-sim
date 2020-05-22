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
	if (P.DoUTM_coords){
		return this->position().distance_squared_to(other->position());
	}
	else
	{
		Vector2<double> delta = (this->position() - other->position()).abs();
		if (P.DoPeriodicBoundaries)
		{
			if (delta.x > P.in_degrees_.width * 0.5) delta.x = P.in_degrees_.width - delta.x;
			if (delta.y > P.in_degrees_.height * 0.5) delta.y = P.in_degrees_.height - delta.y;
		}
		return delta.length_squared();
	}
}
