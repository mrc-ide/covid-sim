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

double Cell::distance_squared_to(Cell *other) const
{
	if (P.DoUTM_coords) {
		return a->position().distance_squared_to(other->position());
	}
	else
	{
		Vector2<double> delta = this->position() - other->position();
		if (P.DoPeriodicBoundaries)
		{
			if (delta.x > P.in_degrees_.width * 0.5) delta.x = P.in_degrees_.width - delta.x;
			if (delta.y > P.in_degrees_.height * 0.5) delta.y = P.in_degrees_.height - delta.y;
		}
		return delta.length_squared();
	}
}


double Cell::distance_squared_to_min(Cell* other) const
{
	int l = (int)(this - Cells);
	int m = (int)(other - Cells);
	int i = l;
	int j = m;

	if (P.DoUTM_coords)
	{
		if (P.in_cells_.width_ * ((double)abs(m / P.nch - l / P.nch)) > PI)
		{
			if (m / P.nch > l / P.nch)
				j += P.nch;
			else if (m / P.nch < l / P.nch)
				i += P.nch;
		}
		else
		{
			if (m / P.nch > l / P.nch)
				i += P.nch;
			else if (m / P.nch < l / P.nch)
				j += P.nch;
		}
		if (m % P.nch > l % P.nch)
			i++;
		else if (m % P.nch < l % P.nch)
			j++;
		return dist2UTM(P.in_cells_.width_ * fabs((double)(i / P.nch)), P.in_cells_.height_ * fabs((double)(i % P.nch)),
		                P.in_cells_.width_ * fabs((double)(j / P.nch)), P.in_cells_.height_ * fabs((double)(j % P.nch)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.in_cells_.width_ * ((double)abs(m / P.nch - l / P.nch)) > P.in_degrees_.width_ * 0.5))
		{
			if (m / P.nch > l / P.nch)
				j += P.nch;
			else if (m / P.nch < l / P.nch)
				i += P.nch;
		}
		else
		{
			if (m / P.nch > l / P.nch)
				i += P.nch;
			else if (m / P.nch < l / P.nch)
				j += P.nch;
		}
		if ((P.DoPeriodicBoundaries) && (P.in_degrees_.height_ * ((double)abs(m % P.nch - l % P.nch)) > P.in_degrees_.height_ * 0.5))
		{
			if (m % P.nch > l % P.nch)
				j++;
			else if (m % P.nch < l % P.nch)
				i++;
		}
		else
		{
			if (m % P.nch > l % P.nch)
				i++;
			else if (m % P.nch < l % P.nch)
				j++;
		}
		double x = P.in_cells_.width_ * fabs((double)(i / P.nch - j / P.nch));
		double y = P.in_cells_.height_ * fabs((double)(i % P.nch - j % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width_ * 0.5) x = P.in_degrees_.width_ - x;
			if (y > P.in_degrees_.height_ * 0.5) y = P.in_degrees_.height_ - y;
		}
		return x * x + y * y;
	}
}