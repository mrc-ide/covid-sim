#include "Cell.h"
#include "Model.h"
#include "Param.h"
#include "Dist.h"
#include <cmath>

using namespace Models;
using namespace Geometry;

int Cell::index() const {
	return (int)(this - Cells);
}

Vector2<double> Cell::position() const {
	int x = std::abs(this->index() / P.number_of_cells.height);
	int y = std::abs(this->index() % P.number_of_cells.height);
	return Vector2<double>(x, y) * (Vector2<double>)P.in_cells_;
}

double Cell::distance_to(Cell *other) const
{
	return std::sqrt(this->distance_to_squared(other));
}

double Cell::distance_to_squared(Cell *other) const
{
	return this->position().distance_to_squared(other->position());
}

double Cell::distance_to_squared_min(Cell* other) const
{
	int l = (int)(this - Cells);
	int m = (int)(other - Cells);
	int i = l;
	int j = m;

	if (P.DoUTM_coords)
	{
		if (P.in_cells_.width * ((double)std::abs(m / P.number_of_cells.height - l / P.number_of_cells.height)) > PI)
		{
			if (m / P.number_of_cells.height > l / P.number_of_cells.height)
				j += P.number_of_cells.height;
			else if (m / P.number_of_cells.height < l / P.number_of_cells.height)
				i += P.number_of_cells.height;
		}
		else
		{
			if (m / P.number_of_cells.height > l / P.number_of_cells.height)
				i += P.number_of_cells.height;
			else if (m / P.number_of_cells.height < l / P.number_of_cells.height)
				j += P.number_of_cells.height;
		}
		if (m % P.number_of_cells.height > l % P.number_of_cells.height)
			i++;
		else if (m % P.number_of_cells.height < l % P.number_of_cells.height)
			j++;
		return dist2UTM(P.in_cells_.width * std::fabs((double)(i / P.number_of_cells.height)),
						P.in_cells_.height * std::fabs((double)(i % P.number_of_cells.height)),
		                P.in_cells_.width * std::fabs((double)(j / P.number_of_cells.height)),
		                P.in_cells_.height * std::fabs((double)(j % P.number_of_cells.height)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.in_cells_.width * ((double)std::abs(m / P.number_of_cells.height - l / P.number_of_cells.height)) > P.in_degrees_.width * 0.5))
		{
			if (m / P.number_of_cells.height > l / P.number_of_cells.height)
				j += P.number_of_cells.height;
			else if (m / P.number_of_cells.height < l / P.number_of_cells.height)
				i += P.number_of_cells.height;
		}
		else
		{
			if (m / P.number_of_cells.height > l / P.number_of_cells.height)
				i += P.number_of_cells.height;
			else if (m / P.number_of_cells.height < l / P.number_of_cells.height)
				j += P.number_of_cells.height;
		}
		if ((P.DoPeriodicBoundaries) && (P.in_degrees_.height * ((double)std::abs(m % P.number_of_cells.height - l % P.number_of_cells.height)) > P.in_degrees_.height * 0.5))
		{
			if (m % P.number_of_cells.height > l % P.number_of_cells.height)
				j++;
			else if (m % P.number_of_cells.height < l % P.number_of_cells.height)
				i++;
		}
		else
		{
			if (m % P.number_of_cells.height > l % P.number_of_cells.height)
				i++;
			else if (m % P.number_of_cells.height < l % P.number_of_cells.height)
				j++;
		}
		double x = P.in_cells_.width * std::fabs((double)(i / P.number_of_cells.height - j / P.number_of_cells.height));
		double y = P.in_cells_.height * std::fabs((double)(i % P.number_of_cells.height - j % P.number_of_cells.height));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.in_degrees_.width * 0.5) x = P.in_degrees_.width - x;
			if (y > P.in_degrees_.height * 0.5) y = P.in_degrees_.height - y;
		}
		return x * x + y * y;
	}
}