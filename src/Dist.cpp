#include "Dist.h"
#include "Param.h"

#include "Model.h"

double dist2(Person* a, Person* b)
{
	return P.distance_->distance_squared(Households[a->hh].loc, Households[b->hh].loc);
}

double dist2_cc(Cell* a, Cell* b)
{
	int l = (int)(a - Cells);
	int m = (int)(b - Cells);

	return P.distance_->distance_squared(l, m, P.nch, P.in_cells_);
}

double dist2_cc_min(Cell* a, Cell* b)
{
	int l = (int)(a - Cells);
	int m = (int)(b - Cells);

	return P.distance_->distance_squared(l, m, P.nch, P.in_cells_, true);
}

double dist2_mm(Microcell* a, Microcell* b)
{
	int l = (int)(a - Mcells);
	int m = (int)(b - Mcells);

	return P.distance_->distance_squared(l, m, P.total_microcells_high_, P.in_microcells_);
}

double dist2_raw(double ax, double ay, double bx, double by)
{
	return P.distance_->distance_squared(CovidSim::Geometry::Vector2d(ax, ay), CovidSim::Geometry::Vector2d(bx, by));
}
